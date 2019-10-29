#ifndef ED_BOUNDED_BUFFER_H_
#define ED_BOUNDED_BUFFER_H_

#include <boost/circular_buffer.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/thread.hpp>
#include <boost/call_traits.hpp>
#include <boost/progress.hpp>
#include <boost/bind.hpp>

//#include <boost/timer/timer.hpp> // for auto_cpu_timer

template <class T>
class bounded_buffer {
public:
   typedef boost::circular_buffer<T> container_type;
   typedef typename container_type::size_type size_type;
   typedef typename container_type::value_type value_type;
   typedef typename boost::call_traits<value_type>::param_type param_type;

   explicit bounded_buffer(size_type capacity) : m_unread(0), m_container(capacity) {}

  void push_front(param_type item) 
  {
      // param_type represents the "best" way to pass a parameter of type value_type to a method
      boost::mutex::scoped_lock lock(m_mutex); // this already locks the mutex
      
      if(!this->is_not_full()) // if buffer full, remove oldest item
      {
              //value_type* pItem;
              //*pItem = 
              m_container[--m_unread]; 
              m_not_full.notify_one();
      }

  //    m_not_full.wait(lock, boost::bind(&bounded_buffer<value_type>::is_not_full, this)); // OLD
      
      m_container.push_front(item);
      ++m_unread;
      lock.unlock();
      m_not_empty.notify_one();
   }

   void pop_back(value_type* pItem) 
   {
      boost::mutex::scoped_lock lock(m_mutex);
      m_not_empty.wait(lock, boost::bind(&bounded_buffer<value_type>::is_not_empty, this));
      *pItem = m_container[--m_unread];
      lock.unlock();
      m_not_full.notify_one();
   }
   
   bool is_not_empty() const { return m_unread > 0; }
   bool is_not_full() const {  return m_unread < m_container.capacity(); }
   
   int capacity() const {return m_container.capacity();}
   int numberUnread() const {return m_unread;}

private:
   bounded_buffer(const bounded_buffer&);              // Disabled copy constructor
   bounded_buffer& operator = (const bounded_buffer&); // Disabled assign operator

   size_type m_unread; // TODO currently, pop_back and m_unread allow for a single consumer only
   container_type m_container;
   boost::mutex m_mutex;
   boost::condition m_not_empty;
   boost::condition m_not_full;
};

#endif