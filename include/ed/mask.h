#ifndef ED_MASK_H_
#define ED_MASK_H_

#include <rgbd/types.h>
#include <opencv2/core/core.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>

namespace ed
{

/** Support for scanning an image as a number of sub-image scans. */
class ImageMask {

public:

    ImageMask() {}

    /**
     * Construct the mask, while setting the mask size.
     * @paran width Width of the mask.
     * @param height Height of the mask.
     */
    ImageMask(int width, int height) : width_(width), height_(height)
    {
    }

    /**
     * Set the size of the image.
     * @paran width Width of the mask.
     * @param height Height of the mask.
     */
    inline void setSize(int width, int height)
    {
        width_ = width;
        height_ = height;
    }

    /** Remove all sub-images. */
    inline void clear()
    {
        points_.clear();
    }

    /**
     * Get the number of sub-images.
     * @return The number of sub-images in the mask.
     */
    inline int getSize() const { return points_.size(); }

    /**
     * Get the width of the mask.
     * @return The width of the mask.
     */
    inline int width() const { return width_; }

    /**
     * Get the height of the mask.
     * @return The height of the mask.
     */
    inline int height() const { return height_; }

    /**
     * Add a sub-image.
     * @param p Base-point of the new sub-image.
     */
    inline void addPoint(const cv::Point2i& p)
    {
         points_.push_back(p);
    }

    /**
     * Add a sub-image.
     * @param x X coordinate of the new sub-image.
     * @param y Y coordinate of the new sub-image.
     */
    inline void addPoint(int x, int y)
    {
        addPoint(cv::Point2i(x, y));
    }

    /**
     * Add a one-pixel sub-image.
     * @param idx Index number of the pixel (scanning horizontally, from top to bottom).
     */
    inline void addPoint(int idx)
    {
        addPoint(idx % width_, idx / width_);
    }

    /**
     * Add a number of sub-images.
     * @param ps Base points of the new sub-images.
     */
    inline void addPoints(const std::vector<cv::Point2i>& ps)
    {
        for(std::vector<cv::Point2i>::const_iterator it = ps.begin(); it != ps.end(); ++it)
        {
            addPoint(*it);
        }
    }

    /**
     * Iterator that scans all sub-images in turn, scanning each sub-image from
     * left to right, from top to bottom (assuming base position of a sub-image
     * is its top-left corner).
     */
    class const_iterator
    {
    public:
        /**
         * Iterator constructor.
         * @param ptr Base points of the sub-images.
         * @param factor Size of a rectangular sub-image.
         */
        const_iterator(const cv::Point2i* ptr, int factor)
            : ptr_(ptr), p_(ptr_->x * factor, ptr_->y * factor), dx_(0), dy_(0), factor_(factor)
        {
        }

        // post increment operator
        inline const_iterator operator++(int)
        {
            const_iterator i(*this);
            ++*this;
            return i;
        }

        /**
         * Pre-increment operator.
         * Increment the X position of the scan line, wrapping to the start of
         * the next scan line, jump to the next sub-image at the end of the
         * current sub-image.
         */
        inline const_iterator& operator++()
        {
            ++dx_;

            // If reached end of the scan line of a sub-image, jump to next scan line.
            if (dx_ == factor_)
            {
                dy_++;
                dx_ = 0;
            }

            // Reached underneath the last scan line, jump back to the start.
            if (dy_ == factor_)
            {
                dx_ = 0;
                dy_ = 0;
                ++ptr_;
            }

            // Compute point of the new position.
            p_ = cv::Point2i(ptr_->x * factor_ + dx_, ptr_->y * factor_ + dy_);

            return *this;
        }

        inline const cv::Point2i& operator*() { return p_; }
        inline const cv::Point2i* operator->() { return &p_; }

        inline bool operator==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
        inline bool operator!=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
    private:

        const cv::Point2i* ptr_;
        cv::Point2i p_; ///< Base position of the current sub-image.
        int dx_, dy_; ///< Variables tracking the x/y position in the current sub-image.
        int factor_; ///< Sub-image X/Y size (sub-image is rectangular).
    };

    const_iterator begin(int width = 0) const
    {
        if (width <= 0)
            return const_iterator(&points_.front(), 1);
        else if (points_.empty())
            return end();
        else
            return const_iterator(&points_.front(), width / width_);
    }

    const_iterator end() const
    {
        return const_iterator(&points_.back() + 1, 0);
    }

private:

    int width_;  ///< Width of the mask.
    int height_; ///< Height of the mask.
    std::vector<cv::Point2i> points_; ///< Base points of the sub-images.
};

} // end namespace ed

#endif
