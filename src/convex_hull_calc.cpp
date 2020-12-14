#include "ed/convex_hull_calc.h"

#include <opencv2/imgproc/imgproc.hpp>

namespace ed
{

namespace tracking
{

FITTINGMETHOD determineCase ( std::vector<geo::Vec2f>& points, unsigned int* cornerIndex, std::vector<geo::Vec2f>::iterator* it_low, std::vector<geo::Vec2f>::iterator* it_high, const geo::Pose3D& sensor_pose )
{
    // Determine is a line or a rectangle should be fitted. In case of a rectangle, the number of elements for both sides should meet the minimum number of points for a line fit
    // for both lines. Otherwise, a line will be fitted on the remaining points.

    *it_low = points.begin();
    *it_high = points.end();

    bool includeCorner = *cornerIndex > 0;

    // In the case of including a corner, check if both segments have enough points to describe a line with. If not, do not use these data.
    if ( includeCorner ) {
        unsigned int nPointsLow = *cornerIndex + 1; // + 1 because the corner can be used for both lines
        unsigned int nPointsHigh = points.size() - nPointsLow + 1; // +1 because the point with max error is considered as the corner point and belongs to both lines
        unsigned int remainingSize = points.size();

        bool fitSingleline = false;
        bool pointsRemoved = false;

        if ( nPointsLow < MIN_POINTS_LINEFIT ) { // Part of section too smal -> remove it from the data which are analyzed and try to fit line again
            *it_low += *cornerIndex;
            remainingSize -= ( nPointsLow - 1 );
            pointsRemoved = true;
        }

        if ( nPointsHigh < MIN_POINTS_LINEFIT ) {
            *it_high -= *cornerIndex;
            remainingSize -= ( nPointsHigh - 1 );
            pointsRemoved = true;
        }

        if ( pointsRemoved && remainingSize < MIN_POINTS_LINEFIT ) {
            *cornerIndex = std::numeric_limits<unsigned int>::quiet_NaN();
            return NONE;
        } else if ( pointsRemoved && remainingSize >= MIN_POINTS_LINEFIT ) {
            *cornerIndex = std::numeric_limits<unsigned int>::quiet_NaN();
            return LINE;
        } else  { // we dit not remove points and a corner is present
            return RECTANGLE;
        }
    } else {
        return LINE;
    }
    return NONE;
}

float fitObject ( std::vector<geo::Vec2f>& points, int FITTINGMETHOD,  unsigned int* cornerIndex, ed::tracking::Rectangle* rectangle, ed::tracking::Circle* circle, std::vector<geo::Vec2f>::iterator* it_low, std::vector<geo::Vec2f>::iterator* it_high, const geo::Pose3D& sensor_pose )
{
    switch ( FITTINGMETHOD )
    {
    case NONE:
    {
        return std::numeric_limits<float>::infinity();
    }
    case LINE:
    {
        return setRectangularParametersForLine ( points,  it_low,  it_high, rectangle, sensor_pose );
    }
    case CIRCLE:
    {
        return fitCircle ( points, circle, sensor_pose );
    }
    case RECTANGLE:
    {
        return fitRectangle ( points, rectangle, sensor_pose , *cornerIndex );
    }
    }
    return false; // end reached without doing something
}

geo::Vec2f avg ( std::vector<geo::Vec2f>& points, std::vector<geo::Vec2f>::const_iterator it_start, std::vector<geo::Vec2f>::const_iterator it_end )
{
    geo::Vec2f avg_point;
    avg_point.x = avg_point.y= 0.0;

    for ( std::vector<geo::Vec2f>::const_iterator it = it_start; it != it_end; ++it ) {
        geo::Vec2f point = *it;
        avg_point.x += point.x;
        avg_point.y += point.y;
    }

    unsigned int nElements = std::distance ( it_start, it_end );
    avg_point.x /= nElements;
    avg_point.y /= nElements;

    return ( avg_point );
}

//Fast Line, Arc/Circle and Leg Detection from Laser Scan Data in a Player Driver: http://miarn.sourceforge.net/pdf/a1738b.pdf
float fitCircle ( std::vector<geo::Vec2f>& points, ed::tracking::Circle* circle, const geo::Pose3D& pose )
{
    // according to https://dtcenter.org/met/users/docs/write_ups/circle_fit.pdf
    float x_avg = 0.0, y_avg = 0.0;
    for ( unsigned int i = 0; i < points.size(); ++i ) {
        x_avg += points[i].x;
        y_avg += points[i].y;
    }

    x_avg /= points.size();
    y_avg /= points.size();

    std::vector<float> ui ( points.size() ), vi ( points.size() );
    float Suu = 0.0, Suv = 0.0, Suuu = 0.0, Suvv = 0.0, Svv = 0.0, Svvv = 0.0, Svuu = 0.0;
    for ( unsigned int i = 0; i < points.size(); ++i ) {
        ui[i] = points[i].x -x_avg;
        vi[i] = points[i].y -y_avg;

        Suu += ui[i]*ui[i];
        Suv += ui[i]*vi[i];
        Suuu += ui[i]*ui[i]*ui[i];
        Suvv += ui[i]*vi[i]*vi[i];

        Svv += vi[i]*vi[i];
        Svvv += vi[i]*vi[i]*vi[i];
        Svuu += vi[i]*ui[i]*ui[i];
    }

    float a = Suu;
    float b = Suv;
    float c = 0.5* ( Suuu+Suvv );
    float d = Suv;
    float e = Svv;
    float f = 0.5* ( Svvv+Svuu );

    float vc = ( f - c*d/a ) / ( e-b*d/a );
    float uc = ( c-vc*b ) /a;

    float xc = uc+x_avg;
    float yc = vc+y_avg;

    float alpha = uc*uc+vc*vc+ ( Suu+Svv ) /points.size();
    float R = std::sqrt ( alpha );

    float sum = 0.0;
    for ( unsigned int i = 0; i < points.size(); ++i ) {
        float error = fabs ( sqrt ( pow ( xc - points[i].x, 2.0 ) + pow ( yc - points[i].y, 2.0 ) ) - R ); // distance between a point and a circle;
        float error2 = pow ( error, 2.0 );
        sum += error2;
    }

    float roll = 0.0, pitch = 0.0, yaw = 0.0;
    circle->setValues ( xc, yc, pose.getOrigin().getZ(),  R, roll, pitch, yaw); // Assumption: object-height identical to sensor-height
    return sum/points.size();
}

Circle::Circle()
{
    float notANumber = 0.0/0.0;
    this->setValues( notANumber, notANumber, notANumber, notANumber, notANumber, notANumber, notANumber  ); // Produces NaN values, meaning that the properties are not initialized yet
}

void Circle::setValues ( float x, float y, float z, float R, float roll, float pitch, float yaw )
{
    x_ = x;
    y_ = y;
    z_ = z;
    R_ = R;
    roll_ = roll;
    pitch_ = pitch;
    yaw_ = yaw;
}

void Circle::printValues ( )
{
    std::cout << "x_ = " << x_;
    std::cout << " y_ = " << y_;
    std::cout << " R_ = " << R_;
    std::cout << " roll_ = " << roll_;
    std::cout << " pitch_ = " << pitch_;
    std::cout << " yaw_ = " << yaw_<< std::endl;
}

void Circle::setMarker ( visualization_msgs::Marker& marker , unsigned int ID )
{
   std_msgs::ColorRGBA color;
   color.a = 0.5;
   color.r = 0.0;
   color.g = 1.0;
   color.b = 0.0;

   this->setMarker ( marker, ID, color );
}

void Circle::setMarker ( visualization_msgs::Marker& marker, unsigned int ID, std_msgs::ColorRGBA color )
{
    marker.header.frame_id = "/map";
    marker.header.stamp = ros::Time::now();
    marker.ns = "default namespace";
    marker.id = ID;
    marker.type = visualization_msgs::Marker::CYLINDER;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = x_;
    marker.pose.position.y = y_;
    marker.pose.position.z = z_;
//     marker.pose.orientation.x = 0.0;
//     marker.pose.orientation.y = 0.0;
//     marker.pose.orientation.z = 0.0;
//     marker.pose.orientation.w = 1.0;
    marker.pose.orientation = tf::createQuaternionMsgFromRollPitchYaw ( roll_, pitch_, yaw_ );
    marker.scale.x = 2*R_;
    marker.scale.y = 2*R_;
    marker.scale.z = 0.1;

    marker.color = color;

    marker.lifetime = ros::Duration ( TIMEOUT_TIME );
}

std::vector< geo::Vec2f > Circle::convexHullPoints(unsigned int nPoints)
{
  std::vector< geo::Vec2f > Points(nPoints);
  float deltaAngle = 2*M_PIl / nPoints;

  for(unsigned int ii = 0; ii < nPoints; ii++)
  {
    float angle = ii*deltaAngle;
    Points[ii].x = x_ + R_*cos(angle);
    Points[ii].y = y_ + R_*sin(angle);
  }

  return Points;
}

float fitRectangle ( std::vector<geo::Vec2f>& points, ed::tracking::Rectangle* rectangle, const geo::Pose3D& pose , unsigned int cornerIndex )
{
    std::vector<geo::Vec2f>::iterator it_start = points.begin();
    std::vector<geo::Vec2f>::iterator it_split = points.begin() + cornerIndex;
    std::vector<geo::Vec2f>::iterator it_end = points.end();

    Eigen::VectorXd beta_hat1 ( 2 ), beta_hat2 ( 2 );
    bool validFit1, validFit2;

    float mean_error1 = fitLine ( points, beta_hat1, &it_start, &it_split );
    float mean_error2 = fitLine ( points, beta_hat2, &it_split, &it_end );

    float x_start1 = points[0].x; // Is this correct in combination with theta?
    float y_start1 = beta_hat1 ( 1 ) * x_start1 + beta_hat1 ( 0 );

    //determine width and height
    float x_end = points[cornerIndex].x;
    float y_end = points[cornerIndex].y;

    float dx = x_start1 - x_end;
    float dy = y_start1 - y_end;
    float width = sqrt ( dx*dx+dy*dy );
    float theta = atan2 ( beta_hat1 ( 1 ), 1 ); // TODO: angle on points low alone?

    float x_start2 = x_end;
    float y_start2 = y_end;

    float x_end2 = points[points.size() - 1].x;
    float y_end2 = points[points.size() - 1].y;

    dx = x_end2 - x_start2;
    dy = y_start2 - y_end2;
    float depth = sqrt ( dx*dx+dy*dy );

    float center_x = 0.5* ( x_start1 + x_end ) + 0.5* ( x_end2-x_start2 );
    float center_y = 0.5* ( y_start1 + y_end ) + 0.5* ( y_end2-y_start2 );

    float roll = 0.0, pitch = 0.0, yaw = theta;
    rectangle->setValues ( center_x, center_y, pose.getOrigin().getZ(), width, depth, ARBITRARY_HEIGHT, roll, pitch, yaw ); // Assumption: object-height identical to sensor-height

    unsigned int low_size = cornerIndex;
    unsigned int high_size = points.size() - cornerIndex + 1;
    return ( mean_error1*low_size+mean_error2*high_size ) / ( low_size + high_size ); // average of error
}

bool findPossibleCorner ( std::vector<geo::Vec2f>& points, std::vector<unsigned int> *IDs, std::vector<geo::Vec2f>::iterator* it_start, std::vector<geo::Vec2f>::iterator* it_end )
{
    float maxDistance = 0.0;
    unsigned int ID = std::numeric_limits<unsigned int>::quiet_NaN();

    geo::Vec2f startPoint = **it_start;
    geo::Vec2f endPoint = * ( *it_end - 1 );

    float a = endPoint.y-startPoint.y;
    float b = endPoint.x-startPoint.x;
    float c = endPoint.x*startPoint.y-endPoint.y*startPoint.x;

    float length = sqrt ( pow ( a,2.0 ) + pow ( b,2.0 ) );

    // See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    for ( std::vector<geo::Vec2f>::iterator it = *it_start + 1; it != *it_end - 1; ++it ) {

        geo::Vec2f point = *it;
        float distance = fabs ( a* ( point.x )-b* ( point.y ) +c ) / length;

        if ( distance > maxDistance ) {
            maxDistance = distance;
            ID = std::distance ( points.begin(), it );
        }
    }

    if ( maxDistance >  MIN_DISTANCE_CORNER_DETECTION ) {
        IDs->push_back ( ID );
        return true;
    } else {
        return false;
    }

}

bool findPossibleCorners ( std::vector<geo::Vec2f>& points, std::vector<unsigned int> *cornerIndices )
{
    // Check in section if a corner is present. If that is the case, split the data at this corner, and check for both parts if another corner is present.
    std::vector<laserSegments> segments;

    std::vector<geo::Vec2f>::iterator it_start = points.begin();
    std::vector<geo::Vec2f>::iterator it_end = points.end();

    if ( findPossibleCorner ( points, cornerIndices, &it_start, &it_end ) ) { // -1 because std::vector::end returns an iterators to one-past-the-end of the container. The element just before is then the last element in the vector.
        laserSegments segmentToAdd;
        segmentToAdd.begin = points.begin();
        segmentToAdd.end = points.begin() + cornerIndices->back();

        if ( segmentToAdd.end - segmentToAdd.begin > MIN_POINTS_LINEFIT ) {
            segments.push_back ( segmentToAdd );
        }

        segmentToAdd.begin = points.begin() + cornerIndices->back();
        segmentToAdd.end = points.end() - 1;

        if ( segmentToAdd.end - segmentToAdd.begin > MIN_POINTS_LINEFIT ) {
            segments.push_back ( segmentToAdd );
        }

        for ( unsigned int ii = 0; ii < segments.size(); ++ii ) {
            laserSegments laserSegment = segments[ii];
            geo::Vec2f pointEnd =  *laserSegment.end;

            bool test = findPossibleCorner ( points, cornerIndices, &laserSegment.begin, &laserSegment.end );

            if ( test ) {

                segmentToAdd.begin = laserSegment.begin;
                segmentToAdd.end = points.begin() + cornerIndices->back();

                if ( segmentToAdd.end - segmentToAdd.begin > MIN_POINTS_LINEFIT ) {
                    segments.push_back ( segmentToAdd );
                }

                segmentToAdd.begin = points.begin() + cornerIndices->back();
                segmentToAdd.end = laserSegment.end;

                if ( segmentToAdd.end - segmentToAdd.begin > MIN_POINTS_LINEFIT ) {
                    segments.push_back ( segmentToAdd );
                }
            }
        }

        std::sort ( cornerIndices->begin(), cornerIndices->end(), greater() );

        return true;
    } else {
        return false;
    }

}

bool checkForSplit ( std::vector<geo::Vec2f>& points, unsigned int &ID,const geo::Pose3D& sensor_pose,  unsigned int cornerIndex )
{
    // check if a split is required: 2 objects close to each other can form a rectangle in the wrong quadrant. Model as 2 separate lines
    geo::Vec2f centerpoint;
    centerpoint.x = 0.5* ( points[0].x + points[points.size() - 1].x );
    centerpoint.y = 0.5* ( points[0].y + points[points.size() - 1].y );

    float centerDist2 = pow ( sensor_pose.getOrigin().getX() - centerpoint.x, 2.0 ) + pow ( sensor_pose.getOrigin().getY() - centerpoint.y, 2.0 );
    float cornerDist2 = pow ( sensor_pose.getOrigin().getX() - points[cornerIndex].x, 2.0 ) + pow ( sensor_pose.getOrigin().getY() - points[cornerIndex].y, 2.0 );

    if ( centerDist2 < cornerDist2 ) {
        return true;
    } else {
        return false;
    }
}

float fitLine ( std::vector<geo::Vec2f>& points, Eigen::VectorXd& beta_hat, std::vector<geo::Vec2f>::iterator* it_start, std::vector<geo::Vec2f>::iterator* it_end )  //, unsigned int& index )
{
    // Least squares method: http://home.isr.uc.pt/~cpremebida/files_cp/Segmentation%20and%20Geometric%20Primitives%20Extraction%20from%202D%20Laser%20Range%20Data%20for%20Mobile%20Robot%20Applications.pdf

    unsigned int size = std::distance ( *it_start, *it_end );;
    Eigen::MatrixXd m ( size, 2 );
    Eigen::VectorXd y ( size );
    unsigned int counter = 0;

    for ( std::vector<geo::Vec2f>::iterator it = *it_start; it != *it_end; ++it ) {
        geo::Vec2f point = *it;
        m ( counter, 0 ) = ( double ) 1.0;
        m ( counter, 1 ) = ( double ) point.x;
        y ( counter ) = ( double ) point.y;
        counter++;
    }

    Eigen::MatrixXd mt ( size, 2 );
    mt = m.transpose();

    beta_hat = ( mt*m ).inverse() * mt * y;

    float error, sum = 0.0;

    counter = 0;
    for ( std::vector<geo::Vec2f>::iterator it = *it_start; it != *it_end; ++it ) {
        // Distance of each point to line
        geo::Vec2f point = *it;
        error = fabs ( -beta_hat ( 1 ) * point.x+point.y - beta_hat ( 0 ) ) /sqrt ( beta_hat ( 1 ) *beta_hat ( 1 ) + 1 ); // distance of a point to a line, see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        float error2 = pow ( error, 2.0 );
        sum += error2;
        counter ++;
    }

    return sum/counter;
}

float setRectangularParametersForLine ( std::vector<geo::Vec2f>& points,  std::vector<geo::Vec2f>::iterator* it_low, std::vector<geo::Vec2f>::iterator* it_high, ed::tracking::Rectangle* rectangle, const geo::Pose3D& sensor_pose )
{
    Eigen::VectorXd beta_hat ( 2 );
    float mean_error2 = fitLine ( points, beta_hat, it_low, it_high ) ;

    float theta = atan2 ( beta_hat ( 1 ), 1 );

    unsigned int ii_start = std::distance ( points.begin(), *it_low );
    float x_start = points[ii_start].x;
    float y_start = points[ii_start].y; // better to rely on the original points: extreme outliers are already filtered out due to segmentation

    unsigned int ii_end = std::distance ( points.begin(), *it_high );
    float x_end = points[ii_end - 1].x; // considered to be rebust to rely on 1 point, as these are filtered out already during the segmantation phase
    float y_end = points[ii_end - 1].y;

    float dx = x_end - x_start;
    float dy = y_start - y_end;
    float width = sqrt ( dx*dx+dy*dy );

    float center_x = 0.5* ( x_start + x_end );
    float center_y = 0.5* ( y_start + y_end );

    float roll = 0.0, pitch = 0.0, yaw = theta;
    rectangle->setValues ( center_x, center_y, sensor_pose.getOrigin().getZ(), width, ARBITRARY_DEPTH, ARBITRARY_HEIGHT, roll, pitch, yaw ); // Assumption: object-height identical to sensor-height

    return mean_error2;
}

Rectangle::Rectangle()
{
    float notANumber = 0.0/0.0;
    P_.setIdentity(2,2);
    this->setValues( notANumber, notANumber, notANumber, notANumber, notANumber, notANumber, notANumber, notANumber, notANumber ); // Produces NaN values, meaning that the properties are not initialized yet
}

void Rectangle::setValues ( float x, float y, float z, float w, float d, float h, float roll, float pitch, float yaw )
{
    x_ = x;
    y_ = y;
    z_ = z;
    w_ = w;
    d_ = d;
    h_ = h;
    roll_ = roll;
    pitch_ = pitch;
    yaw_ = yaw;
}

void Rectangle::printValues ( )
{
    std::cout << "x_ = "      << x_;
    std::cout << " y_ = "     << y_;
    std::cout << " z_ = "     << z_;
    std::cout << " w_ = "     << w_;
    std::cout << " d_ = "     << d_;
    std::cout << " h_ = "     << h_;
    std::cout << " roll_ = "  << roll_;
    std::cout << " pitch_ = " << pitch_;
    std::cout << " yaw_ = "   << yaw_ << std::endl;
}

void Rectangle::setMarker ( visualization_msgs::Marker& marker, unsigned int ID, std_msgs::ColorRGBA color, std::string ns )
{
    marker.header.frame_id = "/map";
    marker.header.stamp = ros::Time::now();
    marker.ns = ns;
    marker.id = ID;
    marker.type = visualization_msgs::Marker::CUBE;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = x_;
    marker.pose.position.y = y_;
    marker.pose.position.z = z_;
    marker.pose.orientation = tf::createQuaternionMsgFromRollPitchYaw ( roll_, pitch_, yaw_ );
    marker.scale.x = w_;
    marker.scale.y = d_;
    marker.scale.z = 0.1;
    marker.color = color;

    marker.lifetime = ros::Duration ( TIMEOUT_TIME );
}

void Rectangle::setMarker ( visualization_msgs::Marker& marker , unsigned int ID )
{
   std_msgs::ColorRGBA color;
   color.a = 0.5;
   color.r = 0.0;
   color.g = 1.0;
   color.b = 0.0;

   this->setMarker ( marker, ID, color );
}

void Rectangle::setMarker ( visualization_msgs::Marker& marker, unsigned int ID, std_msgs::ColorRGBA color)
{
        this->setMarker ( marker, ID, color, "default namespace" );
}

std::vector<geo::Vec2f> Rectangle::determineCorners ( float associationDistance )
{

    float dx = 0.5* ( w_ + associationDistance ); // blow up for associations
    float dy = 0.5* ( d_ + associationDistance );

    float ct = cos ( yaw_ );
    float st = sin ( yaw_ );

    geo::Vec2f originCorner ( x_ + ct*-dx + st* dy, y_ + st*-dx + ct*-dy ); // Rotation matrix
    geo::Vec2f corner1 (      x_ + ct* dx + st* dy, y_ + st* dx + ct*-dy );
    geo::Vec2f corner2 (      x_ + ct* dx - st* dy, y_ + st* dx + ct* dy );
    geo::Vec2f corner3 (      x_ + ct*-dx - st* dy, y_ + st*-dx + ct* dy );

    std::vector< geo::Vec2f > corners;
    corners.push_back ( originCorner );
    corners.push_back ( corner1 );
    corners.push_back ( corner2 );
    corners.push_back ( corner3 );

    return corners;
}

void wrapToInterval ( float* alpha, float lowerBound, float upperBound )
{
    float delta = upperBound - lowerBound;

    if ( *alpha < lowerBound ) {
        while ( *alpha < lowerBound ) {
            *alpha += delta;
        }
    } else if ( *alpha >= lowerBound ) {
        while ( *alpha >= lowerBound ) {
            *alpha -= delta;
        }
    }

}

void FeatureProbabilities::setMeasurementProbabilities ( float errorRectangleSquared, float errorCircleSquared, float circleDiameter, float typicalCorridorWidth )
{
    if ( !std::isinf ( errorRectangleSquared ) && !std::isinf ( errorCircleSquared ) ) {
        float probabilityScaling = 1.0;
        if ( circleDiameter > 0.5*typicalCorridorWidth ) {
            // Circles with a very large radius could have a smaller error compared to fitting a line. For the type of environment, this is very unlikely.
            // Therefore, the probability of large circles (diameter > typicalCorridorWidth) is reduced in an exponential fashion
            probabilityScaling = std::exp ( -1/ ( 0.5*typicalCorridorWidth ) * ( circleDiameter -0.5*typicalCorridorWidth ) );
        }

        float sum = errorRectangleSquared + errorCircleSquared;
        float pCircle = probabilityScaling * errorRectangleSquared/sum;
        float pRectangle =  1.0 - pCircle;  // Only 2 objects now, so the sum of it equals 1

        pmf_.setProbability ( "Rectangle", pRectangle );
        pmf_.setProbability ( "Circle", pCircle );
    } else {
        // Infinity detected on one of the elements: this is the case when the number of points is too small to do a fit. Equal probability set.
        pmf_.setProbability ( "Rectangle", 0.5 );
        pmf_.setProbability ( "Circle", 0.5 );
    }
}

void FeatureProbabilities::update ( float pRectangle_measured, float pCircle_measured )
{
    pbl::PMF pmf_measured = pmf_;

    pmf_measured.setProbability ( "Rectangle", pRectangle_measured );
    pmf_measured.setProbability ( "Circle", pCircle_measured );

    pmf_.update ( pmf_measured );
}

void FeatureProbabilities::update ( FeatureProbabilities& featureProbabilities_in )
{
    this->pmf_.update ( featureProbabilities_in.pmf_ );
}

void FeatureProperties::updateCircleFeatures(float Q_k, float R_k, float z_k) // z = observation
{
  float x_k_1_k_1 = circle_.get_R();
  float x_k_k_1 = x_k_1_k_1; // Predicted state estimate. F = Identity.

  float y_k = z_k - x_k_k_1;// H = Identity

  float P_k_k_1 = circle_.get_P() + Q_k;
  float S_k = P_k_k_1 + R_k;
  float K_k = P_k_k_1/S_k;
  float x_k_k = x_k_1_k_1 + K_k*y_k;
  float P_k_k = (1 - K_k)*P_k_k_1;

  circle_.set_R(x_k_k);
  circle_.set_P(P_k_k);

}

void FeatureProperties::updateRectangleFeatures(Eigen::MatrixXd Q_k, Eigen::MatrixXd R_k, Eigen::VectorXd z_k)
{
//   std::cout << "RectangleFeature update1" << std::endl;
  Eigen::MatrixXd F(2,2), H(2,2), I(2,2);
  I.setIdentity(2,2);
  F = I;
  H = I;

//   std::cout << "RectangleFeature update2" << std::endl;

  Eigen::VectorXd x_k_1_k_1(2);
  x_k_1_k_1 << rectangle_.get_w(),rectangle_.get_d();

//   std::cout << "RectangleFeature update3" << std::endl;

  //Eigen::MatrixXd::Identity(2,2);

  //= circle_.get_R();
  Eigen::VectorXd x_k_k_1 = F*x_k_1_k_1; // Predicted state estimate. F = Identity.
//   std::cout << "RectangleFeature update4" << std::endl;
  Eigen::VectorXd y_k = z_k - H*x_k_k_1;// H = Identity
//   std::cout << "RectangleFeature update5" << std::endl;
//     std::cout << " rectangle_.get_P() = " << rectangle_.get_P() << std::endl;
//   std::cout << " Q_k = " << Q_k << std::endl;
  Eigen::MatrixXd P_k_k_1 = rectangle_.get_P() + Q_k;
//     std::cout << "RectangleFeature update5.1" << std::endl;
  Eigen::MatrixXd S_k = H*P_k_k_1*H + R_k; // H = H transpose
//     std::cout << "RectangleFeature update5.2" << std::endl;
  Eigen::MatrixXd K_k = P_k_k_1*H*S_k.inverse();
//     std::cout << "RectangleFeature update5.3" << std::endl;
  Eigen::VectorXd x_k_k = x_k_1_k_1 + K_k*y_k;
//     std::cout << "RectangleFeature update5.4" << std::endl;
  Eigen::MatrixXd P_k_k = (I - K_k)*P_k_k_1;
//   std::cout << "RectangleFeature update6" << std::endl;
  rectangle_.set_w(x_k_k(0) );
//   std::cout << "RectangleFeature update6.1" << std::endl;
  rectangle_.set_d(x_k_k(1) );
//   std::cout << "RectangleFeature update6.2" << std::endl;

  rectangle_.set_P(P_k_k );
//   std::cout << "RectangleFeature update6.3" << std::endl;
}


}

namespace convex_hull
{
// ----------------------------------------------------------------------------------------------------

/**
 * @brief create fill a ConvexHull and put its origin in the middle of the Convexhull
 * @param points points describing the shape. Points should be in map frame.
 * @param z_min min height
 * @param z_max max height
 * @param c filled ConvexHull
 * @param pose new pose of the convexhull
 */
void create(const std::vector<geo::Vec2f>& points, float z_min, float z_max, ConvexHull& chull, geo::Pose3D& pose)
{
    cv::Mat_<cv::Vec2f> points_2d ( 1, points.size() );
    for ( unsigned int i = 0; i < points.size(); ++i ) {
        points_2d.at<cv::Vec2f> ( i ) = cv::Vec2f ( points[i].x, points[i].y );
    }

    pose = geo::Pose3D::identity();

    pose.t.z = ( z_min + z_max ) / 2;

    std::vector<int> chull_indices;
    cv::convexHull ( points_2d, chull_indices );

    chull.z_min = z_min - pose.t.z;
    chull.z_max = z_max - pose.t.z;

    geo::Vec2f xy_min ( 1e9, 1e9 );
    geo::Vec2f xy_max ( -1e9, -1e9 );

    chull.points.clear();
    for ( unsigned int i = 0; i < chull_indices.size(); ++i ) {
        const cv::Vec2f& p_cv = points_2d.at<cv::Vec2f> ( chull_indices[i] );
        geo::Vec2f p ( p_cv[0], p_cv[1] );

        chull.points.push_back ( p );

        xy_min.x = std::min ( xy_min.x, p.x );
        xy_min.y = std::min ( xy_min.y, p.y );

        xy_max.x = std::max ( xy_max.x, p.x );
        xy_max.y = std::max ( xy_max.y, p.y );
    }

    // Average segment position
    pose.t.x = ( xy_min.x + xy_max.x ) / 2;
    pose.t.y = ( xy_min.y + xy_max.y ) / 2;

    // Move all points to the pose frame
    for ( unsigned int i = 0; i < chull.points.size(); ++i ) {
        geo::Vec2f& p = chull.points[i];
        p.x -= pose.t.x;
        p.y -= pose.t.y;
    }

    // Calculate normals and edges
    convex_hull::calculateEdgesAndNormals ( chull );

    // Calculate area
    calculateArea ( chull );
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief createAbsolute fill a ConvexHull with its origin in map frame.
 * @param points points describing the shape. Points should be in map frame.
 * @param z_min min height
 * @param z_max max height
 * @param c filled ConvexHull
 */
void createAbsolute(const std::vector<geo::Vec2f>& points, float z_min, float z_max, ConvexHull& chull)
{
    cv::Mat_<cv::Vec2f> points_2d ( 1, points.size() );
    for ( unsigned int i = 0; i < points.size(); ++i ) {
        points_2d.at<cv::Vec2f> ( i ) = cv::Vec2f ( points[i].x, points[i].y );
    }

    chull.z_min = z_min;
    chull.z_max = z_max;

    std::vector<int> chull_indices;
    cv::convexHull ( points_2d, chull_indices );

    chull.points.resize(chull_indices.size());
    for(unsigned int i = 0; i < chull_indices.size(); ++i)
    {
        const cv::Vec2f& p_cv = points_2d.at<cv::Vec2f>(chull_indices[i]);
        chull.points[i] = geo::Vec2f(p_cv[0], p_cv[1]);
    }

    // Calculate normals and edges
    convex_hull::calculateEdgesAndNormals(chull);

    // Calculate area
    calculateArea(chull);
}

// ----------------------------------------------------------------------------------------------------

void calculateEdgesAndNormals(ConvexHull& chull)
{
    chull.edges.resize(chull.points.size());
    chull.normals.resize(chull.points.size());

    for(unsigned int i = 0; i < chull.points.size(); ++i)
    {
        unsigned int j = (i + 1) % chull.points.size();

        const geo::Vec2f& p1 = chull.points[i];
        const geo::Vec2f& p2 = chull.points[j];

        // Calculate edge
        geo::Vec2f e = p2 - p1;
        chull.edges[i] = e;

        // Calculate normal
        chull.normals[i] = geo::Vec2f(e.y, -e.x).normalized();
    }
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief collide Check of two ConvexHull collide with padding in xy and in z.
 * @param c1 ConvexHull 1
 * @param pos1 position of ConvexHull 1
 * @param c2 ConvexHull 2
 * @param pos2 position of ConvexHull 2
 * @param xy_padding padding in xy-plane
 * @param z_padding padding in z-plane
 * @return
 */
bool collide(const ConvexHull& c1, const geo::Vector3& pos1,
             const ConvexHull& c2, const geo::Vector3& pos2,
             float xy_padding, float z_padding)
{
    if ( c1.points.size() < 3 || c2.points.size() < 3 ) {
        return false;
    }

    float z_diff = pos2.z - pos1.z;

    if ( c1.z_max < ( c2.z_min + z_diff - 2 * z_padding ) || c2.z_max < ( c1.z_min - z_diff - 2 * z_padding ) ) {
        return false;
    }

    geo::Vec2f pos_diff ( pos2.x - pos1.x, pos2.y - pos1.y );

    for ( unsigned int i = 0; i < c1.points.size(); ++i ) {
        const geo::Vec2f& p1 = c1.points[i];
        const geo::Vec2f& n = c1.normals[i];

        // Calculate min and max projection of c1
        float min1 = n.dot ( c1.points[0] - p1 );
        float max1 = min1;
        for ( unsigned int k = 1; k < c1.points.size(); ++k ) {
            // Calculate projection
            float p = n.dot ( c1.points[k] - p1 );
            min1 = std::min ( min1, p );
            max1 = std::max ( max1, p );
        }

        // Apply padding to both sides
        min1 -= xy_padding;
        max1 += xy_padding;

        // Calculate p1 in c2's frame
        geo::Vec2f p1_c2 = p1 - pos_diff;

        // If this bool stays true, there is definitely no collision
        bool no_collision = true;

        // True if projected points are found below c1's bounds
        bool below = false;

        // True if projected points are found above c1's bounds
        bool above = false;

        // Check if c2's points overlap with c1's bounds
        for ( unsigned int k = 0; k < c2.points.size(); ++k ) {
            // Calculate projection on p1's normal
            float p = n.dot ( c2.points[k] - p1_c2 );

            below = below || ( p < max1 );
            above = above || ( p > min1 );

            if ( below && above ) {
                // There is overlap with c1's bound, so we may have a collision
                no_collision = false;
                break;
            }
        }

        if ( no_collision )
            // definitely no collision
        {
            return false;
        }
    }

    return true;
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief calculateArea calculate the area of the ConvexHull in xy-plane.
 * @param c ConvexHull
 */
void calculateArea(ConvexHull& c)
{
    c.area = 0;
    for ( unsigned int i = 0; i < c.points.size(); ++i )
    {
        unsigned int j = ( i + 1 ) % c.points.size();

        const geo::Vec2f& p1 = c.points[i];
        const geo::Vec2f& p2 = c.points[j];

        c.area += 0.5 * ( p1.x * p2.y - p2.x * p1.y );
    }
}

// ----------------------------------------------------------------------------------------------------

}

}
