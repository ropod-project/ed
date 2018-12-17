#include "shape_loader.h"

#include "xml_shape_parser.h"

#include <tue/filesystem/path.h>

#include <geolib/serialization.h>
#include <geolib/Shape.h>
#include <geolib/CompositeShape.h>
#include <geolib/Importer.h>

// Heightmap generation
#include "polypartition/polypartition.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

// string split
#include <sstream>
#include <iostream>

namespace ed
{
namespace models
{

/*
std::string split implementation by using delimiter as a character. Multiple delimeters are removed.
*/
/**
 * @brief split Implementation by using delimiter as a character. Multiple delimeters are removed.
 * @param strToSplit input string, which is splitted
 * @param delimeter char on which the string is split
 * @return vector of sub-strings
 */
std::vector<std::string> split(std::string& strToSplit, char delimeter)
{
    std::stringstream ss(strToSplit);
    std::string item;
    std::vector<std::string> splittedStrings;
    while (std::getline(ss, item, delimeter))
    {
        if (!item.empty() && item[0] != delimeter)
            splittedStrings.push_back(item);
    }
    return splittedStrings;
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief getFilePath searches GAZEBO_MODEL_PATH and GAZEBO_RESOURCH_PATH for file
 * @param type subpath+filename incl. extension
 * @return full path or empty string in case not found
 */
std::string getFilePath(std::string type)
{
    //ToDo: don't do this for every call. Very ineffecient.
    const char * mpath = ::getenv("GAZEBO_MODEL_PATH");
    const char * rpath = ::getenv("GAZEBO_RESOURCE_PATH");
    if (!mpath && !rpath)
        return "";

    // remove prefix in case of sdf
    std::string str1 = "file://";
    std::string str2 = "model://";

    std::string::size_type i = type.find(str1);
    bool file = true, model = true;
    // Start by both true. So both paths are checked in case of ED yaml.
    // In case of SDF, the other one is set to false.
    if (i != std::string::npos)
    {
       type.erase(i, str1.length());
       model = false;
    }
    i = type.find(str2);
    if (i != std::string::npos)
    {
       type.erase(i, str2.length());
       file = false;
    }

    std::string item;
    std::vector<std::string> model_paths;
    if (model)
    {
        std::stringstream ssm(mpath);
        while (std::getline(ssm, item, ':'))
            model_paths.push_back(item);
    }
    if (file)
    {
        std::stringstream ssr(rpath);
        while (std::getline(ssr, item, ':'))
            model_paths.push_back(item);
    }

    // romove duplicate elements
    std::sort( model_paths.begin(), model_paths.end() );
    model_paths.erase( unique( model_paths.begin(), model_paths.end() ), model_paths.end() );

    for(std::vector<std::string>::const_iterator it = model_paths.begin(); it != model_paths.end(); ++it)
    {
        tue::filesystem::Path file_path(*it + "/" + type);
        if (file_path.exists())
            return file_path.string();
    }

    return "";
}


// ----------------------------------------------------------------------------------------------------

/**
 * @brief findContours
 * @param image Grayscale image
 * @param p_start starting point
 * @param d_start starting direction
 * @param points
 * @param line_starts
 * @param contour_map
 */
void findContours(const cv::Mat& image, const geo::Vec2i& p_start, int d_start, std::vector<geo::Vec2i>& points,
                  std::vector<geo::Vec2i>& line_starts, cv::Mat& contour_map)
{
    static int dx[4] = {1,  0, -1,  0 };
    static int dy[4] = {0,  1,  0, -1 };

    unsigned char v = image.at<unsigned char>(p_start.y, p_start.x);

    int d = d_start; // Current direction
    geo::Vec2i p = p_start;

    while (true)
    {
        if (d == 3) // up
            line_starts.push_back(p);

        if (d == 3)
            contour_map.at<unsigned char>(p.y, p.x - 1) = 1;
        else if (d == 1)
            contour_map.at<unsigned char>(p.y, p.x) = 1;

        if (image.at<unsigned char>(p.y + dy[(d + 3) % 4], p.x + dx[(d + 3) % 4]) == v) // left
        {
            switch (d)
            {
                case 0: points.push_back(p - geo::Vec2i(1, 1)); break;
                case 1: points.push_back(p - geo::Vec2i(0, 1)); break;
                case 2: points.push_back(p); break;
                case 3: points.push_back(p - geo::Vec2i(1, 0)); break;
            }

            d = (d + 3) % 4;
            p.x += dx[d];
            p.y += dy[d];
        }
        else if (image.at<unsigned char>(p.y + dy[d], p.x + dx[d]) == v) // straight ahead
        {
            p.x += dx[d];
            p.y += dy[d];
        }
        else // right
        {
            switch (d)
            {
                case 0: points.push_back(p - geo::Vec2i(0, 1)); break;
                case 1: points.push_back(p); break;
                case 2: points.push_back(p - geo::Vec2i(1, 0)); break;
                case 3: points.push_back(p - geo::Vec2i(1, 1)); break;

            }

            d = (d + 1) % 4;
        }

        if (p.x == p_start.x && p.y == p_start.y && d == d_start)
            return;
    }
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief getHeightMapShape convert grayscale image in a heigtmap mesh
 * @param image_orig grayscale image
 * @param pos position of the origin of the heigtmap
 * @param size dimensions of the final mesh
 * @param inverted false: CV/ROS standard (black = height); true: SDF/GAZEBO (White = height)
 * @param error errorstream
 * @return final mesh; or empty mesh in case of error
 */
geo::ShapePtr getHeightMapShape(cv::Mat& image_orig, const geo::Vec3& pos, const geo::Vec3& size, const bool inverted, std::stringstream& error)
{
    double resolution_x = size.x/image_orig.cols;
    double resolution_y = size.y/image_orig.rows;
    double blockheight = size.z;

    // invert grayscale for SDF
    if (inverted)
    {
        for (int i = 0; i < image_orig.rows; i++)
        {
            for (int j = 0; j < image_orig.cols; j++)
            {
                image_orig.at<uchar>(i, j) = 255 - image_orig.at<uchar>(i, j);
            }
        }
    }

    // Add borders
    cv::Mat image(image_orig.rows + 2, image_orig.cols + 2, CV_8UC1, cv::Scalar(255));
    image_orig.copyTo(image(cv::Rect(cv::Point(1, 1), cv::Size(image_orig.cols, image_orig.rows))));

    cv::Mat vertex_index_map(image.rows, image.cols, CV_32SC1, -1);
    cv::Mat contour_map(image.rows, image.cols, CV_8UC1, cv::Scalar(0));

    geo::CompositeShapePtr shape(new geo::CompositeShape);

    for(int y = 0; y < image.rows; ++y)
    {
        for(int x = 0; x < image.cols; ++x)
        {
            unsigned char v = image.at<unsigned char>(y, x);

            if (v < 255)
            {
                std::vector<geo::Vec2i> points, line_starts;
                findContours(image, geo::Vec2i(x, y), 0, points, line_starts, contour_map);

                unsigned int num_points = (unsigned int) points.size();

                if (num_points > 2)
                {
                    geo::Mesh mesh;

                    double min_z = pos.z;
                    double max_z = pos.z + (double)(255 - v) / 255 * blockheight;

                    std::list<TPPLPoly> testpolys;

                    TPPLPoly poly;
                    poly.Init(num_points);

                    for(unsigned int i = 0; i < num_points; ++i)
                    {
                        poly[i].x = points[i].x;
                        poly[i].y = points[i].y;

                        // Convert to world coordinates
                        double wx = points[i].x * resolution_x + pos.x;
                        double wy = (image.rows - points[i].y - 2) * resolution_y + pos.y;

                        vertex_index_map.at<int>(points[i].y, points[i].x) = mesh.addPoint(geo::Vector3(wx, wy, min_z));
                        mesh.addPoint(geo::Vector3(wx, wy, max_z));
                    }

                    testpolys.push_back(poly);

                    // Calculate side triangles
                    for(unsigned int i = 0; i < num_points; ++i)
                    {
                        int j = (i + 1) % num_points;
                        mesh.addTriangle(i * 2, i * 2 + 1, j * 2);
                        mesh.addTriangle(i * 2 + 1, j * 2 + 1, j * 2);
                    }

                    for(unsigned int i = 0; i < line_starts.size(); ++i)
                    {
                        int x2 = line_starts[i].x;
                        int y2 = line_starts[i].y;

                        while(image.at<unsigned char>(y2, x2) == v)
                            ++x2;

                        if (contour_map.at<unsigned char>(y2, x2 - 1) == 0)
                        {
                            // found a hole, so find the contours of this hole
                            std::vector<geo::Vec2i> hole_points;
                            findContours(image, geo::Vec2i(x2 - 1, y2 + 1), 1, hole_points, line_starts, contour_map);

                            if (hole_points.size() > 2)
                            {
                                TPPLPoly poly_hole;
                                poly_hole.Init(hole_points.size());
                                poly_hole.SetHole(true);

                                for(unsigned int j = 0; j < hole_points.size(); ++j)
                                {
                                    poly_hole[j].x = hole_points[j].x;
                                    poly_hole[j].y = hole_points[j].y;

                                    // Convert to world coordinates
                                    double wx = hole_points[j].x * resolution_x + pos.x;
                                    double wy = (image.rows - hole_points[j].y - 2) * resolution_y + pos.y;

                                    vertex_index_map.at<int>(hole_points[j].y, hole_points[j].x) = mesh.addPoint(geo::Vector3(wx, wy, min_z));
                                    mesh.addPoint(geo::Vector3(wx, wy, max_z));
                                }
                                testpolys.push_back(poly_hole);

                                // Calculate side triangles
                                for(unsigned int j = 0; j < hole_points.size(); ++j)
                                {
                                    const geo::Vec2i& hp1 = hole_points[j];
                                    const geo::Vec2i& hp2 = hole_points[(j + 1) % hole_points.size()];

                                    int i1 = vertex_index_map.at<int>(hp1.y, hp1.x);
                                    int i2 = vertex_index_map.at<int>(hp2.y, hp2.x);

                                    mesh.addTriangle(i1, i1 + 1, i2);
                                    mesh.addTriangle(i2, i1 + 1, i2 + 1);
                                }
                            }
                        }
                    }

                    TPPLPartition pp;
                    std::list<TPPLPoly> result;

                    if (!pp.Triangulate_EC(&testpolys, &result))
                    {
                        error << "[ED::MODELS::LOADSHAPE] Error while creating heightmap: could not triangulate polygon." << std::endl;
                        return geo::ShapePtr();
                    }

                    for(std::list<TPPLPoly>::iterator it = result.begin(); it != result.end(); ++it)
                    {
                        TPPLPoly& cp = *it;

                        int i1 = vertex_index_map.at<int>(cp[0].y, cp[0].x) + 1;
                        int i2 = vertex_index_map.at<int>(cp[1].y, cp[1].x) + 1;
                        int i3 = vertex_index_map.at<int>(cp[2].y, cp[2].x) + 1;

                        mesh.addTriangle(i1, i3, i2);
                    }

                    geo::Shape sub_shape;
                    sub_shape.setMesh(mesh);

                    shape->addShape(sub_shape, geo::Pose3D::identity());

                    cv::floodFill(image, cv::Point(x, y), 255);
                }
            }
        }
    }

    return shape;
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief getHeightMapShape convert grayscale image in a heigtmap mesh
 * @param image_filename full path of grayscale image
 * @param pos position of the origin of the heigtmap
 * @param size dimensions of the final mesh
 * @param inverted false: CV/ROS standard (black = height); true: SDF/GAZEBO (White = height)
 * @param errorerrorstream
 * @return final mesh; or empty mesh in case of error
 */
geo::ShapePtr getHeightMapShape(const std::string& image_filename, const geo::Vec3& pos, const geo::Vec3& size,
                                const bool inverted, std::stringstream& error)
{
    cv::Mat image_orig = cv::imread(image_filename, CV_LOAD_IMAGE_GRAYSCALE);   // Read the file

    if (!image_orig.data)
    {
        error << "[ED::MODELS::LOADSHAPE] Error while loading heightmap '" << image_filename << "'. Image could not be loaded." << std::endl;
        return geo::ShapePtr();
    }

    return getHeightMapShape(image_orig, pos, size, inverted, error);
}


// ----------------------------------------------------------------------------------------------------

/**
 * @brief getHeightMapShape convert grayscale image in a heigtmap mesh
 * @param image_filename full path of grayscale image
 * @param pos position of the origin of the heigtmap
 * @param blockheight height of the heightmap of max grayscale value
 * @param resolution_x resolution in x direction in meters
 * @param resolution_y resolution in y direction in meters
 * @param inverted false: CV/ROS standard (black = height); true: SDF/GAZEBO (White = height)
 * @param errorerrorstream
 * @return final mesh; or empty mesh in case of error
 */
geo::ShapePtr getHeightMapShape(const std::string& image_filename, const geo::Vec3& pos, const double blockheight,
                                const double resolution_x, const double resolution_y, const bool inverted, std::stringstream& error)
{
    cv::Mat image_orig = cv::imread(image_filename, CV_LOAD_IMAGE_GRAYSCALE);   // Read the file

    if (!image_orig.data)
    {
        error << "[ED::MODELS::LOADSHAPE] Error while loading heightmap '" << image_filename << "'. Image could not be loaded." << std::endl;
        return geo::ShapePtr();
    }

    double size_x = resolution_x * image_orig.cols;
    double size_y = resolution_y * image_orig.rows;
    geo::Vec3 size(size_x, size_y, blockheight);

    return getHeightMapShape(image_orig, pos, size, inverted, error);
}


// ----------------------------------------------------------------------------------------------------

/**
 * @brief getHeightMapShape convert grayscale image in a heigtmap mesh
 * @param image_filename image_filename full path of grayscale image
 * @param cfg reader with model/shape information
 * @param error errorstream
 * @return final mesh; or empty mesh in case of error
 */
geo::ShapePtr getHeightMapShape(const std::string& image_filename, tue::config::Reader cfg, std::stringstream& error)
{
    double resolution, origin_x, origin_y, origin_z, blockheight;
    if (!(cfg.value("origin_x", origin_x) &&
            cfg.value("origin_y", origin_y) &&
            cfg.value("origin_z", origin_z) &&
            cfg.value("resolution", resolution) &&
            cfg.value("blockheight", blockheight)))
    {
        error << "[ED::MODELS::LOADSHAPE] Error while loading heightmap parameters at '" << image_filename
              << "'. Required shape parameters: resolution, origin_x, origin_y, origin_z, blockheight" << std::endl;
        return geo::ShapePtr();
    }

    int inverted = 0;
    cfg.value("inverted", inverted);

    return getHeightMapShape(image_filename, geo::Vec3(origin_x, origin_y, origin_z), blockheight, resolution, resolution,
                             (bool) inverted, error);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief createPolygon create polygon mesh from points
 * @param shape filled mesh
 * @param points 2D points which define the mesh
 * @param height height of the mesh
 * @param create_bottom false: open bottom; true: closed bottom
 */
void createPolygon(geo::Shape& shape, const std::vector<geo::Vec2>& points, double height, bool create_bottom)
{
    TPPLPoly poly;
    poly.Init((unsigned int) points.size());

    double min_z = -height / 2;
    double max_z =  height / 2;

    geo::Mesh mesh;

    for(unsigned int i = 0; i < points.size(); ++i)
    {
        poly[i].x = (unsigned int) points[i].x;
        poly[i].y = (unsigned int) points[i].y;

        mesh.addPoint(geo::Vector3(points[i].x, points[i].y, min_z));
        mesh.addPoint(geo::Vector3(points[i].x, points[i].y, max_z));
    }

    // Add side triangles
    for(unsigned int i = 0; i < points.size(); ++i)
    {
        int j = (i + 1) % points.size();
        mesh.addTriangle(i * 2, j * 2, i * 2 + 1);
        mesh.addTriangle(i * 2 + 1, j * 2, j * 2 + 1);
    }

    std::list<TPPLPoly> polys;
    polys.push_back(poly);

    TPPLPartition pp;
    std::list<TPPLPoly> result;

    if (!pp.Triangulate_EC(&polys, &result))
    {
        std::cout << "TRIANGULATION FAILED" << std::endl;
        return;
    }

    for(std::list<TPPLPoly>::iterator it = result.begin(); it != result.end(); ++it)
    {
        TPPLPoly& cp = *it;

        int i1 = mesh.addPoint(cp[0].x, cp[0].y, max_z);
        int i2 = mesh.addPoint(cp[1].x, cp[1].y, max_z);
        int i3 = mesh.addPoint(cp[2].x, cp[2].y, max_z);
        mesh.addTriangle(i1, i2, i3);

        if (create_bottom)
        {
            int i1 = mesh.addPoint(cp[0].x, cp[0].y, min_z);
            int i2 = mesh.addPoint(cp[1].x, cp[1].y, min_z);
            int i3 = mesh.addPoint(cp[2].x, cp[2].y, min_z);
            mesh.addTriangle(i1, i3, i2);
        }
    }

    mesh.filterOverlappingVertices();

    shape.setMesh(mesh);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief createCylinder create a mesh from radius and height
 * @param shape filled mesh
 * @param radius radius of the cylinder
 * @param height height of the cylinder
 * @param num_corners divided the circumference in N points and N+1 lines
 */
void createCylinder(geo::Shape& shape, double radius, double height, int num_corners)
{
    geo::Mesh mesh;

    // Calculate vertices
    for(int i = 0; i < num_corners; ++i)
    {
        double a = 6.283 * i / num_corners;
        double x = sin(a) * radius;
        double y = cos(a) * radius;

        mesh.addPoint(x, y, -height / 2);
        mesh.addPoint(x, y,  height / 2);
    }

    // Calculate top and bottom triangles
    for(int i = 1; i < num_corners - 1; ++i)
    {
        int i2 = 2 * i;

        // bottom
        mesh.addTriangle(0, i2, i2 + 2);

        // top
        mesh.addTriangle(1, i2 + 3, i2 + 1);
    }

    // Calculate side triangles
    for(int i = 0; i < num_corners; ++i)
    {
        int j = (i + 1) % num_corners;
        mesh.addTriangle(i * 2, i * 2 + 1, j * 2);
        mesh.addTriangle(i * 2 + 1, j * 2 + 1, j * 2);
    }

    shape.setMesh(mesh);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief readVec3 read x, y and z into a vector
 * @param cfg reader
 * @param v filled Vec3 vector
 * @param pos_req RequiredOrOptional
 */
void readVec3(tue::config::Reader& cfg, geo::Vec3& v, tue::config::RequiredOrOptional pos_req = tue::config::REQUIRED)
{
    cfg.value("x", v.x, pos_req);
    cfg.value("y", v.y, pos_req);
    cfg.value("z", v.z, pos_req);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief readVec3Group read a config group into a Vec3 group
 * @param cfg reader
 * @param v filled Vec3 vector
 * @param vector_name name of the reader group to be read
 * @param pos_req RequiredOrOptional
 * @return indicates succes
 */
bool readVec3Group(tue::config::Reader& cfg, geo::Vec3& v, const std::string& vector_name, tue::config::RequiredOrOptional pos_req = tue::config::REQUIRED)
{
    std::string vector_string;
    if (cfg.readGroup(vector_name))
    {
        readVec3(cfg, v);
        cfg.endGroup();
    }
    else if (cfg.value(vector_name, vector_string))
    {
        std::vector<std::string> vector_vector = split(vector_string, ' ');
        v.x = std::stod(vector_vector[0]);
        v.y = std::stod(vector_vector[1]);
        v.z = std::stod(vector_vector[2]);
    }
    else
        return false;

    return true;
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief readPose read pose into Pose3D. Both ED yaml and SDF. Also reads pos(position) of SDF.
 * @param cfg reader
 * @param pose filled Pose3D pose
 * @param pos_req position RequiredOrOptional
 * @param rot_req rotation RequiredOrOptional
 * @return indicates succes
 */
bool readPose(tue::config::Reader& cfg, geo::Pose3D& pose, tue::config::RequiredOrOptional pos_req, tue::config::RequiredOrOptional rot_req)
{
    double roll = 0, pitch = 0, yaw = 0;
    std::string pose_string = ""; //sdf pose will be a string
    if (cfg.readGroup("pose"))
    {
        readVec3(cfg, pose.t, pos_req);

        cfg.value("X", roll,  rot_req);
        cfg.value("Y", pitch, rot_req);
        cfg.value("Z", yaw,   rot_req);
        cfg.value("roll",  roll,  rot_req);
        cfg.value("pitch", pitch, rot_req);
        cfg.value("yaw",   yaw,   rot_req);

        cfg.endGroup();
    }
    else if (cfg.value("pose", pose_string, pos_req))
    {
        // pose is in SDF
        std::vector<std::string> pose_vector = split(pose_string, ' ');
        if (pose_vector.size() == 6)
        {
            // ignoring pose, when incorrect/incomplete
            pose.t.x = std::stod(pose_vector[0]);
            pose.t.y = std::stod(pose_vector[1]);
            pose.t.z = std::stod(pose_vector[2]);
            roll = std::stod(pose_vector[3]);
            pitch = std::stod(pose_vector[4]);
            yaw = std::stod(pose_vector[5]);
        }
        else
            return false;
    }
    else if (cfg.value("pos", pose_string, pos_req))
    {
        // position is in SDF
        std::vector<std::string> pose_vector = split(pose_string, ' ');
        if (pose_vector.size() == 3)
        {
            // ignoring pose, when incorrect/incomplete
            pose.t.x = std::stod(pose_vector[0]);
            pose.t.y = std::stod(pose_vector[1]);
            pose.t.z = std::stod(pose_vector[2]);
        }
        else
            return false;
    }

    // Set rotation
    pose.R.setRPY(roll, pitch, yaw);
    return true;
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief loadShape load the shape of a model.
 * @param model_path path of the model
 * @param cfg reader
 * @param shape_cache cache for complex models
 * @param error errorstream
 * @return final mesh; or empty mesh in case of error
 */
geo::ShapePtr loadShape(const std::string& model_path, tue::config::Reader cfg,
                        std::map<std::string, geo::ShapePtr>& shape_cache, std::stringstream& error)
{
    geo::ShapePtr shape;
    geo::Pose3D pose = geo::Pose3D::identity();

    std::string path;
    if (cfg.value("path", path)) // ED YAML ONLY
    {
        if (path.empty())
        {
            error << "[ED::MODELS::LOADSHAPE] Error while loading shape: provided path is empty.";
            return shape;
        }

        tue::filesystem::Path shape_path;

        if (model_path.empty() || path[0] == '/')
            shape_path = path;
        else
            shape_path = model_path + "/" + path;

        // Check cache first
        std::map<std::string, geo::ShapePtr>::const_iterator it = shape_cache.find(shape_path.string());
        if (it != shape_cache.end())
            return it->second;

        if (shape_path.exists())
        {
            std::string xt = shape_path.extension();
            if (xt == ".pgm" || xt == ".png")
            {
                shape = getHeightMapShape(shape_path.string(), cfg, error);
            }
            else if (xt == ".geo")
            {
                geo::serialization::registerDeserializer<geo::Shape>();
                shape = geo::serialization::fromFile(shape_path.string());
            }
            else if (xt == ".3ds" || xt == ".stl" || xt == ".dae")
            {
                shape = geo::Importer::readMeshFile(shape_path.string());
            }
            else if (xt == ".xml")
            {
                std::string error;
                shape = parseXMLShape(shape_path.string(), error);
            }

            if (!shape)
                error << "[ED::MODELS::LOADSHAPE] Error while loading shape at " << shape_path.string() << std::endl;
            else
                // Add to cache
                shape_cache[shape_path.string()] = shape;
        }
        else
        {
            error << "[ED::MODELS::LOADSHAPE] Error while loading shape at " << shape_path.string() << " ; file does not exist" << std::endl;
        }
    }
    else if (cfg.readGroup("box")) // SDF AND ED YAML
    {
        geo::Vec3 min, max, size;
        if (readVec3Group(cfg, min, "min"))
        {
            if (readVec3Group(cfg, max, "max"))
                shape.reset(new geo::Box(min, max));
            else
            {
                error << "[ED::MODELS::LOADSHAPE] Error while loading shape: box must contain 'min' and 'max' (only 'min' specified)";
            }
        }
        else if (readVec3Group(cfg, size, "size"))
            shape.reset(new geo::Box(-0.5 * size, 0.5 * size));
        else
        {
            error << "[ED::MODELS::LOADSHAPE] Error while loading shape: box must contain 'min' and 'max' or 'size'.";
        }

        readPose(cfg, pose);
    }
    else if (cfg.readGroup("cylinder")) // SDF AND ED YAML
    {
        int num_points = 12;
        cfg.value("num_points", num_points, tue::config::OPTIONAL);

        double radius, height;
        if (cfg.value("radius", radius) && (cfg.value("height", height) || cfg.value("length", height))) //length is used in SDF
        {
            shape.reset(new geo::Shape());
            createCylinder(*shape, radius, height, num_points);
        }

        cfg.endGroup();
    }
    else if (cfg.readGroup("polygon")) // ED YAML ONLY
    {
        std::vector<geo::Vec2> points;
        if (cfg.readArray("points", tue::config::REQUIRED) || cfg.readArray("point", tue::config::REQUIRED))
        {
            while(cfg.nextArrayItem())
            {
                points.push_back(geo::Vec2());
                geo::Vec2& p = points.back();
                cfg.value("x", p.x);
                cfg.value("y", p.y);
            }
            cfg.endArray();
        }

        double height;
        if (cfg.value("height", height))
        {
            shape.reset(new geo::Shape());
            createPolygon(*shape, points, height, true);
        }

        cfg.endGroup();
    }
    else if (cfg.readArray("polyline")) // SDF ONLY
    {
        std::vector<geo::Vec2> points;
        double height = 0.0;
        while(cfg.nextArrayItem())
        {
            std::string point_string;
            if (cfg.value("point", point_string))
            {
                points.push_back(geo::Vec2());
                geo::Vec2& p = points.back();
                std::vector<std::string> point_vector = split(point_string, ' ');
                p.x = std::stod(point_vector[0]);
                p.y = std::stod(point_vector[1]);
            }
            else
                cfg.value("height", height);
        }
        cfg.endArray();
        if (height > 0 && std::abs<double>(height) >= 0.01)
        {
            shape.reset(new geo::Shape());
            createPolygon(*shape, points, height, true);
        }
    }
    else if (cfg.readGroup("mesh")) // SDF ONLY
    {
        std::string uri_path;
        if (cfg.value("uri", uri_path))
        {
            // remove prefix in case of sdf
            std::string str1 = "file://";
            std::string str2 = "model://";

            std::string::size_type i = uri_path.find(str1);
            if (i != std::string::npos)
               uri_path.erase(i, str1.length());
            i = uri_path.find(str2);
            if (i != std::string::npos)
               uri_path.erase(i, str2.length());

            tue::filesystem::Path mesh_path;
            if (model_path.empty() || uri_path[0] == '/')
                mesh_path = uri_path;
            else
                mesh_path = model_path + "/" + uri_path;

            if (!mesh_path.exists())
                    mesh_path = getFilePath(uri_path);
            if (mesh_path.exists())
                shape = geo::Importer::readMeshFile(mesh_path.string());
            else
                error << "[ED::MODELS::LOADSHAPE] File: '" << mesh_path.string() << "' doesn't exist." << std::endl;
        }
        std::string dummy;
        if (cfg.value("submesh", dummy))
            error << "[ED::MODELS::LOADSHAPE] 'submesh' of mesh is not supported by ED " << std::endl;

    }
    else if (cfg.readGroup("heightmap")) // SDF AND ED YAML
    {
        std::string image_filename;
        double height, resolution;

        geo::Vec3 size;
        if ((cfg.value("image", image_filename) && !image_filename.empty() && cfg.value("resolution", resolution) && cfg.value("height", height))) // ED YAML ONLY
        {
            std::string image_filename_full = image_filename;
//            if (image_filename[0] == '/')
//                image_filename_full = image_filename;
//            else
//                image_filename_full = model_path + "/" + image_filename;

            shape = getHeightMapShape(image_filename_full, geo::Vec3(0, 0, 0), height, resolution, resolution, false, error);

            readPose(cfg, pose);
        }
        else if(cfg.value("uri", image_filename) && readVec3Group(cfg, size, "size")) // SDF ONLY
        {
            image_filename = getFilePath(image_filename);

            // by default center should be in the middle.
            geo::Vec3 pos = -size/2;
            pos.z = 0;

            shape = getHeightMapShape(image_filename, pos, size, true, error);
            readPose(cfg, pose);
        }
        else
        {
            error << "[ED::MODELS::LOADSHAPE] Error while loading shape: heightmap must contain 'image', 'resolution' and 'height'." << std::endl;
        }
    }
    else if (cfg.readArray("compound") || cfg.readArray("group")) // ED YAML ONLY
    {
        geo::CompositeShapePtr composite(new geo::CompositeShape);
        while(cfg.nextArrayItem())
        {
            std::map<std::string, geo::ShapePtr> dummy_shape_cache;
            geo::ShapePtr sub_shape = loadShape(model_path, cfg, dummy_shape_cache, error);
            composite->addShape(*sub_shape, geo::Pose3D::identity());
        }
        cfg.endArray();

        shape = composite;
    }
    else
    {
        error << "[ED::MODELS::LOADSHAPE] Error while loading shape: must contain one of the following 'path', 'heightmap', 'box', 'compound'." << std::endl;
    }

    readPose(cfg, pose);

    if (shape)
    {
        // Transform shape according to pose
        geo::ShapePtr shape_tr(new geo::Shape());
        shape_tr->setMesh(shape->getMesh().getTransformed(pose));
        shape = shape_tr;
    }

    return shape;
}

}
}
