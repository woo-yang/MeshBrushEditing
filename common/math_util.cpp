#include "mesh_util.h"
#include <igl/triangle_triangle_adjacency.h>

namespace common {

    double cross_2d(
        const Eigen::Vector2d& a,
        const Eigen::Vector2d& b)
    {
        return a[0] * b[1] - b[0] * a[1];
    }

    bool point_in_triangle(
        const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
        const Eigen::Vector2d& p, double& u, double& v)
    {
        Eigen::MatrixX2d tri_vec(3, 2);
        tri_vec.row(0) = p2 - p0;
        tri_vec.row(1) = p1 - p0;
        tri_vec.row(2) = p - p0;

        double dot[2][3];
        for (int i = 0; i < 2; i++)
            for (int j = i; j < 3; j++)
                dot[i][j] = tri_vec.row(i).dot(tri_vec.row(j));

        double inver_deno = 1 / (dot[0][0] * dot[1][1] - dot[0][1] * dot[0][1]);

        u = (dot[1][1] * dot[0][2] - dot[0][1] * dot[1][2]) * inver_deno;
        if (u < 0 || u > 1)return false;

        v = (dot[0][0] * dot[1][2] - dot[0][1] * dot[0][2]) * inver_deno;
        if (v < 0 || v > 1)return false;

        return u + v <= 1;
    }

    double point_line_distance(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& p0,
        const Eigen::Vector3d& p1)
    {
       return ((p - p0).cross(p - p1)).norm() / (p1 - p0).norm();
    }

    Eigen::Vector3d point_line_projection(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& p0,
        const Eigen::Vector3d& p1)
    {
        return p + (p - p0).dot(p1 - p0) / (p1 - p0).dot(p1 - p0) * (p1 - p0);
    }



    Eigen::Vector3d vector_plane_projection(
        const Eigen::Vector3d& n,
        const Eigen::Vector3d& vec)
    {
        Eigen::Vector3d norm = n.normalized();
        return vec - vec.dot(norm) * norm;
    }

    Eigen::Vector3d orthogonal_vector(
        const Eigen::Vector3d& n)
    {
        Eigen::Vector3d vec = Eigen::Vector3d::Constant(1);
        if (n[0] != 0) vec[0] = -(n[1] + n[2]) / n[0];
        else if (n[1] != 0) vec[1] = -(n[0] + n[2]) / n[1];
        else if (n[2] != 0) vec[2] = -(n[0] + n[0]) / n[2];
        return vec.normalized();
    }

    Eigen::Vector2d polar_to_cartesian(
        const Eigen::Vector2d& polar_coord)
    {
        Eigen::Vector2d coord;
        coord << polar_coord[0] * cos(polar_coord[1]), polar_coord[0] * sin(polar_coord[1]);
        return coord;
    }

    Eigen::Vector2d cartesian_to_polar(
        const Eigen::Vector2d& coord)
    {
        Eigen::Vector2d polar_coord;
        double r = sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
        double theta;
        if (coord[0]==0) {
            theta = M_PI;
        }
        else {
            theta = atan(coord[1] / coord[0]);
        }
        if (coord[0] < 0) theta += M_PI;
        else if (coord[1] < 0)theta += 2 * M_PI;

        polar_coord << r, theta;
        return polar_coord;
    }

    Eigen::Vector3d cartesian_to_sphere(
        const Eigen::Vector3d& coord)
    {
        Eigen::Vector3d sphere_coord;
        double r = coord.norm();
        double theta = std::acos(coord[1] / r);
        double phi;
        if (coord[0] == 0) {
            if (coord[2] > 0)phi = 3 * M_PI / 2.;
            else if (coord[2] < 0)phi = M_PI / 2.;
            else phi = 0;
        }
        else {
            phi = std::atan(-coord[2] / coord[0]);
            if (coord[0] < 0) phi += M_PI;
            else if (coord[1] > 0)phi += 2 * M_PI;
        }

        sphere_coord << r, theta, phi;
        return sphere_coord;
    }

    Eigen::Vector3d cartesian_to_cylinder(
        const Eigen::Vector3d& coord)
    {
        Eigen::Vector3d cylinder_coord;
        double r = std::sqrt(coord[0] * coord[0] + coord[2] * coord[2]);
        double h = coord[1];
        double phi;
        if (coord[0] == 0) {
            if (coord[2] > 0)phi = 3 * M_PI / 2.;
            else if (coord[2] < 0)phi = M_PI / 2.;
            else phi = 0;
        }
        else {
            phi = std::atan(-coord[2] / coord[0]);
            if (coord[0] < 0) phi += M_PI;
            else if (coord[1] > 0)phi += 2 * M_PI;
        }
        cylinder_coord << r, phi, h;
        return cylinder_coord;
    }

    bool ray_line_intersection(
        const Eigen::Vector2d& ray_origin,
        const Eigen::Vector2d& ray_dir,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2,
        Eigen::Vector2d& t)
    {

        Eigen::Vector2d v1 = ray_origin - p1;
        Eigen::Vector2d v2 = p2 - p1;
        Eigen::Vector2d v3 (-ray_dir[1], ray_dir[0]);

        float dot = v2.dot(v3);
        if (abs(dot) < 1e-5)
            return false;

        float t1 = cross_2d(v2, v1) / dot;
        float t2 = v1.dot(v3) / dot;

        if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)) 
        {
            t << t1, t2;
            return true;
        }
        return false;
    }
}