#pragma once
#include <math.h>
#include <Eigen/Dense>

namespace common {
	double cross_2d(
		const Eigen::Vector2d& a,
		const Eigen::Vector2d& b);

	double safe_acos(double x)
	{
		if (x < -1.0) return M_PI;
		else if (x > 1.0) return 0;
		return acos(x);
	}

	bool point_in_triangle(
		const Eigen::Vector2d& p0, 
		const Eigen::Vector2d& p1, 
		const Eigen::Vector2d& p2,
		const Eigen::Vector2d& p, 
		double& u, double& v);

	double point_line_distance(
		const Eigen::Vector3d& p,
		const Eigen::Vector3d& p0,
		const Eigen::Vector3d& p1);

	Eigen::Vector3d point_line_projection(
		const Eigen::Vector3d& p,
		const Eigen::Vector3d& p0,
		const Eigen::Vector3d& p1);

	Eigen::Vector3d vector_plane_projection(
		const Eigen::Vector3d& n,
		const Eigen::Vector3d& vec);

	bool ray_line_intersection(
		const Eigen::Vector2d& ray_origin,
		const Eigen::Vector2d& ray_dir,
		const Eigen::Vector2d& p1,
		const Eigen::Vector2d& p2,
		Eigen::Vector2d& t);

	Eigen::Vector3d orthogonal_vector(
		const Eigen::Vector3d& n);

	Eigen::Vector2d polar_to_cartesian(
		const Eigen::Vector2d& polar_coord);

	Eigen::Vector2d cartesian_to_polar(
		const Eigen::Vector2d& coord);

	Eigen::Vector3d cartesian_to_sphere(
		const Eigen::Vector3d& coord);

	Eigen::Vector3d cartesian_to_cylinder(
		const Eigen::Vector3d& coord);

	
}