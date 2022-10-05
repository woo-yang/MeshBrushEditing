#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace common {
	int cal_cot_laplace(
		const Eigen::Matrix3Xd& V,
		const Eigen::Matrix3Xi& F,
		Eigen::SparseMatrix<double>& L);

	void vertex_triangle_adjacency(
		const Eigen::Matrix3Xi& F, 
		std::vector<std::vector<int>>& v_f);

	void vertex_vertex_adjacency(
		const Eigen::Matrix3Xi& F, 
		std::vector<std::vector<int>>& v_v);

	void triangle_triangle_adjacency(
		const Eigen::Matrix3Xi& F,
		std::vector < std::vector<int>>& f_f);

	void get_vertex_normal(
		const Eigen::Matrix3Xd& V,
		const Eigen::Matrix3Xi& F,
		size_t sv,
		Eigen::Vector3d& n);

	void get_mesh_vertex_normal(
		const Eigen::Matrix3Xd& V,
		const Eigen::Matrix3Xi& F,
		Eigen::Matrix3Xd& N );
}