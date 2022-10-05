#pragma once
#include "../common/mesh_type.h"
#include "../common/mesh_util.h"
#include "../common/mesh_io.h"
#include "../common/math_util.h"

#include "remeshing.h"
#include <Eigen/Dense>

namespace alg {
	class Deformer {
	public:
		Deformer(const Eigen::Matrix3Xd& V, const Eigen::Matrix3Xi& F)
			:_v(V), _f(F) {};

		~Deformer() = default;

	public:
		int mesh_deformation(
			const Eigen::Vector3d& dir,
			size_t sv, double r, double ext_radio);

		void get_result_mesh(Eigen::Matrix3Xd& V, Eigen::Matrix3Xi& F) {
			V = _v; F = _f;
		}

	private:
		double get_remesh_length(const Eigen::Vector3d& dir,size_t sv, double r);

		int get_remesh_patch(
			const Eigen::Vector3d& dir,size_t sv, double r, 
			double ext_ratio,std::vector<int>& faces);

		int get_local_patch(size_t sv, double r);

		int get_local_mapping(Eigen::Matrix3Xd& limit_pos, Eigen::Matrix3Xd& limit_lap);

		int get_extend_patch(Eigen::Matrix3Xd& lv, Eigen::Matrix3Xi& lf, std::vector<int>& l2g);

		int get_mesh_from_vertex(
			const Eigen::VectorXi& flag, 
			Eigen::Matrix3Xd& lv, Eigen::Matrix3Xi& lf, 
			std::vector<int>& l2g);

		int Deformer::smooth_buffer_laplace(
			const Eigen::Matrix3Xd& lv,const Eigen::Matrix3Xi& lf,
			const Eigen::SparseMatrix<double>& L,const Eigen::Matrix3Xd& limit_lap,
			const std::vector<int>& l2l,Eigen::Matrix3Xd& laplace);

	private: //data
		Eigen::Matrix3Xd _v;
		Eigen::Matrix3Xi _f;
		Eigen::Matrix3Xd _n;
		std::vector<std::vector<int>> _v_v;

		//cache data
		size_t _source_v;
		double _geo_radius;
		double _ext_ratio;
		double _target_length;

		Eigen::Vector3d _direction;
		Eigen::Vector3d _orth_direction;

		//if _source_v & _geo_radius has not changed,
		//the following data does not have to be updated.
		bool _is_valid = false;
		size_t _geo_r_lsv;
		Eigen::Matrix2Xd _geo_r_lu;
		Eigen::Matrix3Xd _geo_r_lv;
		Eigen::Matrix3Xi _geo_r_lf;
		std::vector<int> _geo_r_l2g;

	private:

		int Deformer::build_linear_system(
			const Eigen::Matrix3Xd& lv,
			const Eigen::Matrix3Xi& lf,
			const Eigen::SparseMatrix<double>& L,
			const Eigen::Matrix3Xd& lap,
			const Eigen::Matrix3Xd& limit_pos,
			const std::vector<int>& l2l,
			Eigen::SparseMatrix<double>& A,
			Eigen::VectorXd& b);

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   _solver;

	};

}