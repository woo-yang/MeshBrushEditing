#include "deformation.h"

#include <set>
#include <unordered_set>
#include <queue>
#include <igl/boundary_loop.h>


namespace alg {

	int Deformer::mesh_deformation(
		const Eigen::Vector3d& dir,
		size_t sv, double r, double ext_ratio)
	{
		double target_length = get_remesh_length(dir, sv, r);

		std::vector<int> range;
		get_remesh_patch(dir, sv, r, ext_ratio, range);

		Eigen::Matrix3Xd re_v;
		Eigen::Matrix3Xi re_f;
		alg::isotropic_remeshing(_v, _f, range, target_length, 5, re_v, re_f);
		_v = re_v; _f = re_f;

		common::get_mesh_vertex_normal(_v, _f, _n);
		common::vertex_vertex_adjacency(_f, _v_v);

		io::save_triangle_obj("../obj/remesh.obj", _v, _f);

		_direction = dir;
		_orth_direction = common::orthogonal_vector(dir);
		_ext_ratio = ext_ratio;

		if (!_is_valid || sv != _source_v || (r - _geo_radius) > 1e-5) {
			_source_v = sv;
			_geo_radius = r;
			get_local_patch(sv, r);
			_is_valid = true;
		}

		Eigen::Matrix3Xd  limit_pos, limit_lap;
		get_local_mapping(limit_pos, limit_lap);

		Eigen::Matrix3Xd lv;
		Eigen::Matrix3Xi lf;
		std::vector<int> l2g;
		get_extend_patch(lv, lf, l2g);

		int num_iv = _geo_r_l2g.size(), num_v = l2g.size();
		std::vector<int> g2l(_v.cols(), -1), l2l(num_iv);
		for (int i = 0; i < num_v; ++i) g2l[l2g[i]] = i;
		for (int i = 0; i < num_iv; ++i) l2l[i] = g2l[_geo_r_l2g[i]];

		Eigen::SparseMatrix<double> L;
		common::cal_cot_laplace(lv, lf, L);
		Eigen::Matrix3Xd lap;
		smooth_buffer_laplace(lv, lf, L, limit_lap, l2l, lap);


		Eigen::SparseMatrix<double> A;
		Eigen::VectorXd b;
		build_linear_system(lv, lf, L, lap, limit_pos, l2l, A, b);

		_solver.analyzePattern(A);
		_solver.factorize(A);
		if (_solver.info() != Eigen::Success)
			std::cout << "Compute failed: " << _solver.info() << std::endl;
		Eigen::VectorXd X = _solver.solve(b);

		for (int i = 0; i < num_v; ++i) {
			for (int j = 0; j < 3; ++j) {
				_v(j, l2g[i]) = X(i + j * num_v);
			}
		}
		return 1;
	}

	double Deformer::get_remesh_length(const Eigen::Vector3d& dir,size_t sv, double r)
	{
		int count = 256 * dir.norm() * r;
		_target_length = sqrt(4 * M_PI * r * r / (count * sqrt(3)));
		return _target_length;
	}

	int  Deformer::get_remesh_patch	(
		const Eigen::Vector3d& dir,
		size_t sv, double r, double ext_ratio,
		std::vector<int>& faces)
	{
		common::vertex_vertex_adjacency(_f, _v_v);
		double avg_length = 0;
		for (const auto v : _v_v[sv]) {
			avg_length += (_v.col(sv) - _v.col(v)).norm();
		}
		avg_length /= _v_v[sv].size();

		double max_d = (1 + ext_ratio) * r;
		int max_level = 2* max_d / avg_length ;

		Eigen::VectorXi flag = Eigen::VectorXi::Constant(_v.cols(), 0);
		std::queue<std::pair<int, int>> que;
		que.push({ sv,0 });
		
		while (!que.empty()) {
			int u = que.front().first;
			int l = que.front().second;
			que.pop();
			flag[u] = 1;
			for (auto v : _v_v[u]) {
				if (flag[v] == 0 && l < max_level &&
					common::point_line_distance(_v.col(v),_v.col(sv),_v.col(sv)+ dir)<max_d) {
					que.push({ v,l +1 });
				}
			}
		}

		for (int i = 0; i < _f.cols(); ++i) {
			if (flag[_f(0, i)] || flag[_f(1, i)] || flag[_f(2, i)]) {
				faces.push_back(i);
			}
		}

		return 1;
	}


	int Deformer::get_local_patch(size_t sv, double r)
	{
		int num_v = _v_v.size();
		Eigen::VectorXd min_distance = Eigen::VectorXd::Constant(num_v, 1, std::numeric_limits<double>::infinity());
		Eigen::VectorXi previous = Eigen::VectorXi::Constant(num_v, 1, -1);
		Eigen::VectorXi visited = Eigen::VectorXi::Constant(num_v, 1, -1);
		Eigen::Matrix2Xd uv = Eigen::Matrix2Xd::Constant(2, num_v, std::numeric_limits<double>::infinity());
		Eigen::VectorXi flag = Eigen::VectorXi::Constant(num_v, 1, 0);

		min_distance[sv] = 0;
		visited[sv] = 0;

		std::set<std::pair<double, int>> vertex_queue;
		vertex_queue.insert({ min_distance[sv],sv });

		while (!vertex_queue.empty())
		{
			int u = vertex_queue.begin()->second;
			double dist = vertex_queue.begin()->first;
			vertex_queue.erase(vertex_queue.begin());
			visited[u] = 1;
			//calc uv 
			Eigen::Vector3d sv_pos = _v.col(sv);
			Eigen::Vector3d sv_e1 = _n.col(sv);
			Eigen::Vector3d sv_e2 = common::vector_plane_projection(sv_e1, _orth_direction).normalized();

			if (u == sv) {
				uv.col(u).setZero();
			}
			else {
				Eigen::Vector3d q_pos = _v.col(previous[u]), r_pos = _v.col(u);
				Eigen::Vector3d q_e1 = _n.col(previous[u]);
				Eigen::Vector3d q_e2 = common::vector_plane_projection(q_e1, _orth_direction).normalized();
				Eigen::Vector3d proj_r = common::vector_plane_projection(q_e1, q_pos - r_pos).normalized();

				double polar_angle = common::safe_acos(q_e2.dot(proj_r));
				if (q_e2.cross(proj_r).dot(q_e1) < 0) {
					polar_angle = 2 * M_PI - polar_angle;
				}
				double geodesic_dis = (q_pos - r_pos).norm();
				Eigen::Vector2d Urq(geodesic_dis * cos(polar_angle), geodesic_dis * sin(polar_angle));
				Eigen::Matrix2d M;
				
				if (previous[u] == sv) {
					M.setIdentity();
				}
				else {
					Eigen::AngleAxisd aa(common::safe_acos(q_e1.dot(sv_e1)), (q_e1.cross(sv_e1)).normalized());
					q_e2 = (aa.toRotationMatrix() * q_e2).normalized();

					double theta = common::safe_acos(q_e2.dot(sv_e2));

					if (sv_e2.cross(q_e2).dot(sv_e1) < 0) {
						theta = 2 * M_PI - theta;
					}
					M << cos(theta), -sin(theta), sin(theta), cos(theta);
				}
				uv.col(u) = uv.col(previous[u]) + M * Urq;
			}

			const std::vector<int>& neighbors = _v_v[u];
			for (auto neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
			{
				int v = *neighbor_iter;
				double distance_through_u = dist + (_v.col(u) - _v.col(v)).norm();
				if (visited[v] != 1 && distance_through_u < min_distance[v]) {
					if (visited[v] == 0)vertex_queue.erase({ min_distance[v], v });
					min_distance[v] = distance_through_u;
					previous[v] = u;
					if (visited[v] == 0)vertex_queue.insert({ min_distance[v], v });
				}
			}
			if (uv.col(u).squaredNorm() < r * r)
			{
				flag[u] = 1;
				for (auto neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
				{
					int v = *neighbor_iter;
					if (visited[v] == -1) {
						vertex_queue.insert({ min_distance[v],v });
						visited[v] = 0;
					}
				}
			}
		}

		get_mesh_from_vertex(flag, _geo_r_lv, _geo_r_lf, _geo_r_l2g);

		_geo_r_lu.resize(2, _geo_r_l2g.size());
		for (int i = 0; i < _geo_r_l2g.size(); ++i) {
			if (_geo_r_l2g[i] == sv)_geo_r_lsv = i;
			_geo_r_lu.col(i) = uv.col(_geo_r_l2g[i]);
		}
		io::save_triangle_obj("../obj/patch.obj", _geo_r_lv, _geo_r_lf, _geo_r_lu);
		return 1;
	}

	int Deformer::get_local_mapping(
		Eigen::Matrix3Xd& limit_pos,
		Eigen::Matrix3Xd& limit_lap)
	{
		int num_v = _geo_r_lv.cols();
		limit_pos.resize(3, num_v);
		limit_lap.resize(3, num_v);

		double H = _direction.norm();
		Eigen::Vector3d e1 = _orth_direction.normalized();
		Eigen::Vector3d e2 = -_direction.normalized();

		Eigen::Vector2d ori_uv = _geo_r_lu.col(_geo_r_lsv);
		Eigen::Vector3d origin = _geo_r_lv.col(_geo_r_lsv);

		Eigen::VectorXi bnd;
		igl::boundary_loop(_geo_r_lf.transpose(), bnd);

		std::map<double, int> rb_tree;
		for (auto i : bnd) {
			double theta = common::cartesian_to_polar(_geo_r_lu.col(i))(1);
			rb_tree.insert({ theta,i });
		}
		for (int i = 0; i < num_v; ++i) {
			if (i == _geo_r_lsv) {
				limit_pos.col(i) << 0, -_geo_radius, 0;
				continue;
			}
			Eigen::Vector2d uv = _geo_r_lu.col(i);
			double theta = common::cartesian_to_polar(uv)(1);

			int p1, p2;
			auto iter = rb_tree.upper_bound(theta);
			if (iter == rb_tree.begin()) {
				p2 = iter->second;
				p1 = (--rb_tree.end())->second;
			}
			else if (iter == rb_tree.end()) {
				p2 = rb_tree.begin()->second;
				p1 = (--rb_tree.end())->second;
			}
			else {
				p2 = iter->second;
				p1 = (--iter)->second;
			}

			Eigen::Vector2d uv1 = _geo_r_lu.col(p1), uv2 = _geo_r_lu.col(p2);
			Eigen::Vector2d para;
			common::ray_line_intersection(ori_uv, uv - ori_uv, uv1, uv2, para);
			double k = 1. / para[0];
			Eigen::Vector3d inter_3d = para[1] * _geo_r_lv.col(p1) + (1 - para[1]) * _geo_r_lv.col(p2);
			double h = H - _geo_radius - (inter_3d - origin).dot(_direction.normalized());
			double t = k * (M_PI * _geo_radius / 2 + h);

			Eigen::Vector3d pos_2d;
			if (t <= M_PI * _geo_radius / 2) {
				t += 3 * M_PI * _geo_radius / 2;
				pos_2d << _geo_radius * cos(t / _geo_radius), _geo_radius* sin(t / _geo_radius), 0;
			}
			else {
				t -= M_PI * _geo_radius / 2;
				pos_2d << _geo_radius, t, 0;
			}
			limit_pos.col(i) = Eigen::AngleAxisd(theta - M_PI, Eigen::Vector3d(0, 1, 0)).
				toRotationMatrix() * pos_2d;
		}

		Eigen::Matrix3d rotate;
		rotate << e1, e2, e1.cross(e2).normalized();
		limit_pos = rotate * limit_pos;
		Eigen::Vector3d translate = origin + (H - _geo_radius) * _direction.normalized();
		limit_pos.colwise() += translate;

		Eigen::SparseMatrix<double> L;
		common::cal_cot_laplace(limit_pos, _geo_r_lf, L);
		limit_lap = (L * limit_pos.transpose()).transpose();

		io::save_triangle_obj("../obj/local.obj", limit_pos, _geo_r_lf);
		return 1;
	}

	int Deformer::get_extend_patch(
		Eigen::Matrix3Xd& lv,
		Eigen::Matrix3Xi& lf,
		std::vector<int>& l2g)
	{
		Eigen::VectorXi flag = Eigen::VectorXi::Constant(_v.cols(), 0);

		for (int i = 0; i < _geo_r_l2g.size(); ++i) {
			flag[_geo_r_l2g[i]] = 1;
		}

		std::queue<std::pair<int, int>> que;
		Eigen::VectorXi inner_bnd;
		igl::boundary_loop(_geo_r_lf.transpose(), inner_bnd);

		for (const auto& i : inner_bnd) {
			que.push({ _geo_r_l2g[i],0 });
		}

		int max_level = 2 * _ext_ratio * _geo_radius / _target_length;
		while (!que.empty()) {
			int i = que.front().first;
			int level = que.front().second;
			que.pop();
			flag[i] = 1;
			for (const auto& j : _v_v[i]) {
				if (flag[j] == 0 && level < max_level &&
					common::point_line_distance(_v.col(j), _v.col(_source_v),
						_v.col(_source_v) + _direction) < (1 + _ext_ratio) * _geo_radius)
				{
					que.push({ j,level + 1 });
				}
			}
		}
		get_mesh_from_vertex(flag, lv, lf, l2g);
		io::save_triangle_obj("../obj/patch2.obj", lv, lf);
		return 1;
	}

	int Deformer::get_mesh_from_vertex(
		const Eigen::VectorXi& flag,
		Eigen::Matrix3Xd& lv,
		Eigen::Matrix3Xi& lf,
		std::vector<int>& l2g)
	{
		std::vector<int> faces;
		faces.reserve(lv.cols() / 3);

		for (int i = 0; i < _f.cols(); ++i) {
			if (flag[_f(0, i)] && flag[_f(1, i)] && flag[_f(2, i)]) {
				faces.push_back(_f(0, i));
				faces.push_back(_f(1, i));
				faces.push_back(_f(2, i));
			}
		}

		l2g = faces;
		std::sort(l2g.begin(), l2g.end());
		l2g.erase(std::unique(l2g.begin(), l2g.end()), l2g.end());

		std::vector<int> g2l(_v.cols(), -1);
		for (int i = 0; i < l2g.size(); ++i) {
			g2l[l2g[i]] = i;
		}
		std::for_each(faces.begin(), faces.end(), [&](int& f) {f = g2l[f]; });
		lf = Eigen::Map<Eigen::Matrix3Xi>(faces.data(), 3, faces.size() / 3);
		lv.resize(3, l2g.size());

		for (int i = 0; i < lv.cols(); ++i) {
			lv.col(i) = _v.col(l2g[i]);
		}

		return 1;
	}


	int Deformer::smooth_buffer_laplace(
		const Eigen::Matrix3Xd& lv,
		const Eigen::Matrix3Xi& lf,
		const Eigen::SparseMatrix<double>& L,
		const Eigen::Matrix3Xd& limit_lap,
		const std::vector<int>& l2l,
		Eigen::Matrix3Xd& laplace)
	{
		laplace = (L * lv.transpose()).transpose();

		Eigen::VectorXi bnd;
		igl::boundary_loop(lf.transpose(), bnd);

		int num_v = lv.cols(), num_iv = _geo_r_lv.cols();
		std::vector<bool> fix(num_v, false);
		for (int i = 0; i < num_iv; ++i) {
			laplace.col(l2l[i]) = limit_lap.col(i);
			fix[l2l[i]] = true;
		}
		for (const auto& i : bnd) {
			fix[i] = true;
		}

		std::vector<std::vector<int>> v_v;
		common::vertex_vertex_adjacency(lf, v_v);

		int nb_iter = 10;
		for (int m = 0; m < nb_iter; ++m) {
			for (int i = 0; i < num_v; ++i) {
				if (fix[i])continue;
				Eigen::Vector3d nerb_lap = Eigen::Vector3d::Zero();
				for (auto j : v_v[i]) {
					nerb_lap += laplace.col(j);
				}
				nerb_lap /= v_v[i].size();
				laplace.col(i) = 1. / 2 * laplace.col(i) + 1. / 2 * nerb_lap;
			}
		}
		return 1;
	}

	int Deformer::build_linear_system(
		const Eigen::Matrix3Xd& lv,
		const Eigen::Matrix3Xi& lf,
		const Eigen::SparseMatrix<double>& L,
		const Eigen::Matrix3Xd& lap,
		const Eigen::Matrix3Xd& limit_pos,
		const std::vector<int>& l2l,
		Eigen::SparseMatrix<double>& A,
		Eigen::VectorXd& b)
	{
		Eigen::SparseMatrix<double> A1, A2, A3;
		Eigen::VectorXd b1, b2, b3;

		int num_v = lv.cols(), num_iv = _geo_r_lv.cols();
		A1.resize(num_v * 3, num_v * 3);
		b1 = Eigen::VectorXd::Zero(num_v * 3);

		std::vector<Eigen::Triplet<double>> triple;

		Eigen::VectorXi bnd;
		igl::boundary_loop(lf.transpose(), bnd);
		for (const auto& i : bnd) {
			for (int j = 0; j < 3; ++j) {
				triple.push_back({ i + j * num_v,i + j * num_v,1 });
				b1(i + j * num_v) = lv(j, i);
			}
		}
		A1.setFromTriplets(triple.begin(), triple.end());

		triple.clear();
		A2.resize(num_v * 3, num_v * 3);
		for (int i = 0; i < L.outerSize(); ++i)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				int row = it.row(), col = it.col();
				triple.push_back({ row,col,it.value() });
				triple.push_back({ row + num_v,col + num_v,it.value() });
				triple.push_back({ row + 2 * num_v,col + 2 * num_v,it.value() });
			}
		}
		A2.setFromTriplets(triple.begin(), triple.end());
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m = lap;
		b2 = Eigen::Map<Eigen::VectorXd>(m.data(), m.size());

		triple.clear();
		A3.resize(num_v, num_v * 3);
		b3 = Eigen::VectorXd::Zero(num_v);

		Eigen::Matrix3Xd mapping_ln;
		common::get_mesh_vertex_normal(limit_pos, _geo_r_lf, mapping_ln);
		for (int i = 0; i < num_iv; ++i) {
			int l_id = l2l[i];
			Eigen::Vector3d n = mapping_ln.col(i);
			triple.push_back({ l_id,l_id,n[0] });
			triple.push_back({ l_id,l_id + num_v,n[1] });
			triple.push_back({ l_id,l_id + 2 * num_v,n[2] });
			b3.row(l_id) = mapping_ln.col(i).transpose() * limit_pos.col(i);
		}
		A3.setFromTriplets(triple.begin(), triple.end());

		A =  A1.transpose() * A1 +  A2.transpose() * A2 +  A3.transpose() * A3;
		b =  A1.transpose() * b1 +  A2.transpose() * b2 +  A3.transpose() * b3;
		return 1;
	}
}