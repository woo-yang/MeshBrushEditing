#pragma once
#include "../common/mesh_type.h"

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <boost/foreach.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace alg {
	class Remesher {
	public:
		Remesher() = default;
        void Remesher::isotropic_remeshing(
            const Eigen::Matrix3Xd& V,
            const Eigen::Matrix3Xi& F,
            const std::vector<int>& range,
            double target_edge_length,
            unsigned int nb_iter,
            Eigen::Matrix3Xd& V_new,
            Eigen::Matrix3Xi& F_new);

	private:
		common::Mesh  mesh;
	};

}