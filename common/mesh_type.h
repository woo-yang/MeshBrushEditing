#pragma once

#include <Eigen/Dense>

#include <gmp/gmp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

namespace common {
	namespace PMP = CGAL::Polygon_mesh_processing;

	using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
	using vertex_descriptor = Mesh::Vertex_index;
	using face_descriptor =  Mesh::Face_index;
	using edge_descriptor = Mesh::Edge_index;
	using halfedge_descriptor = Mesh::Halfedge_index;
}