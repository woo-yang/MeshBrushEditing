#include "remeshing.h"
#include "../common/mesh_io.h"

#include <CGAL/Surface_mesh/IO/OFF.h>

namespace alg {

    void Remesher::isotropic_remeshing(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        const std::vector<int>& face_range,
        double target_edge_length,
        unsigned int nb_iter,
        Eigen::Matrix3Xd& V_new,
        Eigen::Matrix3Xi& F_new)
    {
        using namespace common;

        mesh.reserve(V.cols(), V.cols() + F.cols() - 2, F.cols());
        for (int i = 0; i < V.cols(); ++i) {
            mesh.add_vertex(CGAL::Epick::Point_3(V(0, i), V(1, i), V(2, i)));
        }
        for (auto i: face_range) {
            mesh.add_face(CGAL::SM_Vertex_index(F(0, i)), CGAL::SM_Vertex_index(F(1, i)), CGAL::SM_Vertex_index(F(2, i)));
        }
        for (int i = 0; i < F.cols(); ++i) {
            mesh.add_face(CGAL::SM_Vertex_index(F(0, i)), CGAL::SM_Vertex_index(F(1, i)), CGAL::SM_Vertex_index(F(2, i)));
        }
        
        CGAL::IO::write_OFF("../obj/1.off", mesh);

        // give each vertex a name, the default is empty
        Mesh::Property_map<edge_descriptor, bool> is_constrained =
            mesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;

        //detect sharp features
        BOOST_FOREACH(edge_descriptor e, mesh.edges())
        {
            halfedge_descriptor hd = halfedge(e, mesh);
            if (is_border(e, mesh))  continue;

            double angle = CGAL::Mesh_3::dihedral_angle(
                mesh.point(source(hd, mesh)),
                mesh.point(target(hd, mesh)),
                mesh.point(target(next(hd, mesh), mesh)),
                mesh.point(target(next(opposite(hd, mesh), mesh), mesh)));
            if (CGAL::abs(angle) < 100)
                is_constrained[e] = true;
        }

        //remesh

        Mesh::Face_range range(mesh.faces_begin(), mesh.faces_begin() + face_range.size());
        PMP::isotropic_remeshing(
            range,
            target_edge_length,
            mesh,
            PMP::parameters::number_of_iterations(nb_iter)
            .edge_is_constrained_map(is_constrained));


        CGAL::IO::write_OFF("../obj/2.off", mesh);
        mesh.collect_garbage();

        std::vector<double> verts;
        std::vector<int> indices;

        verts.reserve(3 * mesh.number_of_vertices());
        indices.reserve(3 * mesh.number_of_faces());

        for (Mesh::Vertex_index vi : mesh.vertices()) {
            Kernel::Point_3 pt = mesh.point(vi);    
            verts.push_back(pt.x());
            verts.push_back(pt.y());
            verts.push_back(pt.z());
   
        }

        for (Mesh::Face_index fi : mesh.faces()) {
            Mesh::Halfedge_index hfi = mesh.halfedge(fi), end = hfi;
            do {
                Mesh::Vertex_index vi = target(hfi, mesh);
                indices.push_back(vi);
                hfi = mesh.next(hfi);
            } while (hfi != end);
        }

        V_new = Eigen::Map<Eigen::Matrix3Xd>(verts.data(), 3, mesh.number_of_vertices());
        F_new = Eigen::Map<Eigen::Matrix3Xi>(indices.data(), 3, mesh.number_of_faces());

        io::save_triangle_obj("../obj/3.obj", V_new, F_new);
    }



}