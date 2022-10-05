#include "remeshing.h"
#include "../common/mesh_io.h"

#include <CGAL/Surface_mesh/IO/OFF.h>

namespace alg {

    void isotropic_remeshing(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        const std::vector<int>& face_range,
        double target_edge_length,
        unsigned int nb_iter,
        Eigen::Matrix3Xd& V_new,
        Eigen::Matrix3Xi& F_new)
    {
        using namespace common;

        Mesh mesh;
        mesh.reserve(V.cols(), V.cols() + F.cols() - 2, F.cols());
        for (int i = 0; i < V.cols(); ++i) {
            mesh.add_vertex(CGAL::Epick::Point_3(V(0, i), V(1, i), V(2, i)));
        }
        for (int i = 0; i < F.cols(); ++i) {
            mesh.add_face(CGAL::SM_Vertex_index(F(0, i)), CGAL::SM_Vertex_index(F(1, i)), CGAL::SM_Vertex_index(F(2, i)));
        }
        
        CGAL::IO::write_OFF("../obj/1.off", mesh);

        std::vector<Mesh::Face_index> range;
        for (const auto& i : face_range) {
            range.push_back(Mesh::Face_index(i));
        }
        PMP::isotropic_remeshing(range,target_edge_length,mesh);

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