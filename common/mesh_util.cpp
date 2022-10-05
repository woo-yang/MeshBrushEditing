#include "mesh_util.h"
#include <igl/triangle_triangle_adjacency.h>

namespace common {

    int cal_cot_angles(
        const Eigen::MatrixXd& V,
        const Eigen::Matrix3Xi& F,
        Eigen::Matrix3Xd& cot_angles);

    int cal_cot_laplace(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        Eigen::SparseMatrix<double>& L)
    {
        Eigen::Matrix3Xd cot_angles;
        cal_cot_angles(V, F, cot_angles);
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(F.cols() * 9);
        for (size_t j = 0; j < F.cols(); ++j) {
            const Eigen::Vector3i& fv = F.col(j);
            const Eigen::Vector3d& ca = cot_angles.col(j);
            for (size_t vi = 0; vi < 3; ++vi) {
                const size_t j1 = (vi + 1) % 3;
                const size_t j2 = (vi + 2) % 3;
                const int fv0 = fv[vi];
                const int fv1 = fv[j1];
                const int fv2 = fv[j2];
                triple.push_back(Eigen::Triplet<double>(fv0, fv0, ca[j1] + ca[j2]));
                triple.push_back(Eigen::Triplet<double>(fv0, fv1, -ca[j2]));
                triple.push_back(Eigen::Triplet<double>(fv0, fv2, -ca[j1]));
            }
        }
        L.resize(V.cols(), V.cols());
        L.setFromTriplets(triple.begin(), triple.end());
        return 1;
    }

    int cal_cot_angles(
        const Eigen::MatrixXd& V,
        const Eigen::Matrix3Xi& F,
        Eigen::Matrix3Xd& cot_angles) 
    {
        cot_angles.resize(3, F.cols());
        for (size_t j = 0; j < F.cols(); ++j) {
            const Eigen::Vector3i& fv = F.col(j);
            for (size_t vi = 0; vi < 3; ++vi) {
                const Eigen::VectorXd& p0 = V.col(fv[vi]);
                const Eigen::VectorXd& p1 = V.col(fv[(vi + 1) % 3]);
                const Eigen::VectorXd& p2 = V.col(fv[(vi + 2) % 3]);
                const double angle = std::acos(std::max(-1.0,
                    std::min(1.0, (p1 - p0).normalized().dot((p2 - p0).normalized()))));
                cot_angles(vi, j) = 1.0 / std::tan(angle);
            }
        }
        return 1;
    }

    void vertex_triangle_adjacency(
        const Eigen::Matrix3Xi& F, 
        std::vector<std::vector<int>>& v_f)
    {
        v_f.resize(F.maxCoeff() + 1);
        for (int i = 0; i < F.cols(); ++i)
            for (int j = 0; j < 3; ++j)
                v_f[F(j, i)].push_back(i);
    }

    void vertex_vertex_adjacency(
        const Eigen::Matrix3Xi& F,
        std::vector<std::vector<int>>& v_v)
    {
        v_v.clear();
        v_v.resize(F.maxCoeff() + 1);
        for (int i = 0; i < F.cols(); i++) {
            for (int j = 0; j < F.rows(); j++) {
                int s = F(j, i);
                int d = F((j + 1) % F.rows(), i);
                v_v.at(s).push_back(d);
                v_v.at(d).push_back(s);
            }
        }
        // Remove duplicates
        for (int i = 0; i < (int)v_v.size(); ++i) {
            std::sort(v_v[i].begin(), v_v[i].end());
            v_v[i].erase(std::unique(v_v[i].begin(), v_v[i].end()), v_v[i].end());
        }

    }

    void triangle_triangle_adjacency(
        const Eigen::Matrix3Xi& F,
        std::vector < std::vector<int>>& f_f)
    {
        f_f.resize(F.cols());
        std::vector<std::vector<std::vector<int>>> tt;
        igl::triangle_triangle_adjacency(F.transpose(), tt);

        for (int i = 0; i < tt.size(); ++i) {
            for (int j = 0; j < tt[i].size(); ++j) {
                for (auto f : tt[i][j]) {
                    f_f[i].push_back(f);
                }
            }
        }
    }

    void get_vertex_normal(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        size_t sv,
        Eigen::Vector3d& n)
    {
        Eigen::Matrix3Xd N;
        get_mesh_vertex_normal(V, F, N);
        n = N.col(sv);
    }


    void get_mesh_vertex_normal(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        Eigen::Matrix3Xd& N)
    {
        N.setZero(3, V.cols());
        for (int i = 0; i < F.cols(); ++i) {
            auto v1 = V.col(F(1, i)) - V.col(F(0, i));
            auto v2 = V.col(F(2, i)) - V.col(F(0, i));
            Eigen::Vector3d n = v1.cross(v2);
            for (int j = 0; j < 3; ++j) {
                N.col(F(j, i)) += n;
            }
        }
        N.colwise().normalize();
    }

    void get_cylinder_laplace(double r, Eigen::Vector3d& lap)
    {
        Eigen::Matrix3Xd V(3, 7);
        Eigen::Matrix3Xi F(3, 6);
        double sin_10 = std::sin(10 * M_PI / 180);
        double cos_10 = std::cos(10 * M_PI / 180);
        double sin_350 = std::sin(350 * M_PI / 180);
        double cos_350 = std::cos(350 * M_PI / 180);
        double h = r * sin_10;
        V.col(0) << r, 0, 0;
        V.col(1) << r, h, 0;
        V.col(2) << r * cos_350, h / 2., -r * sin_350;
        V.col(3) << r * cos_350, -h / 2., -r * sin_350;
        V.col(4) << r, -h, 0;
        V.col(5) << r * cos_10, -h / 2., -r * sin_10;
        V.col(6) << r * cos_10, h / 2., -r * sin_10;

        for (int i = 0; i < F.cols(); ++i) {
            F.col(i) << 0, i + 1, i + 2;
        }
        F(2, 5) = 1;

        Eigen::SparseMatrix<double> L;
        common::cal_cot_laplace(V, F, L);

        auto m = L * V.transpose();
        lap = m.row(0).transpose();
    }

    void get_sphere_top_laplace(double r, Eigen::Vector3d& lap)
    {
        Eigen::Matrix3Xd V(3, 7);
        Eigen::Matrix3Xi F(3, 6);
        double sin_10 = std::sin(10 * M_PI / 180);
        double cos_10 = std::cos(10 * M_PI / 180);
        V.col(0) << 0, r, 0;
        for (int i = 1; i < V.cols(); ++i) {
            double phi = (i - 1) * 60 * M_PI / 180;
            V.col(i) <<
                r * sin_10 * std::cos(phi),
                r* cos_10,
                -r * sin_10 * std::sin(phi);
        }
        for (int i = 0; i < F.cols(); ++i) {
            F.col(i) << 0, i + 1, i + 2;
        }
        F(2, 5) = 1;
        //io::save_triangle_obj("../obj/test.obj",V, F);
        Eigen::SparseMatrix<double> L;
        common::cal_cot_laplace(V, F, L);

        auto m = L * V.transpose();
        lap = m.row(0).transpose();
    }
}