#pragma once
#include "../common/mesh_type.h"

namespace alg {

    void isotropic_remeshing(
        const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        const std::vector<int>& range,
        double target_edge_length,
        unsigned int nb_iter,
        Eigen::Matrix3Xd& V_new,
        Eigen::Matrix3Xi& F_new);
}