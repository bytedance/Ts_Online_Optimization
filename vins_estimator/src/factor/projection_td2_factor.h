#pragma once

#include <ros/assert.h>
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../parameters.h"

class ProjectionTd2Factor : public ceres::SizedCostFunction<2, 7, 7, 7, 1, 1>
{
  public:
    ProjectionTd2Factor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j,
    				   const double _td_i, const double _td_j, const Eigen::Vector3d &_veli, const Eigen::Vector3d &_velj, const Eigen::Vector3d &_un_gyr_i, const Eigen::Vector3d &_un_gyr_j);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Vector3d pts_i, pts_j;
    double td_i, td_j;
    Eigen::Matrix<double, 2, 3> tangent_base;
    double row_i, row_j;
    static Eigen::Matrix2d sqrt_info;
    static double sum_t;
    Eigen::Vector3d veli, velj;
    Eigen::Vector3d un_gyr_i;
    Eigen::Vector3d un_gyr_j;
};
