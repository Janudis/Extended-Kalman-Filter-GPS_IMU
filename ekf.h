#ifndef EKF_HPP
#define EKF_HPP

#pragma once
#include <Eigen/Dense>

class ExtendedKalmanFilter {
public:
    ExtendedKalmanFilter(const Eigen::Vector3d& x, const Eigen::Matrix3d& P);

    void update(const Eigen::Vector2d& z, const Eigen::Matrix2d& Q);
    void propagate(const Eigen::Vector2d& u, double dt, const Eigen::Matrix3d& R);

//    const Eigen::Vector3d& getState() const;
//    const Eigen::Matrix3d& getCovariance() const;
    Eigen::Vector3d x_;
    Eigen::Matrix3d P_;
};
#endif // EKF_HPP