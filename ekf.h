#ifndef EKF_HPP
#define EKF_HPP

#pragma once
#include <Eigen/Dense>

class ExtendedKalmanFilter {
public:
    ExtendedKalmanFilter() = default;
    void initialize(const Eigen::Vector3d& x, double xy_obs_noise_std, double yaw_rate_noise_std, double forward_velocity_noise_std, double initial_yaw_std);
    void update(const Eigen::Vector2d& z);
    void propagate(const Eigen::Vector2d& u, double dt);

//    const Eigen::Vector3d& getState() const;
//    const Eigen::Matrix3d& getCovariance() const;
    Eigen::Vector3d x_;
    Eigen::Matrix3d P_;
    Eigen::Matrix3d P;

    // Measurement error covariance Q
    Eigen::Matrix2d Q;

    // state transition noise covariance R
    Eigen::Matrix3d R;
};
#endif // EKF_HPP