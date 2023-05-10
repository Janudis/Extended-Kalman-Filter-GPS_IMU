#include "ekf.h"
#include <iostream>
#include "geo_ned.h"


ExtendedKalmanFilter::ExtendedKalmanFilter(const Eigen::Vector3d& x, const Eigen::Matrix3d& P)
        : x_(x), P_(P) {}

void ExtendedKalmanFilter::update(const Eigen::Vector2d& z, const Eigen::Matrix2d& Q) {
    Eigen::Matrix<double, 2, 3> H;
    H << 1.0, 0.0, 0.0,
            0.0, 1.0, 0.0;

    Eigen::Matrix<double, 3, 2> K = P_ * H.transpose() * (H * P_ * H.transpose() + Q).inverse();

    Eigen::Vector2d z_ = x_.head<2>();
    x_ = x_ + K * (z - z_);

    P_ = P_ - K * H * P_;
}

void ExtendedKalmanFilter::propagate(const Eigen::Vector2d& u, double dt, const Eigen::Matrix3d& R) {
    double v = u[0];
    double omega = u[1];
    //std::cout<<"omega"<<omega<<std::endl;
    double theta = x_[2];
    double dtheta = omega * dt; //swsto
    double r = v / omega;
    //std::cout<<"theta  "<<theta<<std::endl;
    //std::cout<<"theta"<<theta<<std::endl;
    double dx = -r * sin(theta) + r * sin(theta + dtheta);
    double dy = r * cos(theta) - r * cos(theta + dtheta);
    //std::cout<<"dx"<<dx<<std::endl;
    x_[0] += dx;
    x_[1] += dy;
    x_[2] += dtheta;

    Eigen::Matrix3d G;
    G << 1.0, 0.0, -r * cos(theta) + r * cos(theta + dtheta),
            0.0, 1.0, -r * sin(theta) + r * sin(theta + dtheta),
            0.0, 0.0, 1.0;

    P_ = G * P_ * G.transpose() + R;
}


//const Eigen::Vector3d& ExtendedKalmanFilter::getState() const {
//    return x_;
//}
//
//const Eigen::Matrix3d& ExtendedKalmanFilter::getCovariance() const {
//    return P_;
//}