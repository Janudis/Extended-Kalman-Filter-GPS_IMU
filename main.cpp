#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include "ekf.h"
#include "geo_ned.h"
using namespace std;
using namespace Eigen;

int main() {
    std::vector<std::array<double, 3>> gt_trajectory_lla; // [longitude(deg), latitude(deg), altitude(meter)] x N
    std::vector<double> gt_yaws; // [yaw_angle(rad),] x N
    std::vector<double> obs_yaw_rates; // [vehicle_yaw_rate(rad/s),] x N
    std::vector<double> obs_forward_velocities; // [vehicle_forward_velocity(m/s),] x N
    std::vector<double> ts;

    std::ifstream raw_data("localization_log_drivepark_3_python.csv");
    std::string line, value;

// Read CSV file and store values in vectors
    std::getline(raw_data, line); // Skip header line

    long long int first_timestamp = 0;
    bool first_line = true;

    while (std::getline(raw_data, line)) {
        std::stringstream line_stream(line);
        std::vector<double> row;

        std::getline(line_stream, value, ',');
        long long int timestamp = std::stoll(value);

        if (first_line) {
            first_timestamp = timestamp;
            first_line = false;
        }
        double elapsed_time = (timestamp - first_timestamp) / 1000000.0;
        ts.push_back(elapsed_time); // Store the elapsed time in seconds
        while (std::getline(line_stream, value, ',')) {
            row.push_back(std::stold(value));
        }
        gt_trajectory_lla.push_back({row[1], row[0], row[2]});
        gt_yaws.push_back(deg2rad(row[9]));
        obs_yaw_rates.push_back(deg2rad(row[12]));
        obs_forward_velocities.push_back(row[3]);
    }
//    while (std::getline(raw_data, line)) { //gia drivepark_dim
//        std::stringstream line_stream(line);
//        std::vector<double> row;
//        while (std::getline(line_stream, value, ',')) {
//            if (value.find("E") != std::string::npos || value.find("e") != std::string::npos) {
//                std::istringstream iss(value);
//                double temp;
//                iss >> temp;
//                row.push_back(temp);
//            } else {
//                row.push_back(std::stod(value));
//            }
//        }
//        ts.push_back(row[0]); // Store the timestamp
//
//        while (std::getline(line_stream, value, ',')) {
//            row.push_back(std::stod(value));
//        }
//
//        gt_trajectory_lla.push_back({row[1], row[2], row[3]});
//        gt_yaws.push_back(deg2rad(row[10]));
//        obs_yaw_rates.push_back(deg2rad(row[13]));
//        obs_forward_velocities.push_back(row[4]);
//    }

    std::array<double, 3> origin = gt_trajectory_lla[0]; // Set the initial position to the origin
    std::vector<std::array<double, 3>> obs_trajectory_xyz = lla_to_enu(gt_trajectory_lla, origin);
//    std::cout<<gt_trajectory_lla[199][0]<<std::endl;
//    std::cout<<gt_trajectory_lla[200][0]<<std::endl;
    size_t N = ts.size(); // Number of data points
//    std::cout<<obs_trajectory_xyz[1000][0]<<std::endl;
//    std::cout<<obs_trajectory_xyz[1000][1]<<std::endl;
//    std::cout<<obs_trajectory_xyz[1000][2]<<std::endl;
    double xy_obs_noise_std = 5.0; // Standard deviation of observation noise of x and y in meters
    double yaw_rate_noise_std = 0.02; // Standard deviation of yaw rate in rad/s
    double forward_velocity_noise_std = 0.3; // Standard deviation of forward velocity in m/s

    // Prepare initial estimate and its error covariance
    double initial_yaw_std = M_PI;
    //double initial_yaw = gt_yaws[0] + sample_normal_distribution(0, initial_yaw_std);
    double initial_yaw = gt_yaws[0];

    Eigen::Vector3d x(obs_trajectory_xyz[0][0], obs_trajectory_xyz[0][1], initial_yaw);
    //std::cout<<x<<std::endl;
    Eigen::Matrix3d P;
    P << xy_obs_noise_std * xy_obs_noise_std, 0, 0,
            0, xy_obs_noise_std * xy_obs_noise_std, 0,
            0, 0, initial_yaw_std * initial_yaw_std;

    // Prepare measurement error covariance Q
    Eigen::Matrix2d Q;
    Q << xy_obs_noise_std * xy_obs_noise_std, 0,
            0, xy_obs_noise_std * xy_obs_noise_std;

    // Prepare state transition noise covariance R
    Eigen::Matrix3d R;
    R << forward_velocity_noise_std * forward_velocity_noise_std, 0, 0,
            0, forward_velocity_noise_std * forward_velocity_noise_std, 0,
            0, 0, yaw_rate_noise_std * yaw_rate_noise_std;

    // Initialize Kalman filter
    ExtendedKalmanFilter kf(x, P);

    std::vector<double> mu_x = {x[0]};
    std::vector<double> mu_y = {x[1]};
    std::vector<double> mu_theta = {x[2]};
    std::vector<double> var_x = {P(0, 0)};
    std::vector<double> var_y = {P(1, 1)};
    std::vector<double> var_theta = {P(2, 2)};
//
    double t_last = 0.0;
    for (size_t t_idx = 1; t_idx < N; ++t_idx) {
        double t = ts[t_idx];
        double dt = t - t_last;
        //std::cout<<"t" <<t<<std::endl;
        //std::cout<<"dt" <<dt<<std::endl;
        // Get control input `u = [v, omega] + noise`
        Eigen::Vector2d u(obs_forward_velocities[t_idx], obs_yaw_rates[t_idx]);
        // Because velocity and yaw rate are multiplied with `dt` in the state transition function,
        // its noise covariance must be multiplied with `dt**2.`
        Eigen::Matrix3d R_ = R * (dt * dt);
        // Propagate!
        kf.propagate(u, dt, R);
        //cout<<"predict :" <<kf.x_[0]<<endl;
        // Get measurement `z = [x, y] + noise`
        Eigen::Vector2d z(obs_trajectory_xyz[t_idx][0], obs_trajectory_xyz[t_idx][1]);
        //cout<<"z"<<z<<endl;
        // Update!
        kf.update(z, Q);
//        cout<<"update :" <<kf.x_[0]<<endl;
        // Save estimated state to analyze later
        mu_x.push_back(kf.x_[0]);
        mu_y.push_back(kf.x_[1]);
        mu_theta.push_back(normalize_angles(kf.x_[2]));

        // Save estimated variance to analyze later
        var_x.push_back(kf.P_(0, 0));
        var_y.push_back(kf.P_(1, 1));
        var_theta.push_back(kf.P_(2, 2));
        t_last = t;
    }
    // mu_x, mu_y, and mu_theta are the estimated 2D pose [x, y, theta]
    // var_x, var_y, and var_theta are the estimated error variances of 2D pose

    std::ofstream output_file("output.csv");
    output_file << "longitude,latitude,altitude,state_x,state_y,state_yaw" << std::endl;

    for (size_t i = 0; i < obs_trajectory_xyz.size(); ++i) {
        double lon = obs_trajectory_xyz[i][0];
        double lat = obs_trajectory_xyz[i][1];
        double alt = obs_trajectory_xyz[i][2];
        double state_x = mu_x[i];
        double state_y = mu_y[i];
        double state_yaw = mu_theta[i];
        output_file << lon << "," << lat << "," << alt << "," << state_x << "," << state_y << "," << state_yaw << std::endl;
//        if (i < mu_x.size()) {
//            double state_x = mu_x[i];
//            double state_y = mu_y[i];
//            double state_yaw = mu_theta[i];
//            output_file << lat << "," << lon << "," << alt << "," << state_x << "," << state_y << "," << state_yaw << std::endl;
//        } else {
//            output_file << lat << "," << lon << "," << alt << ",," << ",," << std::endl;
//        }
    }

    output_file.close();

    return 0;
}