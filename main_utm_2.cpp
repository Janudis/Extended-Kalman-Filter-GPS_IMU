#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include "ekf.h"
#include "geo_ned.h"
//#include "utm.cpp"
#include "utm.h"
using namespace std;
using namespace Eigen;
#include <iomanip>

int main() {
    // Create an instance of utmconv
    utmconv::utm_coords coords{};

    std::vector<std::array<double, 3>> obs_trajectory_xyz;
    std::vector<std::array<double, 3>> gt_trajectory_xyz;
    std::vector<double> gt_yaws; // [yaw_angle(rad),] x N
    std::vector<double> obs_yaw_rates; // [vehicle_yaw_rate(rad/s),] x N
    std::vector<double> obs_forward_velocities; // [vehicle_forward_velocity(m/s),] x N
    std::vector<double> ts;
    double last_non_empty_speed = 0.0;
    double last_non_empty_yaw = 0.0;
    double last_non_empty_yaw_rate = 0.0;
    double last_non_empty_x = 0.0;
    double last_non_empty_y = 0.0;

    std::ifstream raw_data("localization_log2.csv");
    std::string line, value;

// Read CSV file and store values in vectors
    std::getline(raw_data, line); // Skip header line

    long long int first_timestamp = 0;
    bool first_line = true;

    while (std::getline(raw_data, line)) {
        std::stringstream line_stream(line);
        std::vector<std::string> string_row;
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
            string_row.push_back(value);
            if (!value.empty()) {
                row.push_back(std::stold(value));
            }
            else {
                row.push_back(std::numeric_limits<double>::quiet_NaN()); // Push NaN to the row for empty value
            }
        }
        if (string_row[0].empty()) {  // DEN EXOUME GPS
            utmconv::geodetic_to_utm(last_non_empty_x, last_non_empty_y, coords);
            gt_trajectory_xyz.push_back({coords.easting,coords.northing});
            obs_trajectory_xyz.push_back({row[0], row[1]}); //push nan nan
            obs_forward_velocities.push_back(last_non_empty_speed);
            last_non_empty_yaw_rate = row[12];
            last_non_empty_yaw = row[9];
            gt_yaws.push_back(deg2rad(last_non_empty_yaw));
            obs_yaw_rates.push_back(last_non_empty_yaw_rate);
        } else { //exoume GPS
            last_non_empty_x = row[0];
            last_non_empty_y = row[1];
            utmconv::geodetic_to_utm(row[0], row[1], coords);
            obs_trajectory_xyz.push_back({coords.easting,coords.northing});
            gt_trajectory_xyz.push_back({coords.easting,coords.northing});
            last_non_empty_speed = row[3];
            obs_forward_velocities.push_back(last_non_empty_speed);
            gt_yaws.push_back(deg2rad(last_non_empty_yaw));
            obs_yaw_rates.push_back(last_non_empty_yaw_rate);
        }
    }

    size_t N = ts.size(); // Number of data points
    double xy_obs_noise_std = 5.0; // Standard deviation of observation noise of x and y in meters
    double yaw_rate_noise_std = 0.02; // Standard deviation of yaw rate in rad/s
    double forward_velocity_noise_std = 0.3; // Standard deviation of forward velocity in m/s

    // Prepare initial estimate and its error covariance
    double initial_yaw_std = M_PI;
    double initial_yaw = gt_yaws[0] + sample_normal_distribution(0, initial_yaw_std);

    Eigen::Vector3d x(obs_trajectory_xyz[0][0], obs_trajectory_xyz[0][1], initial_yaw);
    // Initialize Kalman filter
    ExtendedKalmanFilter kf;
    kf.initialize(x, xy_obs_noise_std, yaw_rate_noise_std, forward_velocity_noise_std, initial_yaw_std);

    std::vector<double> mu_x = {x[0]};
    std::vector<double> mu_y = {x[1]};
    std::vector<double> mu_theta = {x[2]};
    std::vector<double> var_x = {kf.P(0, 0)};
    std::vector<double> var_y = {kf.P(1, 1)};
    std::vector<double> var_theta = {kf.P(2, 2)};

//
    double t_last = 0.0;
    for (size_t t_idx = 1; t_idx < N; ++t_idx) {
        double t = ts[t_idx];
        double dt = t - t_last;
        // Get control input `u = [v, omega] + noise`
        Eigen::Vector2d u(obs_forward_velocities[t_idx], obs_yaw_rates[t_idx]);
        // Because velocity and yaw rate are multiplied with `dt` in the state transition function,
        // its noise covariance must be multiplied with `dt**2.`
        //Eigen::Matrix3d R_ = R * (dt * dt);
        // Propagate!
        kf.propagate(u, dt);
        // Get measurement `z = [x, y] + noise`
        if (!std::isnan(obs_trajectory_xyz[t_idx][0])) {
            Eigen::Vector2d z(obs_trajectory_xyz[t_idx][0], obs_trajectory_xyz[t_idx][1]);
            // Update!
            kf.update(z);
        }
//        Eigen::Vector2d z(obs_trajectory_xyz[t_idx][0], obs_trajectory_xyz[t_idx][1]);
//        // Update!
//        kf.update(z);
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

    std::ofstream output_file("output_utm.csv");
    output_file << std::setprecision(12) << "easting,northing,yaw,state_x,state_y,state_yaw" << std::endl;
    for (size_t i = 0; i < obs_trajectory_xyz.size(); ++i) {
        double lon = gt_trajectory_xyz[i][0];
        double lat = gt_trajectory_xyz[i][1];
        double yaw = gt_yaws[i];
        double state_x = mu_x[i];
        double state_y = mu_y[i];
        double state_yaw = mu_theta[i];
        output_file << lon << "," << lat << "," << yaw << "," << state_x << "," << state_y << "," << state_yaw << std::endl;
    }
    output_file.close();
    return 0;
}