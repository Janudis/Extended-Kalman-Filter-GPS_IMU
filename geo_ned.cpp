#include "geo_ned.h"
using namespace Eigen;

double a = 6378137.0;
double f = 1.0 / 298.257223563;
double b = (1.0 - f) * a;
double e = sqrt(a * a - b * b) / a;

double deg2rad(double deg){
    return deg/180*M_PI;
}

double sample_normal_distribution(double mean, double std_dev) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(mean, std_dev);
    return distribution(gen);
}

double normalize_angles(double angle) {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    //if (angle < 0) angle += 2.0 * M_PI;
    return angle - M_PI;
}

Eigen::Matrix3d Rx(double theta){
    Eigen::Matrix3d mat;
    mat << 1, 0, 0,
            0, cos(theta), sin(theta),
            0, -sin(theta), cos(theta);
    return mat;
}

Eigen::Matrix3d Ry(double theta){
    Eigen::Matrix3d mat;
    mat << cos(theta), 0, -sin(theta),
            0, 1, 0,
            sin(theta), 0, cos(theta);
    return mat;
}

Eigen::Matrix3d Rz(double theta){
    Eigen::Matrix3d mat;
    mat << cos(theta), sin(theta), 0,
            -sin(theta), cos(theta), 0,
            0, 0, 1;
    return mat;
}
std::vector<std::array<double, 3>> lla_to_ecef(const std::vector<std::array<double, 3>>& points_lla) {
    std::vector<std::array<double, 3>> points_ecef(points_lla.size());

    for (int i = 0; i < points_lla.size(); i++) {
        double lon = points_lla[i][0] * M_PI / 180.0;
        double lat = points_lla[i][1] * M_PI / 180.0;
        double alt = points_lla[i][2];
        double N = a / sqrt(1.0 - (e * sin(lat)) * (e * sin(lat)));
        double x = (N + alt) * cos(lat) * cos(lon);
        double y = (N + alt) * cos(lat) * sin(lon);
        double z = (N * (1.0 - e * e) + alt) * sin(lat);
        points_ecef[i] = { x, y, z };
    }
    return points_ecef;
}

std::vector<std::array<double, 3>> ecef_to_enu(const std::vector<std::array<double, 3>>& points_ecef, const std::array<double, 3>& ref_lla) {
    std::vector<std::array<double, 3>> ref_ecef_vec = lla_to_ecef({ ref_lla });
    std::array<double, 3> ref_ecef = ref_ecef_vec[0];

    MatrixXd dX(points_ecef.size(), 3);
    for (int i = 0; i < points_ecef.size(); i++) {
        dX.row(i) = Vector3d(points_ecef[i][0], points_ecef[i][1], points_ecef[i][2]).transpose() - Vector3d(ref_ecef[0], ref_ecef[1], ref_ecef[2]).transpose();
    }

    double ref_lon = ref_lla[0] * M_PI / 180.0;
    double ref_lat = ref_lla[1] * M_PI / 180.0;

    Matrix3d Renu;
    Renu << -sin(ref_lon), cos(ref_lon), 0,
            -sin(ref_lat) * cos(ref_lon), -sin(ref_lat) * sin(ref_lon), cos(ref_lat),
            cos(ref_lat) * cos(ref_lon), cos(ref_lat) * sin(ref_lon), sin(ref_lat);

    MatrixXd dXenu = Renu * dX.transpose();

    std::vector<std::array<double, 3>> points_enu(points_ecef.size());
    for (int i = 0; i < points_ecef.size(); i++) {
        points_enu[i] = { dXenu(0, i), dXenu(1, i), dXenu(2, i) };
    }
    return points_enu;
}

std::vector<std::array<double, 3>> lla_to_enu(const std::vector<std::array<double, 3>>& points_lla, const std::array<double, 3>& ref_lla) {
    std::vector<std::array<double, 3>> points_ecef = lla_to_ecef(points_lla);
    return ecef_to_enu(points_ecef, ref_lla);
}