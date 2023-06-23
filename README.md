# GPS, Velocity and IMU fusion with Extended Kalman Filter
Use main.cpp for ENU coordinate system 
Use main_utm for UTM coordinate system
# Algorithm
![alt text](https://github.com/Janudis/EKF_GPS_IMU/blob/master/Extended-Kalman-Filter-Step.png)

This code implements an Extended Kalman Filter (EKF) for fusing Global Positioning System (GPS) and Inertial Measurement Unit (IMU) measurements. The goal is to estimate the state (position and orientation) of a vehicle using both GPS and IMU data.

In our case, IMU provide data more frequently than GPS. Here is a step-by-step description of the process:
1) Initialization: Firstly, initialize your EKF with initial estimates for the position, velocity, and orientation, typically using your first GPS reading for position and velocity. Your orientation could come from IMU sensor data. The covariance matrix (P) should also be initialized, usually with relatively high values on the diagonals to reflect initial uncertainty.
 
2) Prediction step (also known as Time Update): In the prediction step, you make a prediction of the current state using the process model and the previously estimated state. In our case, the state can be represented as [position, velocity, orientation] or [x, y, z, dx/dt, dy/dt, dz/dt, roll, pitch, yaw] in 3D. Using the IMU readings, we can integrate acceleration to estimate change in velocity and integrate again to estimate change in position. Similarly, we can integrate angular velocity readings to estimate changes in orientation. We also need to predict the state covariance matrix (P) at this point.
 
3) Update step (also known as Measurement Update): In the update step, we correct the predicted state using the measurement data. Here, we would use the GPS readings to correct the position and velocity predictions from the IMU. We compute the Kalman gain, which balances the weights given to the prediction and the new measurement based on their respective uncertainties. If GPS readings are available (say, every 25th step), we use them for this correction.
 
4) Estimate update: With the Kalman gain, we can compute the updated (a posteriori) state estimate by combining the predicted state estimate and the weighted difference between the actual measurement and the measurement predicted by the a priori estimate (also known as the measurement residual).
   
5) Covariance update: We also compute the updated (a posteriori) estimate covariance (P). The a posteriori state and covariance estimates at the current time become the a priori estimates for the next time step.

6) Repeat steps 2-5: This process is then repeated for each time step, using the a posteriori estimates from the previous time step as the a priori estimates for the current step.
   
7) Handling of GPS updates: When the GPS update is available, it is incorporated in the update step, otherwise the process continues with the prediction step using only the IMU data.

# Dependencies
C++ compiler supporting C++11 or higher
Eigen library (for linear algebra operations)

# Usage
Change the CMakeLists.txt

# Code Structure
ekf.h: Header file containing the declaration of the ExtendedKalmanFilter class, which implements the EKF algorithm.

geo_ned.h: Header file containing helper functions for converting between geodetic (WGS84) and East-North-Down (ENU) coordinate systems.

utm.h: Header file containing the definition of the utm_coords struct and the utmconv namespace, which provides functions for converting between geodetic and UTM coordinates.

main_utm.cpp: The main C++ file that reads input data from a CSV file, performs the GPS and IMU fusion using the EKF, and outputs the estimated position and orientation. (Use main.cpp for ENU coordinate system).

# Input Data Format
The input data is expected to be in a CSV file (localization_log2.csv) with the following columns:

Timestamp (in nanoseconds)
Latitude (in degrees, WGS84)
Longitude (in degrees, WGS84)
Altitude (in meters, WGS84)
Forward velocity (in meters per second)
Yaw rate (in radians per second)

# Output
The output of the code is a CSV file (output_utm.csv) containing the estimated position and orientation in UTM coordinates. The file has the following columns:

Easting (UTM coordinate, in meters)
Northing (UTM coordinate, in meters)
Yaw (orientation angle, in radians)
Estimated X position (in meters)
Estimated Y position (in meters)
Estimated yaw (in radians)

# Adjusting Parameters
The code provides options for adjusting the standard deviation of the observation noise for x and y coordinates (xy_obs_noise_std), the yaw rate (yaw_rate_noise_std), and the forward velocity (forward_velocity_noise_std). These parameters can be modified in the code to suit your specific scenario and sensor characteristics.


