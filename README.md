# GPS, Velocity and IMU fusion with Extended Kalman Filter
Use main.cpp for ENU coordinate system and change the CMakeLists.txt
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


