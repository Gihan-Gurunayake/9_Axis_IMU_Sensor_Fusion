# 9 Axis IMU Sensor Fusion

THis repository contains necessary files to interface IMU sensors in a mobile phones ( Accelerometer, gyroscope and Magnetometer) and fuse the sensor data in real time using kalman filters to obtain attitude and heading thereby orientation as the end.

# Note
All the scripts are written to accept sensor data in the following order.
- Accelerometer (ax, ay, az)
- Magnetometer (mx, my, mz)
- Gyroscope (gx, gy, gz)

# Instructions
## Lab 0
Install Hyper IMU Android Application and follow the guide to configure.

## Lab 1
Enable Accelerometer and Gyroscope in the sensor list in the app to do the activity

## Lab 2
When running this script, run each section individually, do not run the whole script at once.
- Section 1: Configuration Settings ( you may change the calibration time)
- Section 2: Magnetometer calibration (Rotate the mobile phone in all orientations slowly during the calibration period)
- Section 3: Real time vizualization
## Lab 3
When running this script, run each section individually, do not run the whole script at once.
- Section 1: Configuration Settings ( you may change the calibration time)
- Section 2: Magnetometer calibration (Rotate the mobile phone in all orientations slowly during the calibration period)
- Section 3: Still intialization ( Keep the phone still on a table and run this section)
- Section 4: Real time vizualization
