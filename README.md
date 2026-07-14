# 9 Axis IMU Sensor Fusion

THis repository contains necessary files to interface IMU sensors in a mobile phones ( Accelerometer, gyroscope and Magnetometer) and fuse the sensor data in real time using kalman filters to obtain attitude and heading thereby orientation at the end.

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

Apple users may use Phyphox App
https://apps.apple.com/lk/app/phyphox/id1127319693
1. Install phyphox
Download phyphox on your iPhone. Connect the iPhone and laptop to the same WinFi network or
personal hotspot.
2. Select the sensor
Open phyphox and choose Accelerometer, Gyroscope, Magnetometer, Attitude, or another
required sensor.
3. Start data collection
Press the Play (n) button to begin recording live sensor data.
4. Enable Remote Access
Open the menu, enable Remote Access, and note the displayed local web address (for example
http://192.168.x.x:8080).
5. Open on your laptop
Enter the displayed address in a web browser on your laptop to view live graphs and controls.
6. Get the sensor data
View live values, export data as CSV/Excel, or access the data through the phyphox REST API.
7. Stop and save
Press Stop when finished and export or save the recorded data
