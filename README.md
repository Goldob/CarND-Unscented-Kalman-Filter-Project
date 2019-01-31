# Unscented Kalman Filter Project

The following is a solution for [UKF Project](https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project) from Term 2 of Udacity's Self Driving Car Nanodegree. Refer to the main repository for compilation and testing instructions.

The following tables show results (RMSE) achieved on both datasets provided in [the simulator](https://github.com/udacity/self-driving-car-sim/releases/tag/v1.45).

 #### Dataset 1

|        | Sensor Fusion | Radar Only | Lidar Only    |
| -------| ------------- | ---------- | ------------- |
| **X**  | 0.0817        | 0.1466     | 0.2115        |
| **Y**  | 0.0978        | 0.2004     | 0.1058        |
| **VX** | 0.2878        | 0.3057     | 0.6815        |
| **VY** | 0.1804        | 0.2219     | 0.2459        |

 #### Dataset 2

|        | Sensor Fusion | Radar Only | Lidar Only    |
| -------| ------------- | ---------- | ------------- |
| **X**  | 0.2581        | 0.2866     | 0.1932        |
| **Y**  | 0.0940        | 0.2025     | 0.0906        |
| **VX** | 0.7963        | 0.8267     | 0.7331        |
| **VY** | 0.2454        | 0.2839     | 0.2664        |

The performance of the algorithm was tested in three operation modes:

* process both measurement types (radar & lidar)
* process radar measurements only, discard lidar
* process lidar measurements only, discard radar
