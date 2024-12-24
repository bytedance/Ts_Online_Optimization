# Universal Online Temporal Calibration for Optimization-based Visual-Inertial Navigation Systems
 6-Degree of Freedom (6DoF) motion estimation with a combination of visual and inertial sensors is a growing area with numerous real-world applications. However, precise calibration of the time offset between these two sensor types is a prerequisite for accurate and robust tracking. To address this, we propose a universal online temporal calibration strategy for optimization-based visual-inertial navigation systems. Technically, we incorporate the time offset t d as a state parameter in the optimization residual model to align the IMU state to the corresponding image timestamp using td, angular velocity and translational velocity. This allows the temporal misalignment t d to be optimized alongside other tracking states during the process. As our method only modifies the structure of the residual model, it can be applied to various optimization-based frameworks with different tracking frontends. We evaluate our calibration method with both EuRoC and simulation data and extensive experiments demonstrate that our approach provides more accurate time offset estimation and faster convergence, particularly in the presence of noisy sensor data.

The main contributions of this study are concluded as follows:
- A universal online temporal calibration strategy for optimization-based VINS algorithms is proposed.
- The proposed approach is integrated into the existing VINS framework to validate its feasibility.
- Comprehensive experiments are conducted to assess its impact on accuracy.

# 1. License

The code is licensed under GPLv3.

The Experimental code is developed on VINS-Mono(**[VINS-Mono](https://github.com/HKUST-Aerial-Robotics/VINS-Mono)**), and thus its license is retained at the beginning of the related files.

**Related Publication:**  
The related PrePrint is released soon.

# 2. Prerequisites
1.1 **Ubuntu** and **ROS**
Ubuntu 16.04. or Ubuntu 18.04.

ROS Kinetic. [Installation](https://wiki.ros.org/kinetic/Installation/Ubuntu)
or 
ROS Melodic. [Installation](https://wiki.ros.org/melodic/Installation/Ubuntu)


1.2. **Ceres Solver**
Follow [Ceres Installation](http://ceres-solver.org/installation.html).

Our testing environment: Ubuntu 18.04, ROS Melodic, OpenCV 3.4.1, Eigen 3.3.4, Ceres-Solver 1.14.0rc1

# 2. Building on ROS
Clone the repository and catkin_make:
```
    cd ~/catkin_ws/src
    git clone project
    cd ../
    catkin_make
    source ~/catkin_ws/devel/setup.bash
```

# 3. Run
```
roslaunch vins_estimator euroc.launch
roslaunch vins_estimator vins_rviz.launch
rosbag play EUROC_DATASET_DIR/MH_01_easy.bag
```

# 4. Comments

Some parameters in euroc_config.yaml need to be set at the beginning.
If the two parameters with "estimate_td2: 1" and "estimate: 0" are set, the proposed method in this paper is enabled. Otherwise, the default approach of timestamp offset calibration integrated on VINS-Mono will be enabled when two parameters are set with "estimate_td2: 0" and "estimate: 1".

The parameter td_perturbation in euroc_config.yaml is the channel to inject the timestamp misalignment manually when the dataset is absolutely without temporal misalignment.

# 5. Security

If you discover a potential security issue in this project, or think you may
have discovered a security issue, we ask that you notify Bytedance Security via our [security center](https://security.bytedance.com/src) or [vulnerability reporting email](sec@bytedance.com).

# 6. Acknowledgement
The Experimental code is developed on [VINS-Mono](https://github.com/HKUST-Aerial-Robotics/VINS-Mono). We extend our gratitude to the authors of the software.


# 7. Citation
TBD

# 8. We are Hiring!
Our team is hiring FTEs with background in Deep Learning, SLAM, and 3D Vision. We are based in Beijing and Shanghai. If you are interested, please send your resume to frank.01[AT]bytedance[DOT]com.