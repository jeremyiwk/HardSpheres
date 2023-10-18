#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>


struct Trajectory_Data
{
	Eigen::MatrixXd pos_traj_x;
	Eigen::MatrixXd pos_traj_y;
	Eigen::MatrixXd vel_traj_x;
	Eigen::MatrixXd vel_traj_y;
	std::vector<double> pressure;
};

Trajectory_Data integrate(Initial_Config initial_config, Simulation simulation)
{
	return 0;
};


