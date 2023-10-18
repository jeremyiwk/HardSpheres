#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>

struct Initial_Config
{
	Eigen::VectorXd pos_x;
	Eigen::VectorXd vel_x;
	Eigen::VectorXd pos_y;
	Eigen::VectorXd vel_y;
};


Initial_Config gen_config(int &N_particles) //, int dimension)
{

	int max_particles = 10000;

	if (N_particles > max_particles)
	{
		printf("Max particle count exceeded.\n");
		printf("Particle count set to maximum: %d \n", max_particles);
		N_particles = max_particles;	
	}

	Initial_Config initial_config;

	initial_config.pos_x = Eigen::VectorXd::Random(N_particles);
	initial_config.vel_x = Eigen::VectorXd::Random(N_particles);
	initial_config.pos_y = Eigen::VectorXd::Random(N_particles);
	initial_config.vel_y = Eigen::VectorXd::Random(N_particles);
	
	return initial_config;
}
