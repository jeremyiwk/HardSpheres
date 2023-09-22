#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>

Eigen::VectorXf gen_config(int &N_particles)
{

	int max_particles = 10;

	if (N_particles > max_particles)
	{
		printf("Max particle count exceeded.\n");
		printf("Particle count set to maximum: %d \n", max_particles);
		N_particles = max_particles;	
	}


	Eigen::VectorXf vec = Eigen::VectorXf::Random(6*N_particles);;
	
	return vec;
}
