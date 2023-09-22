#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include "input_reader.h"
#include "integrator.h"
#include "gen_config.h"


int main()
{

  Simulation simulation;
	integrate();

	Eigen::VectorXf vec = gen_config(simulation.particles);
//	std::cout << vec << std::endl;
	return 0;
}
