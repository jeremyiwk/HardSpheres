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


  printf("Initiating Hard Sphere simulation.\n");
  printf("Simulation Parameters:\n");

  printf("dimension      = %d\n",   simulation.dimension );
  printf("particles      = %d\n",   simulation.particles );
  printf("time_steps     = %d\n",   simulation.time_steps );
  printf("particle radii = %f\n",   simulation.particle_radii);
  printf("time_step      = %f\n",   simulation.time_step );
  printf("sim_box_x      = %f\n",   simulation.sim_box_x );
  printf("sim_box_y      = %f\n",   simulation.sim_box_y );
  printf("sim_box_z      = %f\n",   simulation.sim_box_z );


	Eigen::VectorXf position = gen_config(simulation.particles);
  Eigen::VectorXf velocity = gen_config(simulation.particles);
  Eigen::MatrixXf position_traj = Eigen::MatrixXf::Zero(simulation.time_steps,3*simulation.particles);

  std::cout << position_traj << std::endl;

  std::cout << "Position trajectory\n" << position_traj(0,Eigen::all) << std::endl;

  position_traj(0,Eigen::all) = position;

  std::cout << position_traj << std::endl;

  std::cout << "Position trajectory\n" << position_traj(0,Eigen::all) << std::endl;

  // Eigen::VectorXf position2 = integrate(position, velocity, simulation.time_step, simulation.particle_radii);

  // std::cout << simulation << std::endl;

	return 0;
}
