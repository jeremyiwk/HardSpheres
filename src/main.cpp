#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "integrator.h"


struct Simulation {

  int dimension = 2;
  int time_steps = 250;
  int particles = 10;
  float particle_radii = 0.01;
  float time_step = 0.1;
  float sim_x_min = -1.0;
  float sim_y_min = -1.0;
  float sim_x_max = 1.0;
  float sim_y_max = 1.0;
};

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
  printf("sim_x_min      = %f\n",   simulation.sim_x_min );
  printf("sim_y_min      = %f\n",   simulation.sim_y_min );
  printf("sim_x_max      = %f\n",   simulation.sim_x_max );
  printf("sim_y_max      = %f\n",   simulation.sim_y_max );

  std::string sim_output_name_x = "sim_traj_x.dat";
  std::string sim_output_name_y = "sim_traj_y.dat";
	Eigen::VectorXf pos_x = Eigen::VectorXf::Random(simulation.particles);
  Eigen::VectorXf vel_x = Eigen::VectorXf::Random(simulation.particles);
  Eigen::VectorXf pos_y = Eigen::VectorXf::Random(simulation.particles);
  Eigen::VectorXf vel_y = Eigen::VectorXf::Random(simulation.particles);
  Eigen::MatrixXf pos_traj_x = Eigen::MatrixXf::Zero(simulation.time_steps, simulation.particles);
  Eigen::MatrixXf pos_traj_y = Eigen::MatrixXf::Zero(simulation.time_steps, simulation.particles);

  printf("Running hard sphere simulation for %d steps\n", simulation.time_steps);

  float sim_time = 0.0;

  for (int it=0; it<simulation.time_steps; ++it)
  {
    printf("Step: %d / %d \t Time: %f\n", it, simulation.time_steps, sim_time);
    sim_time += simulation.time_step;

    pos_x += simulation.time_step*vel_x;
    pos_y += simulation.time_step*vel_y;

    // Resolve Overlaps

    for (int i=0; i<simulation.particles; ++i)
    {
      for (int j=i+1; j<simulation.particles; ++j)
      {

        Eigen::VectorXf r_ij(2);

        r_ij(0) = pos_x(i) - pos_x(j);
        r_ij(1) = pos_y(i) - pos_y(j);

        double dist = sqrt(r_ij.dot(r_ij));

        r_ij /= dist;

        if (dist < 2*simulation.particle_radii)
        {
          pos_x(i) += (2*simulation.particle_radii - dist)*r_ij(0);
          pos_y(i) += (2*simulation.particle_radii - dist)*r_ij(1);
          pos_x(j) -= (2*simulation.particle_radii - dist)*r_ij(0);
          pos_y(j) -= (2*simulation.particle_radii - dist)*r_ij(1);

          vel_x(i) = sqrt(vel_x(j)*vel_x(i) + vel_y(j)*vel_y(i))*r_ij(0);
          vel_y(i) = sqrt(vel_x(j)*vel_x(i) + vel_y(j)*vel_y(i))*r_ij(1);
          vel_x(j) = -sqrt(vel_x(i)*vel_x(j) + vel_y(i)*vel_y(j))*r_ij(0);
          vel_y(j) = -sqrt(vel_x(i)*vel_x(j) + vel_y(i)*vel_y(j))*r_ij(1);
        }
      }

      // Resolve Boundary Conditions

      if (pos_x(i) - simulation.particle_radii < simulation.sim_x_min)
      {
        pos_x(i) = simulation.sim_x_min + simulation.particle_radii;
        vel_x(i) *= -1;
      }
      else if (pos_y(i) - simulation.particle_radii < simulation.sim_y_min)
      {
        pos_y(i) = simulation.sim_y_min + simulation.particle_radii;
        vel_y(i) *= -1;
      }
      else if (pos_x(i) + simulation.particle_radii > simulation.sim_x_max)
      {
        pos_x(i) = simulation.sim_x_max - simulation.particle_radii;
        vel_x(i) *= -1;
      }
      else if (pos_y(i) + simulation.particle_radii > simulation.sim_y_max)
      {
        pos_y(i) = simulation.sim_y_max - simulation.particle_radii;
        vel_y(i) *= -1;
      }

    }

    // Copy coordinates to trajectory vector

    pos_traj_x.row(it) = pos_x;
    pos_traj_y.row(it) = pos_y;
  }


  // Write arrays to file

  std::ofstream filex(sim_output_name_x);
  if (filex.is_open())
  {
    printf("Writing simulation result to file\n");
    filex << pos_traj_x << '\n';
  }
  filex.close();
  std::ofstream filey(sim_output_name_y);
  if (filey.is_open())
  {
    printf("Writing simulation result to file\n");
    filey << pos_traj_y << '\n';
  }
  filey.close();

	return 0;
}
