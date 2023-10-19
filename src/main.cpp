#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "input_params.h"
#include "gen_config.h"
#include "integrator.h"
#include "file_writer.h"



int main()
{

    printf("Initiating Hard Sphere simulation.\n");

    // Simulation simulation; // defined in input_params.h
    // display_simulation_params(simulation); // defined in input_params.h
    // Initial_Config initial_config = gen_config(simulation.particles); // defined in gen_config.h
    // Trajectory_Data trajectory_data = integrate(simulation, initial_config); // defined in integrator.h
    // save_trajectory(simulation, trajectory_data); // defined in file_writer.h



    Eigen::VectorXd particle_radii = Eigen::VectorXd::LinSpaced(100, 0.001, 0.07);

    for (int i; i < particle_radii.rows(); ++i)
    {
    	// std::cout << "particle radii" << g << std::endl;

    	double radius = particle_radii(i);

    	Simulation simulation; // defined in input_params.h

    	simulation.particle_radii = radius;

    	simulation.pos_output_name_x = "position_traj_x_" + std::to_string(i) + ".dat";
  		simulation.pos_output_name_y = "position_traj_y_" + std::to_string(i) + ".dat";
  		simulation.vel_output_name_x = "velocity_traj_x_" + std::to_string(i) + ".dat";
  		simulation.vel_output_name_y = "velocity_traj_y_" + std::to_string(i) + ".dat";
  		simulation.pressure_output_name = "pressure_" + std::to_string(i) + ".dat";

	    display_simulation_params(simulation); // defined in input_params.h
	    Initial_Config initial_config = gen_config(simulation.particles); // defined in gen_config.h
	    Trajectory_Data trajectory_data = integrate(simulation, initial_config); // defined in integrator.h
	    save_trajectory(simulation, trajectory_data); // defined in file_writer.h

    }

    std::ofstream filer("radii.dat");
    if (filer.is_open())
    {
        filer << particle_radii << '\n';
    }
    filer.close();



    return 0;
}
