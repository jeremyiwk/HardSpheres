#include <iostream>
#include <stdio.h>
#include <fstream>
#include <map>
#include <string>
#include <any>


struct Simulation {

  int dimension = 2;
  int time_steps = 10000;
  int particles = 100;
  double particle_radii = 0.005;
  double time_step = 0.001;
  double sim_x_min = -1.0;
  double sim_y_min = -1.0;
  double sim_x_max = 1.0;
  double sim_y_max = 1.0;
  std::string pos_output_name_x = "position_traj_x.dat";
  std::string pos_output_name_y = "position_traj_y.dat";
  std::string vel_output_name_x = "velocity_traj_x.dat";
  std::string vel_output_name_y = "velocity_traj_y.dat";
  std::string pressure_output_name = "pressure.dat";
};


void display_simulation_params(Simulation simulation)
{

  printf("Simulation Parameters:\n");
  printf("dimension = %d\n", simulation.dimension);
  printf("particles = %d\n", simulation.particles);
  printf("time_steps = %d\n", simulation.time_steps);
  printf("particle radii = %f\n", simulation.particle_radii);
  printf("time_step = %f\n", simulation.time_step);
  printf("sim_x_min = %f\n", simulation.sim_x_min);
  printf("sim_y_min = %f\n", simulation.sim_y_min);
  printf("sim_x_max = %f\n", simulation.sim_x_max);
  printf("sim_y_max = %f\n", simulation.sim_y_max);

  std::cout << "pos_output_name_x = " << simulation.pos_output_name_x << std::endl;
  std::cout << "pos_output_name_y = " << simulation.pos_output_name_y << std::endl;
  std::cout << "vel_output_name_x = " << simulation.pos_output_name_x << std::endl;
  std::cout << "vel_output_name_y = " << simulation.pos_output_name_y << std::endl;
  std::cout << "pressure_output_name = " << simulation.pressure_output_name << std::endl;
};


/*

Below is some code to hopefully someday make this simulation run on input files
that are read into a map object or something similar.
*/

// std::map<std::string, std::string> read_input(std::string filename) 
// {
//     std::ifstream inputFile(filename);
//     // if (!inputFile.is_open()) {
//     //     std::cerr << "Could not open the input file." << std::endl;
//     //     return 1;
//     // }

//     std::map<std::string, std::string> inputMap;
//     std::string key, value;

//     while (inputFile >> key >> value) {
//         inputMap[key] = value;
//     }

//     inputFile.close();

//     // Print the contents of the map
//     for (const auto& entry : inputMap) {
//         std::cout << "Key: " << entry.first << ", Value: " << entry.second << std::endl;
//     }

//     return inputMap;
// }