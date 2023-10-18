#include <iostream>
#include <stdio.h>
#include <fstream>
#include <map>
#include <string>
#include <any>


struct Simulation {

  int dimension = 2;
  int time_steps = 50;
  int particles = 50;
  double particle_radii = 0.005;
  double time_step = 0.01;
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