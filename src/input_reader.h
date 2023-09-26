#include <iostream>
#include <stdio.h>
#include <map>
#include <string>
#include <any>

struct Simulation {

  int dimension = 3;
  int time_steps = 10;
  int particles = 1;
  float particle_radii = 0.001;
  float time_step = 0.01;
  float sim_box_x = 1.0;
  float sim_box_y = 1.0;
  float sim_box_z = 1.0;
};

// std::map<std::string, std::any> sim_params;


// sim_params["dimension"]  = 3;
// sim_params["particles"]  = 11;
// sim_params["sim_box_x"]  = 1.0;
// sim_params["sim_box_y"]  = 1.0;
// sim_params["sim_box_z"]  = 1.0;
// sim_params["total_time"] = 1.0;
// sim_params["time_step"]  = 0.01;
