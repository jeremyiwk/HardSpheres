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

    Simulation simulation; // defined in input_params.h

    display_simulation_params(simulation); // defined in input_params.h

    Initial_Config initial_config = gen_config(simulation.particles); // defined in gen_config.h

    Trajectory_Data trajectory_data = integrate(simulation, initial_config); // defined in integrator.h

    save_trajectory(simulation, trajectory_data); // defined in file_writer.h


    return 0;
}
