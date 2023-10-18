#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>


void save_trajectory(Simulation simulation, Trajectory_Data trajectory_data)
{
	// Write arrays to file

    std::ofstream filex(simulation.pos_output_name_x);
    if (filex.is_open())
    {
        printf("Writing x-position trajectory to file: %s\n", simulation.pos_output_name_x.c_str());
        filex << trajectory_data.pos_traj_x << '\n';
    }
    filex.close();
    std::ofstream filey(simulation.pos_output_name_y);
    if (filey.is_open())
    {
        printf("Writing y-position trajectory to file: %s\n", simulation.pos_output_name_y.c_str());
        filey << trajectory_data.pos_traj_y << '\n';
    }
    filey.close();

    std::ofstream filevx(simulation.vel_output_name_x);
    if (filevx.is_open())
    {
        printf("Writing x-velocity trajectory to file: %s\n", simulation.vel_output_name_x.c_str());
        filevx << trajectory_data.vel_traj_x << '\n';
    }
    filevx.close();
    std::ofstream filevy(simulation.vel_output_name_y);
    if (filevy.is_open())
    {
        printf("Writing y-velocity trajectory to file: %s\n", simulation.vel_output_name_y.c_str());
        filevy << trajectory_data.vel_traj_y << '\n';
    }
    filevy.close();

    printf("Writing pressure to file: %s\n", simulation.pressure_output_name.c_str());

    std::ofstream filep(simulation.pressure_output_name);
    if (filep.is_open())
    {
        for (const auto& item : trajectory_data.pressure)
        {
            filep << item << std::endl;
        }
    }
    filep.close();
}