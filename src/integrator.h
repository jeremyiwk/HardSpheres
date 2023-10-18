#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>


struct Trajectory_Data
{
	Eigen::MatrixXd pos_traj_x;
	Eigen::MatrixXd pos_traj_y;
	Eigen::MatrixXd vel_traj_x;
	Eigen::MatrixXd vel_traj_y;
	std::vector<double> pressure;
};

Trajectory_Data integrate(Simulation simulation, Initial_Config initial_config)
{
	Trajectory_Data trajectory_data;

	trajectory_data.pos_traj_x = Eigen::MatrixXd::Zero(simulation.time_steps, simulation.particles);
	trajectory_data.pos_traj_y = Eigen::MatrixXd::Zero(simulation.time_steps, simulation.particles);
	trajectory_data.vel_traj_x = Eigen::MatrixXd::Zero(simulation.time_steps, simulation.particles);
	trajectory_data.vel_traj_y = Eigen::MatrixXd::Zero(simulation.time_steps, simulation.particles);

	// printf("Running hard sphere simulation for %d steps\n", simulation.time_steps);

    float sim_time = 0.0;

    for (int it=0; it<simulation.time_steps; ++it)
    {
        // printf("Step: %d / %d \t Time: %f\n", it, simulation.time_steps, sim_time);

        sim_time += simulation.time_step;

        initial_config.pos_x += simulation.time_step*initial_config.vel_x;
        initial_config.pos_y += simulation.time_step*initial_config.vel_y;

        // Resolve Overlaps

        for (int i=0; i<simulation.particles; ++i)
        {
            for (int j=i+1; j<simulation.particles; ++j)
            {

                Eigen::VectorXd r_ij(2);
                Eigen::VectorXd v_ij(2);

                r_ij(0) = initial_config.pos_x(i) - initial_config.pos_x(j);
                r_ij(1) = initial_config.pos_y(i) - initial_config.pos_y(j);

                double dist = sqrt(r_ij.dot(r_ij));

                r_ij /= dist;

                if (dist < 2*simulation.particle_radii)
                {
                    initial_config.pos_x(i) += (2*simulation.particle_radii - dist)*r_ij(0)/2;
                    initial_config.pos_y(i) += (2*simulation.particle_radii - dist)*r_ij(1)/2;
                    initial_config.pos_x(j) -= (2*simulation.particle_radii - dist)*r_ij(0)/2;
                    initial_config.pos_y(j) -= (2*simulation.particle_radii - dist)*r_ij(1)/2;

                    r_ij(0) = initial_config.pos_x(i) - initial_config.pos_x(j);
                    r_ij(1) = initial_config.pos_y(i) - initial_config.pos_y(j);

                    v_ij(0) = initial_config.vel_x(i) - initial_config.vel_x(j);
                    v_ij(1) = initial_config.vel_y(i) - initial_config.vel_y(j);

                    double b = r_ij.dot(v_ij);

                    dist = sqrt(r_ij.dot(r_ij));

                    r_ij /= dist;

                    initial_config.vel_x(i) -= b*r_ij(0)/(2*simulation.particle_radii);
                    initial_config.vel_y(i) -= b*r_ij(1)/(2*simulation.particle_radii);
                    initial_config.vel_x(j) += b*r_ij(0)/(2*simulation.particle_radii);
                    initial_config.vel_y(j) += b*r_ij(1)/(2*simulation.particle_radii);
                }
            }

            // Resolve Boundary Conditions

            if (initial_config.pos_x(i) - simulation.particle_radii < simulation.sim_x_min)
            {
                trajectory_data.pressure.push_back(abs(initial_config.vel_x(i)));
                initial_config.pos_x(i) = simulation.sim_x_min + simulation.particle_radii;
                initial_config.vel_x(i) *= -1;
            }
            else if (initial_config.pos_y(i) - simulation.particle_radii < simulation.sim_y_min)
            {
                trajectory_data.pressure.push_back(abs(initial_config.vel_y(i)));
                initial_config.pos_y(i) = simulation.sim_y_min + simulation.particle_radii;
                initial_config.vel_y(i) *= -1;
            }
            else if (initial_config.pos_x(i) + simulation.particle_radii > simulation.sim_x_max)
            {
                trajectory_data.pressure.push_back(abs(initial_config.vel_x(i)));
                initial_config.pos_x(i) = simulation.sim_x_max - simulation.particle_radii;
                initial_config.vel_x(i) *= -1;
            }
            else if (initial_config.pos_y(i) + simulation.particle_radii > simulation.sim_y_max)
            {
                trajectory_data.pressure.push_back(abs(initial_config.vel_y(i)));
                initial_config.pos_y(i) = simulation.sim_y_max - simulation.particle_radii;
                initial_config.vel_y(i) *= -1;
            }

        }

        // Copy coordinates to trajectory vector

        trajectory_data.pos_traj_x.row(it) = initial_config.pos_x;
        trajectory_data.pos_traj_y.row(it) = initial_config.pos_y;
        trajectory_data.vel_traj_x.row(it) = initial_config.vel_x;
        trajectory_data.vel_traj_y.row(it) = initial_config.vel_y;

    }

	return trajectory_data;
};

