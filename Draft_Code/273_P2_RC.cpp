

#include <iostream>
#include <cmath>
#include <vector>
#include <random>


std::vector<std::vector<double>> CoordsToPotential(const std::vector<double>& Xinit, const std::vector<double>& Yinit, double a, double b, int Particles) {

    std::vector<std::vector<double>> Phi(Particles, std::vector<double>(Particles, 0.0));
    
    for (int i = 0; i < Particles; ++i) {
        for (int j = i + 1; j < Particles; ++j) {
            if (i == j) continue;

            double Dx = Xinit[i] - Xinit[j];
            double Dy = Yinit[i] - Yinit[j];
            double r = std::sqrt(Dx*Dx + Dy*Dy);

            if (r > 0) {
                Phi[i][j] = (a / std::pow(r, 12)) - (b / std::pow(r, 6));
            } else {
                Phi[i][j] = std::pow(10, -16);
            }
            Phi[j][i] = Phi[i][j];
        }
    }
    return Phi;
}

std::vector<double> Utot_Calc(const std::vector<std::vector<double>>& Phi) {
    int Particles = Phi.size();
    std::vector<double> Utot(Particles, 0.0);
    
    for (int j = 0; j < Particles; ++j) {
        for (int i = 0; i < Particles; ++i) {
            Utot[j] += Phi[i][j];
        }
    }
    return Utot;
}


int main(){
    //set variables here
    int tsteps = 500;
    int particles = 500;
    double const_a = .01;
    double const_b = .01;
    double const_k = .2;
    double temp = 298.15;
    double plot_size = 100;
    double max_mass = 20;
    double max_move = 1;

    //create matrices
    std::vector<std::vector<double>> x_positions(
        particles,
        std::vector<double>(tsteps,0.0)
    );
    std::vector<std::vector<double>> y_positions(
        particles,
        std::vector<double>(tsteps,0.0)
    );
    //mass
    std::vector<double> masses(particles, 0.0);
    std::vector<double> mass_effect(particles, 0.0);
    std::vector<double> bm(particles, 0.0);

    //generators
    std::random_device rd;  //seed
    std::mt19937 gen(rd()); //Mersenne Twister engine for random values

    std::uniform_real_distribution<double> rand_positions(-plot_size, plot_size); //starting position randomizer
    std::uniform_real_distribution<double> rand_mass(1, max_mass); //starting mass randomizer
    std::uniform_real_distribution<double> rand_move(-max_move, max_move); //mover
    std::uniform_real_distribution<double> bm_odds(0, 1.0); //brownian motion odds

    //fill with random positions for t=0 and masses
    for (int i = 0; i < particles; i++) {
        x_positions[i][0] = rand_positions(gen);
        y_positions[i][0] = rand_positions(gen);
        masses[i] = rand_mass(gen);
        mass_effect[i] = 1/masses[i];
        bm[i] = bm_odds(gen);
    }

    //potential calculations for t=0
    std::vector<double> x_t0(particles), y_t0(particles);
    for (int p = 0; p < particles; ++p) {
        x_t0[p] = x_positions[p][0];
        y_t0[p] = y_positions[p][0];
    }

    //generate all random moves to use later


    std::vector<std::vector<double>> Phi_new = CoordsToPotential(x_t0, y_t0, const_a, const_b, particles);
    std::vector<double> Utot = Utot_Calc(Phi_new);

    //time step loop
    for (std::size_t n = 0; n < tsteps - 1; ++n) {
        std::vector<double> Utot_new(particles, 0.0);
        // Loop over all particles
        for (std::size_t i = 0; i < particles; ++i) {
            //random step for x and y directions
            double dx = rand_move(gen);
            double dy = rand_move(gen);
            //new positions
            double X_new = x_positions[i][n] + dx;
            double Y_new = y_positions[i][n] + dy;

            //get positions of particles at

            //potential energy
            Phi_new = CoordsToPotential(x_positions[n], y_positions[n], const_a, const_b, particles);
            Utot_new = Utot_Calc(Phi_new);

            //energy difference
            double Utot_diff = Utot[i] - Utot_new[i];

            //acceptance probability
            double boltz_prob = std::exp(-Utot_new[i] / const_k);

            //brownian motion probability
            double bm_prob = bm_odds(gen);

            //move particle if accepted, otherwise set equal to previous
            if (boltz_prob > bm_prob && Utot_diff > 0) {
                x_positions[i][n+1] = X_new;
                y_positions[i][n+1] = Y_new;
            } else {
                x_positions[i][n+1] = x_positions[i][n];
                y_positions[i][n+1] = y_positions[i][n];
            }
        }

        // Update Utot for next iteration
        Utot = Utot_new;
    }


return 0;

}
