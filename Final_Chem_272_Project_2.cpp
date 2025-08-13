#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//Function for calculating particles potential based on their cordinates
std::vector<double> CoordsToPotential(const std::vector<double>& Xinit, const std::vector<double>& Yinit, double a, double b, int Particles) {

    std::vector<double> Utot(Particles, 0.0); //Initializes vector for Utot values that has single dimension of the number of Particles. Utot values are zero to be added to later
    
    double sigma = std::pow((a/b), (1.0/6.0)); //Point at which minimum potential occurs.
    double r_cutoff = 3.0 * sigma; //Establishes a cutoff distance at which computation isn't performed to save memory based on the minimum potential point.

    for (int i = 0; i < Particles; ++i) { // Iterates through particles i
        for (int j = i+1; j < Particles; ++j) { //Iterates through particles j avoiding self interaction through j=i+1 for each i before moving to next particle

            double Dx = Xinit[i] - Xinit[j]; //Calcualtes difference of x coordinates between particle i and particle j
            double Dy = Yinit[i] - Yinit[j]; //Calcualtes difference of y coordinates between particle i and particle j
            double r = std::sqrt(Dx*Dx + Dy*Dy); //Calculates distance between particles i and j

            if (r > r_cutoff) { //Skips potential calculations for particles that have a distance farther than the cutoff
                continue;
            }

            else {
            double Utot_part = (a / std::pow(r, 12)) - (b / std::pow(r, 6)); //Calculates Lennard-Jones Potential

            Utot[i] += Utot_part; //Adds Lennard Jones potential to particle i.
            Utot[j] += Utot_part; //Adds same Lennard Jones potential to particle j. Helps with computation as we only need to calculate potentials once.
            }
        }
    }
    return Utot; //Returns the total Lennard-Jones Potential of each Particle
}

int main() {

int Niter = 10000; //Defines the number of time steps
int Num_Particles = 1000; //Defines the number of particles
double a = 0.01; //Defines value of a in Lennard-Jones equation
double b = 0.01; //Defines value of b in Lennard-Jones equation
double k = 0.2; //Defines Boltzmann constant in Boltzmann Factor Calculation (arbitrary value not actual Boltzmann constant so particles actually move)
double Temp = 298.15; //Defines Temperature for Boltzmann Factor Calculation (In this case room temperature on a Kelvin scale. Can be 0 to inf depending on how much particles should move due to brownian motion)
double TempScale = k*Temp; //Defines factor of k and temp as one in Boltzmann Factor Calculation for easier access to variable

std::random_device rd; // Provides a non-deterministic seed
std::mt19937 gen(rd()); // Mersenne Twister engine seeded by random_devic

std::vector<double> Xinit(Num_Particles, 0); //Initializes vector with Particle number of values filled with 0 for initial x positions
std::vector<double> Yinit(Num_Particles, 0); //Initializes vector with Particle number of values filled with 0 for initial y positions
std::vector<double> Mass(Num_Particles, 0); //Initializes vector with Particle number of values filled with 0 for initial masses of particles
std::vector<double> No_BM_Probz(Num_Particles, 0); //Initializes vector with Particle number of values filled with 0 for probability of no Brownian motion

std::uniform_real_distribution<double> positions(-100.0, 100.0); //Produces a random double value normally distributed between -100 and 100 to be used for initial x and y positions
std::uniform_real_distribution<double> distru(0.0, 1.0); //Produces a random double value normally distributed between 0 and 1 to be used for Brownian motion probability
std::uniform_real_distribution<double> mass_val(1.0, 20.0); //Produces a random double value normally distributed between 1 and 20 to be used as the mass of the particle
int moves[] = {-1, 1}; //Creates moves of -1 or 1 for dist_move to index
std::uniform_int_distribution<int> dist_move(0, 1); //Produces int -1 or 1 at random for movement of Particles

for (auto &x : Xinit) { //Auto determines variable type in Yinit container and & references actual Xinit memory location not a copy
    x = positions(gen); //Generates a random double value for each particles initial position x
}
   
for (auto &y : Yinit) { //Determines value and accesses memory location
    y = positions(gen); //Generates a random double value for each particles initial position y
}

for (auto &m : Mass) { //Determines value and accesses memory location
    m = mass_val(gen); //Generates a random double value for each particles mass
}

std::vector<double> Mass_Effect(Num_Particles, 0.0); //Initializes vector for the effect mass has on each particle.
for (int i = 0; i < Num_Particles; ++i) { //Iterates through each particle i
        Mass_Effect[i] = 1.0 / Mass[i]; //Calculates the effect of mass for each particle. 1/Mass ensures larger particles move less.
}

std::vector<double> Utot = CoordsToPotential(Xinit, Yinit, a, b, Num_Particles); //Calculates Lennard-Jones Potential based on initial positions of each particle

std::vector<double> dx_move(Num_Particles, 0.0); //Initializes vector with Particle number of values for movement of particles x coordinates
std::vector<double> dy_move(Num_Particles, 0.0); //Initializes vector with Particle number of values for movement of particles y coordinates
std::vector<double> Xinit_new(Num_Particles, 0.0); //Initializes vector with Particle number of values for new x coordinates
std::vector<double> Yinit_new(Num_Particles, 0.0); //Initializes vector with Particle number of values for new y coordinates
std::vector<double> Utot_new (Num_Particles, 0.0); //Initializes vector with Particle number of values  for new Lennard=Jones Potential
std::vector<double> Utot_diff(Num_Particles, 0.0); //Initializes vector with Particle number of values for difference of Lennard-Jones potential before and after moving the particle
std::vector<double> Boltzmann_Factor(Num_Particles, 0.0); //Initializes vector with Particle number of values for Boltzmann_Factor calculation

for (int n = 0; n < Niter; ++n){ //Iterates through all time steps
    for (int i = 0; i < Num_Particles; ++i) {
        dx_move[i] = moves[dist_move(gen)]*Mass_Effect[i]; //Moves particle in x plane based on random chance of -1 or 1 value. Indexs into moves for value to move. Mass limits the move throug Mass_Effect
        dy_move[i] = moves[dist_move(gen)]*Mass_Effect[i]; //Moves particle in y plane based on random chance of -1 or 1 value. Indexs into moves for value to move. Mass limits the move throug Mass_Effect
    }

    for (int i = 0; i < Num_Particles; ++i) { //Iterates through particles to add moves to current positions
         Xinit_new[i] = Xinit[i] + dx_move[i];
         Yinit_new[i] = Yinit[i] + dy_move[i];
    }

    Utot_new = CoordsToPotential(Xinit_new, Yinit_new, a, b, Num_Particles); //Calculates potential based on new positions

    for (int i = 0; i < Num_Particles; ++i) { //Iterates through particles to determine a potential difference between original positions and new positions.
        Utot_diff[i] = Utot[i] - Utot_new[i];
        Boltzmann_Factor[i] = std::exp(-1 * Utot_diff[i]/ TempScale); //Creates a Boltzmann Facotr for each particle at time step based on temperature scale
    }

    for (auto &p : No_BM_Probz){ //Iterates through vector No_BM_Probz referencing the vector to modify it directly with random probabilities to compare to Boltzmann Factor
        p = distru(gen);
        }

    for (size_t i=0; i < Num_Particles; ++i) { //Iterates through particles to determine move acceptance
        if (Utot_diff[i]>0) { //Determines if the potential decreased to accept or reject move
            Xinit[i] = Xinit_new[i];
            Yinit[i] = Yinit_new[i];
        }
        else if (Boltzmann_Factor[i] > No_BM_Probz[i] && Utot_diff[i]<0) { //Determines if Brownian motion occurred based on if random probability is less than Boltzmann Factor
            Xinit[i] = Xinit_new[i];
            Yinit[i] = Yinit_new[i];
        }
    }

    Utot = CoordsToPotential(Xinit, Yinit, a, b, Num_Particles); //Calculates potential of each particle at new positions

    if (n%100 == 0) { //Plots every 100th time step
        plt::title("After N = " + std::to_string(n) + " iterations");
        plt::scatter(Xinit, Yinit);
        plt::xlim(-100, 100);
        plt::ylim(-100, 100);
        plt::show();
    }
}
return 0;
}