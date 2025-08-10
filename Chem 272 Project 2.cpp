#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::vector<std::vector<double>> CoordsToPotential(const std::vector<double>& Xinit, const std::vector<double>& Yinit, double a, double b, int Particles) {

    std::vector<std::vector<double>> Phi(Particles, std::vector<double>(Particles, 0.0));
    
    double sigma = std::pow((a/b), (1.0/6.0));
    double r_cutoff = 3 * sigma;

    for (int i = 0; i < Particles; ++i) {
        for (int j = 0; j < Particles; ++j) {
                Phi[i][j] = std::pow(10, -16);
            }
        }

    for (int i = 0; i < Particles; ++i) {
        for (int j = i+1; j < Particles; ++j) {

            double Dx = Xinit[i] - Xinit[j];
            double Dy = Yinit[i] - Yinit[j];
            double r = std::sqrt(Dx*Dx + Dy*Dy);

            if (r > r_cutoff) {
                Phi[i][j] = Phi[j][i] = std::pow(10, -16);
                continue;
            }

            double Utot = (a / std::pow(r, 12)) - (b / std::pow(r, 6));
            Phi[i][j] = Phi[j][i] = Utot;
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

int main() {

int Niter = 10000;
int Num_Particles = 1000;
double a = 0.01;
double b = 0.01;
double k = 0.2;
double Temp = 298.15;

std::random_device rd; // Provides a non-deterministic seed
std::mt19937 gen(rd()); // Mersenne Twister engine seeded by random_devic

std::vector<double> Xinit(Num_Particles, 0);
std::vector<double> Yinit(Num_Particles, 0);
std::vector<double> Mass(Num_Particles, 0);
std::vector<double> No_BM_Probz(Num_Particles, 0);

std::uniform_real_distribution<double> positions(-100.0, 100.0);
std::uniform_real_distribution<double> distru(0.0, 1.0);
std::uniform_real_distribution<double> mass_val(1.0, 20.0);
std::uniform_real_distribution<double> dist_move(-1.0, 1.0);

for (auto &x : Xinit) {
    x = positions(gen);
}
   
for (auto &y : Yinit) {
    y = positions(gen);
}

for (auto &m : Mass) {
    m = mass_val(gen);
}

std::vector<double> Mass_Effect(Num_Particles, 0.0);
for (int i = 0; i < Num_Particles; ++i) {
        Mass_Effect[i] = 1.0 / Mass[i];
}

std::vector<std::vector<double>> Phi = CoordsToPotential(Xinit, Yinit, a, b, Num_Particles);
std::vector<double> Utot = Utot_Calc(Phi);

std::vector<double> dx_move(Num_Particles, 0.0);
std::vector<double> dy_move(Num_Particles, 0.0);
std::vector<double> Xinit_new(Num_Particles, 0.0);
std::vector<double> Yinit_new(Num_Particles, 0.0);
std::vector<std::vector<double>> Phi_new(Particles, std::vector<double>(Particles, 0.0));
std::vector<double> Utot_new (Num_Particles, 0.0);
std::vector<double> Utot_diff(Num_Particles, 0.0);
std::vector<double> Boltzmann_Dist(Num_Particles, 0.0);

for (int n = 0; n < Niter; ++n){
    for (int i = 0; i < Num_Particles; ++i) {
        dx_move[i] = dist_move(gen)*Mass_Effect[i];
    }

    for (int i = 0; i < Num_Particles; ++i) {
         dy_move[i] = dist_move(gen)*Mass_Effect[i];
    }

    for (int i = 0; i < Num_Particles; ++i) {
         Xinit_new[i] = Xinit[i] + dx_move[i];
    }

    for (int i = 0; i < Num_Particles; ++i) {
        Yinit_new[i] = Yinit[i] + dy_move[i];
    }

    Phi_new = CoordsToPotential(Xinit_new, Yinit_new, a, b, Num_Particles);
    Utot_new = Utot_Calc(Phi_new);

    for (int i = 0; i < Num_Particles; ++i) {
        Utot_diff[i] = Utot[i] - Utot_new[i];
    }

    for (int i = 0; i < Num_Particles; ++i) {
        Boltzmann_Dist[i] = std::exp(-1 * Utot_diff[i]/ (k*Temp));
    }

    for (auto &p : No_BM_Probz){
        p = distru(gen);
        }

    for (size_t i=0; i < Utot_diff.size(); ++i) {
        if (Utot_diff[i]>0) {
            Xinit[i] = Xinit_new[i];
            Yinit[i] = Yinit_new[i];
        }
    }

    for (size_t i=0; i < Boltzmann_Dist.size(); ++i) {
        if (Boltzmann_Dist[i] > No_BM_Probz[i] && Utot_diff[i]<0) {
            Xinit[i] = Xinit_new[i];
            Yinit[i] = Yinit_new[i];
        }
    }

    Phi = CoordsToPotential(Xinit, Yinit, a, b, Num_Particles);
    Utot = Utot_Calc(Phi);

    if (n%1000 == 0) {
        plt::title("Particle Movement");
        plt::scatter(Xinit, Yinit);
        plt::xlim(-100, 100);
        plt::ylim(-100, 100);
        plt::show();
    }
    }
    return 0;
}