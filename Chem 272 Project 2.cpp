#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::vector<std::vector<double>> CoordsToPotential(const std::vector<double>& Xinit, const std::vector<double>& Yinit, double a, double b, int Particles) {

    int Particles = Xinit.size()

    std::vector<std::vector<double>> Phi(Particles, std::vector<double>(Particles, 0.0));
    
    for (int i = 0; i < Particles; ++i) {
        for (int j = 0; j < Particles; ++j) {
            if (i == j) continue;

            double Dx = Xinit[i] - Xinit[j];
            double Dy = Yinit[i] - Yinit[j];
            double r = std::sqrt(Dx*Dx + Dy*Dy);

            if (r > 0) {
                Phi[i][j] = (a / std::pow(r, 12)) - (b / std::pow(r, 6));
            } else {
                Phi[i][j] = std::pow(10, -16);
            }
        }
    }
    return Phi;

    std::vector<double> Utot_Calc(const std::vector<std::vector<double>>& Phi); {
    int Particles = Phi.size()
    std::vector<double> Utot(Particle, 0.0)
    
        for (int j = 0; j < Particles; ++j) {
            for (int i = 0; i < Particles; ++i)
                Utot[j] += Phi[i][j];
        }
    }
    return Utot;
}

int main() {

int Num_Particles = 1000
double a = 0.01
double b = 0.01

std::random_device rd; // Provides a non-deterministic seed
std::mt19937 gen(rd()); // Mersenne Twister engine seeded by random_devic

std::vector<std::vector<double>> Xinit(1, std::vector<double>(Num_Particles));
std::vector<std::vector<double>> Yinit(1, std::vector<double>(Num_Particles));
std::vector<std::vector<double>> Mass(1, std::vector<double>(Num_Particles));

std::uniform_real_distribution<double> dist(0.0, 1.0);
std::uniform_real_distribution<double> dist_mass(1, 20.0);

for (auto &row : Xinit){
    for (auto &x : Xinit) {
        x = dist(gen);
    }
}    
for (auto &row : Yinit){
    for (auto &y : Yinit) {
        y = dist(gen);
    }
} 
for (auto &row : Mass){
    for (auto &m : Mass) {
        m = dist(gen);
    }
}

std::vector<std::vector<double>> Phi = CoordsToPotential(Xinit, Yinit, a, b);
std::vector<double> Utot = Utot_Calc(Phi)

return 0;

}




