#include <vector>
#include <iostream>
#include <random>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

CoordsToPotential( std::vector<double> Xinit, std::vector<double> Yinit, double a, double b, double Particles) {

std::vector<std::vector<double>> Eye(double n) {
    std::vector<std::vector<double>> Id_Matrix(n, std::vector<double>(n, 0.0));
    for (double i = 0; i < n; ++i) {
        mat[i][i] = 1.0;
    }
    return Id_Matrix;
}

std::vector<std::vector<double>> Ones(Particles, std::vector<double>(Particles, 0.0));

std::vector<std::vector<double>> (std::vector<std::vector<double>> Dist) {
    for (<std::vector<double>> i = 0; i <n; ++i) {
        sqrt()
    }
    return 
}



}





int main() {
std::random_device rd; // Provides a non-deterministic seed
std::mt19937 gen(rd()); // Mersenne Twister engine seeded by random_devic

std::vector<double> Xinit(Particles);
std::vector<double> Yinit(Particles);
std::vector<double> Mass(Particles);

std::uniform_real_distribution<double> dist(0.0, 1.0);
std::uniform_real_distribution<double> dist_mass(1, 20.0);


for (auto &x : Xinit) {
    x = dist(gen);
}
for (auto &y : Yinit) {
    y = dist(gen);
}
for (auto &m : Mass) {
    m = dist_mass(gen);
}

Utot_Move_Particel

return 0;

}




