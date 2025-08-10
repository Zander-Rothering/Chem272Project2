#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 16:12:10 2025

@author: anika
"""
import numpy as np
import matplotlib.pyplot as plt

def PlotLocations(Xinit, Yinit, Niter, mass):
    plt.title(f'After N = {Niter} iterations')
    plt.scatter(Xinit, Yinit, marker='o', s=mass, facecolors='grey', edgecolors='black', alpha=0.5)
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    plt.show()
    
    #routine that plots current locations of the particles
    
def Potential(r, a, b, Particles):
    
    Ones = np.ones((Particles,)) #1s
    Phi = (a / r**12) - (b / r**6) #r with large diagonal

    Utot = np.dot(Phi, Ones) #sums rows up
    
    return Utot
    
    
def DistToPotential(Xinit, Yinit, Particles):
    X = np.tile(Xinit, (Particles, 1))
    Y = np.tile(Yinit, (Particles, 1))
    Dx = X - X.transpose()
    Dy = Y - Y.transpose()

    Eye = np.eye(Particles) #1s on diagonal
    
    #Dist += 1e16*np.eye(N)
    
    r = np.sqrt((Dx**2) + (Dy**2))
    r += Eye*(10**16) #large diagonal
    
    return r
    
    
def MoveParticle(Xinit, Yinit, Niter, Particles, a, b, k, Temp, mass):
    
    r = DistToPotential(Xinit, Yinit, Particles)
    Utot = Potential(r, a, b, Particles)

    for n in range(Niter):
        
        Mass_Effect = 1/mass
        dx = Mass_Effect * np.random.choice([-1, 1], (Particles,))
        dy = Mass_Effect * np.random.choice([-1, 1], (Particles,))

        Xinit_new = Xinit + dx
        Yinit_new = Yinit + dy

        r_new = DistToPotential(Xinit_new, Yinit_new, Particles)
        Utotnew = Potential(r_new, a, b, Particles)
        Utotdiff = Utotnew - Utot

        Utot_change_index = np.argwhere(Utotdiff < 0).flatten()
        Xinit[Utot_change_index] = Xinit_new[Utot_change_index]
        Yinit[Utot_change_index] = Yinit_new[Utot_change_index]
        Utot[Utot_change_index] = Utotnew[Utot_change_index]

        Prob_No_BM = np.random.rand(Particles) #Probability no brownian motion occurs
            
        Boltzmann_Dist = np.exp(-Utotdiff/(k*Temp))
        Brown_Index = np.argwhere((Utotdiff > 0) & (Boltzmann_Dist > Prob_No_BM)).flatten()

        Xinit[Brown_Index] = Xinit_new[Brown_Index]
        Yinit[Brown_Index] = Yinit_new[Brown_Index]
        Utot[Brown_Index] = Utotnew[Brown_Index]

        r = DistToPotential(Xinit, Yinit, Particles)
        Utot = Potential(r, a, b, Particles)

        if n % 100 == 0:
            PlotLocations(Xinit, Yinit, n, mass)

    return Xinit, Yinit, Utot



class SimulateParticles():
    def __init__(self, Niter=10000, Particles=1000, a=0.01, b=0.01, Pmin=-100, Pmax= 100, k=0.2, Temp=298.15):
        self.Niter = Niter
        self.Particles = Particles
        self.a = a
        self.b = b
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.k = k
        self.Temp = Temp

        self.Xinit = np.random.uniform(Pmin, Pmax, (Particles,))
        self.Yinit = np.random.uniform(Pmin, Pmax, (Particles,))
        self.mass = np.random.uniform(1, 20, (Particles,))
        
    def run(self):
        MoveParticle(self.Xinit, self.Yinit, self.Niter, self.Particles, self.a, self.b, self.k, self.Temp, self.mass)
        

SimulateParticles().run()
