import numpy as np
import matplotlib.pyplot as plt

#Lennard-Jones potential function
def lj_Potential(r, a, b, epsilon= 1e-12):       #calculate the lj_potential between two particles
    r = np.maximum(r, epsilon)
    return (a/ r**12)- (b/r**6)   

def Plot_Particles(x, y, Niter, mass, title= "Particle Locations"):
    plt.figure(figsize=(5,5))
    plt.scatter(x, y, marker='o', s=mass, facecolors='grey', edgecolors='black', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable = 'box')
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    plt.title(f'{title} After N = {Niter} iterations')
    plt.show() 

#Function to calculate the total potential energy
def CoordsToPotential(Xinit, Yinit, a, b, Particles, rcut2=100): 
    #calculate pairwise potential for all particle pairs
    dx = Xinit[:, None] - Xinit  #pairwise difference in x
    dy = Yinit[:, None] - Yinit  #pairwise difference in y 
    r2 = dx**2 + dy**2   #calculated squared distance
    
    mask = r2 < rcut2  #only consider interactions within cutoff distance
    r2[~mask] = 0
    r = np.sqrt(r2)    #Pairwise distances 
    
    #calculate Lennard-Jones Potential for all pairs 
    potential = lj_Potential(r, a, b)

    #Sum the potential only for the upper triangle
    Utot = np.sum(np.triu(potential, k=1))
    return Utot

#Metropolis Monte Carlo move function to update particle positions 
def Utot_Move_Particle(Xinit, Yinit, Niter, Particles, a, b, k, Temp, mass):  
    Utot = CoordsToPotential(Xinit, Yinit, a, b, Particles)

    for n in range(Niter):
        
        dx = np.random.choice([-1, 1], (Particles,)) * (10/mass)
        dy = np.random.choice([-1, 1], (Particles,)) * (10/mass)

        Xinit_new = Xinit + dx
        Yinit_new = Yinit + dy

        #calculate new potential after move
        Utotnew = CoordsToPotential(Xinit_new, Yinit_new, a, b, Particles) 
        Utotdiff = Utotnew - Utot

        #accept the move if energy decreases with probability
        Prob_No_BM = np.random.rand(Particles) #Probability no Brownian motion occurs
        exponent = -Utotdiff / (k * Temp)
        exponent = np.clip(exponent, -500, 500)
        Boltzmann_Dist = np.exp(exponent) #Boltzmann factor
        Brown_Index = np.argwhere((Utotdiff > 0) & (Boltzmann_Dist > Prob_No_BM)).flatten()
    
        Xinit[Brown_Index] = Xinit_new[Brown_Index]
        Yinit[Brown_Index] = Yinit_new[Brown_Index]
        
        #update total potential energy after the move
        Utot = CoordsToPotential(Xinit, Yinit, a, b, Particles)

        if n % 1000 == 0:                    #plot particles every 100 steps
            Plot_Particles(Xinit, Yinit, n, mass)
            plt.pause(0.1)

    return Xinit, Yinit, Utot

class ParticleMotion:

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

    def Run_Particle_Movement(self):
        #plot the initial particle locations
        Plot_Particles(self.Xinit, self. Yinit, 0, self.mass, title= "Initial Particles Locations")

        #Perform Monte Carlo simulation to move particles
        Xfinal, Yfinal, Utot_final = Utot_Move_Particle(self.Xinit, self.Yinit, self.Niter, self.Particles, self.a, self.b, self.k, self.Temp, self.mass)

        #Final state of particles 
        Plot_Particles(Xfinal, Yfinal, self.Niter, self.mass, title= "Final Particle Locations")

#run the particle motion simulation
ParticleMotion().Run_Particle_Movement()