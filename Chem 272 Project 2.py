import numpy as np
import matplotlib.pyplot as plt

def Plot_Particles(Xinit, Yinit, Niter, mass):
    plt.title(f'After N = {Niter} iterations')
    plt.scatter(Xinit, Yinit, marker='o', s=mass, facecolors='grey', edgecolors='black', alpha=0.5)
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    plt.show()

def CoordsToPotential(Xinit, Yinit, a, b, Particles):
    
    X = np.tile(Xinit, (Particles, 1)) #tile 1D array of x coords into matrix
    Y = np.tile(Yinit, (Particles, 1)) #tile 1D array of y coords into matrix
    
    Dx = X - X.transpose() #matrix of distances between xi and xj
    Dy = Y - Y.transpose() #matrix of distances between yi and yj

    Eye = np.eye(Particles) #square matrix with diagonal of 1s with size of Particles
    Ones = np.ones((Particles,)) #1D array of 1s

    r = np.sqrt((Dx**2) + (Dy**2)) #matrix of distances between particles
    r += Eye*(10**16) #diagonal of large infinitely large
    Phi = (a / r**12) - (b / r**6) #calculate potential with lennard jones

    Utot = np.dot(Phi, Ones) #sum of potential each particle feels into one row
    
    return Utot

def Utot_Move_Particle(Xinit, Yinit, Niter, Particles, a, b, k, Temp, mass):
    
    Utot = CoordsToPotential(Xinit, Yinit, a, b, Particles) #starting potential

    for n in range(Niter):
        
        Mass_Effect = 1/mass #----------------------------
        dx = Mass_Effect * np.random.choice([-1, 1], (Particles,)) #choose -1 OR 1 --> 1D array size of Particles for dx
        dy = Mass_Effect * np.random.choice([-1, 1], (Particles,)) #choose -1 OR 1 --> 1D array size of Particles for dy

        Xinit_new = Xinit + dx #move x direction
        Yinit_new = Yinit + dy #move y direction

        Utotnew = CoordsToPotential(Xinit_new, Yinit_new, a, b, Particles) #potential after move
        Utotdiff = Utotnew - Utot

        Utot_change_index = np.argwhere(Utotdiff < 0).flatten() #indices where potential decreases
        
        Xinit[Utot_change_index] = Xinit_new[Utot_change_index] #move particles (by taking indices where potential decreases) to new x coord
        Yinit[Utot_change_index] = Yinit_new[Utot_change_index] #move particles (by taking indices where potential decreases) to new x coord

        Prob_No_BM = np.random.rand(Particles) #Probability no brownian motion occurs
            
        Boltzmann_Dist = np.exp(-Utotdiff/(k*Temp)) #probability of particle moving
        Brown_Index = np.argwhere((Utotdiff > 0) & (Boltzmann_Dist > Prob_No_BM)).flatten() #particles at indices where potential increases AND boltz prob > probability no motion

        Xinit[Brown_Index] = Xinit_new[Brown_Index] #particles who's energy will increase can move if boltzman probability is high
        Yinit[Brown_Index] = Yinit_new[Brown_Index] 
        
        Utot = CoordsToPotential(Xinit, Yinit, a, b, Particles) #recalculate Utot


        if n % 100 == 0:
            Plot_Particles(Xinit, Yinit, n, mass)

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

        self.Xinit = np.random.uniform(Pmin, Pmax, (Particles,)) #initialize x coordinates of particles (1D np array)
        self.Yinit = np.random.uniform(Pmin, Pmax, (Particles,)) #initialize y coordinates of particles (1D np array)
        self.mass = np.random.uniform(1, 20, (Particles,)) #initialize random masses of particles (1D np array)

    def Run_Particle_Movement(self):
        Utot_Move_Particle(self.Xinit, self.Yinit, self.Niter, self.Particles, self.a, self.b, self.k, self.Temp, self.mass)

ParticleMotion().Run_Particle_Movement()
