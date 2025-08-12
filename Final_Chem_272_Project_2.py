import numpy as np
import matplotlib.pyplot as plt

def Plot_Particles(x, y, Niter, mass, title= "Particle Locations"):
    plt.figure(figsize=(6,6))
    plt.title(f'{title} After N = {Niter} iterations')
    plt.scatter(x, y, marker='o', s=mass, facecolors='grey', edgecolors='black', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable = 'box')
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    plt.show()
    
#Lennard-Jones potential function
def LJ_Potential(r, a, b, epsilon= 1e-12):       #calculate the lj_potential between two particles
    r = np.maximum(r, epsilon)
    return (a/ r**12)- (b/r**6)   

#energy constribution of single particle
def _Ei(i, X, Y, a, b, rcut=1000):
    dx = X[i] -X #Vector difference between particle and all other particles in x direction
    dy = Y[i] -Y #Vector difference between particle and all other particles in y direction
    r2= dx**2 + dy**2 #Calculates difference between particle and all other particles
    r2[i] = np.inf #To avoid self interation sets distance with self to an infinite distance
    
    mask= r2 < (rcut *rcut) #Creates a mask where an particle distance beyond the cut off is not calculated
    
    r= np.sqrt(r2[mask]) #Calcualtes r based on mask to save compuational time.
    
    return np.sum(LJ_Potential(r, a, b)) #returns energy of a single particle

#Function to calculate the total potential energy
def CoordsToPotential(X, Y, a, b, rcut=1000): 
    #calculate pairwise potential for all particle pairs
    dx = X[:, None] - X[None, :] #pairwise difference in x
    dy = Y[:, None] - Y[None, :] #pairwise difference in y 
    r2 = dx**2 + dy**2   #Pairwise distances 
    np.fill_diagonal(r2, np.inf)  # exclude self interaction
    
    #creates mask to ignore distances beyond a cut off point
    mask = r2 < rcut**2
    
    #calculate Lennard-Jones Potential for all pairs 
    potential = np.zeros_like(r2)
    r = np.sqrt(r2[mask])
    potential[mask] = LJ_Potential(r, a, b)

    #Sum the potential energy per particle
    Utot = np.sum(np.triu(potential, k=1))
    return Utot

#Metropolis Monte Carlo move function to update particle positions 
def Utot_Move_Particle(X, Y, Niter, a, b, k, Temp, mass, Particles, rcut=20.0, step=1.0):
    
    kBT = k*Temp #Combines Botlzmann constant (arbitrary value) and Temp into one variable
    X = X.copy() #Copies the position arrays so the original values are not changed outside of Utot_Move_Particle
    Y = Y.copy()
    
    Utot = CoordsToPotential(X, Y, a, b, rcut) #Calculates initial Lennard-Jones Potential

    for n in range(Niter): #Iterates over time steps Niter
        for i in range(Particles): #Iterates over each particle at each time step
            dx = step * np.random.choice([-1, 1]) /  mass[i] #Determines step size in random direction x with movement limitted by the mass
            dy = step * np.random.choice([-1, 1]) / mass[i] #Determines step size in random direction y with movement limitted by the mass
        
            X_new = X[i] + dx #Calculates new position after moving in x direction
            Y_new = Y[i] + dy #Calculates new position after moving in y direction
           
            Ei_old =_Ei(i,X,Y,a, b, rcut=rcut) #Calculates energy of particle with current positions
            X_temp = X.copy() #Copies current x position
            Y_temp = Y.copy() #Copies current y position
            X_temp[i] = X_new #Changes X position for particle to new position
            Y_temp[i] = Y_new #Changes Y position for particle to new position
            
            Ei_new = _Ei( i, X_temp, Y_temp, a, b, rcut) #Calculates energy of particle with new positions
            dE = Ei_new - Ei_old #Calculates difference between new and old energy of particle
        
            if dE < 0: #If energy difference is less than 0 allows move and adds energy to total energy
                X[i], Y[i] = X_new, Y_new
                Utot += dE
            else:
                Prob_No_BM = np.random.rand() 
                if np.exp(-dE / kBT) > Prob_No_BM:  #If energy difference is larger allows Brownian Motion to occur with random probability
                    X[i], Y[i] = X_new, Y_new
                    Utot += dE
        
        if n == 0:                          #plots particles at initial locations
            Plot_Particles(X, Y, n, mass, title= "Initial Particles Locations")
        
        elif n % 100 == 0:                    #plot particles every 100 steps
            Plot_Particles(X, Y, n, mass)
            
        if n == Niter -1:                 #plots particles at final locations
            Plot_Particles(X, Y, n, mass, title= "Final Particle Locations")
            

    return X, Y, Utot

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

        #Perform Monte Carlo simulation to move particles and plot them
        X, Y, Utot = Utot_Move_Particle(self.Xinit, self.Yinit, self.Niter, self.a, self.b, self.k, self.Temp, self.mass, self.Particles)

#run the particle motion simulation
ParticleMotion().Run_Particle_Movement()