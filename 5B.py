# This is the .py file for Section 5B

# Task 1
# For Z = grand partition function
# mu = chemical potential
# eta(i) = index specific energy level
# Pi = giant pi = product thingy
# beta = 1/kBT

# omega = Pi[i=1, M](1 + exp(-beta(eta(i) - mu)))

# Task 2

# Part a 

# There can be i particles at the E=0 energy state and N-k particles at the E=eta energy state. 
# This means that in total, there are N + 1 microstates for the system.

# Part b

# Let Zc = classical partition function under canonical ensemble
# n = energy specific microstate
# Zc = sum[n=0, N](exp(-beta*n*eta))

# In order to evanluate the propability of a specific microstate, divide the expression for the microstate
# by the partition function. This makes sense because you are weighing the specific state against the 
# probabilistic authority of all possible states.

# Part c

import numpy as np
import matplotlib.pyplot as plt

def calc_Z(N, epsilon, beta):
    return (1 - np.exp(-beta * epsilon * (N + 1))) / (1 - np.exp(-beta * epsilon))

def avg_particles(N, epsilon, T):
    kB = 1.0  # Boltzmann constant, set to 1 to properly scale the plot and display the relationship
    beta = 1 / (kB * T)  
    
    # Partition function
    Z_C = calc_Z(N, epsilon, beta)
    
    dlnZ_depsilon = -N * (N + 1) * np.exp(-beta * epsilon) / (
        (1 - np.exp(-beta * epsilon * (N + 1))) * (1 - np.exp(-beta * epsilon))) * epsilon

    n_epsilon_avg = -1 / beta * dlnZ_depsilon
    
    n_0_avg = N - n_epsilon_avg
    
    return n_0_avg, n_epsilon_avg

def plot_avg_particles(N, epsilon, Tmin, Tmax):
    T = np.linspace(Tmin, Tmax, 100) 
    n_0_avg, n_epsilon_avg = avg_particles(N, epsilon, T)
    
    plt.figure(figsize=(10, 5))
    plt.plot(T, n_0_avg, label=r'$\langle n_0 \rangle_C$', color='blue')
    plt.plot(T, n_epsilon_avg, label=r'$\langle n_\epsilon \rangle_C$', color='red')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Average Number of Particles')
    plt.title('Average Number of Particles in the Ground and Excited States')
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig('5B_c_plot.png')


N = 100  
epsilon = 1.0 

# Generate the plot
plot_avg_particles(N, epsilon, 0.1, 10)

# The plot conveys the interdependence the average values have on each other. The average values themselves
# are also embedded in the code above.

# Part d
# I am far less familiar with quantum boseman statistics so this is where the work is likely gpoing to
# taper off a little bit but I will keep plugging away. Don't mind my commit time! I am actually going to
# wrap this up a bit more tomorrow.

