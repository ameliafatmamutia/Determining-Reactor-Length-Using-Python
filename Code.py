# Case: Determining Reactor Length

# Import library
import numpy as np
import matplotlib.pyplot as plt


# Data
CA0 = 0.2  # Concentration of A at inlet reactor (kmol/m3)
CB0 = 0    # Concentration of B at inlet reactor (kmol/m3)
U = 7.5    # Gas velocity (m/s)
k1 = 8     # constants
k2 = 3     # constants
k3 = 0.01  # constants
x0 = 0     # initial conversion at z=0
xf = 0.90  # target conversion


# Function for reaction rate
def reactionRate(xi):
    CA = CA0*(1-xi)
    CB = CB0*(1-xi)
    _RA = k1*CA/(1+k2*CA+k3*CB)
    fi = 1/_RA
    return fi


# Function for reactor length
def reactorLength(xfj):
    #Integration is carried out using Simpson method
    n = 20 # number of integration segment
    dx = (xfj-x0)/n
    xj = np.linspace(x0,xfj,n+1)
    f = np.zeros(n+1)
    fsim = np.zeros(n+1)
    for j in range(n+1):
        f[j] = reactionRate(xj[j])
        if j == 0:
            fsim[j] = f[j]
        elif j == n:
            fsim[j] = f[j]
        elif (-1)**j < 0:
            fsim[j] = 4*f[j]
        else:
            fsim[j] = 2*f[j]
    z = U*CA0*dx/3*sum(fsim)       
    return z


# Calling out the results
nx = 21
x = np.linspace(x0,xf,nx)
z = np.zeros(nx)
for i in range(0,nx):
    z[i] = reactorLength(x[i])
    
    
# Showing the Table
print('Calculation Results:')
garis = '-'*53
tabel = np.zeros([nx,2])
tabel[:,0] = z
tabel[:,1] = x*100
header = ['Reactor Length, m','Reaction Conversion, %']
print(garis)
print('{:^25s} {:^25s}'.format(*header))
print(garis)
for baris in tabel:
    print('{:^25.2f} {:^25.2f}'.format(*baris))
print(garis)


# Showing the Graph
plt.plot(z,x*100,'o:r')
plt.xlabel('Reactor Length, m')
plt.ylabel('Reaction Conversion, %')


