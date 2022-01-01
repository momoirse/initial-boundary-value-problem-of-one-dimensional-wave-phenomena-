import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

LengthS = input("please input the LengthS in meter? (Press enter for default is 1)") or 1
max_iter_time = input("please input the max itteraation type? (Press enter for default is 1000)") or 1000
x0 = input("please input the X0 The place of initial displacement as number between 0 and 1? (Press enter for default is 0.8)") or 0.8
A = input("please input the A as the amount of initial isplacemen in unit of meter? (Press enter for default is 0.01) ") or 0.01
F = input("please input the Frequency in Hertz? (Press enter for default is 400)") or 400

delta_x = 0.04
Landa = 2*LengthS
C = F*Landa
sizeDiscrete=int(LengthS/delta_x)
print (f"the plate brokedown is: {sizeDiscrete} Pcs")

delta_t = 0.1*delta_x/C
LandaCoff = ((C**2*delta_t**2)/delta_x**2)/2

print(delta_t, LandaCoff)
# Initialize solution: the grid of u(k, i, j)
u = np.empty((max_iter_time, sizeDiscrete))

# Initial condition everywhere inside the grid
u_initial = 0

# Set the initial condition
u.fill(u_initial)

for i in range(1, sizeDiscrete-1,1):
    if i < x0*sizeDiscrete:
        u[0, i] = (A*i/sizeDiscrete)/x0
    else:
        u[0, i] = (A*(LengthS - i/sizeDiscrete))/(LengthS-x0)

# Boundary conditions
u_0_t = 0
u_L_t = 0

# Set the boundary conditions
u[:, (sizeDiscrete-1):] = u_L_t
u[:, :1] = u_0_t

def calculate0(u):
    for k in range(0, 2, 1):
        for i in range(1, sizeDiscrete-1, 1):
            u[k + 1, i] = LandaCoff * (u[k][i+1] - 2 * u[k][i] + u[k][i-1]) + u[k][i]
    return u

def calculate(u):
    for k in range(2, max_iter_time-1, 1):
        for i in range(1, sizeDiscrete-1, 1):
            u[k + 1, i] = 2*LandaCoff * (u[k][i+1] - 2 * u[k][i] + u[k][i-1]) +2* u[k][i]-u[k-1][i]
    return u
# Do the calculation here
u = calculate0(u)
u = calculate(u)

#Plot
xplot= np.linspace(0,LengthS,sizeDiscrete) 
fig, ax = plt.subplots()

# def animate(i): 
#     cont=plt.plot(xplot, u[i,:], linewidth=4.0,color="black")
#     return cont
k=0
def animate(i): 
    fig.clear()
    plt.title("Assignment 1 PMNEC - 1D Wave equation CTCS")
    plt.xlabel("Lenght Unit=Meter")
    plt.ylabel("Displacement Unit=Meter")
    plt.grid(True)
    plt.ylim([-A,A])
    plt.xlim([0,LengthS])
    cont=plt.plot(xplot, u[i,:], linewidth=2.0,color="blue")
    return cont
ani = animation.FuncAnimation(fig, animate,frames=max_iter_time, blit=True,interval=2)
plt.show()
# ani.save("Assignment1.mp4")
print("Done!")
