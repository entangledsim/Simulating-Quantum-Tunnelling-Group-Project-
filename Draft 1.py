import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
"""
simulating infinite square well; must have real boundary coniditions where 
psi is 0
"""

"""
Constants
"""
m = 1
hbar = 1
xmin = -6.5 #boundaries (xaxis size)
xmax = 6.5
N = 1000 #pieces it breaks into
x = np.linspace(xmin,xmax,N+1) #x values; adding 1 extra value to N (1001 pieces)
dx = x[1]-x[0]
p = 40
V0 = p**2/(2*m)
sig = 0.15 #width?
x0 = 2 #where its going to start;not to interfere with wall as +-6.5
V=0*x
#print(len(V)); prints length of V


"""
Rectangle Barrier
"""
for i in range(len(V)):
    if x[i]>0 and x[i]<1: #width of well
        V[i]=V0

#plt.plot(x,V) #plotting potential

"""
Initial wavefunction: 
        - must be 2 less than the number of x values (x=N+1) so this must be N-1
            x[1:-1] means x omitting first and last value
"""
Psi0 = np.exp( -(x[1:-1]+x0)**2/sig**2)*np.exp(1j*p*(x[1:-1]+x0)) #1j = i
#normalising Psi0:
A = np.sum(np.abs(Psi0)**2*dx) #this is the integral
Psi0 = Psi0/np.sqrt(A) #this makes the integral above equal 1 (normalised)
#plt.plot(x[1:-1],np.abs(Psi0)) #need the np.abs() on Psi0 because its complex. again x[1:-1] because do not want to plot it on edges



"""
Example matrix code:
    Nn = 6
    H = np.diag(np.ones(Nn)) -> diagonal matrix with ones
    
    H = np.diag(np.ones(Nn)) + 3*np.diag(np.ones(Nn-1),1)
       - Nn-1 because diagonal above the first one has 1 fewer elements
       - ,1 moves it up
    
    print(H)
"""


"""
Making matrix of Hamiltonian to find the eigenvalues
"""
H = (hbar**2/(m*dx**2))*np.diag(np.ones(N-1))+V[1:-1]*np.diag(np.ones(N-1))+ (-hbar**2/(2*m*dx**2))*np.diag(np.ones(N-2),1)+ (-hbar**2/(2*m*dx**2))*np.diag(np.ones(N-2),-1)
#print(H.shape)

E,psi = np.linalg.eigh(H) #returns eigenvalues and eigenvectors; vectors are in columns but needed in rows
psi = psi.T #transposes so theyre rows

#normalise them
A = np.sum(np.abs(psi[0])**2*dx) #integral over all space
psi = psi/np.sqrt(A)
#some of the solutions:
#plt.plot(x[1:-1],psi[0])
#plt.plot(x[1:-1],psi[1])
#plt.plot(x[1:-1],psi[2])

"""
Finding the C coefficients
"""
#plt.plot(x,0.003*V)

#---------

c = 0*Psi0 #creates an array of complex numbers; timesing 0 by complex numbers is still complex
for i in range(len(c)):
    c[i] = np.sum(np.conj(psi[i])*Psi0*dx)   
#print(c[0],c[1],c[2])

"""
Calculating Psi(x,t)
    multiple plots in a row with a bigger change in t showing the wavefunction as its travelling and hitting the potential barrier (steady state snapshots)
"""

t = 0
dt = 0.001
Psi = 0*psi[0] #set psi to 0
for i in range(len(c)):
    Psi = Psi + c[i]*psi[i]*np.exp(-1j*E[i]*t/hbar)

#plt.plot(x[1:-1],np.abs(Psi))
#plt.show()


#t = t + 15*dt
#Psi = 0*psi[0] #set psi to 0
#for i in range(len(c)):
#    Psi = Psi + c[i]*psi[i]*np.exp(-1j*E[i]*t/hbar)

#plt.plot(x[1:-1],np.abs(Psi))
#plt.show()

#t = t + 100*dt
#Psi = 0*psi[0] #set psi to 0
#for i in range(len(c)):
#    Psi = Psi + c[i]*psi[i]*np.exp(-1j*E[i]*t/hbar)

#plt.plot(x[1:-1],np.abs(Psi))
#plt.show()



fig, ax = plt.subplots(figsize=(10,6), dpi=300)



#wave packet line
line, = ax.plot([], [], 'b-', linewidth=1.5, label='|Ψ(x,t)|', color='magenta', zorder=1) #creates empty line object that updates each time. [] means start with no data
# Potential line (thinner, on top)
potential_line, = ax.plot(x, 0.003*V, 'b-', linewidth=1.5, label='V(x) \times 0.003', alpha=0.7, zorder=2) #plots potential barrier
#setting plot boundaries
ax.set_xlim(xmin, xmax) #-6.5 to 6.5 for x
ax.set_ylim(-3, 7) #0 to 1.2 for y
ax.axhline(y=0, color='gray', linewidth=0.5,alpha=0.5)
ax.set_xlabel('Position (x)', fontsize=12, fontweight='bold')
ax.set_ylabel('Wavefunction', fontsize=12, fontweight='bold')
ax.set_title('Wave Packet Through Potential Barrier ', fontsize=14, fontweight='bold', pad=20)
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

t = 0
dt = 0.001
frames = 300
pause_frames = 30


def init():
    line.set_data([], [])
    return line, 

def animate(frame):
    global t
    
    if frame >= frames:
        frame = frames - 1
        
    t = frame * dt
    #print(f"Frame: {frame}, Time: {t}")
    
    #original physics from before calculating Psi(x,t)
    Psi = 0*psi[0]
    for i in range(len(c)):
        Psi = Psi + c[i]*psi[i]*np.exp(-1j*E[i]*t/hbar)
    
    line.set_data(x[1:-1], np.abs(Psi)**2) #updates with new data. prob density (blue envelope)
    ax.set_title(f'Wave Packet Through Potential Barrier (t= {t:.3f})', fontsize=14) #updates title to show current time
    return line, 

anim = FuncAnimation(fig, animate, init_func=init, frames=frames + pause_frames, 
                     interval=20, blit=False, repeat=True) #faster with blit=True, but doesnt update time in the title

plt.tight_layout()
plt.show()