# -*- coding: utf-8 -*-
"""
Barrier Library
"""

#Triangle
V = 0*x
for i in range(len(V)):
    if x[i] > -0.5 and x[i] < 0.25:
        
        V[i] = V0 * (x[i] - (-0.5)) / (0.25 - (-0.5))
    elif x[i] >= 0.25 and x[i] < 1.0:
        
        V[i] = V0 * (1.0 - x[i]) / (1.0 - 0.25)
"""
starts at x=-0.5 with V=0, linearly increases to V=V0 at x=0.25 (peak), 
linearly decreases back to V=0 at x=1.0. total width =1.5 units, peak at 0.25

    -> shows how smooth, gradual potential affects transmission. less
    reflection than sharp rectangular barrier because the wave has time to
    "adjust" to the changing potential
"""

#Gaussian
V = 0*x
centre = 0.25
width = 0.3
for i in range(len(V)):
    V[i] = V0 * np.exp(-(x[i] - centre)**2 / (2*width**2))