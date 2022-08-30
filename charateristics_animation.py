# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:59:28 2022

@author: valen
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from charateristics_functions_constants import *
from charateristics_scheme import *

V = scheme(dx,dt,V,Lambd)

U = V[0] + V[1]
H = ((V[1] - V[0])/2)**2 / g


plt.plot(H[int(0.7/dt)])
plt.show()

'''
fig = plt.figure()
plt.ylim(np.min(H),np.max(H))

for i in range(mesh_time.size):
    plt.plot(mesh_space,H[i])
    plt.ylim(np.min(H),np.max(H))
    plt.pause(0.05)
    plt.clf()    
'''

'''
# Animation

fig = plt.figure()  # initialise la figure
line, = plt.plot([], []) 

def animate(i):
    line.set_data(mesh_space, U[i])
    return line,

ani = animation.FuncAnimation(fig, animate, frames=mesh_time.size)

plt.show()
'''