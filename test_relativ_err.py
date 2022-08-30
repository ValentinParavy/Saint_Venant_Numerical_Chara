# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 22:23:28 2022

@author: valen
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:59:50 2022

@author: valen
"""

import numpy as np
import matplotlib.pyplot as plt
import charateristics_functions_constants as cf
import charateristics_scheme as scheme

'''
ad hoc test
'''

mesh_container = cf.MeshContainer(dx = 0.005, dt = 0.005 ,length = 1, final_time = 0.1)

physical_container = cf.PhysicalContainer(length = 1, final_time = 0.1, g = 9.81, \
                                          beta = 229674, air_p = 0.01 ,\
                                              air_rho = 1.06e3, vessel_area = 1, viscosity = 0)

err = []

for k in range(30):
    v = np.zeros((2, mesh_container.mesh_time_size, mesh_container.mesh_space_size))
    v[:,0,:] = 1 + 0.5 *  np.sin(mesh_container.mesh_space * np.pi / physical_container.length)
    w = v.copy()
    c_0 = (1 + np.exp(-(mesh_container.mesh_space - (physical_container.length/2))**2))**(1/4)
    v = scheme.chara_scheme_explicit(v, c_0, physical_container, mesh_container)
    w = scheme.chara_scheme_implicit(w, c_0, k, physical_container, mesh_container)
    err.append(np.linalg.norm(v - w) / np.linalg.norm(v))

plt.plot(err)