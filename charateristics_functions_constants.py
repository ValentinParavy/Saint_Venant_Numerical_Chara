# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:27:54 2022

@author: valen
"""

import numpy as np

class PhysicalContainer:
    '''
    Contains the physical variables of the prob
    lem.
    
    Assume that both the time and space origins are zero
    '''
    def __init__(self,length,final_time,g,beta,air_p,air_rho,vessel_area,viscosity):
        self.length = length
        self.final_time = final_time
        self.g = g
        self.beta = beta
        self.air_p = air_p
        self.air_rho = air_rho
        self.vessel_area = vessel_area
        self.mu = viscosity

class MeshContainer:
    '''
    Contains the variables which relate to the meshing
    
    '''
    def __init__(self,dx,dt,length,final_time):
        self.dt = dt
        self.dx = dx
        self.mesh_space_size = int(length / dx) + 1
        self.mesh_time_size = int(final_time / dt) + 1
        self.mesh_space = np.linspace(0,length, self.mesh_space_size)
        self.mesh_time = np.linspace(0,final_time, self.mesh_time_size)
        
        
def compute_pressure(physical_constants,a):
    a = a * (a > 0)
    return physical_constants.air_p + physical_constants.beta*(np.sqrt(a) - np.sqrt(physical_constants.vessel_area))

def compute_wave_speed(physical_constants,a):
    a = a * (a > 0)
    return np.sqrt((physical_constants.beta * np.sqrt(a)) / (2 * physical_constants.air_rho))

def wave_speed_to_area(physical_constants,c):
    return (((c**2) * 2 * physical_constants.air_rho) / physical_constants.beta)**2
      
def riemann_to_physical(v, c_0, physical_constants, converter_speed = wave_speed_to_area):
    return np.array(((v[0] + v[1])/2, converter_speed(physical_constants, c_0 + (v[0] - v[1]) / 8)))

def physical_to_riemann(physical_constants, u, a, compute_wave_speed = compute_wave_speed):
    a_0 = a[0,:]
    c_0 = compute_wave_speed(physical_constants, a_0)
    c = compute_wave_speed(physical_constants,a)
    return np.array(((u + 4*(c - c_0)), u - 4*(c - c_0)))

def eigen_values(v_1, v_2, c_0):
    return np.array(((5/8)*v_1 + (3/8)*v_2 + c_0, (3/8)*v_1 + (5/8)*v_2 - c_0))

def q_eigen_values(v_1, v_2, c_0, dt, rectangle_rule = True):
    return  int(rectangle_rule) * dt * np.array(((5/8)* v_1 + (3/8)* v_2 + c_0, (3/8)* v_1 + (5/8)* v_2 - c_0)) 

def r_terms(v, c_0, dc_0, physical_constants, dx, riemann_to_physical, compute_wave_speed):
    physical_variables = riemann_to_physical(v, c_0, physical_constants)
    u = physical_variables[0]
    a = physical_variables[1] * (physical_variables[1] > 0)
    c = compute_wave_speed(physical_constants, a)
    return  np.array((-8 * np.pi * physical_constants.mu * (u / a) - 4 * (dc_0) *  (u + c - c_0), -8 * np.pi * physical_constants.mu * (u / a) + 4 * (dc_0) *  (u - c + c_0)))

def r_terms_simplified(v, c_0, dc_0, physical_constants):
    u = (v[0] + v[1]) / 2
    delta_c = (v[0] - v[1]) / 8
    a = (4 * (physical_constants.air_rho**2) / physical_constants.beta) * (delta_c + c_0)**4 
    return 4 * dc_0 * np.array([-u - delta_c, u - delta_c]) - 8 * np.pi * physical_constants.mu * (u/a)

def r_terms_simplified_no_viscosity(v, dc0):
    return 4 * dc0 * np.array([(5/8) * v[0] + (3/8) * v[1], (3/8) * v[0] + (5/8) * v[1]])
            
def interpo(v, x, n, i, meshcontainer):
    '''
    
    Compute the interpolated value of v_1 and v_2 at the instant t_n and point (x_1,x_2).
    
    '''
    dx = meshcontainer.dx
    M = meshcontainer.mesh_space_size
    y = int(int(1) * (x <= dx) + int((M-1)) * (x > (M-1) * dx) + int(x/dx) * (x > dx and x < (M-1) * dx))
    return v[i, n, y] + (v[i, n, y-1] - v[i, n, y]) * (y - (x/dx))

'''
def contra(g, dt, v, n, i, c_0, meshcontainer, eigen_values = eigen_values, interpo = interpo):
    v  = interpo(v, g, n, i, meshcontainer)
    return g - dt * eigen_values(v[0,n,i], v_2, c_0)[i] 
'''        
        
                 