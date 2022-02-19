from audioop import add
from dis import dis
from xmlrpc.server import DocXMLRPCRequestHandler
import pygame
import random
import numpy as np
import matplotlib.pyplot as plt
import math
import time

#import torch
#import torchvision
#import torch.nn as nn
#import torch.optim as opti
#import torch.nn.functional as funct
#from torch.utils.data import DataLoader
#import torchvision.datasets as datasets
#import torchvision.trsndforms as transforms

background_colour = (255,255,255)
width, height = 500, 500
aafac = 2 # anti-aliasing factor screen to off-screen image
K_b = 1.380649*10**(-23)
tau = 0.1

l=28
file = open('start.txt','r')
lines = file.readlines()
for i in range(l):
    lines[i] = lines[i].rstrip('\n')
    print(f"Вывод строки: {i}) - {lines}")


number_x_entered = lines[1]
number_y_entered = lines[3]
sigma = lines[5]
sigma_stop = lines[7]
r = lines[9]
a = lines[11]
m = lines[13]
bx = lines[15]
by = lines[17]
epsilon = lines[19] 
additional_particles = lines[21] 
vacancy = lines[23] 
termo_temp = lines[25]
iterations = lines[27]

number_x_entered = int(number_x_entered)
number_y_entered = int(number_y_entered)
sigma = float(sigma)
sigma_stop = float(sigma_stop)
r = float(r)
a = float(a)
m = float(m)
bx = float(bx)
by = float(by)
epsilon = float(epsilon)
additional_particles = int(additional_particles)
vacancy = int(vacancy)
termo_temp = float(termo_temp)
iterations = int(iterations)
file.close()

sigma2 = sigma*sigma
dt = 0.1 # simulation time interval between frames

print('Setting the grid')
asq3 = a*(3)**(1/2)
coords = []

for i in range(number_y_entered): #odd
    y_odd = (sigma+i*asq3)
    for j in range(number_x_entered):
       x_odd = (sigma+j*a)
       coords.append([x_odd, y_odd])

for i in range(number_y_entered): #even
    y_even = sigma+(asq3/2)+i*asq3
    for j in range(number_x_entered-1):
        x_even = sigma+(a/2)+j*a
        coords.append([x_even, y_even])

#if there are additional particles, it enters them into a system with random coordinates
if additional_particles > 0: 
    for i in range(additional_particles):
        rand_x = random.randint(0, 5)
        rand_y = random.randint(0, 5)
        coords.append([rand_x, rand_y])

#if there is a need to remove particles, removes them from a random lattice node
if vacancy > 0:
    for i in range(vacancy):
        coords.pop(random.randrange(len(coords)))

#counting the total number of particles
total_particles = len(coords) - vacancy + additional_particles
print(total_particles)

#--------- distances ----------
def dist(coords):
    distances = []
    for i in range(total_particles):
        temp_dist = []
        x_i = coords[i][0]
        y_i = coords[i][1]
        for j in range(total_particles):
            x_j = coords[j][0]
            y_j = coords[j][1]
            temp_dist.append(((x_j-x_i)**2+(y_j-y_i)**2)**(1/2))
        distances.append(temp_dist)
    return distances

#xy
def dist_xy(coords,xy):
    xy_distances = []
    for i in range(total_particles):
        temp_dist = []
        xy_i = coords[i][xy]
        for j in range(total_particles):
            xy_j = coords[j][xy]
            temp_dist.append(abs(xy_j-xy_i))            
        xy_distances.append(temp_dist)
    return xy_distances


#--------- angles ---------- 
def angl(coords):
    angles_x =[]
    angles_y =[]
    for i in range(total_particles):
        temp_x_angl = []
        temp_y_angl = []
        x_i = coords[i][0]
        y_i = coords[i][1]
        for j in range(total_particles):
            x_j = coords[j][0]
            y_j = coords[j][1]
            angle_rad = math.atan2(y_i-y_j,x_i-x_j) #x
            angle_degrees = math.degrees(angle_rad)
            temp_x_angl.append(angle_degrees)
            temp_y_angl.append(angle_degrees-90)
        angles_x.append(temp_x_angl)
        angles_y.append(temp_y_angl)
            #if i==j:
            #    angles_x[i][j] = 0
            #    angles_y[i][j] = 0
    return angles_x,  angles_y
#x
#y

#---------potential_energy---------- 
def LJ_potential_energy(r):
    U = []
    for i in range(total_particles):
        temp_u = []
        for j in range(total_particles):
            if (r[i][j] > sigma_stop or i==j):
                temp_u.append(0)
            else:
                temp_u.append(abs(4*epsilon*((sigma/r[i][j])**12-(sigma/r[i][j])**6)))
        U.append(temp_u)
    return U
    
#--------- interaction force ---------    
def LJ_force(r):
    F = []
    for i in range(total_particles):
        temp_f = []
        for j in range(total_particles):
            if (r[i][j] > sigma_stop or i==j):
                temp_f.append(0)
            else:
                temp_f.append(abs((24/sigma)*epsilon*((sigma/r[i][j])**13-(sigma/r[i][j])**7)))
        F.append(temp_f)
    return F 

#--------- acceleration ----------
def acc(F,angles_x,angles_y):
    acceleration = []
    x_acceleration = []
    y_acceleration = []
    for i in range(total_particles):
        temp_acceleration = []
        temp_x_acceleration = []
        temp_y_acceleration = []
        for j in range(total_particles):
            temp_acceleration.append(F[i][j]/m)
            temp_x_acceleration.append(math.cos(angles_x[i][j])*(F[i][j]/m))
            temp_y_acceleration.append(math.cos(angles_y[i][j])*(F[i][j]/m))
            acceleration.append(temp_acceleration)
            x_acceleration.append(temp_x_acceleration)
            y_acceleration.append(temp_y_acceleration)
    return acceleration, x_acceleration, y_acceleration

def add_of_acc(x_acceleration, y_acceleration):
        x_acc_sum = []
        y_acc_sum = []
        for i in range(total_particles):
            all_accelerations_for_ix = 0
            all_accelerations_for_iy = 0
            for j in range(total_particles):
                all_accelerations_for_ix += x_acceleration[i][j]
                all_accelerations_for_iy += y_acceleration[i][j]
            x_acc_sum.append(all_accelerations_for_ix)
            y_acc_sum.append(all_accelerations_for_iy)
        return x_acc_sum, y_acc_sum

  #  for i in range(total_particles):
  #      for j in range(total_particles):
  #          x_acc_sum += x_acceleration[i][j]
  #          y_acc_sum += y_acceleration[i][j]

  #  return x_acc_sum, y_acc_sum 

#--------- berendsen thermostat ----------
def stop(x_speed, y_speed):
    for i  in range(total_particles):
        x_speed.append(0.0)
        y_speed.append(0.0)
    return x_speed, y_speed

def fs_energy(x_speed, y_speed):
    full_sys_energy = 0
    for i in range(total_particles):
       full_sys_energy += m*(x_speed[i]**2 + y_speed[i]**2)*10**(-20)/2
    return full_sys_energy

def sys_temp(x, y):
    temperature = []
    for i in range(total_particles):
        temp_temp = []
        for j in range(total_particles):
            energy = (m*(x**2+y**2))/2
            temperature[i][j] = (2*energy[i][j])/(3*K_b)
    return temperature

def temp(energy):
    temperature = (2*energy)/(3*K_b*total_particles)

def termo_temp(sys_temp):
    tau = 0.01
    math.sqrt(1+(dt/tau)*((termo_temp/sys_temp)-1))



def scatterplot(x_data, y_data, x_label="", y_label="", title="", color = "r"):

    # Create the plot object
    _, ax = plt.subplots()

    # Plot the data, set the size (s), color and transparency (alpha)
    # of the points
    ax.scatter(x_data, y_data, s = 10, color = color, alpha = 0.75)


    # Label the axes and provide a title
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.clear() 

#------------ end class particle ------------

#------------ start main program ------------

#the first cycle will be calculated manually
x_speed = []
y_speed = []

x_speed, y_speed = stop(x_speed, y_speed)
#full_sys_energy = 0
#full_sys_energy = fs_energy(x_speed, y_speed)

#for i in range(total_particles):
#    x_speed = x_speed*termo_temp
#    x_speed = x_speed*termo_temp

#temperature = temp(x_speed, y_speed)
distances = dist(coords)
x_distances = dist_xy(coords, 0)
y_distances = dist_xy(coords, 1)
U = LJ_potential_energy(distances)
F = LJ_force(distances)
angles_x, angles_y = angl(coords)
acceleration, x_acceleration, y_acceleration = acc(F,angles_x,angles_y)
#here it is calculated how each particle acts on the other

x_acceleration, y_acceleration = add_of_acc(x_acceleration, y_acceleration)
#and here is the total acceleration for each particle, that is, the sum of the lines

#for plot
x_label="X-axis"
y_label="Y-axis"
title="Dynamics"
color = "r"

# Create the plot object
_, ax = plt.subplots()

ax.set_title(title)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)

Time = 0

#-----------------------------------------iterations-----------------------------------------

#old parameters are written to the corresponding variables
for i in range(iterations):
    old_coords = coords
    #old_acceleration = acceleration
    old_x_acceleration = x_acceleration
    old_y_acceleration = y_acceleration
    old_x_speed = x_speed
    old_y_speed = y_speed
    

    for i in range(total_particles):
        ax.clear()
        coords[i][0] = old_coords[i][0] + (old_x_speed[i]*dt) + (old_x_acceleration[i][0])*dt**2
        coords[i][1] = old_coords[i][1] + (old_y_speed[i]*dt) + (old_y_acceleration[i][0])*dt**2
        if coords[i][0] > bx:
            coords[i][0] = coords[i][0] - bx
        if coords[i][0] < bx:
            coords[i][0] = coords[i][0] + bx
        if coords[i][1] > by:
            coords[i][1] = coords[i][1] - by
        if coords[i][1] < by:
            coords[i][1] = coords[i][1] + by
    
    distances = dist(coords)
    x_distances, x_distances = dist_xy(coords)
    U = LJ_potential_energy(r)
    F = LJ_force(r)
    angles_x, angles_y = angl(coords)

    acceleration, x_acceleration, y_acceleration = acc(F,angles_x,angles_y)
    

    for i in range(total_particles):
        x_speed[i] = old_x_speed[i] + ((old_x_acceleration[i] + x_acceleration[i])/2)*dt
        y_speed[i] = old_y_speed[i] + ((old_y_acceleration[i] + y_acceleration[i])/2)*dt
    
    #temperature = temp(x_speed, y_speed)
    time += dt

    ax.scatter(coords, s = 10, color = color, alpha = 0.75)
    #--------- pygame event loop ----------
    
    
    # Plot the data, set the size (s), color and transparency (alpha)
    
    # of the points
# Label the axes and provide a title


