from audioop import add
from dis import dis
from multiprocessing.context import ForkContext
from xmlrpc.server import DocXMLRPCRequestHandler
import pygame
import random
import numpy as np
import matplotlib.pyplot as pyplot
import math
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

file = open('start.txt','r')
number_x_entered = file.readline(2) 
number_y_entered = file.readline(4)
sigma = file.readline(4)
sigma_stop = file.readline(4) 
r = file.readline(6)
a = file.readline(6)
m = file.readline(6)
bx = file.readline(8) 
by = file.readline(10) 
epsilon = file.readline(12) 
additional_particles = file.readline(12) 
vacancy = file.readline(12) 
#D = 1
file.close()
#number_entered = input("Enter the number of nodes horizontally")
#ax = "Enter the boundary conditions X"
#ay = "Enter the boundary conditions Y"
#number_of_particles = number_entered*number_entered
#my_particles = []

#sigma = 10
sigma2 = sigma*sigma
#e = 5
dt = 0.1 # simulation time interval between frames
timesteps = 10 # intermediate invisible steps of length dt/timesteps 

print('Setting the grid')
asq3 = a*(3)**(1/2)
coords = []
for i in range(number_y_entered): #odd
    y_odd = []
    y_odd[i] = sigma+i*asq3
    print(y_odd)
    for j in range(number_x_entered):
       x_odd = []
       x_odd[j] = sigma+j*a
       print(x_odd)
       coords.append((x_odd, y_odd))
       print(coords)

for i in range(number_y_entered): #even
   y_even = sigma+asq3/2+y_even*asq3
   for j in range(number_x_entered-1):
       x_even = sigma+(a/2)+j*a
       print(x_even)
       coords.append((x_even, y_even))
       print(coords)

#if there are additional particles, it enters them into a system with random coordinates
if additional_particles > 0: 
    for i in range(additional_particles):
        rand_x = random.randint(0, 5)
        rand_y = random.randint(0, 5)
        coords.append((rand_x, rand_y))

#if there is a need to remove particles, removes them from a random lattice node
if vacancy > 0:
    for i in range(vacancy):
        coords.pop(random.randrange(len(coords)))

#counting the total number of particles
total_particles = len(coords) - vacancy + additional_particles
print(total_particles)

#функции это конечно здорово, но зачем, если можно сразу считать? Какой тогда у меня ввод?
#если они одноразовые, то смысл их определять?
#--------- distances ----------
#full
def dist(coords):
    distances = []
    for i in range(total_particles):
        x_i = coords[i][0]
        y_i = coords[i][1]
        for j in range(total_particles):
            x_j = coords[j][0]
            y_j = coords[j][1]            
            distances[i][j] = ((x_j-x_i)**2+(y_j-y_i)**2)**(1/2)
            if i==j:
                distances[i][j] = 0
    return distances

#x
def dist_x(x):
    x_distances = []
    for i in range(total_particles):
        x_i = coords[i][0]
        for j in range(total_particles):
            x_j = coords[j][0]            
            x_distances[i][j] = (x_j-x_i)
            if i==j:
                x_distances[i][j] = 0
    return x_distances

#y
def dist_y(y):
    y_distances = []
    for i in range(total_particles):
        y_i = coords[i][1]
        for j in range(total_particles):
            y_j = coords[j][1]            
            y_distances[i][j] = (y_j-y_i)
            if i==j:
                y_distances[i][j] = 0
    return y_distances

#--------- angles ---------- 
#full
def angl(x,y):
    angles_x =[]
    angles_y =[]
    for i in range(total_particles):
        x_i = coords[i][0]
        y_i = coords[i][1]
        for j in range(total_particles):
            x_j = coords[j][0]
            y_j = coords[j][1]
            angle_rad = math.atan2(y_i-y_j,x_i-x_j) #x
            angle_degrees = math.degrees(angle_rad)
            angles_x[i][j] = angle_degrees
            angles_y[i][j] = angle_degrees-90 #y
            if i==j:
                angles_x[i][j] = 0
                angles_y[i][j] = 0
    return angles_x,  angles_y

#x
#y
#---------potential_energy---------- 
def LJ_potential_energy():
    U = []
    for i in range(total_particles):
        for j in range(total_particles):
            r = dist[i]
            U[i][j] = 4*epsilon((sigma/r)**12-(sigma/r)**6)
            if i==j:
                U[i][j] = 0
    return U

    
#--------- interaction force ---------    
def LJ_force(distances):
    F = []
    for i in range(total_particles):
        for j in range(total_particles):
            r = distances[i][j] #надо вынести функции за определения
            F[i][j] = (24/sigma)*epsilon*((sigma/r)**13-(sigma/r)**7) # interaction force
            if i==j:
                F[i][j] = 0
    return F 

#    r2s = r2/sigma2+1
#    r6s = r2s*r2s*r2s
#    f = 24*e*( 2/(r6s*r6s) - 1/(r6s) )
#    f = (12*D)/sigma*((sigma/a)^7-(sigma/a)^13) 

#--------- acceleration ----------
def acc(F,m,angles_x,angles_y):
    acceleration = []
    x_acceleration = []
    y_acceleration = []
    for i in range(total_particles):
        for j in range(total_particles):
            acceleration[i][j] = F[i][j]/m
            x_acceleration[i][j] = math.cos(angles_x[i][j])*acceleration[i][j]
            y_acceleration[i][j] = math.cos(angles_y[i][j])*acceleration[i][j]
    return acceleration, x_acceleration, y_acceleration

def full_x_acc():
    full_x_acceleration = [] #это будет основная матрица ускорений
    for i in range(total_particles): #начинаем просматривать основную матрицу
        for j in range(total_particles):
            full_x_acceleration_temp = 0
            for k in range(total_particles): #подцикл, который ищет подходящие частицы
                for l in range(total_particles):
                    if distances < sigma_stop:
                        full_x_acceleration_temp += x_acceleration[k][l] #складываем ускорение во временную переменную
            full_x_acceleration[i][j] = full_x_acceleration_temp #выгружаем результат в матрицу
            return full_x_acceleration
            #в теории работает, но я явно что то сделал не так

def full_y_acc():
    full_y_acceleration = []
    for i in range(total_particles): 
        for j in range(total_particles):
            full_y_acceleration_temp = 0
            for k in range(total_particles):
                for l in range(total_particles):
                    if distances < sigma_stop:
                        full_y_acceleration_temp += y_acceleration[k][l]
            full_y_acceleration[i][j] = full_y_acceleration_temp
            return full_y_acceleration

#--------- berendsen thermostat ----------
def temp(x, y):
    temperature = []
    for i in range(total_particles):
        for j in range(total_particles):
            energy = (m*(x**2+y**2))/2
            temperature[i][j] = (2*energy[i][j])/(3*K_b)
    return temperature


def borders(x,y):
    if x > bx:
        x = 0

#--------- rewriting parameters ----------
#def update_param(x,y):
#    for i in range(total_particles):
#        x

#class Particle():
    def __init__(self, x, y, vx, vy, size):
        #задание координатa=2
#borders x - 2*sigma+(number_x_entered-1)*a 
#borders y - 2*sigma+number_y_entered*asq3 
#не совсем понял как работают массивы. поэтому разделю координаты в обоих циклах по разным массивам
#а затем, объединю их в один координатный массив, чтобы избедать перезаписи и конфликтов


        self.x = x
        self.y = y
        #задание скоростей
        self.vx = vx
        self.vy = vy
        #для отображения
        self.size = size
        self.colour = (0, 0, 255)
        self.thickness = 2
        self.ax = 0
        self.ay = 0

    def display(self,screen, aa):
        pygame.draw.circle(screen, self.colour, (int(aa*self.x+0.5), int(aa*self.y+0.5)), aa*self.size, aa*self.thickness)

#   def interaction(self):
          

#------------ end class particle ------------
#------------ start main program ------------

#the first cycle will be calculated manually
#так как температура зависит только от скорости, в первом цикле она не считается
#
distances = dist
x_distances = dist_x
y_distances = dist_y
U = LJ_potential_energy
F = LJ_force
acceleration = acc
x_acceleration = acc
y_acceleration = acc
x_velocity = 0
y_velocity = 0
#iterations
#old parameters are written to the corresponding variables
old_coords = coords
old_acceleration = acceleration
old_x_acceleration = x_acceleration
old_y_acceleration = y_acceleration
old_x_velocity = x_velocity
old_y_velocity = y_velocity
x_velocity = old_x_velocity + old_x_acceleration*dt
y_velocity = old_y_velocity + old_y_acceleration*dt
#как аккуратно пройтись по всем матрицам и ничго не сломать?
for i in range(total_particles):
    coords[i][0] = old_coords[i][0] + (x_velocity*dt)
    coords[i][1] = old_coords[i][1] + (y_velocity*dt)
    if coords[i][0] > bx:
        coords[i][0] = coords[i][0] - bx
    if coords[i][0] < bx:
        coords[i][0] = coords[i][0] + bx
    if coords[i][1] > by:
        coords[i][1] = coords[i][1] - by
    if coords[i][1] < by:
        coords[i][1] = coords[i][1] + by

#функуии видимо не работают, проще сразу вставлять их в основной код
distances = dist
x_distances = dist_x
y_distances = dist_y
U = LJ_potential_energy
F = LJ_force

#--------- pygame event loop ----------
screen = pygame.display.set_mode((width, height))
offscreen = pygame.Surface((aafac*width, aafac*height))

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    offscreen.fill(background_colour)

    for k in range(timesteps):
        V_step(my_particles, dt/timesteps)

    for particle in my_particles:
        particle.display(offscreen, aafac)

    pygame.transform.smoothscale(offscreen, (width,height), screen)
    pygame.display.flip()
pygame.quit()

file2 = open('xyz.txt','w')
file2.write('Hello \n World') 
file2.close() 