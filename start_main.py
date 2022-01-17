from dis import dis
from multiprocessing.context import ForkContext
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

file1 = open('start.txt','r')
number_x_entered = readline(2) 
number_y_entered = readline(4)
sigma = readline(4)
sigma_stop = readline(4) 
r = readline(6)
a = readline(6)
m = readline(6)
bx = readline(8) 
by = readline(10) 
epsilon = readline(12) 
#D = 1
file1.close()
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
while i in range(number_y_entered): #odd
    y_odd = []
    y_odd[i] = sigma+i*asq3
    print(y_odd)
   for j in range(number_x_entered):
       x_odd = []
       x_odd[j] = sigma+j*a
       print(x_odd)
       coords.append((x_odd, y_odd))
       print(coords)

while i in range(number_y_entered): #even
   y_even = sigma+asq3/2+y*asq3
   for j in range(number_x_entered-1):
       x_even = asq3/2+i*a
       print(x_even)
       coords.append((x_even, y_even))
       print(coords)
total_particles = len(coords)
print(total_particles)

#функции это конечно здорово, но зачем, если можно сразу считать? Какой тогда у меня ввод?
#если они одноразовые, то смысл их определять?
#--------- distances ----------
def dist(x,y):
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
    print(distances)
#--------- angles ---------- 
def angl(x,y):
    angles =[]
    for i in range(total_particles):
        x_i = coords[i][0]
        y_i = coords[i][1]
        for j in range(total_particles):
            x_j = coords[j][0]
            y_j = coords[j][1]
            angle_rad = math.atan(y_i-y_j,x_i-x_j)
            angle_degrees = math.degrees(angle_rad)
            angles[i][j] = angle_degrees
            if i==j:
                angles[i][j] = 0

    print(angles)
       # x_proj = F*math.cos(angle_degrees) #?
#---------potential_energy---------- 
def LJ_potential_energy():
    U = []
    for i in range(total_particles):
        for j in range(total_particles):
            r = dist[i]
            U[i][j] = 4*epsilon((sigma/r)**12-(sigma/r)**6)
            if i==j:
                U[i][j] = 0
    print(U)

    
#--------- interaction force ---------    
def LJ_force(p1,p2):
    F = []
    for i in range(total_particles):
        for j in range(total_particles):
            r = distances[i][j] #надо вынести функции за определения
            F[i][j] = (24/sigma)*epsilon*((sigma/r)**13-(sigma/r)**7) # interaction force
            if i==j:
                F[i][j] = 0
    print(F)

    p1.ax += f*(rx/r2)
    p1.ay += f*(ry/r2)
    p2.ax -= f*(rx/r2)
    p2.ay -= f*(ry/r2)

#    r2s = r2/sigma2+1
#    r6s = r2s*r2s*r2s
#    f = 24*e*( 2/(r6s*r6s) - 1/(r6s) )
#    f = (12*D)/sigma*((sigma/a)^7-(sigma/a)^13) 

#--------- acceleration ----------



class Particle():
    def __init__(self, x, y, vx, vy, size):
        #задание координатa=2
#borders x - 2*sigma+(number_x_entered-1)*a 
#borders y - 2*sigma+number_y_entered*asq3 
#не совсем понял как работают массивы. поэтому разделю координаты в обоих циклах по разным массивам
#а затем, объединю их в один координатный массив, чтобы избедать перезаписи и конфликтов
while i in range(number_y_entered): #odd
    y_odd = []
    y_odd[i] = sigma+i*asq3
    print(y_odd)
   for j in range(number_x_entered):
       x_odd = []
       x_odd[j] = sigma+j*a
       print(x_odd)
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

    def interaction(self):
          

#------------ end class particle ------------
#------------ start main program ------------









#--------- berendsen thermostat ----------



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
        Verlet_step(my_particles, dt/timesteps)

    for particle in my_particles:
        particle.display(offscreen, aafac)

    pygame.transform.smoothscale(offscreen, (width,height), screen)
    pygame.display.flip()
pygame.quit()

file2 = open('xyz.txt','w')
file2.write('Hello \n World') 
file2.close() 