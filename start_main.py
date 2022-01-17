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
r = readline(6)
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
def LJ_potential_energy():
    U = 4*epsilon((sigma/r)**12-(sigma/r)**6)
def LJ_force(p1,p2):
    rx = p1.x - p2.x
    ry = p1.y - p2.y
    r = 

    f = (24/sigma)*epsilon*((sigma/r)**13-(sigma/r)**7) # interaction force

    p1.ax += f*(rx/r2)
    p1.ay += f*(ry/r2)
    p2.ax -= f*(rx/r2)
    p2.ay -= f*(ry/r2)

#    r2s = r2/sigma2+1
#    r6s = r2s*r2s*r2s
#    f = 24*e*( 2/(r6s*r6s) - 1/(r6s) )
#    f = (12*D)/sigma*((sigma/a)^7-(sigma/a)^13) 
def distances(x,y):
    while i < count:


#--------- distances ----------  

if i > sigma_stop
Force
angle_rad = math.atan(y-y,x-x)
angle_degrees = math.degrees(angle_rad)


class Particle():
    def __init__(self, x, y, vx, vy, size):
        #задание координат
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
print('Setting the grid')
#number_entered = 10 #test!
#asq3 = 1.732 #test!
#a = 1 #test!
#by = 10 #test!

asq3 = a*(3)**(1/2)
x = 0
y = 0
i = 0
j = 0
sigma = 0
#borders x - 2*sigma+(number_x_entered-1)*a 
#borders y - 2*sigma+number_y_entered*asq3 
#не совсем понял как работают массивы. поэтому разделю координаты в обоих циклах по разным массивам
#а затем, объединю их в один координатный массив, чтобы избедать перезаписи и конфликтов
while i in range(number_y_entered): #odd
   y_odd = sigma+i*asq3
   for j in range(number_x_entered):
       x_odd = sigma+j*a
   vx, vy = 0., 0.
   particle = Particle((x, y),(vx,vy), 10) #red to end

while i in range(number_y_entered): #even
   y_even = sigma+asq3/2+y*asq3
   for j in range(number_x_entered-1):
       x_even = asq3/2+i*a
   vx, vy = 0., 0.
   particle = Particle((x, y),(vx,vy), 10)



#--------- acceleration ----------




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