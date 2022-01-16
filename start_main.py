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
number_entered = readline(2) 
sigma = readline(4) 
a = readline(6)
bx = readline(8) 
by = readline(10) 
e = readline(12) 
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

def LJ_force(p1,p2):
    rx = p1.x - p2.x
    ry = p1.y - p2.y
    r2 = rx*rx+ry*ry

    f = -(24)*(2*(sigma/a)^14-(sigma/a)^8) # interaction force
    p1.ax += f*(rx/r2)
    p1.ay += f*(ry/r2)
    p2.ax -= f*(rx/r2)
    p2.ay -= f*(ry/r2)

#    r2s = r2/sigma2+1
#    r6s = r2s*r2s*r2s
#    f = 24*e*( 2/(r6s*r6s) - 1/(r6s) )
#    f = (12*D)/sigma*((sigma/a)^7-(sigma/a)^13) 

def Verlet_step(particles, h):
    for p in particles:
        p.verlet1_update_vx(h); 
        p.bounce(self)
    #t += h;
    for i, p1 in enumerate(particles):
        for p2 in particles[i+1:]:
            LJ_force(p1,p2)
    for p in particles:
        p.verlet2_update_v(h)

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

    def verlet1_update_vx(self,h):
        self.vx += self.ax*h/2
        self.vy += self.ay*h/2
        self.x += self.vx*h
        self.y += self.vy*h
        self.ax = 0
        self.ay = 0

    def verlet2_update_v(self,h):
        self.vx += self.ax*h/2
        self.vy += self.ay*h/2

    def display(self,screen, aa):
        pygame.draw.circle(screen, self.colour, (int(aa*self.x+0.5), int(aa*self.y+0.5)), aa*self.size, aa*self.thickness)

    def bounce(self):
        if self.x > width - self.size:
            self.x = 2*(width - self.size) - self.x
            self.vx = - self.vx

        elif self.x < self.size:
            self.x = 2*self.size - self.x
            self.vx = - self.vx

        if self.y > height - self.size:
            self.y = 2*(height - self.size) - self.y
            self.vy = - self.vy

        elif self.y < self.size:
            self.y = 2*self.size - self.y
            self.vy = - self.vy            

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

while y < by: #odd
   y += y*asq3
   for j in range(number_entered):
       x += j*a
   vx, vy = 0., 0.
   particle = Particle((x, y),(vx,vy), 10)

while y < by: #even
   y += asq3/2+y*asq3
   for j in range(number_entered-1):
       x += asq3/2+i*a
   vx, vy = 0., 0.
   particle = Particle((x, y),(vx,vy), 10)


#for n in range(number_of_particles):
#    x = 1.0*random.randint(15, width-15)
#    y = 1.0*random.randint(15, height-15)
#    vx, vy = 0., 0.
#    for k in range(6):
#        vx += random.randint(-10, 10)/2.
#        vy += random.randint(-10, 10)/2.
#
#    particle = Particle((x, y),(vx,vy), 10)
#
#    my_particles.append(particle)

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