import matplotlib.pyplot as plt
import numpy as np
import math as math



u = np.loadtxt('trial1.txt',skiprows=1 , unpack=True)
u2 = np.loadtxt('trial2.txt',skiprows=1 , unpack=True)
u3 = np.loadtxt('trial3.txt',skiprows=1 , unpack=True)

radius = []
r1 = 0.0
for i in range(len(u)):
	r2 = r1+ i
	radius.append(r2)

# create a figure window
fig1 = plt.figure(1, figsize=(10,6))
# add extra white space between subplots
# (needed for axes labels)
fig1.subplots_adjust(wspace=0.3, hspace=0.3)
ax1 = fig1.add_subplot(1,1,1)


ax1.plot(radius, u, '--b')
ax1.plot(radius, u2, '--g')
ax1.plot(radius, u3, '--r')



##ax1.legend(loc = 'best')
ax1.grid()
ax1.set_xlabel('radius')
ax1.set_ylabel('X')
ax1.set_title('Functions')



plt.savefig('trial')


plt.show() 
