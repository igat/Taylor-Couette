import matplotlib.pyplot as plt
import numpy as np
import math as math



##u = np.loadtxt('trial1.txt',skiprows=1 , unpack=True)
##u2 = np.loadtxt('trial2.txt',skiprows=1 , unpack=True)
uphi_shooting = np.loadtxt('Uphi2.txt',skiprows=0 , unpack=True)
uphi_sparse = np.loadtxt('uphi_new.txt',skiprows=0 , unpack=True)
pressure = np.loadtxt('pressure.txt',skiprows=0 , unpack=True)
pressure2 = np.loadtxt('pressured1.txt',skiprows=0 , unpack=True)


radius = []
radius2 = []
x = 1.0-0.01
radius2.append(x)
r1 = 1.0
for i in range(len(uphi_shooting)):
	r2 = r1+ (i*0.01)
	radius.append(r2)
	radius2.append(r2)
##y = 2.0+0.01	
##radius2.append(y)

# create a figure window
fig1 = plt.figure(1, figsize=(10,8))
# add extra white space between subplots
# (needed for axes labels)
fig1.subplots_adjust(wspace=0.3, hspace=0.3)
ax1 = fig1.add_subplot(2,1,1)


##ax1.plot(radius, u, '--b')
##ax1.plot(radius, u2, '--g')
ax1.plot(radius, uphi_shooting, '--r', label = 'shooting')
ax1.plot(radius, uphi_sparse, '--b', label = 'sparse')



ax1.legend(loc = 'best')
ax1.set_xlim(0.8, 2.2)
ax1.grid()
ax1.set_xlabel('radius')
ax1.set_ylabel(r'$u_{\phi}(r)$')
ax1.set_title(r'$u_{\phi}(r)$ vs radius')


ax2 = fig1.add_subplot(2,1,2)
ax2.set_xlim(0.8, 2.2)
ax2.plot(radius2, pressure, 'b', label = '2nd deriv')
ax2.plot(radius2, pressure2, 'g', label = 'first deriv')
ax2.legend(loc = 'best')
ax2.grid()
ax2.set_xlabel('radius')
ax2.set_ylabel('Pressure')
ax2.set_title('Pressure versus radius')

##plt.savefig('shooting_and_sparce')
plt.savefig('pressure_and_velocity')


plt.show() 
