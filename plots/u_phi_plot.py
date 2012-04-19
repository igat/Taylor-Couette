import matplotlib.pyplot as plt
import numpy as np
import math as math



radius, u = np.loadtxt('Uphi3.txt',skiprows=0 , unpack=True)
radius2, u2 = np.loadtxt('u_phi.txt',skiprows=0 , unpack=True)

# create a figure window
fig1 = plt.figure(1, figsize=(10,6))
# add extra white space between subplots
# (needed for axes labels)
fig1.subplots_adjust(wspace=0.3, hspace=0.3)
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(radius, u, '--b', label = '1d')
ax1.plot(radius2, u2, '--g', label = '2d slice')

ax1.legend(loc = 'best')
ax1.grid()
ax1.set_xlabel('radius')
ax1.set_ylabel(r'velocity, $U_{\phi}$')
ax1.set_title(r'$U_{\phi}$ vs $r$')



plt.savefig('u_phi_graphs')


plt.show() 
