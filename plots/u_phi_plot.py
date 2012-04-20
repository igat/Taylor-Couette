import matplotlib.pyplot as plt
import numpy as np
import math as math



radius256, u256 = np.loadtxt('Uphi256.txt',skiprows=0 , unpack=True)
radius2256, u2256 = np.loadtxt('u_phi256.txt',skiprows=0 , unpack=True)
radius128, u128 = np.loadtxt('Uphi128.txt',skiprows=0 , unpack=True)
radius2128, u2128 = np.loadtxt('u_phi128.txt',skiprows=0 , unpack=True)
radius64, u64 = np.loadtxt('Uphi64.txt',skiprows=0 , unpack=True)
radius264, u264 = np.loadtxt('u_phi64.txt',skiprows=0 , unpack=True)
radius32, u32 = np.loadtxt('Uphi32.txt',skiprows=0 , unpack=True)
radius232, u232 = np.loadtxt('u_phi32.txt',skiprows=0 , unpack=True)

# create a figure window
fig1 = plt.figure(1, figsize=(10,6))
# add extra white space between subplots
# (needed for axes labels)
fig1.subplots_adjust(wspace=0.3, hspace=0.3)
ax1 = fig1.add_subplot(1,1,1)
'''ax1.plot(radius32, u32, '-.b', label = '1d-32')
ax1.plot(radius232, u232, '--b', label = '2d slice-32')
ax1.plot(radius64, u64, '-.g', label = '1d-64')
ax1.plot(radius264, u264, '--g', label = '2d slice-64')
ax1.plot(radius128, u128, '-.r', label = '1d-128')
ax1.plot(radius2128, u2128, '--r', label = '2d slice-128')
ax1.plot(radius256, u256, '-.c', label = '1d-256')
ax1.plot(radius2256, u2256, '--c', label = '2d slice-256')'''

ax1.plot(radius32, u32, '--b', label = '1d-32')
ax1.plot(radius232, u232, 'b', label = '2d slice-32')
ax1.plot(radius64, u64, '--g', label = '1d-64')
ax1.plot(radius264, u264, 'g', label = '2d slice-64')
ax1.plot(radius128, u128, '--r', label = '1d-128')
ax1.plot(radius2128, u2128, 'r', label = '2d slice-128')
ax1.plot(radius256, u256, '--k', label = '1d-256')
ax1.plot(radius2256, u2256, 'k', label = '2d slice-256')



ax1.legend(loc = 'best')
ax1.grid()
ax1.set_xlabel('radius')
ax1.set_ylabel(r'velocity, $U_{\phi}$')
ax1.set_title(r'$U_{\phi}$ vs $r$')



plt.savefig('u_phi_all')


plt.show() 
