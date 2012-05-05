import matplotlib.pyplot as plt
import numpy as np
import math as math
'''
for m in range(100):
	string1 = str(m+52)
	dim = len(string1)
	zeros = 6-dim
	velx_filename = 'helmholtz_velX_per' + zeros*'0' + string1 + '.txt'
	vely_filename = 'helmholtz_velY_per' + zeros*'0' + string1+ '.txt'
	pres_filename = 'helmholtz_pressure_per' + zeros*'0' + string1+ '.txt'
	dens_filename = 'helmholtz_density_per' + zeros*'0' + string1+ '.txt'
	print velx_filename
	d = np.loadtxt(dens_filename,skiprows=0 , unpack=True)
	p = np.loadtxt(pres_filename,skiprows=0 , unpack=True)
	vx = np.loadtxt(velx_filename,skiprows=0 , unpack=True)
	vy = np.loadtxt(vely_filename,skiprows=0 , unpack=True)



	# create a figure window
	fig1 = plt.figure(1, figsize=(15,8))
	# add extra white space between subplots
	# (needed for axes labels)
	fig1.subplots_adjust(wspace=0.1, hspace=0.3)
	ax1 = fig1.add_subplot(2,2,1)
	im = ax1.imshow(d, origin = 'lower')
	ax1.set_title("Density")
	plt.colorbar(im, shrink = 0.8, extend = 'both')

	ax2 = fig1.add_subplot(2, 2, 2)
	im2 = ax2.imshow(p, origin = 'lower')
	plt.colorbar(im2, shrink = 0.8, extend = 'both')
	ax2.set_title("Pressure")

	ax3 = fig1.add_subplot(2, 2, 3)
	im3 = ax3.imshow(vx, origin = 'lower')
	plt.colorbar(im3, shrink = 0.8, extend = 'both')
	ax3.set_title("Velocity in x")

	ax4 = fig1.add_subplot(2, 2, 4)
	im4 = ax4.imshow(vy, origin = 'lower')
	plt.colorbar(im4, shrink = 0.8, extend = 'both')
	ax4.set_title("Velocity in y")

	savename = 'helmholtz_contour_periodic_2test' + zeros*'0' + string1 + '.png'
	plt.savefig(savename)
	plt.close()

##plt.show()
'''
p = np.loadtxt('Pressure_timedep.txt',skiprows=0 , unpack=True)
vr = np.loadtxt('ur_timedep.txt',skiprows=0 , unpack=True)
vphi = np.loadtxt('uphi_timedep.txt',skiprows=0 , unpack=True)



# create a figure window
fig1 = plt.figure(1, figsize=(15,8))
# add extra white space between subplots
# (needed for axes labels)
fig1.subplots_adjust(wspace=0.1, hspace=0.3)
ax1 = fig1.add_subplot(2,2,1)
im = ax1.imshow(p, origin = 'lower')
ax1.set_title("Pressure")
plt.colorbar(im, shrink = 0.8, extend = 'both')

ax2 = fig1.add_subplot(2, 2, 2)
im2 = ax2.imshow(vr, origin = 'lower')
plt.colorbar(im2, shrink = 0.8, extend = 'both')
ax2.set_title("Vr")

ax3 = fig1.add_subplot(2, 2, 3)
im3 = ax3.imshow(vphi, origin = 'lower')
plt.colorbar(im3, shrink = 0.8, extend = 'both')
ax3.set_title("Vphi")
'''
ax4 = fig1.add_subplot(2, 2, 4)
im4 = ax4.imshow(vy, origin = 'lower')
plt.colorbar(im4, shrink = 0.8, extend = 'both')
ax4.set_title("Velocity in y")

'''
##ax1.legend(loc = 'best')
##ax1.axhline(color='gray')
##ax1.axvline(color='gray')
##ax1.grid()

plt.savefig("test")

plt.show()

