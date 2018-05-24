import math
import numpy as np
from pylab import *
from astropy.io import fits

hdul=fits.open('check_files/resi_modified_psfex.fits')
#hdul=fits.open('check_files/chi_run2-r0.fits')

image_chi = hdul[0].data


#samp_tot=391
nsize = 61
xnstamps = image_chi.shape[0]/nsize
ynstamps = image_chi.shape[1]/nsize

chi = np.zeros(xnstamps*ynstamps)



x = np.arange(xnstamps)
y = np.arange(ynstamps)
x = x*nsize
y = y*nsize

xy = meshgrid(x,y)
x = xy[0].flatten()
y = xy[1].flatten()



for j in range(xnstamps*ynstamps):
	#print x[j],x[j]+nsize,y[j],y[j]+nsize
	chi[j] = np.mean(image_chi[x[j]:x[j]+nsize,y[j]:y[j]+nsize])
	
mask = (chi != 0.)*~np.isnan(chi)

print 'total de estrellas ', len(chi[mask])

plt.hist(chi[mask],100)
