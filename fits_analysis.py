import numpy as np
from pylab import *
import os
import sys
from astropy.wcs import WCS
import time
from astropy.io import fits
import psfex


# Vemos los fits de chequeo:
# chi.fits, proto.fits, samp.fits, resi.fits, snap.fits

path='/home/martin/Documentos/Doctorado/Lentes/PSFEx/DES_LN/'
check_path = path + 'check_files/'

chi_f = check_path+'chi_run2-r1.fits'
proto_f = check_path+'proto_run2-r1.fits'
samp_f = check_path+'samp_run2-r1.fits'
resi_f = check_path+'resi_run2-r1.fits'
snap_f = check_path+'snap_run2-r1.fits'
psf_file = path+'psfex_files/run2-r1.psf'
VIG=35

plt.ion()


def extract_vignettes_means(data, mask=True, plot=False):
	x_split = data.shape[0]/VIG
	y_split = data.shape[1]/VIG
	print 'Vignetas: ', x_split*y_split
	arr = np.split(data, x_split, axis=0)
	arr = [np.split(M, y_split, axis=1) for M in arr]
	means=[]
	if mask:
		for rows in arr:
			for squares in rows:
				if squares.mean()==0:
					#means.append(0)
					continue
				mask = squares>0
				means.append(squares[mask].mean())
				if plot:
					plt.figure()
					plt.title('mean: '+str(squares[mask].mean()))
					plt.imshow(squares, interpolation='none')
					plt.show()
	else:
		for rows in arr:
			for squares in rows:
				means.append(squares.mean())
				if plot:
					plt.figure()
					plt.title('mean: '+str(squares.mean()))
					plt.imshow(squares, interpolation='none')
					plt.show()				
		#[[means.append(squares.mean()) for squares in rows] for rows in arr]		

	return np.array(means)

# CHI.FITS -----------------------------------------------------
hdul = fits.open(chi_f)
chi_data = hdul[0].data

chi_means = extract_vignettes_means(chi_data, mask=True, plot=False)

plt.figure()
plt.axvline(chi_means.mean(), color='r', linestyle='-', label='Mean: '+str(round(chi_means.mean(),3)))
plt.hist(chi_means,30, color='k',histtype='step', label='Std: '+str(round(chi_means.std(),3)))
plt.title('Chi^2')
plt.ylabel('N')
plt.xlabel('chi')
plt.legend()
plt.show()
del hdul

# RESI.FITS -----------------------------------------------------
hdul = fits.open(resi_f)
resi_data = hdul[0].data

resi_means = extract_vignettes_means(resi_data, mask=False, plot=False)

plt.figure()
plt.axvline(resi_means.mean(), color='r', linestyle='-', label='Mean: '+str(round(resi_means.mean(),3)))
plt.hist(resi_means,30, color='k',histtype='step', label='Std: '+str(round(resi_means.std(),3)))
plt.title('Resi')
plt.ylabel('N')
plt.xlabel('resi')
plt.legend()
plt.show()

# PROTO.FITS (componentes de la PSF)-----------------------------------------------------
hdul = fits.open(proto_f)
proto_data = hdul[0].data
PSF_SIZE=25
x_split = proto_data.shape[0]/PSF_SIZE
y_split = proto_data.shape[1]/PSF_SIZE

arr = np.split(proto_data, x_split, axis=0)
arr = [np.split(M, y_split, axis=1) for M in arr]

plt.figure(figsize=(16,3))
cont=1
PSFbasis = []
for i in range(x_split):
	for j in range(y_split):
		plt.subplot(x_split,y_split,cont)
		plt.imshow(arr[i][j], interpolation='bilinear', origin='low')
		PSFbasis.append(arr[i][j])
		cont+=1
plt.suptitle('PSF components', fontsize=16)
plt.tight_layout()
plt.show()


plt.figure()
psf1 = psfex.PSF(psf_file)
psfmodel = psf1.get_psf(col=1000,row=1000)
plt.imshow(psfmodel, origin='low')
plt.show()



plt.ioff()
