import os
import sys
import numpy as np
from astropy.io import fits
from seeing import *
from pipeio import *


def correlate_magnitudes(sex_salida,RA0,DEC0):
	
	'''
	Uses KDTree to corralate the magnitudes of the output SExtractor catalog
	
	Input:
	
	sex_salida =  (str) Name of SExtractor output
	RA0        =  (float array) RA coordinate to be correlate
	DEC0       =  (float array) DEC coordinate to be correlate
	
	Output
	MAG        =  (float array) Array of correlated magnitudes
	'''
	
	
	MAG = np.zeros(len(RA0))
	MAG.fill(-99.)

	sexcat = np.loadtxt(sex_salida)
	
	ra = sexcat[:,3]
	dec = sexcat[:,4]
	mag = sexcat[:,5]


	tree=spatial.KDTree(np.array([RA0,DEC0]).T)
	dist,ind=tree.query(np.array([ra,dec]).T)
	mdist=dist<1./3600.

	MAG[ind[mdist]]=mag[mdist]
	
	return MAG


def runsex_other_images(image,hdul_gx,filtros,images,simult,pixsize,corrida):
	
	
	'''
	Run SExtractor in all the image filters and correlate magnitudes with the galaxy catalog
	
	Input:
	image      =   (str) Main image where the source classification is performed
	hdul_gx    =   (Hdulist object) Catalog of sources classified as galaxies
	filtros    =   (str array) Filter of the main image and where the filters of other images are going
	               to be recorded
	images     =   (str array) array containing filter name, image name, and zero point of 
	                all the available images but the main.
	simult     =   (str) Parameter that indicates if SExtractor is going to be run using the main image
	               to do the detection ('si')
	pixsize    =   (float) Size of the pixel in arcsec
	corrida    =   (int) Core number
	
	Output:
	MAGNITUDES =   (float array) Catalog of magnitudes corralted with hdul_gx
	FILTROS    =   (str array) All the filters used including the main image
	'''
	
	nimages=images.shape[0]/3
	MAGNITUDES=np.zeros((len(hdul_gx[2].data),nimages))

	RA0=hdul_gx[2].data['ALPHA_J2000']
	DEC0=hdul_gx[2].data['DELTA_J2000']



	for j in range(nimages):
		FILTRO,IMAGE,ZERO=images[3*j],images[3*j+1],images[3*j+2]
		filtros=np.append(filtros,FILTRO)
		hdu = fits.open(IMAGE)
		GAIN = str(hdu[0].header['GAIN'])
		
		print 'CORRIENDO SExtractor en el filtro ',FILTRO, IMAGE
		
		print seeing_func
		
		SEEING, SATUR=seeing_func(IMAGE,pixsize,ZERO,GAIN,corrida,FILTRO,23.,12.,3.,plot)
		print 'SEEING ',SEEING
		sex_conf, sex_salida = sex_config_file('third', FILTRO, corrida, pixsize, ZERO, GAIN, SEEING, SATUR)
		if simult in ('s', 'S', 'si', 'Si', 'SI'):
			call_sex = 'sextractor '+image+','+IMAGE+' -c '+sex_conf#+' > sex_output'
			os.system(call_sex)
			MAGNITUDES[:,j] = correlate_magnitudes(sex_salida,RA0,DEC0)
		else:
			call_sex = 'sextractor '+IMAGE+' -c '+sex_conf#+' > sex_output'
			os.system(call_sex)
			MAGNITUDES[:,j] = correlate_magnitudes(sex_salida,RA0,DEC0)
			
	
	return MAGNITUDES,filtros
