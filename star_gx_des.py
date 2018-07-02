import math
import numpy as np
import os
from astropy.io import fits
import sys
from scipy.optimize import curve_fit


def star_gx(sex_salida, fwhm, plot, PSFEx_manual):
	'''
	Clasifies sources in stars, galaxies and false detections. Also makes the PSFEx catalog.
	 input:
	 	sex_salida:			String. SExtractor output filename. The file needs to be in FITS_LDAC format
	 	fwhm:				Seeing of the image in pixel units
	 	plot:				Flag. String with 'si' to make the plots
	 	PSFEx_manual:		Flag. If True we modify SExtractor output to remove bright stars

	 output:
	 	hdul_stars:			Hdulist object with the selected stars
	 	hdul_gx:			Hdulist object with the selected galaxies
	 	salida_sex_mod:		String. Modified SExtractor output filename with '_mod' added at the end.
	 						If PSFEx_manual = False, the file is not modified and returns salida_sex
	'''
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		from pylab import *
	
	sex_salida_mod = sex_salida + '_mod'
	os.system('rm '+sex_salida_mod)

	# Lee el catalogo de sextractor -------------------------------------------------
	hdul = fits.open(sex_salida)
	data = hdul[2].data
	del hdul

	FWHM 	= data['FWHM_IMAGE']
	MUMAX 	= data['MU_MAX']
	CLASS 	= data['CLASS_STAR']
	FLAG 	= data['FLAGS']
	MAGBEST 	= data['MAG_BEST']
	MAGAUTO 	= data['MAG_AUTO']
	FLUXRADIUS 	= data['FLUX_RADIUS']

	# Ajusta la relacion lineal mu vs mag para las estrellas ------------------------
	m_candidatas = (MAGAUTO < 23.)*(FWHM>(fwhm-0.5))*(FWHM < (fwhm+0.5))

	X = MAGAUTO[m_candidatas]
	Y = MUMAX[m_candidatas]
	stars_fit = lambda x,m,n: x*m+n
	popt, pcov = curve_fit(stars_fit, X, Y)
	m,n = popt	
	varx = np.array([X.min(), X.max()])
	vary = m*varx+n

	zero=MAGAUTO.min()
	x=np.arange(15)+zero
	y=m*x+n
			
	# Clasifica los objetos ----------------------------------------------------------
	ancho=0.4  #width in magnitudes 
	mumin=11.5
	mumax=17.
	mu=m*MAGBEST+n
	mag = (MUMAX-n)
	mask_good = (MUMAX > mumin) * (FWHM > fwhm-0.8) * (FLAG < 4.0) # Quitamos falsas detecciones
	mask_stars = (MUMAX < mu+ancho) * (MUMAX < mumax) * (FWHM < fwhm+1.0) * (MUMAX > mu-ancho) * mask_good	# Seleccionamos las estrellas
	mask_gx = (~mask_stars) * (CLASS < 0.8) * ((MUMAX > mumax)+(MUMAX > mu-ancho)) * mask_good	# Seleccionamos las galaxias

	print 'cantidad de estrellas seleccionadas: ', mask_stars.sum()
	print 'cantidad de galaxias seleccionadas: ', mask_gx.sum()

	# Seleccionamos las estrellas a mano para PSFEx -----------------------------------
	if PSFEx_manual:
		mu_brightest = 18.										# Descartamos las estrellas mas brillantes
		mask_bright = MUMAX<mu_brightest
		mask_psfex = (~mask_bright)*mask_stars #+ mask_gx		# incluimos las galaxias ???

		hdul_psfex = fits.open(sex_salida)						# Abro una copia y guardo los objetos para PSFEx
		hdul_psfex[2].data = hdul_psfex[2].data[mask_psfex]
		#sex_salida_mod = './sex_files/modified_psfex.cat'
		os.system('rm '+sex_salida_mod)
		hdul_psfex.writeto(sex_salida_mod)
		del hdul_psfex
	else:
		sex_salida_mod = sex_salida		

	# Separo los hdulist de estrellas y galaxias ---------------------------------------
	hdul_gx 	= fits.open(sex_salida)			# Abro una copia y modifico
	hdul_stars 	= fits.open(sex_salida)			# Abro una copia y modifico
	hdul_gx[2].data 	= data[mask_gx]
	hdul_stars[2].data 	= data[mask_stars]

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.figure()
		plt.plot(MAGAUTO,MUMAX,'k.')
		plt.plot(MAGAUTO[mask_stars],MUMAX[mask_stars],'b.',x,y,'r',varx,vary,'b')
		plt.plot(MAGAUTO[mask_gx],MUMAX[mask_gx],'g.')
		#plt.plot(MAGAUTO[m_candidatas],MUMAX[m_candidatas],'g.')
		if PSFEx_manual: plt.scatter(MAGAUTO[mask_psfex], MUMAX[mask_psfex], s=20,edgecolor='r',facecolor='None')
		plt.xlabel('MAG_AUTO')
		plt.ylabel('MU_MAX')
		plt.show()

		plt.figure()
		plt.plot(FLUXRADIUS, MAGAUTO,'k.',markersize=2)
		plt.plot(FLUXRADIUS[mask_stars], MAGAUTO[mask_stars],'b.', FLUXRADIUS[mask_gx], MAGAUTO[mask_gx], 'g.',markersize=2)
		#plt.plot(FLUXRADIUS[m_candidatas],MUMAX[m_candidatas],'g.')
		if PSFEx_manual: plt.scatter(FLUXRADIUS[mask_psfex], MAGAUTO[mask_psfex], s=20,edgecolor='r',facecolor='None')
		plt.xlabel('FLUX_RADIUS')
		plt.ylabel('MAG_AUTO')
		plt.show()
		
		plt.figure()
		plt.plot(FWHM, MAGAUTO,'k.')
		plt.plot(FWHM[mask_stars], MAGAUTO[mask_stars],'b.', FWHM[mask_gx], MAGAUTO[mask_gx], 'g.',markersize=2)
		#plt.plot(FWHM[m_candidatas],MAGAUTO[m_candidatas],'g.')
		if PSFEx_manual: plt.scatter(FWHM[mask_psfex], MAGAUTO[mask_psfex], s=20,edgecolor='r',facecolor='None')
		plt.xlabel('FWHM')
		plt.ylabel('MAG_AUTO')
		plt.show()


		#Ns = len(hdulist[2].data)		# esto es para ver todas las estrellas juntas, muy lento
		#xg = int(Ns**0.5)+1
		#plt.figure()
		#for i_s in xrange(Ns):
		#	plt.subplot(xg, xg, i_s+1)
		#	plt.imshow(hdulist[2].data[i_s]['VIGNET'], origin='low')
		#plt.show()


	return hdul_stars, hdul_gx, sex_salida_mod

def background_gx(hdul_gx, MAGMIN, MAGMAX, plot):
	'''
	Selects background galaxies. (Por ahora cortamos en magnitud. Podemos agregar cortes en color y z_phot)
	 input:
	 	hdul_gx:	Hdulist object with the galaxies
	 	MAGMIN:		Lower magnitud bound: mag > MAGMIN
	 	MAGMAX:		Upper magnitud bound: mag < MAGMAX
	 	plot:		Flag. String with 'si' to make the plots

	 output:
		hdul_gx:	Hdulist object with the selected background galaxies
	'''

	MAGAUTO = hdul_gx[2].data['MAG_AUTO']
	FWHM 	= hdul_gx[2].data['FWHM_IMAGE']

	mask_back = (MAGAUTO<MAGMAX) * (MAGAUTO>MAGMIN) * (FWHM>4.)

	hdul_gx[2].data = hdul_gx[2].data[mask_back]

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.figure()
		plt.plot(FWHM, MAGAUTO, 'k.', markersize=2)
		plt.scatter(FWHM[mask_back], MAGAUTO[mask_back], s=20, edgecolor='r', facecolor='None')
		plt.xlabel('FWHM')
		plt.ylabel('MAG_AUTO')
		plt.show()		

	return hdul_gx
