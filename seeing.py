import os
#import pyfits
import sys
#from pyraf import iraf
from pipeio import *
	
	
def seeing_func(imagen,pixsize,zeropoint,gain,corrida,filtro,magmax,magmin,fwhmmax,plot):

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		from pylab import *

	#print '----------------------------------------------------------------'
	#print '                RUNNING SExtractor (pre-PSFEx)                  '
	#print '================================================================'

	sexfile, salida_sex = sex_config_file('first', filtro, corrida, pixsize, zeropoint, gain, seeing=1., satur=50000)
	
	call_sex = imagen+' -c '+sexfile+' > sex_output'		
	out_sex_code = os.system('sextractor '+call_sex)
	if out_sex_code != 0: os.system('sex '+call_sex)
	#-------------------------------------------------------------------------

	cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor

	#compute the SATURATION LEVEL and SEEING
	FLUXMAX = cat[:,8]
	SATUR=0.8*FLUXMAX.max()
		
	#now, starts to compute the seeing
	FWHM = cat[:,5]
	FLAG = cat[:,11]
	MAGBEST = cat[:,6]
	FLUXRADIUS = cat[:,7]
		
	if plot in ('s', 'S', 'si', 'Si', 'SI'):

		plt.figure()	
		plt.plot(MAGBEST,FLUXMAX, 'b.')
		plt.axhline(SATUR)
		plt.ylabel('FLUXMAX')
		plt.xlabel('MAG_BEST')
		plt.show()
		plt.figure()	
		plt.plot(FWHM,MAGBEST, 'b.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		plt.show()
		plt.figure()	
		plt.plot(FLUXRADIUS,MAGBEST, 'm.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FLUX_RADIUS')
		plt.show()

	#make the cuts in magnitude to get mostly the stars
	mask2=(MAGBEST < magmax)*(MAGBEST > magmin)*(FLUXRADIUS < fwhmmax)
	fwhm=FWHM[mask2]
	mag=MAGBEST[mask2]
	fluxradius=FLUXRADIUS[mask2]
	j=len(fluxradius)

	#and get the maximun value of the FWHM distribution	
	bin_height, bin_edges = np.histogram(fwhm, bins=np.arange(fwhm.min(),fwhm.max(),0.05))
	ind=np.argmax(bin_height)
	moda=(bin_edges[ind]+bin_edges[ind+1])/2.
			
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.figure()	
		plt.plot(FLUXRADIUS,MAGBEST,'k.')
		plt.plot(fluxradius,mag,'r.')
		plt.show()
		print 'Seeing in pix',moda
		plt.figure()	
		plt.hist(fwhm,15)
		plt.show()
	
	#print 'SEEING in arcsec',moda*pixsize
	
	seeing = moda*pixsize
	print ' '
	print ' ------------ seeing: ', seeing, ' arcsec'
	print ' ------------ imagen: ', imagen
	print ' '

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.figure()
		plt.plot(FWHM,MAGBEST, 'b.',moda,18.,'ro')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		#~ plt.axis([0,30,16,26])
		plt.show()
	#---------------------------------------------------------------------------------------------	
	
	
	return seeing, SATUR
	


