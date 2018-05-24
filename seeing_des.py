from pylab import *
import os
#import pyfits
import sys
#from pyraf import iraf
from config_file import *
	
	
def seeing(imagen,pixsize,corrida,filtro,magmax,magmin,fwhmmax,plot):


	#print '----------------------------------------------------------------'
	#print '                RUNNING SExtractor (pre-PSFEx)                  '
	#print '================================================================'

	sexfile, salida_sex = sex_config_file('first', filtro, corrida, pixsize, seeing=1., satur=50000)
	
	callsex='sextractor '+imagen+' -c '+sexfile+' > sex_output'		
	os.system(callsex)
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
		plt.plot(FWHM,MAGBEST, 'b.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		plt.show()
		plt.plot(FLUXRADIUS,MAGBEST, 'm.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FLUX_RADIUS')
		plt.show()

	#make the cuts in magnitude to get mostly the stars
	mask2=(MAGBEST < magmax)*(MAGBEST > magmin)*(FWHM < fwhmmax)
	fwhm=FWHM[mask2]
	mag=MAGBEST[mask2]
	j=len(fwhm)

	#and get the maximun value of the FWHM distribution	
	bin_height, bin_edges = np.histogram(fwhm, bins=np.arange(fwhm.min(),fwhm.max(),0.05))
	ind=np.argmax(bin_height)
	moda=(bin_edges[ind]+bin_edges[ind+1])/2.
			
	if plot in ('s', 'S', 'si', 'Si', 'SI'):	
		plt.plot(FWHM,MAGBEST,'k.')
		plt.plot(fwhm,mag,'r.')
		plt.show()
		print 'SEEING in pix',moda
		plt.hist(fwhm,15)
		plt.show()
	
	#print 'SEEING in arcsec',moda*pixsize
	seeing=moda*pixsize
	print ' '
	print ' ------------ seeing: ',seeing
	print ' ------------ imagen: ',imagen
	print ' '

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.plot(FWHM,MAGBEST, 'b.',moda,18.,'ro')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		#~ plt.axis([0,30,16,26])
		plt.show()
	#---------------------------------------------------------------------------------------------	
	
	
	return seeing, SATUR
	
