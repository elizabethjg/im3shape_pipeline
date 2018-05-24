import os
import sys
from pylab import *


def config_file(param):

	#print '----------------------------------------------------------'
	#print '           WRITING SExtractor CONFIGURATION FILE          '
	#print '=========================================================='
	
	#param = {'run': 'first', 'sexfile': 'nombre', 'seeing': 0.9, 'pixscale': pixscale, 'satu': satur}
	
	sexrun 	= param['run'] # 'first', 'second', 'third'
	sexfile = param['sexfile']	# file name
	seeing 	= str(param['seeing'])
	pixscale = str(param['pixscale'])
	satur 	= str(param['satur'])
	
	if sexrun == 'first':
		output_format = 'ASCII_HEAD'
		thresh = '5'		
	elif sexrun == 'second':
		output_format = 'FITS_LDAC'
		thresh = '5'		
	elif sexrun == 'third':
		output_format = 'ASCII_HEAD'
		thresh = '1.5'

	param_file = sexrun+'.param'

	os.system('rm '+sexfile)
	os.system('rm '+salida_sex)	
	f1=open(sexfile,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+salida_sex+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     '+output_format+'     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  '+param_file+'     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE		CCD			# "CCD" or "PHOTO" (*)\n')
	f1.write('FLAG_IMAGE		flag.fits		# filename for an input FLAG-image\n')
	f1.write('DETECT_MINAREA	5			# minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH		'+thresh+' 			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH	2.			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER		Y			# apply filter for detection ("Y" or "N")?\n')
	f1.write('FILTER_NAME	  	default.conv	# name of the file containing the filter\n')
	f1.write('DEBLEND_NTHRESH	32			# Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT	0.005			# Minimum contrast parameter for deblending\n')
	f1.write('CLEAN			Y			# Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM		1.0			# Cleaning efficiency\n')
	f1.write('MASK_TYPE		CORRECT		# type of detection MASKing: can be one of\n')
	f1.write('					# "NONE", "BLANK" or "CORRECT"\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES	10			# MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS	2.5, 3.5		# MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('SATUR_LEVEL		'+satur+'		# level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT		26.73			# magnitude zero-point\n')
	f1.write('MAG_GAMMA		4.0			# gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN			2.00			# detector gain in e-/ADU.\n')
	f1.write('PIXEL_SCALE		'+pixsize+'	# size of pixel in arcsec (0=use FITS WCS info).\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM		'+seeing+'			# stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME	default.nnw		# Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE		64			# Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE	3			# Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)\n')
	f1.write('BACKPHOTO_THICK	24			# thickness of the background LOCAL annulus (*)\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE	NONE			# can be one of "NONE", "BACKGROUND",\n')
	f1.write('						# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",\n')
	f1.write('						# "-OBJECTS", "SEGMENTATION", "APERTURES",\n')
	f1.write('						# or "FILTERED" (*)\n')
	f1.write('CHECKIMAGE_NAME	apertures.fits	# Filename for the check-image (*)\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK	8000			# number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK	400000		# number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE	1024			# number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)\n')
	f1.close()

def sex_run(image, sexrun):
	
	
	sexfile='first'+filtro+str(corrida)+'.sex' #name of configuration file of SExtractor
	salida_sex='run1-'+filtro+str(corrida)+'.cat' #exit catalogue from SExtractor

	#callsex='sextractor '+imagen+' -c '+sexfile+' > sex_output'
