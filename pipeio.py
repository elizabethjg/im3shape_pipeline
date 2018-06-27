import os
import sys
import numpy as np
from astropy.io import fits
import collections
from scipy import spatial


# Directories where the outputs are going to be written. If they dont exists they are created
hostname = os.uname()[1]


if 'clemente' in hostname: 
	sex_conf_file_path = '/opt/external/gcc/7.2.0/sextractor/2.19.5/share/sextractor/'                #clemente
else:
	sex_conf_file_path = '/usr/share/sextractor/'

# Directories where the outputs are going to be written. If they dont exists they are created
sex_path           = './sex_files/'
psfex_path           = './psfex_files/'
check_path         = './check_files/'
im3_path           = './im3_files/'



if not os.path.isdir(sex_path): os.makedirs(sex_path)
if not os.path.isdir(psfex_path): os.makedirs(psfex_path)
if not os.path.isdir(check_path): os.makedirs(check_path)
if not os.path.isdir(im3_path): os.makedirs(im3_path)

class clr:
	'''
	Shell output colors. The ENDC sets the shell back to white. Use: print clr.WARNING+'warning text'+clr.ENDC. 
	'''
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

def sex_config_file(sexrun, filtro, corrida, pixsize, zeropoint, gain, seeing, satur):
	'''
	SExtractor configuration file
	Input:
		sexrun:		(str) 'first', 'second', 'third'
		filtro: 	(str) 'r','i', ...
		corrida: 	(int) core number
		pixsize: 	(flt) pixel size in arcsec
		seeing: 	(flt) seeing in arcsec
		satur: 		(flt) saturation level in ADUs
	Outpu:
		sex_conf: 	(str) config file name
		sex_salida:	(str) output file of sextractor
	'''
	
	Seeing_file = np.round(seeing/pixsize)
	
	corrida   = str(corrida)
	seeing 	  = str(seeing)
	pixsize   = str(pixsize)
	satur 	  = str(satur)
	gain      = str(gain)
	zeropoint = str(zeropoint)

	

	if Seeing_file < 3.:
		filter_file = 'gauss_2.0_5x5.conv'
	elif Seeing_file == 3.:
		filter_file = 'gauss_3.0_7x7.conv'
	elif Seeing_file == 4.:
		filter_file = 'gauss_4.0_7x7.conv'
	else:
		filter_file = 'gauss_5.0_9x9.conv'


	sex_conf 	= sex_path + sexrun+filtro+corrida+'.sex' #name of configuration file of SExtractor
	param_file	= sexrun+'.param'	# parameters to compute
	
	if sexrun == 'first':
		sex_salida		= sex_path + 'run1-'+filtro+corrida+'.cat'		
		output_format	= 'ASCII_HEAD'
		thresh			= '5'
		check_type      = 'NONE'
	elif sexrun == 'second':
		sex_salida		= sex_path + 'run2-'+filtro+corrida+'.cat'
		output_format	= 'FITS_LDAC'
		thresh			= '1.5'	
		check_type      = 'SEGMENTATION'
	elif sexrun == 'third':
		sex_salida		= sex_path + 'run2-'+filtro+corrida+'.cat'
		output_format	= 'ASCII_HEAD'
		thresh			= '1.5'	
		check_type      = 'NONE'


	os.system('rm '+sex_conf)
	#os.system('rm '+sex_salida)	
	f1=open(sex_conf,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+sex_salida+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     '+output_format+'     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  '+param_file+'     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE		CCD			# "CCD" or "PHOTO" (*)\n')
	f1.write('FLAG_IMAGE		flag.fits		# filename for an input FLAG-image\n')
	f1.write('DETECT_MINAREA	1.5			# minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH		'+thresh+' 			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH	1.5			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER		Y			# apply filter for detection ("Y" or "N")?\n')
	f1.write('FILTER_NAME	  	'+sex_conf_file_path+filter_file+'	# name of the file containing the filter\n') #gauss_3.0_5x5.conv
	f1.write('DEBLEND_NTHRESH	32			# Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT	0.005			# Minimum contrast parameter for deblending\n')
	f1.write('CLEAN			Y			# Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM		1.0			# Cleaning efficiency\n')
	f1.write('MASK_TYPE		CORRECT		# type of detection MASKing: can be one of\n')
	f1.write('					# "NONE", "BLANK" or "CORRECT"\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES	5			# MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS	2.5, 3.5		# MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>\n')
	f1.write('SATUR_LEVEL		'+satur+'		# level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT		'+zeropoint+'			# magnitude zero-point\n') #26.73
	f1.write('MAG_GAMMA		4.0			# gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN			'+gain+'			# detector gain in e-/ADU.\n')
	f1.write('PIXEL_SCALE		'+pixsize+'	# size of pixel in arcsec (0=use FITS WCS info).\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM		'+seeing+'			# stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME		'+sex_conf_file_path+'default.nnw	# Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE		64			# Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE	3			# Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)\n')
	f1.write('BACKPHOTO_THICK	24			# thickness of the background LOCAL annulus (*)\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE	'+check_type+'			# can be one of "NONE", "BACKGROUND",\n')
	f1.write('						# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",\n')
	f1.write('						# "-OBJECTS", "SEGMENTATION", "APERTURES",\n')
	f1.write('						# or "FILTERED" (*)\n')
	f1.write('CHECKIMAGE_NAME	'+sex_path+'seg_'+corrida+'.fits	# Filename for the check-image (*)\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK	8000			# number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK	400000		# number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE	1024			# number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)\n')
	f1.close()

	return sex_conf, sex_salida



def psfex_config_file(sex_salida, filtro, corrida):
	'''
	PSFEx configuration file
	Input:
		sex_salida: (str) sextractor output catalog name.
		filtro: 	(str) 'r','i', ...
		corrida: 	(int) core number
	Outpu:
		psfex_conf: 	(str) config file name
		psfex_salida:	(str) output file of psfex

	Notes: PSFEx doesn't let the user specify a name for the PSF file. It uses the same name as sextractor
	output catalog but with the '.psf' extension. The user can only specify the output directory path.
	'''	
	
	corrida = str(corrida)

	psfex_conf	 = psfex_path + filtro+corrida+'.psfex'
	sex_name     = sex_salida[sex_salida.rfind('/')+1 : sex_salida.rfind('.')]		# Removes path and extension
	psfex_salida = psfex_path + sex_name+'.psf'
	
	f1=open(psfex_conf,'w')
	f1.write('#-------------------------------- PSF model ----------------------------------\n')
	f1.write(' \n')
	f1.write('BASIS_TYPE      PIXEL_AUTO      # NONE, PIXEL, GAUSS-LAGUERRE or FILE\n')
	f1.write('BASIS_NUMBER    20              # Basis number or parameter\n')
	f1.write('PSF_SAMPLING    1.0             # Sampling step in pixel units (0.0 = auto)\n')
	f1.write('PSF_ACCURACY    0.01            # Accuracy to expect from PSF "pixel" values\n')
	f1.write('PSF_SIZE        120,120           # Image size of the PSF model\n')
	f1.write(' \n')
	f1.write('#------------------------- Point source measurements -------------------------\n')
	f1.write(' \n')
	f1.write('CENTER_KEYS     X_IMAGE,Y_IMAGE # Catalogue parameters for source pre-centering\n')
	f1.write('PHOTFLUX_KEY    FLUX_APER(1)    # Catalogue parameter for photometric norm.\n')
	f1.write('PHOTFLUXERR_KEY FLUXERR_APER(1) # Catalogue parameter for photometric error\n')
	f1.write(' \n')
	f1.write('#----------------------------- PSF variability -------------------------------\n')
	f1.write(' \n')
	f1.write('PSFVAR_KEYS     X_IMAGE,Y_IMAGE # Catalogue or FITS (preceded by :) params\n')
	f1.write('PSFVAR_GROUPS   1,1             # Group tag for each context key\n')
	f1.write('PSFVAR_DEGREES  2               # Polynom degree for each group\n')
	f1.write(' \n')
	f1.write('#----------------------------- Sample selection ------------------------------\n')
	f1.write(' \n')
	f1.write('SAMPLE_AUTOSELECT  Y            # Automatically select the FWHM (Y/N) ?\n')
	f1.write('SAMPLEVAR_TYPE     SEEING       # File-to-file PSF variability: NONE or SEEING\n')
	f1.write('SAMPLE_FWHMRANGE   2, 10     # Allowed FWHM range\n')	#2,10
	f1.write('SAMPLE_VARIABILITY 0.5          # 0.5 Allowed FWHM variability (1.0 = 100%)\n')
	f1.write('SAMPLE_MINSN       20.           # 20. Minimum S/N for a source to be used\n')
	f1.write('SAMPLE_MAXELLIP    0.2          #0.2 Maximum (A-B)/(A+B) for a source to be used\n')
	f1.write(' \n')
	f1.write('#------------------------------- Check-plots ---------------------------------- \n')
	f1.write('CHECKPLOT_ANTIALIAS   Y \n')
	f1.write('  \n')
	f1.write('CHECKPLOT_DEV       SVG         # NULL, XWIN, TK, PS, PSC, XFIG, PNG, \n')
	f1.write('								# JPEG, AQT, PDF or SVG \n')
	f1.write('CHECKPLOT_TYPE      FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS #FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS \n')
	f1.write('								# or NONE \n')
	f1.write('CHECKPLOT_NAME      fwhm, ellipticity, counts, countfrac, chi2, resi \n')
	f1.write('  \n')
	f1.write('#------------------------------ Check-Images --------------------------------- \n')
	f1.write('  \n')
	f1.write('CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS \n')
	f1.write('								# or MOFFAT,-MOFFAT,-SYMMETRICAL \n')
	f1.write('CHECKIMAGE_NAME '+check_path+'chi.fits, '+check_path+'proto.fits, ' \
			 +check_path+'samp.fits, '+check_path+'resi.fits, '+check_path+'snap.fits \n')
	f1.write('								# Check-image filenames \n')
	f1.write('  \n')
	f1.write('#----------------------------- Miscellaneous --------------------------------- \n')
	f1.write('  \n')
	f1.write('PSF_DIR       '+psfex_path+'    # Where to write PSFs (empty=same as input) \n')
	f1.write('VERBOSE_TYPE    NORMAL          # can be QUIET,NORMAL,LOG or FULL \n')
	f1.write('WRITE_XML       N               # Write XML file (Y/N)? \n')
	f1.write('XML_NAME        psfex.xml       # Filename for XML output \n')
	f1.write('XSL_URL       http://www.w3.org/1999/XSL/Transform       # XSLT translation url \n')
	f1.write('NTHREADS        1               # 1 single thread \n')
	f1.close()

	return psfex_conf, psfex_salida


def im2_config_file(image, im2_entrada, im2_psf, corrida, chains=2, niter=200):
	'''
	Im2Shape configuration file
	Input:
		image: 			(str) Name of the .fits image
		im2_entrada: 	(str) Name of the source catalog. Needs at least sex_id, sex_x, sex_y for each source
		im2_psf: 		(str) Name of the psf catalog. PSF computed at each source position
		corrida:		(int) Core number
		chains:			(int) Number of MCMC chains to run
		niter:			(int) Number of iteration of the MCMC

	Outpu:
		im2_conf: 		(str) config file name
		im2_salida:		(str) output file of im2shape
	'''	
	
	corrida  = str(corrida)
	chains 	 = str(chains)
	niter    = str(niter)

	im2_conf 	= im2_path + 'conf_'+corrida+'.in'
	im2_salida 	= im2_path + 'output_im2_'+corrida+'.dat'

	#os.system('rm '+im2_salida)		# delete any previous file
	f1=open(im2_conf,'w')
	f1.write('0 \n')
	f1.write(im2_salida+' \n')
	f1.write(image+'\n')
	f1.write(im2_entrada+'\n')
	f1.write('16\n')
	f1.write('2 0.000000 20.000000\n')
	f1.write('1 0.000000 1000.000000\n')
	f1.write('1\n')
	f1.write('1 0.000000 16.000000\n')
	f1.write('1 0.000000 16.000000\n')
	f1.write('3 0.000000 1.000000\n')
	f1.write('1 0.000000 3.141600\n')
	f1.write('2 -2.000000 8.000000\n')
	f1.write('2 0.000000 20.000000\n')
	f1.write(im2_psf+'\n')
	f1.write(chains+'\n')
	f1.write(niter+'\n')
	f1.write('model.txt\n')
	f1.write('2.000000\n')
	f1.write('0\n')
	f1.close()	

	return im2_conf, im2_salida
	
	
def im3_config_file(stampsize,niter,psf_input,cores,corrida):
	
	'''
	Im3Shape configuration file
	Input:
		stampsize: 		(str) The size of a postage-stamp patch with a single object in to analyze.
		niter: 	        (str) The maximum number of iterations of the minimizer
		psf_input:      (str) Model of the PSF to use
		corrida:		(int) Core number

	Output:
		im3_conf        (str) Name of the configuration catalog
	'''	
	
	corrida = str(corrida)
	cores = str(cores)
	im3_conf 	= im3_path + 'conf_'+cores+'.ini'


	#os.system('rm '+im2_salida)		# delete any previous file
	f1=open(im3_conf,'w')
	f1.write('model_name = sersics \n')
	f1.write('psf_input = '+psf_input+' \n')
	f1.write('use_segmentation_mask = NO \n')
	f1.write('segmentation_mask_filename = '+sex_path+'seg_'+corrida+'.fits \n')
	f1.write('noise_sigma = 1. \n')
	f1.write('rescale_stamp = Y \n')
	f1.write('sersics_x0_min = '+str((int(stampsize)/2.)-2)+' \n')
	f1.write('sersics_x0_max = '+str((int(stampsize)/2.)+2)+' \n')
	f1.write('sersics_y0_min = '+str((int(stampsize)/2.)-2)+' \n')
	f1.write('sersics_y0_max = '+str((int(stampsize)/2.)+2)+' \n')
	f1.write('sersics_bulge_A_max = 10000.0 \n')
	f1.write('sersics_disc_A_max = 10000.0 \n')
	f1.write('sersics_bulge_A_start = 0.5 \n')
	f1.write('sersics_disc_A_start = 0.5 \n')
	f1.write('sersics_disc_A_min = 0.0 \n')
	f1.write('upsampling = 5 \n')
	f1.write('n_central_pixel_upsampling = 9 \n')
	f1.write('n_central_pixels_to_upsample = 5 \n')
	f1.write('padding = 2 \n')
	f1.write('stamp_size = '+stampsize+' \n')
	f1.write('minimizer_max_iterations = '+niter+' \n')
	f1.write('levmar_eps1 = -1 \n')
	f1.write('levmar_eps2 = 1e-20 \n')
	f1.write('levmar_eps3 = -1 \n')
	f1.write('levmar_eps4 = 1e-07 \n')
	f1.write('levmar_tau = 1e-10 \n')
	f1.write('levmar_LM_INIT_MU = 1e-08 \n')
	f1.write('minimizer_loops = 1 \n')
	f1.write('sersics_bulge_A_min =  0.0 \n')
	f1.write('sersics_radius_start = 2.7 \n')
	f1.write('save_images = y \n')
	f1.write('sersics_beermat_amplitudes_positive = NO \n')
	f1.write('background_subtract = y \n')
	f1.write('psf_truncation_pixels = 10.0 \n')
	f1.write('minimizer_verbosity = -1 \n')
	f1.write('verbosity=4 \n')
	f1.close()	

	return im3_conf




	

def gx_catalog_header(FILTROS):
	'''
	Defines the columns that will be recorded in the galaxy catalog.
	'''
	# Header dictionary. The Key will be the field name in the catalog
	# ordered_dict  =  (  Key , [flag, SexFieldName, Format, Comment]  )

	sex_H = [
		('ID_SEX'		,	[True,	'NUMBER',		'1J', 'SExtractor run id']),
		('X_SEX'		,	[True,	'X_IMAGE',		'1E', 'Object position along x. SExtractor']),
		('Y_SEX'		,	[True,	'Y_IMAGE',		'1E', 'Object position along y. SExtractor']),
		('XWIN_SEX'		,	[True,	'XWIN_IMAGE',		'1E', 'Object position along x. SExtractor']),
		('YWIN_SEX'		,	[True,	'YWIN_IMAGE',		'1E', 'Object position along y. SExtractor']),
		('ALPHA'		,	[True,	'ALPHA_J2000',	'1D', 'Right ascension (J2000). SExtractor']),
		('DELTA'		,	[True,	'DELTA_J2000',	'1D', 'Declination (J2000). SExtractor']),
		('MAG_BEST'		,	[False,	'MAG_BEST',		'1E', 'Best of MAG_AUTO and MAG_ISOCOR']),
		('MAGERR_BEST'	,	[False,	'MAGERR_BEST',	'1E', 'RMS error for MAG_BEST']),
		('MAG_APER'		,	[False,	'MAG_APER',		'1E', 'Fixed aperture magnitude vector']),
		('MAGERR_APER'	,	[False,	'MAGERR_APER',	'1E', 'RMS error vector for fixed aperture mag']),
		('MAG_AUTO'		,	[True,	'MAG_AUTO',		'1E', 'Kron-like elliptical aperture magnitude in filter '+FILTROS[0]]),
		('MAGERR_AUTO'	,	[True,	'MAGERR_AUTO',	'1E', 'RMS error for AUTO magnitude']),
		('FWHM'			,	[True,	'FWHM_IMAGE',	'1D', 'FWHM assuming a gaussian core']),
		('FLUX_RADIUS'	,	[False,	'FLUX_RADIUS',	'1D', 'Fraction-of-light radii']),
		('FLUX_MAX'		,	[False,	'FLUX_MAX',		'1D', 'Peak flux above background']),
		('MU_MAX'		,	[False,	'MU_MAX',		'1D', 'Peak surface brightness above background']),
		('CLASS_STAR'	,	[False,	'CLASS_STAR',	'1D', 'S/G classifier output']),
		('FLAGS'		,	[False,	'FLAGS',		'1I', 'Extraction flags. SExtractor']),
		('SNR_WIN'		,	[False,	'SNR_WIN',		'1E', 'Gaussian-weighted SNR. SExtractor']),
		('FLUX_APER'	,	[False,	'FLUX_APER',	'1E', 'Flux vector within fixed circular aperture(s)'])
		]
	# python output
	im3_H = [
		('x0                  '	,	[True, 	     0	,    '1E', 'IM3SHAPE parameter']),
		('y0                  '	,	[True, 	     1	,    '1E', 'IM3SHAPE parameter']),
		('e1                  '	,	[True, 	     2	,    '1D', 'IM3SHAPE parameter']),
		('e2                  '	,	[True, 		 3	,    '1D', 'IM3SHAPE parameter']),
		('radius              '	,	[True, 		 4	,    '1E', 'IM3SHAPE parameter']),
		('radius_ratio        '	,	[True, 	     5	,    '1E', 'IM3SHAPE parameter']),
		('bulge_A             '	,	[True, 	     6	,    '1E', 'IM3SHAPE parameter']),
		('disc_A              '	,	[True, 	     7	,    '1E', 'IM3SHAPE parameter']),
		('bulge_index         '	,	[True, 	     8	,    '1E', 'IM3SHAPE parameter']),
		('disc_index          '	,	[True, 	     9	,	'1E', 'IM3SHAPE parameter']),
		('delta_e_bulge       '	,	[True, 	    10	,	'1E', 'IM3SHAPE parameter']),
		('delta_theta_bulge   '	,	[True, 		11	,	'1E', 'IM3SHAPE parameter']),
		('identifier          '	,	[True, 		12	,	'1E', 'IM3SHAPE parameter']),
		('time                '	,	[True, 		13	,	'1E', 'IM3SHAPE parameter']),
		('bulge_flux          '	,	[True, 		14	,	'1E', 'IM3SHAPE parameter']),
		('disc_flux           '	,	[True, 		15	,	'1E', 'IM3SHAPE parameter']),
		('flux_ratio          '	,	[True, 		16	,	'1E', 'IM3SHAPE parameter']),
		('snr                 '	,	[True,		17	,	'1E', 'IM3SHAPE parameter']),
		('old_snr             '	,	[True,		18	,	'1E', 'IM3SHAPE parameter']),
		('min_residuals       '	,	[True,		19	,	'1E', 'IM3SHAPE parameter']),
		('max_residuals       '	,	[True,		20	,	'1E', 'IM3SHAPE parameter']),
		('model_min           ' ,   [True,      21  ,    '1E', 'IM3SHAPE parameter']),
		('model_max           ' ,   [True,      22  ,    '1E', 'IM3SHAPE parameter']),
		('likelihood          ' ,   [True,      23  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_start_error  ' ,   [True,      24  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_end_error    ' ,   [True,      25  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_resid_grad   ' ,   [True,      26  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_vector_diff  ' ,   [True,      27  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_error_diff   ' ,   [True,      28  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_comp_grad    ' ,   [True,      29  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_iterations   ' ,   [True,      30  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_reason       ' ,   [True,      31  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_like_evals   ' ,   [True,      32  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_grad_evals   ' ,   [True,      33  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_sys_evals    ' ,   [True,      34  ,    '1E', 'IM3SHAPE parameter']),
		('mean_flux           ' ,   [True,      35  ,    '1E', 'IM3SHAPE parameter']),
		('number_varied_params' ,   [True,      36  ,    '1E', 'IM3SHAPE parameter']),
		('covariance matrix   ' ,   [True,      np.arange(37,86)  ,    '49E', 'IM3SHAPE parameter']),
		('psf_fwhm            ' ,   [False,      86  ,    '1E', 'IM3SHAPE parameter']),
		('psf_beta            ' ,   [False,      87  ,    '1E', 'IM3SHAPE parameter']),
		('psf_e2              ' ,   [False,      88  ,    '1E', 'IM3SHAPE parameter']),
		('psf_e1              ' ,   [False,      89  ,    '1E', 'IM3SHAPE parameter'])
		]
	'''
		
	im3_H = [
		('ID'	,	        [True, 	     0	,    '1J', 'IM3SHAPE parameter']),
		('catalog_x'	,	[True, 	     1	,    '1E', 'IM3SHAPE parameter']),
		('catalog_y'	,	[True, 	     2	,    '1D', 'IM3SHAPE parameter']),
		('Likelihood'	,	[True, 		 3	,    '1D', 'IM3SHAPE parameter']),
		('x0'	,	        [True, 		 4	,    '1E', 'IM3SHAPE parameter']),
		('y0'	,	        [True, 	     5	,    '1E', 'IM3SHAPE parameter']),
		('e1'	,	        [True, 	     6	,    '1E', 'IM3SHAPE parameter']),
		('e2'	,	        [True, 	     7	,    '1E', 'IM3SHAPE parameter']),
		('radius'	,	    [True, 	     8	,    '1E', 'IM3SHAPE parameter']),
		('bulge_A'	,	    [True, 	     9	,	'1E', 'IM3SHAPE parameter']),
		('disc_A'	,	    [True, 	    10	,	'1E', 'IM3SHAPE parameter']),
		('flux_ratio'	,	[True, 		11	,	'1E', 'IM3SHAPE parameter']),
		('fwhm'	,	        [True, 		12	,	'1E', 'IM3SHAPE parameter']),
		('snr'	,	        [True, 		13	,	'1E', 'IM3SHAPE parameter']),
		('min_residuals',	[True, 		14	,	'1E', 'IM3SHAPE parameter']),
		('max_residuals',	[True, 		15	,	'1E', 'IM3SHAPE parameter']),
		('model_min'	,	[True, 		16	,	'1E', 'IM3SHAPE parameter']),
		('model_max'	,	[True,		17	,	'1E', 'IM3SHAPE parameter']),
		('psf_e1'	,	    [True,		18	,	'1E', 'IM3SHAPE parameter']),
		('psf_e2'	,	    [True,		19	,	'1E', 'IM3SHAPE parameter']),
		('psf_fwhm'	,	    [True,		20	,	'1E', 'IM3SHAPE parameter']),
		('psf_beta' ,       [True,      21  ,    '1E', 'IM3SHAPE parameter']),
		('Rgpp_Rp' ,        [True,      22  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_e0' ,      [True,      23  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_e' ,       [True,      24  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_J' ,       [True,      25  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_Dp' ,      [True,      26  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_mu' ,      [True,      27  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_De' ,      [True,      28  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_nit' ,     [True,      29  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_reason' ,  [True,      30  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_neval' ,   [True,      31  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_njac' ,    [True,      32  ,    '1E', 'IM3SHAPE parameter']),
		('levmar_nlin' ,    [True,      33  ,    '1E', 'IM3SHAPE parameter']),
		('time_taken' ,     [True,      34  ,    '1E', 'IM3SHAPE parameter'])
		]
	'''
	
	dic=[]
	
	
	
	for i,j in enumerate(FILTROS[1:]):
		d = ('MAG_'+j, [True, i, '1E', 'MAG_AUTO in filter '+j])
		dic.append(d)

	if len(FILTROS[1:]) == 0:
		mag_H = []
	else:
		mag_H = [d]
		
	rot_H = [
		('e1_sky                  '	,	[True, 	     0	,    '1E', 'e1 wcs component']),
		('e2_sky                 '	,	[True, 	     1	,    '1E', 'e2 wcs component'])]

	sex_H = collections.OrderedDict(sex_H)
	im3_H = collections.OrderedDict(im3_H) 
	mag_H = collections.OrderedDict(mag_H) 
	rot_H = collections.OrderedDict(rot_H) 

	return sex_H, im3_H, mag_H, rot_H

def write_header(hdul, full_H):
	'''
	Writes the header
	'''
	for key, field in hdul.header.items():
		try:
			hdul.header.comments[key] = full_H[field][3]		# Embeded in a try to avoid extra header keys
		except:
			pass
	return hdul

#tbhdu = fits.BinTableHDU.from_columns(



def merge_gx_catalog(hdul_gx, MAGNITUDES, FILTROS, im3_txt, rot_coords):
	'''
	Juntamos las salidas de sextractor e im2shape en un catalogo binario
	'''
	#hdul_gx = fits.open('./sex_files/run2-r1.cat')

	#im2_salida = './im2_files/gx01.out'
	#im3_txt = np.loadtxt(im3_salida, skiprows=1)


	sex_H, im3_H, mag_H, rot_H = gx_catalog_header(FILTROS)

	# SEXTRACTOR ---------------------------------------------------------------------
	sex_cols = []
	for key in sex_H.keys():
		if not sex_H[key][0]: continue							# if False this field wont be recorded
		try:
			data = hdul_gx[2].data[sex_H[key][1]]				# extracts the column for sextractor catalog
			fmt  = sex_H[key][2]								# format of that column
			c = fits.Column(name=key, format=fmt, array=data)
			sex_cols.append(c)
		except:
			print clr.WARNING+'The column '+sex_H[key][1]+' does not exist in '+hdul_gx.fileinfo(0)['filename']+clr.ENDC

	#sex_cols = fits.ColDefs(sex_cols)							# Defines columns
	del hdul_gx

	magnitude_cols = []
	
	# Other filter magnitudes
	
	for key in mag_H:
		if not mag_H[key][0]: continue
		c = fits.Column(name=key, format=mag_H[key][2], array=MAGNITUDES[:, mag_H[key][1]])
		magnitude_cols.append(c)
	

	# IM3SHAPE -----------------------------------------------------------------------
	im3_cols = []
	for key in im3_H.keys():
		if not im3_H[key][0]: continue							# if False this field wont be recorded
		try:
			len(im3_H[key][1])

			data = im3_txt[:, im3_H[key][1]]
			data = data.reshape((-1,7,7))
			fmt  = im3_H[key][2]									# format of that column
			c = fits.Column(name=key, format=fmt, array=data, dim = '(7,7)')
			
		except:	
			data = im3_txt[:, im3_H[key][1]]						# extracts the column for im2shape catalog
			fmt  = im3_H[key][2]									# format of that column
			c = fits.Column(name=key, format=fmt, array=data)

		im3_cols.append(c)

	#im3_cols = fits.ColDefs(im3_cols)							# Defines columns
	del im3_txt
	
	
	# ROT COORDS ----------------------------------------------------------------------
	
	rot_cols = []
	
	for key in rot_H:
		if not rot_H[key][0]: continue
		c = fits.Column(name=key, format=rot_H[key][2], array=rot_coords[:, rot_H[key][1]])
		rot_cols.append(c)

	# MERGE ---------------------------------------------------------------------------
	all_cols = sex_cols + magnitude_cols + im3_cols + rot_cols
	hdul = fits.BinTableHDU.from_columns(all_cols)

	full_H = sex_H.copy()
	full_H.update(mag_H)
	full_H.update(im3_H)
	full_H.update(rot_H)


	hdul = write_header(hdul, full_H)
	
	return hdul

def write_2final_catalog(salida_file, hdul_gx):
    '''
    Writes the output to the final catlog
    '''

    # We first check if the file already exists. If False it means it is the first loop to write and has to create it.
    # If True it just appends to the end
    exists = os.path.isfile(salida_file)

    sex_H, im3_H, mag_H, rot_H = gx_catalog_header('r')
    full_H = sex_H.copy()
    full_H.update(mag_H)
    full_H.update(im3_H)
    full_H.update(rot_H)                # Merge both dict

    if exists:
        hdul_final     = fits.open(salida_file)
        nrows_old     = len(hdul_final[1].data)
        nrows_gx     = len(hdul_gx.data)
        hdul_aux     = fits.BinTableHDU.from_columns(hdul_final[1].columns,
                                                    nrows=nrows_old+nrows_gx,
                                                    header=hdul_final[1].header)
        
        for colname in hdul_final[1].columns.names:
            hdul_aux.data[colname][nrows_old:] = hdul_gx.data[colname]

        # Writes the header
        hdul_aux = write_header(hdul_aux, full_H)

        hdul_final.close()
        hdul_aux.writeto(salida_file, overwrite=True)
        del hdul_final, hdul_aux
    else:
        # Writes the header
        hdul_gx = write_header(hdul_gx, full_H)

        hdul_gx.writeto(salida_file)
        del hdul_gx

    return None


def merge_hdul_parallel(salida):
    '''
    Juntar las tablas de todos los procesos en una sola
     input:
         salida:            List with multiple threads output. Each element corresponds to a single thread.
                         It should contain the hdul and number of galaxies of each thread.

     output:
         hdul_merged:    Hdul object with the data of every thread
    '''
    # Specify the index in wich these variables appear in salida
    i_hdul = 1        # hdul index in each element of salida
    i_ngx  = 2        # N_gx index in each element of salida

    N_tot = sum([corrida[i_ngx] for corrida in salida])        # Sums the number of gx of each process
    hdul_0 = salida[0][i_hdul]
    hdul_merged    = fits.BinTableHDU.from_columns(hdul_0.columns,
                                                nrows=N_tot,
                                                header=hdul_0.header)
    prev_rows = 0
    for i in xrange(1,len(salida)):
        hdul = salida[i][i_hdul]
        i_rows = salida[i][i_ngx]
        prev_rows += salida[i-1][i_ngx]
        for colname in hdul_0.columns.names:
            hdul_merged.data[colname][prev_rows : prev_rows+i_rows] = hdul.data[colname]

    del salida, hdul_0, hdul

    return hdul_merged

