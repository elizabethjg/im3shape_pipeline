import numpy as np
from pylab import *
from astropy.io import fits
import shape_models
from ctypes import *
import psfex

class PSF:
	'''
	Class object that reconstructs the PSF in the CCD
	 __init__: Creates the object
		input: 
			psf_file:	(str)	PSFEx .psf file containing the psf basis
	
	 get_psf:	Reconstructs the PSF
		input:
			col, row:	(flt) 	Positions X,Y in the CCD

		output:
			psf_model:	(array2D) PSF vignette centred at X,Y
	'''

	def __init__(self, psf_file):
		hdul = fits.open(psf_file)

		self.sources = hdul[1].header['ACCEPTED'] #Number of accepted sources
		self.degree = hdul[1].header['POLDEG1']
		self.dim 	= hdul[1].header['PSFAXIS3']
		self.size 	= hdul[1].header['PSFAXIS1']
		if hdul[1].header['PSFAXIS1'] != hdul[1].header['PSFAXIS2']: print 'WARNING: PSF vignete is not square.'
		
		self.basis 	= hdul[1].data['PSF_MASK'][0]
		self.zero_x	 = hdul[1].header['POLZERO1']
		self.scale_x = hdul[1].header['POLSCAL1']
		self.zero_y	 = hdul[1].header['POLZERO2']
		self.scale_y = hdul[1].header['POLSCAL2']
		#self.samp 	= hdul[1].header['PSF_SAMP']		# que pasa cuando no es 1
		del hdul

	def get_psf(self, col, row):
		x_s = (col - self.zero_x) / self.scale_x
		y_s = (row - self.zero_y) / self.scale_y
		poly = [[x_s**i * y_s**j for i in xrange(self.degree+1-j)] for j in xrange(self.degree+1)]
		poly = np.array(sum(poly))	# fletten list and convert to array

		terms = np.multiply(self.basis, poly[:, np.newaxis, np.newaxis])
		psf_model = terms.sum(axis=0)
		del poly,terms
		return psf_model


def psf_map(psf_file, CCD_shape):
	'''
	Creates a map of the PSF across the CCD
	 input:
	 	psf_file:	(str) 	PSFEx .psf file containing the psf basis
		CCD_shape:	(tuple)	Shape of the CCD (X, Y)
	 output:
	 	Nothing yet... it only shows the plot. Maybe it could save a figure and return the file name, or just the plt object
	'''

	#psf_file = './psfex_files/run2-r2.psf'
	#CCD_shape = (2048,4096)

	#psf = PSF(psf_file)
	pex = psfex.PSFEx(psf_file)
	
	x = CCD_shape[0] * np.linspace(0.05,0.95,10)
	y = CCD_shape[1] * np.linspace(0.05,0.95,10)
	xx, yy = np.meshgrid(x,y, indexing='ij')

	ellip = np.zeros(xx.shape)
	theta = np.zeros(xx.shape)

	e1 = np.zeros((len(x),len(y)))
	e2 = np.zeros((len(x),len(y)))
	
	out = np.zeros(4)
	
	for j in xrange(len(y)):
		for i in xrange(len(x)):
			#psf_stamp = psf.get_psf(col=x[i], row=y[j])
			psf_stamp = pex.get_rec(x[i], y[j])
			p = shape_models.FitMoffat2D(psf_stamp/psf_stamp.max())
			out = np.vstack((out,np.array([p['beta'],p['fwhm'],p['e1'],p['e2']])))
			e1[i,j] = p['e1']
			e2[i,j] = p['e2']

	out = out[1:,:]
	#theta = np.pi/2.-np.arctan2(e2,e1)/2.0
	ellip = np.sqrt(e1**2+e2**2)
	theta = (np.arctan2(e2,e1)/2.)
	theta = np.pi/2. - theta
	
	ex = np.cos(theta)
	ey = np.sin(theta)

	plt.figure()
	plt.quiver(xx-ex/2, yy-ey/2, ex, ey, headwidth=1,headlength=0,units='xy')
	plt.axes().set_aspect('equal')
	plt.xlim([0,CCD_shape[0]])
	plt.ylim([0,CCD_shape[1]])
	plt.show()

	return out

def compute_psf_4im2shape(psf_file, hdu, corrida):
	'''
	Computes the PSF at the location of each source and creates the input files for im2shape
	 input:
	 	psf_file:	(str)    PSFEx .psf file containing the psf basis
	 	hdu:		(object) hdulist object from the fits files
	 	corrida:	(int)    core number	

	 output:
		im2_entrada:	(str) File name of the input catalog of sources for im2shape 
		im2_psf:		(str) File name of the psf input file for im2shape
	'''

	im2_path 	= './im2_files/'
	corrida 	= str(corrida)
	im2_psf 	= im2_path + 'psf_im2_'+corrida+'.dat'
	im2_entrada = im2_path + 'input_im2_'+corrida+'.dat'
	
	# Computes psf for im2shape-----------------------------------------------------------
	psf = PSF(psf_file)
	ids = hdu[2].data['NUMBER']
	x   = hdu[2].data['X_IMAGE']
	y   = hdu[2].data['Y_IMAGE']

	N 	= len(x)								# Number of sources
	psf_src = np.zeros((N, 8), dtype=float)		# Array to save the psf parameters

	for i in xrange(N):
		psf_stamp = psf.get_psf(col=x[i], row=y[i])
		p = shape_models.FitGaussian2D(psf_stamp)
		psf_src[i,:] = [x[i], y[i], 0,0, p['ellip'], p['theta'], p['ab'], p['amp']]		# deberiamos usar x,y ajustados ??

	# Writes psf file for im2shape-----------------------------------------------------------
	f1=open(im2_psf,'w')
	f1.write('1\n')			# Number of gaussians
	f1.write(str(N)+'\n')	# Number of psf locations
	f1.write('6\n')			# Number of gaussian parameters
	np.savetxt(f1, psf_src, fmt=['%15.10f']*4+['%13.10f']*4)
	f1.close()

	# Writes input catalog file for im2shape------------------------------------------------
	cat = np.vstack((ids, x, y)).T 		# This are the columns required by im2shape

	f1=open(im2_entrada,'w')
	f1.write('3 '+str(N)+'\n')			# Number of columns, number of rows
	np.savetxt(f1, cat, fmt=['%5i']+['%15.8f']*2)
	f1.close()
	
	return im2_entrada, im2_psf
	
	
def compute_psf_4im3shape(psf_file, options, ids, x, y, corrida):

	'''
	Computes the PSF at the location of each source and creates the input files for im32shape
	 input:
	 	psf_file:	(str)    PSFEx .psf file containing the psf basis
	 	hdu:		(object) hdulist object from the fits files
	 	corrida:	(int)    core number	
	 	corrida:	(obj)    im3shape options

	 output:
		im3_entrada:	(str) File name of the input catalog of sources for im3shape 
		im3_psf:		(str) File name of the psf input file for im3shape
	'''
	import galsim, galsim.des
	from py3shape.utils import *

	im3_path 	= './im3_files/'
	
	corrida 	= str(corrida)
	im3_entrada = im3_path + 'input_im3_'+corrida+'.dat'
	
	# Computes psf for im3shape-----------------------------------------------------------
	#psf = PSF(psf_file)
	
	N 	= len(x)								# Number of sources
	


	if options.psf_input == 'psf_image_cube':
		
		im3_psf = im3_path + 'psf_im3_cube_'+corrida+'.fits'
		u = options.upsampling
		pad = options.padding
		stamp_size = options.stamp_size
		Nx=(stamp_size+pad)*u
		psfex_galsim = galsim.des.DES_PSFEx(psf_file)
		
		stamps = []
		
		for i in xrange(N):

			psfex_stamp = getPSFExarray(psfex_galsim, x[i], y[i], nsidex=Nx, nsidey=Nx, upsampling=u)
			
			stamps.append(galsim.ImageF(psfex_stamp))
			
		galsim.fits.writeCube(stamps, im3_psf)
			
				
	elif options.psf_input == 'moffat_catalog':
		
		pex = psfex.PSFEx(psf_file)
		
		psf_src = np.zeros((N, 4), dtype=float)		# Array to save the psf parameters
		im3_psf 	= im3_path + 'psf_im3_'+corrida+'.dat'

		for i in xrange(N):
			
			#psf_stamp = psf.get_psf(col = x[i], row = y[i])
			psf_stamp = pex.get_rec(x[i], y[i])
			p = shape_models.FitMoffat2D(psf_stamp/psf_stamp.max())
			e = (p['e1']**2+p['e2']**2)**0.5
			th = (np.arctan2(p['e2'],p['e1'])/2.)
			th = np.pi/2.0 - th
			#if th > np.pi/2.:
				#th = th - np.pi
			
			psf_src[i,:] = [p['beta'], p['fwhm'], e*np.cos(2.*th), e*np.sin(2*th)]		# deberiamos usar x,y ajustados ?? no
	
		# Writes psf file for im2shape-----------------------------------------------------------
		f1 = open(im3_psf,'w')
		f1.write('# beta fwhm e1 e2 \n')			# Number of gaussians
		np.savetxt(f1, psf_src, fmt=['%15.10f']*4)
		f1.close()
		
		del pex,psf_src
	elif options.psf_input == 'no_psf':
		
		im3_psf = im3_path + 'psf_im3_none_'+corrida+'.dat'
		f1 = open(im3_psf,'w')
		f1.write('# beta fwhm e1 e2 \n')
		f1.close()

	# Writes input catalog file for im3shape------------------------------------------------
	cat = np.vstack((ids, x, y)).T 		# This are the columns required by im3shape

	f1 = open(im3_entrada,'w')
	f1.write(' #ID   x    y \n')			# Number of columns, number of rows
	np.savetxt(f1, cat, fmt=['%5i']+['%15.8f']*2)
	f1.close()
	
	
	
	return im3_entrada, im3_psf
