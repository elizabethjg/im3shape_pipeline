import numpy as np
from astropy.modeling import fitting, models
import random
import os
from astropy.io import fits
#from pylab import *
from scipy import signal

###############################################################################################################
###############################################################################################################
#
# GENERAL TOOLS
#

def shear_matrix(e1,e2):
	'''
	Shear matrix as used by im3shape
	'''
	e = (e1*e1 + e2*e2)**0.5

	# By default, the shear is clockwise from y axis.
	e1 *= -1		# Multiply e1 by -1 to change to x axis and make it counter-clockwise
	#e2 *= -1		# Multiply e2 by -1 to to keep y axis but counter-clockwise

	m1 = (1 + e - 2*e1)/(1-e)
	m2 = -2*e2/(1-e)
	m3 = (1 + e + 2*e1)/(1-e)

	return m1,m2,m3

def convolution(A,B):
	'''
	Computes the convolution between two image stamps
	'''
	if A.shape != B.shape: print 'WANRING: Input arrays have different shapes.'
	C = signal.fftconvolve(A, B, mode='same')

	return C


###############################################################################################################
###############################################################################################################
#
# GAUSSIAN MODEL SECTION
#

def Gauss2D(x,y,amp=1,ab=0.1,e1=0.0,e2=0.0,x0=0,y0=0):
	'''
	Evaluates a 2D gaussian using 6 parameters.
	 input:
	 	x,y:	Position to evaluate the gaussian.
	 			Can be scalar values or grids: x, y = np.meshgrid(range(Nx), range(Ny), indexing='ij')
		amp:	Amplitude of the central peak [counts]
		ab:		Product of semi-axes or area scale [counts**2]
		x0,y0:	Position of the peak [pixels]
		e1:		Ellipticity component along x,y axis
		e2:		Ellipticity component along x=y axis (i.e. 45 degrees of e1)
	
	 output:
		G:		Value of 2D gaussian at x,y [counts]
	'''

	dx = x - x0
	dy = y - y0
	m1,m2,m3 = shear_matrix(e1,e2)
	r2 = m1*dx*dx + 2*m2*dx*dy + m3*dy*dy

	G = amp * np.exp(-r2/ab)		# Gaussian 2D
	
	return G

def Gauss2D_im2shape(x,y,amp=1,ab=0.1,ellip=0.0,x0=0,y0=0,theta=0):
	'''
	Evaluates a 2D gaussian using 6 parameters just like Im2Shape.
	 input:
	 	x,y:	Position to evaluate the gaussian.
	 			Can be scalar values or grids: x, y = np.meshgrid(range(Nx), range(Ny), indexing='ij')
		amp:	Amplitude as used by im2shape. This is not the height of the central peak [counts]
		ab:		Product of semi-axes or area scale [counts**2]
		ellip:	Ellipticity. Defined as (a-b)/(a+b)
		x0,y0:	Position of the peak [pixels]
		theta:	Rotation angle ranging from 0 to pi counterclockwise from the positive y axis [radians]
	
	 output:
		G:		Value of 2D gaussian at x,y [counts]
	'''

	b_a = (1.-ellip)/(1.+ellip)				# This is = b/a
	a = np.sqrt(ab/b_a)						
	b = np.sqrt(ab*b_a)						

	dx = x - x0
	dy = y - y0
	# Covariance matrix. C12 = C21
	C11 = 2 * ((np.cos(theta)/a)**2 + (np.sin(theta)/b)**2)
	C12 = (-1/(a*a) + 1/(b*b)) * np.sin(2*theta)
	C22 = 2 * ((np.sin(theta)/a)**2 + (np.cos(theta)/b)**2)
	detC = C11 * C22 - C12 * C12

	G = amp/(2*np.pi*detC) * np.exp(-0.5*(C11*dx*dx + 2*C12*dx*dy + C22*dy*dy))		# Gaussian 2D
	
	return G


def Gauss2D_stamp(amp=1,ab=20,e1=0.,e2=0.,x0=25,y0=25,shape=(51,51),noise=0,zerolevel=0):
	'''
	Contructs a 2D gaussian vignette with background noise. Uses 6 parameters just like Im2Shape.
	 input:
		amp:	Amplitude as used by im2shape. This is not the height of the central peak [counts]
		ab:		Product of semi-axes or area scale [counts**2]
		ellip:	Ellipticity. Defined as (a-b)/(a+b)
		x0,y0:	Position of the peak [pixels]
		theta:	Rotation angle ranging from 0 to pi counterclockwise from the positive y axis [radians]
		shape:	Tuple indicating the shape of the vignette. Make it odd for a central peak (51, 51)
		noise:	Fraction of amp to be used as the gaussian dispersion of the background noise N(0,noise*amp)
		zerolevel:	Constant background floor [counts]
	
	 output:
		V:		Vignette with the 2D gaussian
	'''
	
	x, y = np.meshgrid(range(shape[0]), range(shape[1]), indexing='ij')
	V = Gauss2D(x,y,amp=amp,ab=ab,e1=e1,e2=e2,x0=x0,y0=y0)

	if noise!=0: 		V += np.random.normal(loc=0, scale=noise*amp, size=shape)		# Adds normal noise
	if zerolevel!=0: 	V += zerolevel													# Adds constant background

	return V

def FitGaussian2D(B, reconstruct_model=False):
	'''
	Fits a 2D gaussian model to a vignette using Levenberg-Marquardt algorithm and least squares statistic.
	 input:
	 	B:		Vignette with the object
	 	reconstruct_model:	Boolean flag. Indicates if the reconstructed model and its residuals should be computed
	
	 output:
	 	params:		Dictionary with the 6 fitted parameters. Keys = ['amp', 'ab', 'ellip', 'theta', 'x0', 'y0']
		*rec_model:	Vignette with the fitted model. Only if reconstruct_model=True
		*rec_resi:	Vignette with the residuals, ie: rec_resi=B-rec_model. Only if reconstruct_model=True

	Notes: Given the polar nature of the ellipticity definition, the fit can be confusing. For ex., given a
	gaussian with ellip = 0.1 and theta = 0 (where theta is measured from the y axis), it is posible that the
	fitter returns ellip = -0.1 and theta = pi/2. This time the minus sign in ellip indicates that it is meassuring
	it on the oposite axis, thus theta is being measured from the x axis.
	We fix the position x0,y0 of the model at the center of the stamp. There is no need to fit them
	Notes 2: Add the posibility of fitting a sum of gaussians instead of just one.
	'''	

	p = model_Gaussian2D(amp=B.max()-B.min(), ab=5, e1=0., e2=0.2,x0=B.shape[0]/2,y0=B.shape[1]/2) + \
		models.Const2D(amplitude=B.min())

	# Constraints on some parameters
	p.x0_0.fixed = True
	p.y0_0.fixed = True

	x, y = np.meshgrid(range(B.shape[0]), range(B.shape[1]), indexing='ij')
	out_gauss, out_const= fitter(p,x,y,B,maxiter=10000,acc=1e-16)		# fitter() is defined globally

	# retrieves the 6 gaussian parameters
	params = {'amp': out_gauss.amp.value,
			  'ab': out_gauss.ab.value,
			  'x0': out_gauss.x0.value,
			  'y0': out_gauss.y0.value,
			  'e1': out_gauss.e1.value, 
			  'e2':out_gauss.e2.value}

	# generates a vignette with the model and computes residuals
	if reconstruct_model:
		rec_model = Gauss2D_stamp(amp=params['amp'],ab=params['ab'],e1=params['e1'],
									e2=params['e2'],x0=params['x0'],y0=params['y0'],
									noise=0,zerolevel=out_const.amplitude, shape=B.shape)
		resi = B - rec_model
		return params, rec_model, resi
	else:
		return params	

###############################################################################################################
###############################################################################################################
#
# MOFFAT MODEL SECTION
#

def Moffat2D(x,y,x0=0,y0=0,e1=0.3,e2=0.0,fwhm=5,beta=1):
	'''
	Evaluates a 2D Moffat profile of amplitude 1 just like im3shape.
	 input:
	 	x,y:	Position to evaluate the moffat.
	 			Can be scalar values or grids: x, y = np.meshgrid(range(Nx), range(Ny), indexing='ij')
		x0,y0:	Position of the peak [pixels]
		e1:		Ellipticity component along x,y axis
		e2:		Ellipticity component along x=y axis (i.e. 45 degrees of e1)
		fwhm:	Full-widht at half maximum [counts]
		beta:	Profile power index
	
	 output:
		M:		Value of 2D moffat at x,y [counts]

	Notes: The Moffat profile is given by f(r) = (1+r^2/d^2)^(-beta).  We shear by changing what r means.
	The equations were taken from im3shape source code. We used these functions:
	i3_unit_shear_matrix()			- tools/i3_image.c
	i3_image_add_moffat()			- tools/i3_psf.c
	i3_image_add_truncated_moffat()	- tools/i3_psf.c
	i3_moffat_fwhm()				- tools/i3_psf.c
	i3_moffat_radius()				- tools/i3_psf.c
	
	We compute the radius d with the FWHM. f(0)=1, so to compute the FWHM we want to solve f(r)=1/2
	This gives: r = d (2^(1/beta)-1)^(1/2). Then FWHM = 2r = 2 d (2^(1/beta)-1)^(1/2)
	'''

	radius2 = 0.25 * fwhm**2 / (2**(1./beta)-1)		# radius squared, d**2
	ellip=np.sqrt(e1**2+e2**2)

	dx = x - x0
	dy = y - y0
	# Shear matrix as used by im3shape
	m1 = (1 + ellip - 2*e1)/(1-ellip)
	m2 = -2*e2/(1-ellip)
	m3 = (1 + ellip + 2*e1)/(1-ellip)
	sheared_r2 = m1*dx*dx + 2*m2*dx*dy + m3*dy*dy

	M = (1 + sheared_r2/radius2)**(-beta)			# Moffat 2D

	return M

def Moffat2D_stamp(x0=25,y0=25,e1=0.3,e2=0.0,fwhm=5,beta=1,shape=(51,51),noise=0,zerolevel=0):
	'''
	Contructs a 2D moffat vignette of amplitude 1 with background noise. Just like Im3Shape.
	 input:
		x0,y0:	Position of the peak [pixels]
		e1:		Ellipticity component along x,y axis
		e2:		Ellipticity component along x=y axis (i.e. 45 degrees of e1)
		fwhm:	Full-widht at half maximum [counts]
		beta:	Profile power index
		shape:	Tuple indicating the shape of the vignette. Make it odd for a central peak (51, 51)
		noise:	Gaussian dispersion of the background noise N(0,noise). Value should be 0 < noise < 1
		zerolevel:	Constant background floor [counts]
	
	 output:
		V:		Vignette with the 2D Moffat
	'''
	
	x, y = np.meshgrid(range(shape[0]), range(shape[1]), indexing='ij')
	V = Moffat2D(x,y,x0=x0,y0=y0,e1=e1,e2=e2,fwhm=fwhm,beta=beta)

	if noise!=0: 		V += np.random.normal(loc=0, scale=noise, size=shape)		# Adds normal noise
	if zerolevel!=0: 	V += zerolevel												# Adds constant background

	return V

def FitMoffat2D(B, reconstruct_model=False):
	'''
	Fits a 2D moffat model to a vignette using Levenberg-Marquardt algorithm and least squares statistic.
	 input:
	 	B:		Vignette with the object
	 	reconstruct_model:	Boolean flag. Indicates if the reconstructed model and its residuals should be computed
	
	 output:
	 	params:		Dictionary with the 6 fitted parameters. Keys = ['beta', 'fwhm', 'e1', 'e2', 'x0', 'y0']
		*rec_model:	Vignette with the fitted model. Only if reconstruct_model=True
		*rec_resi:	Vignette with the residuals, ie: rec_resi=B-rec_model. Only if reconstruct_model=True

	Notes: We fix the position x0,y0 of the model at the center of the stamp. There is no need to fit them
	'''	

	p = model_Moffat2D(fwhm=5, beta=1, x0=B.shape[0]/2,y0=B.shape[1]/2,e1=0.1, e2=0.1) + \
		models.Const2D(amplitude=B.min())
	
	# Constraints on some parameters
	p.x0_0.fixed = True
	p.y0_0.fixed = True
	
	x, y = np.meshgrid(range(B.shape[0]), range(B.shape[1]), indexing='ij')
	out_moffat, out_const= fitter(p,x,y,B,maxiter=10000,acc=1e-16)			# fitter() is defined globally

	# retrieves the 6 moffat parameters
	params = {'beta': out_moffat.beta.value,
			  'fwhm': out_moffat.fwhm.value,
			  'x0': out_moffat.x0.value,
			  'y0': out_moffat.y0.value,
			  'e1': out_moffat.e1.value, 
			  'e2': out_moffat.e2.value} 

	# generates a vignette with the model and computes residuals
	if reconstruct_model:
		rec_model = Moffat2D_stamp(beta=out_moffat.beta,fwhm=out_moffat.fwhm,x0=out_moffat.x0,
									y0=out_moffat.y0,e1=out_moffat.e1,e2=out_moffat.e2,noise=0,
									zerolevel=out_const.amplitude, shape=B.shape)
		resi = B - rec_model
		return params, rec_model, resi
	else:
		return params	

###############################################################################################################
###############################################################################################################
#
# SERSIC MODEL SECTION
#

def Sersic2D(x,y,amp=1,radius=5,e1=0,e2=0,x0=0,y0=0,sersic_n=4):
	'''
	Evaluates a 2D Sersic profile.
	 input:
	 	x,y:		Position to evaluate the profile.
	 				Can be scalar values or grids: x, y = np.meshgrid(range(Nx), range(Ny), indexing='ij')
		x0,y0:		Position of the peak [pixels]
		e1:			Ellipticity component along x,y axis
		e2:			Ellipticity component along x=y axis (i.e. 45 degrees of e1)
		amp:		Amplitud of the central peak [counts]
		radius: 	Radius of half-light [pixels]
		sersic_n:	Sersic index. n=4 gives the de Vaucouleurs profile wich is a good fit for elliptical
					galaxies. n=1 gives the normal exponential profile wich is a good fit for spirals.
	
	 output:
		S:			Value of 2D sersic at x,y [counts]
	'''

	kappa = 1.9992*sersic_n - 0.3271 	# Voigt and Bridle 2010: https://arxiv.org/pdf/0905.4801.pdf

	dx = x - x0
	dy = y - y0
	m1,m2,m3 = shear_matrix(e1,e2)
	r2 = m1*dx*dx + 2*m2*dx*dy + m3*dy*dy
	radius2 = radius*radius

	S = amp * np.exp(-kappa * (r2/radius2)**(0.5/sersic_n))		# Sersic 2D

	return S

def Sersic2D_stamp(amp=1,radius=5,e1=0,e2=0,x0=25,y0=25,sersic_n=4,shape=(51,51),noise=0,zerolevel=0):
	'''
	Contructs a 2D sersic vignette with background noise.
	 input:
		x0,y0:		Position of the peak [pixels]
		e1:			Ellipticity component along x,y axis
		e2:			Ellipticity component along x=y axis (i.e. 45 degrees of e1)
		amp:		Amplitud of the central peak [counts]
		radius: 	Radius of half-light [pixels]
		sersic_n:	Sersic index. n=4 gives the de Vaucouleurs profile wich is a good fit for elliptical
					galaxies. n=1 gives the normal exponential profile wich is a good fit for spirals.
		shape:		Tuple indicating the shape of the vignette. Make it odd for a central peak (51, 51)
		noise:		Gaussian dispersion of the background noise N(0,noise). Value should be 0 < noise < 1
		zerolevel:	Constant background floor [counts]
	
	 output:
		V:			Vignette with the 2D Sersic profile
	'''
	
	x, y = np.meshgrid(range(shape[0]), range(shape[1]), indexing='ij')
	V = Sersic2D(x,y,amp=amp,radius=radius,x0=x0,y0=y0,e1=e1,e2=e2,sersic_n=sersic_n)

	if noise!=0: 		V += np.random.normal(loc=0, scale=noise, size=shape)		# Adds normal noise
	if zerolevel!=0: 	V += zerolevel												# Adds constant background

	return V

def FitSersic2D(B, reconstruct_model=False):
	'''
	Fits a 2D sersic model to a vignette using Levenberg-Marquardt algorithm and least squares statistic.
	 input:
	 	B:		Vignette with the object
	 	reconstruct_model:	Boolean flag. Indicates if the reconstructed model and its residuals should be computed
	
	 output:
	 	params:		Dictionary with the 7 fitted parameters. Keys = ['amp', 'radius', 'sersic_n', 'e1', 'e2', 'x0', 'y0']
		*rec_model:	Vignette with the fitted model. Only if reconstruct_model=True
		*rec_resi:	Vignette with the residuals, ie: rec_resi=B-rec_model. Only if reconstruct_model=True

	Notes: We fix the position x0,y0 of the model at the center of the stamp. There is no need to fit them
	'''	

	p = model_Sersic2D(amp=B.max()-B.min(), radius=10, sersic_n=4, x0=B.shape[0]/2,y0=B.shape[1]/2,e1=0.1, e2=0.1) + \
		models.Const2D(amplitude=B.min())
	
	# Constraints on some parameters
	p.x0_0.fixed = True
	p.y0_0.fixed = True
	
	x, y = np.meshgrid(range(B.shape[0]), range(B.shape[1]), indexing='ij')
	out_sersic, out_const= fitter(p,x,y,B,maxiter=10000,acc=1e-16)			# fitter() is defined globally

	# retrieves sersic parameters
	params = {'amp': out_sersic.amp.value,
			  'sersic_n': out_sersic.sersic_n.value,
			  'radius': out_sersic.radius.value,
			  'x0': out_sersic.x0.value,
			  'y0': out_sersic.y0.value,
			  'e1': out_sersic.e1.value, 
			  'e2': out_sersic.e2.value} 

	# generates a vignette with the model and computes residuals
	if reconstruct_model:
		rec_model = Sersic2D_stamp(amp=out_sersic.amp,radius=out_sersic.radius,sersic_n=out_sersic.sersic_n,
									x0=out_sersic.x0,y0=out_sersic.y0,e1=out_sersic.e1,e2=out_sersic.e2,
									noise=0,zerolevel=out_const.amplitude, shape=B.shape)
		resi = B - rec_model
		return params, rec_model, resi
	else:
		return params	



##################################################################################################################
##################################################################################################################

# We define the models globally. Inside their functions were eating up the memory with each call
fitter = fitting.LevMarLSQFitter()						# Defines the fitter()
model_Gaussian2D = models.custom_model(Gauss2D)			# Creates model object with a gaussian
model_Moffat2D   = models.custom_model(Moffat2D)		# Creates model object with a moffat
model_Sersic2D   = models.custom_model(Sersic2D)		# Creates model object with a sersic

##################################################################################################################
##################################################################################################################

# Ejemplo ------
'''
plt.ion()
Nx = 31
Ny = 31

e=6.6e-2
theta = -1.20689592#np.deg2rad(120)
e1 = e*np.cos(2*theta)
e2 = e*np.sin(2*theta)

B=Gauss2D_stamp(x0=Nx/2,y0=Ny/2,e1=e1,e2=e2,amp=1.0,ab=6.5,noise=0.,zerolevel=100,shape=(Nx,Ny))
#p, mod, res = FitGaussian2D(B, reconstruct_model=True)
#print p
#print 'th',np.rad2deg(np.arctan2(p['e2'],p['e1'])/2.),',	e ', (p['e1']**2+p['e2']**2)**0.5

#B=Moffat2D_stamp(x0=Nx/2,y0=Ny/2,e1=e1,e2=e2,fwhm=4.1284675291,beta=6.1991757010,noise=0.,zerolevel=0,shape=(Nx,Ny))
p, mod, res = FitMoffat2D(B, reconstruct_model=True)
print p
th=(np.arctan2(p['e2'],p['e1'])/2.)

th = np.pi/2.0 - th
print e*np.cos(2.*th), e*np.sin(2*th)

print 'th',np.rad2deg(th),',	e ', (p['e1']**2+p['e2']**2)**0.5

#B=Sersic2D_stamp(amp=1,x0=Nx/2,y0=Ny/2,e1=e1,e2=e2,radius=20,sersic_n=4,noise=0.,zerolevel=100,shape=(Nx,Ny))
#p, mod, res = FitSersic2D(B, reconstruct_model=True)
#print p
#print 'th',np.rad2deg(np.arctan2(p['e2'],p['e1'])/2.),',	e ', (p['e1']**2+p['e2']**2)**0.5

# Make the plot-----------------------------------------------------------
plt.ion()
fig = plt.figure(figsize=(12,3))
ax = fig.add_subplot(1, 3, 1)
ims = ax.imshow(B, origin='low') ; fig.colorbar(ims, ax=ax)
#ims = ax.imshow(psf_stamp, origin='low') ; fig.colorbar(ims, ax=ax)
ax.set_title('input')
ax = fig.add_subplot(1, 3, 2)
ims = ax.imshow(mod, origin='low') ; fig.colorbar(ims, ax=ax)
ax.set_title('model')
ax = fig.add_subplot(1, 3, 3)
ims = ax.imshow(res, origin='low') ; fig.colorbar(ims, ax=ax)
ax.set_title(r'residuals')
plt.show()

# Residuals distribution-----------------------------------------------
m = res.mean()
s = res.std()
plt.figure()
plt.title('residuals')
plt.hist(res.flatten(), 15)
plt.axvline(m, label='Mean = '+str(round(m,3)), linestyle='solid', c='k')
plt.axvline(m+s, label='Std = '+str(round(s,3)), linestyle='dashed', c='k')
plt.axvline(m-s, linestyle='dashed', c='k')
plt.legend()
plt.show()
'''
