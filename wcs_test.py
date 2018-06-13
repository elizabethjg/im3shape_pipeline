import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from pylab import *


class my_WCS:

	def __init__(self, image):
		
		hdul = fits.open(image)
		header = hdul[0].header
		self.r1 = header['CRPIX1']
		self.r2 = header['CRPIX2']
		self.CD11 = header['CD1_1']
		self.CD12 = header['CD1_2']
		self.CD21 = header['CD2_1']
		self.CD22 = header['CD2_2']
		self.ra_p = np.deg2rad(header['CRVAL1'])
		self.dec_p = np.deg2rad(header['CRVAL2'])
		self.phi_p = np.deg2rad(180)					#Native longitude of celestial pole (0 or 180)

	def pix2world(self,px,py):

		# Projection plane coordinates
		x = np.deg2rad(self.CD11 * (px - self.r1) + self.CD12 * (py - self.r2))
		y = np.deg2rad(self.CD21 * (px - self.r1) + self.CD22 * (py - self.r2))

		# Native spherical coordinates
		phi = np.arctan2(x,-y)
		R_th = np.sqrt(x*x + y*y)
		theta = np.arctan2(1, R_th)

		# Celestial spherical coordinates
		sin_th = np.sin(theta)
		cos_th = np.cos(theta)
		sin_dec_p = np.sin(self.dec_p)
		cos_dec_p = np.cos(self.dec_p)
		sin_dphi = np.sin(phi-self.phi_p)
		cos_dphi = np.cos(phi-self.phi_p)

		ra  = self.ra_p + np.arctan2(-cos_th*sin_dphi , sin_th*cos_dec_p - cos_th*sin_dec_p*cos_dphi)
		dec = np.arcsin(sin_th*sin_dec_p + cos_th*cos_dec_p*cos_dphi)
		
		return np.rad2deg(ra), np.rad2deg(dec)

	def get_omega(self,px,py):
		'''
		Returns the angle omega between the x-axis and ra-axis at the position px,py... maybe
		'''
		# Projection plane coordinates
		x = np.deg2rad(self.CD11 * (px - self.r1) + self.CD12 * (py - self.r2))
		y = np.deg2rad(self.CD21 * (px - self.r1) + self.CD22 * (py - self.r2))

		# Native spherical coordinates
		phi = np.arctan2(x,-y)
		R_th = np.sqrt(x*x + y*y)
		theta = np.arctan2(1, R_th)

		# Celestial spherical coordinates
		sin_th = np.sin(theta)
		cos_th = np.cos(theta)
		sin_dec_p = np.sin(self.dec_p)
		cos_dec_p = np.cos(self.dec_p)
		sin_dphi = np.sin(phi-self.phi_p)
		cos_dphi = np.cos(phi-self.phi_p)

		# Derivatives		
		x_px = self.CD11
		y_px = self.CD21
		x_py = self.CD12
		y_py = self.CD22

		theta_xy_denominator = (R_th + R_th**3)
		theta_x = -x / theta_xy_denominator
		theta_y = -y / theta_xy_denominator

		phi_xy_denominator = R_th**2
		phi_x = -y / phi_xy_denominator
		phi_y =  x / phi_xy_denominator

		ra_thetaphi_denominator = (sin_th*cos_dec_p - cos_th*sin_dec_p*cos_dphi)**2 + (cos_th*sin_dphi)**2
		ra_theta = -cos_dec_p*sin_dphi / ra_thetaphi_denominator
		ra_phi   = (cos_th*sin_th*cos_dec_p*cos_dphi - sin_dphi*cos_th**2) / ra_thetaphi_denominator

		term_x = (ra_phi*phi_x + ra_theta*theta_x)
		term_y = (ra_phi*phi_y + ra_theta*theta_y)
		ra_px = term_x * x_px + term_y * y_px
		ra_py = term_x * x_py + term_y * y_py

		#return np.rad2deg(ra_px), np.rad2deg(ra_py)
		omega = np.arctan2(ra_py,ra_px)
		#3er y 4to cuadrante
		m_omega = (ra_py<0.)
		omega[m_omega] = np.pi - abs(omega[m_omega])

		return omega


'''
#image = 'cfht_01.fits'
image = 'sdss_01.fits' ; shape = (2048,1489)
posx = np.array([100,500,1500, 1020.25724925, 369])
posy = np.array([200,650,1500., 982.67858193626, 1340])

w = WCS(image)
w_ra,w_dec=w.all_pix2world(posx,posy,1)

my_w = my_WCS(image)
my_ra,my_dec=my_w.pix2world(posx,posy)

print '	Ref		My'
for i in xrange(len(posx)):
	print w_ra[i],'	', my_ra[i],'	|	',w_dec[i],'	', my_dec[i]

xx,yy = np.meshgrid(np.arange(0,shape[0],10),np.arange(0,shape[1],10), indexing='ij')

ra,dec = my_w.pix2world(xx.flatten(), yy.flatten())

#ra_px, ra_py = my_w.get_omega(xx, yy)
omega = my_w.get_omega(xx, yy)

#plt.quiver(xx,yy,ra_px,ra_py,cmap=cm.Reds,headwidth=1,headlength=0,units='xy')#,scale=10)
#plt.show()

plt.imshow(omega,origin='low')
plt.colorbar()
plt.show()
'''
