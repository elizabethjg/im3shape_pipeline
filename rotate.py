import numpy as np
from astropy.io import fits
from wcs_test import my_WCS 

def rotate(im3salida, hdul_gx, image):
	
	posx   = hdul_gx[2].data['X_IMAGE']
	posy   = hdul_gx[2].data['Y_IMAGE']
	
	e1 = im3salida[:,2]
	e2 = im3salida[:,3]
	
	e=(e1**2+e2**2)**0.5
	theta=(np.arctan2(e2,e1)/2.0)
	theta[(theta<0.)] = np.pi + theta[(theta<0.)]
	
	my_w = my_WCS(image)
	
	omega = my_w.get_omega(posx,posy)	# angulo del eje de ascencion recta en cada x,y
	e1 = e*np.cos(2.*(theta+(np.pi/2.)))
	e2 = e*np.sin(2.*(theta+(np.pi/2.)))
	
	e1_rot = (e1*np.cos(2*omega)+e2*np.sin(2*omega))
	e2_rot = (-e1*np.sin(2*omega)+e2*np.cos(2*omega))
	
	return np.array([e1_rot, e2_rot]).T
