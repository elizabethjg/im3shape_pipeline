import math
import numpy as np
import os
#import pyfits
import sys

def psffile(goodstars,objcat,vecinas,psfout):

		idgood=goodstars[:,1]
		e1=goodstars[:,18]
		e2=goodstars[:,19]
		amp=goodstars[:,11]
		ab=goodstars[:,10]
		posx= goodstars[:,2]+goodstars[:,6]
		posy= goodstars[:,3]+goodstars[:,7]
		e=(e1**2+e2**2)**0.5
		theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
				
		a=(ab*((1.0+e)/(1.0-e)))**0.5
		a1=a*cos(theta)
		a2=a*sin(theta)
		
		nobj=objcat.shape[0]
		psfobj=np.zeros((nobj,8),float)
		
		#print 'n gx',nobjcat	
		
		if nobj == 1:
			gxposx= np.array([objcat[1]])
			gxposy= np.array([objcat[2]])
		else:
			gxposx= objcat[:,1]
			gxposy= objcat[:,2]

		for n in range(nobj):
			x=gxposx[n]
			y=gxposy[n]
			difpos=((posx-x)**2+(posy-y)**2)**0.5
			cercanas=argsort(difpos)[0:vecinas]

			prome1=e1[cercanas].mean()
			prome2=e2[cercanas].mean()
			promab=ab[cercanas].mean()
			
			thetaprom=((np.arctan2(prome2,prome1)/2.0))+np.pi
			eprom=(prome1**2+prome2**2)**0.5
			psfobj[n,:]=[x,y,0.0,0.0,eprom,thetaprom,promab,1.0]
		
		#print 'Writing the psf file for galaxies in ',psfobjfile
		f1=open(psfout,'w')
		f1.write('1\n')
		f1.write(str(nobj))
		f1.write('\n')
		f1.write('6\n')
		np.savetxt(f1, psfobj, fmt=['%15.10f']*4+['%13.10f']*4)
		f1.close()
		
		return psfout
