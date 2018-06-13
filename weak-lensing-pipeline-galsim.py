import math
import numpy as np
from pylab import *
import cosmolopy.distance as cd
from multiprocessing import Pool
from multiprocessing import Process
import os
import sys
from star_gx_des import *
from rot import *
from psffile import *
#from selectstars import *
#from astropy.wcs import WCS
from seeing import *
from pipeio import *
import time
from psfex_pipe import *
from analisis import *
from astropy.io import fits
from other_filters import *
from im3_paralelo import *

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}



ok_flag   	 = 'no'	                #raw_input('Desea continuar ejecucion (S/N) ')
plot 		 = 'no'	                #raw_input('Graficar (S/N) ')
archivo_in	 = 'data_obj.cat'		#int(raw_input('Emepzar en imagen (0 a n) '))
proc  		 = '001'	                #raw_input('Maquina ')
simult 		 = 'si'	                #raw_input('Correr filtros en simultaneo ')
PSFEx_manual = True
corrida      = 0


entrada=np.loadtxt(archivo_in,comments='#',dtype='str')
print entrada

print 'entrada para la funcion analisis'
print entrada
print '-------------------------------'


ALFA0		= np.float(entrada[1])
DELTA0		= np.float(entrada[2])
z_nobj		= np.float(entrada[3])
pixsize     = np.float(entrada[4])

# IMAGEN QUE VA A USAR

filtro    = entrada[5]
filtros   = filtro
image     = entrada[6]
zeropoint = np.float(entrada[7])

#hdu = fits.open(image)
gain = 1.0#str(hdu[0].header['GAIN'])


#print '----------------------------------------------------------------'
#print '                    DETERMINING SEEING                          ' 

SEEING, SATUR=seeing_func(image,pixsize,zeropoint,gain,corrida,filtro,100.,-100.,10.,plot)


#print '----------------------------------------------------------------'
#print '             RUNNING SExtractor for stars                       ' 

sex_conf, sex_salida = sex_config_file('second', filtro, corrida, pixsize, zeropoint, gain, SEEING, SATUR)
call_sex = 'sextractor '+image+' -c '+sex_conf#+' > sex_output'
print clr.OKBLUE + call_sex	+ clr.ENDC
os.system(call_sex)


#merge_sex_multiband() # aca juntamos los dos filtros de sextractor

#print '----------------------------------------------------------------'
#print '                     SOURCE SELECTION                           '

# It also selects objects for PSFEx
hdul_stars, hdul_gx, sex_salida_mod = star_gx(sex_salida, SEEING/pixsize, plot, PSFEx_manual)


#print '----------------------------------------------------------------'
#print '                      COMPUTING PSF                             '

psfex_conf, psfex_salida = psfex_config_file(sex_salida_mod, filtro, corrida)
call_psfex = 'psfex '+sex_salida_mod+' -c '+psfex_conf#+' > psfex_output'
print clr.OKBLUE + call_psfex + clr.ENDC
os.system(call_psfex)
#psfex_salida = './psfex_files/modified_psfex.psf'

#print '----------------------------------------------------------------'
#print '             RUNNING SExtractor for galaxies                    ' 

#sex_conf, sex_salida = sex_config_file('second', filtro, corrida, pixsize, zeropoint, gain, SEEING, SATUR)
#call_sex = 'sextractor gal.004.fits -c '+sex_conf#+' > sex_output'
#print clr.OKBLUE + call_sex	+ clr.ENDC
#os.system(call_sex)
#hdul_gx 	= fits.open(sex_salida)

#print '----------------------------------------------------------------'
#print '                     MEASURING SHAPES                           '


input_cat=np.loadtxt('gs_gal_004.cat')

ids = input_cat[:,0]
x   = input_cat[:,1]
y   = input_cat[:,2]


im3_salida = im3call(['gal.004.fits', psfex_salida, ids, x, y, corrida])

	
print 'END OF THE PROGRAM :)'
