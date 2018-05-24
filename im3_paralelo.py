import os
import sys
import numpy as np
from astropy.io import fits
from seeing import *
from pipeio import *
from psfex import *
from multiprocessing import Pool
from multiprocessing import Process
from psfex_pipe import *
from py3shape.options import Options

def im3call(entrada):

	image, psfex_salida, ids, x, y, corrida = entrada
	
	print ids[:10]

	#psf_input = 'moffat_catalog'
	#psf_input = 'no_psf'
	psf_input = 'psf_image_cube'
	im3_conf = im3_config_file('37', '500', psf_input, corrida)
	
	options = Options(im3_conf)
	
	im3_entrada, im3_psf = compute_psf_4im3shape(psfex_salida, options, ids, x, y, corrida)
	
	
	im3_salida = 'im3_'+str(corrida)+'.out'

	os.system('rm '+im3_salida)

	#call_im3 = 'python ~/IM3SHAPE/bin/im3shape.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	#call_im3 = '~/IM3SHAPE/bin/im3shape '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	#call_im3 = 'python ~/IM3SHAPE/bin/im3shape.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 5'#+str(len(ids)-1)
	call_im3 = 'python ~/IM3SHAPE/bin/im3shape.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	print clr.OKBLUE + call_im3 + clr.ENDC
	os.system(call_im3)
	
	
	im3out = np.loadtxt(im3_salida)
	
	return im3out

def im3paralelo(image, psfex_salida, hdul_gx, cores):
	

	
	ids = hdul_gx[2].data['NUMBER']
	x   = hdul_gx[2].data['X_IMAGE']
	y   = hdul_gx[2].data['Y_IMAGE']


	ngxcat = len(ids)
		
	step = 0
	step0=int(ngxcat/cores)+1
	
	entrada=[]
	
	for j in range(cores):

		if j==(cores-1):
			arreglo = [image, psfex_salida, ids[step:], x[step:], y[step:], j]
			entrada.append(arreglo)
		else: 
			arreglo = [image, psfex_salida, ids[step:step+step0], x[step:step+step0], y[step:step+step0], j]
			entrada.append(arreglo)
		
		step=step+step0

	
	
	
	pool = Pool(processes=(cores))
	salida=np.array(pool.map(im3call, entrada))
	pool.terminate()
	
	mask=np.array(map(len,salida))>0

	
	salida = salida[mask]
		
	im3_out=salida[0]	
	
	for j in salida[1:]:
		im3_out=np.concatenate((im3_out,j),axis = 0)

	

	return im3_out
