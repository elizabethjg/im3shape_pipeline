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
from py3shape_main import *

def im3call(entrada):

	image, psfex_salida, ids, x, y, cores, corrida = entrada
	
	print ids[:10]

	#psf_input = 'moffat_catalog'
	#psf_input = 'no_psf'
	psf_input = 'psf_image_cube'
	im3_conf = im3_config_file('37', '500', psf_input, cores, corrida)
	
	options = Options(im3_conf)
	
	im3_entrada, im3_psf = compute_psf_4im3shape(psfex_salida, options, ids, x, y, cores)
	
	im3_salida = 'im3_'+str(cores)+'.out'

	os.system('rm '+im3_salida)

	#call_im3 = 'python ~/IM3SHAPE/bin/im3shape.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	#call_im3 = '~/IM3SHAPE/bin/im3shape '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	#call_im3 = 'python ~/IM3SHAPE/bin/im3shape_203.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 2'#+str(len(ids)-1)
	
	entrada_py3shape = [im3_conf,image,im3_entrada,im3_psf,im3_salida,'0','2']
	
	main_203(entrada_py3shape)
	
	#call_im3 = 'python ~/IM3SHAPE/bin/im3shape_203.py '+im3_conf+' '+image+' '+im3_entrada+' '+im3_psf+' '+im3_salida+' 0 '+str(len(ids)-1)
	#print clr.OKBLUE + call_im3 + clr.ENDC
	#os.system(call_im3)
	
	
	im3_tmp = np.loadtxt(im3_salida)
	
	
	if len(ids) != len(im3_tmp):
		print clr.FAIL+'WARNING: SExtractor and Im3Shape catalogs dont have the same number of objects!'+clr.ENDC

	
	
		# CORRELATING ARRAYS -------------------------------------------------------------
				
		#im3id = im3_salida[:,0] # PARA C
		im3id = im3_tmp[:,12] # PARA PYTHON
		nrows = len(ids)
		
		
		mask = np.in1d(ids,im3id)
		
		im3_out = np.zeros((nrows,im3_tmp.shape[1]))
		im3_out[mask,:] = im3_tmp
	
	else: 
		im3_out = im3_tmp
	
	return im3_out

def im3paralelo(image, psfex_salida, hdul_gx, cores,corrida):
	

	
	ids = hdul_gx[2].data['NUMBER']
	x   = hdul_gx[2].data['X_IMAGE']
	y   = hdul_gx[2].data['Y_IMAGE']


	ngxcat = len(ids)
		
	step = 0
	step0=int(round(ngxcat/cores, 0))
	
	entrada=[]
	
	for j in range(cores):

		if j==(cores-1):
			arreglo = [image, psfex_salida, ids[step:], x[step:], y[step:], j,corrida]
			entrada.append(arreglo)
			print len(x[step:])
		else: 
			arreglo = [image, psfex_salida, ids[step:step+step0], x[step:step+step0], y[step:step+step0], j,corrida]
			entrada.append(arreglo)
			print len(x[step:step+step0])
		
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
