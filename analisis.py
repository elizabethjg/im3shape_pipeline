def analisis(entrada):

	print 'entrada para la funcion analisis'
	print entrada
	print '-------------------------------'
	iden		= int(entrada[0])
	beta		= entrada[1]
	z_nobj		= entrada[2]
	ALFA0		= entrada[3]
	DELTA0		= entrada[4]
	clase		= entrada[5]	
	corrida		= int(entrada[6])
	seeing		= entrada[7]
	extincion	= entrada[8]
	MAGMIN		= entrada[9]
	MAGMAX		= entrada[10]
	xy0			= [entrada[11],entrada[12]]
	zback		= entrada[13]
	satur		= entrada[14]
	
	#IMAGEN QUE VA A USAR
	image=im_path+'im_04_%02d.fits' % (iden)
	t0 = time.time()
	#print '----------------------------------------------------------------'
	#print '                    RUNNING SExtractor                          ' 

	sex_conf, sex_salida = sex_config_file('second', 'r', corrida, pixsize, seeing, satur)
	call_sex = 'sextractor '+image+' -c '+sex_conf#+' > sex_output'
	print clr.OKBLUE + call_sex	+ clr.ENDC
	#os.system(call_sex)

	#merge_sex_multiband() # aca juntamos los dos filtros de sextractor
	
	#print '----------------------------------------------------------------'
	#print '                     SOURCE SELECTION                           '

	# It also selects objects for PSFEx
	hdul_stars, hdul_gx, sex_salida_mod = star_gx(sex_salida, seeing/pixsize, plot, PSFEx_manual)
	hdul_gx = background_gx(hdul_gx, MAGMIN, MAGMAX, plot) 

	#print '----------------------------------------------------------------'
	#print '                      COMPUTING PSF                             '

	psfex_conf, psfex_salida = psfex_config_file(sex_salida_mod, 'r', corrida)
	call_psfex = 'psfex '+sex_salida_mod+' -c '+psfex_conf#+' > psfex_output'
	print clr.OKBLUE + call_psfex + clr.ENDC
	#os.system(call_psfex)

	#print '----------------------------------------------------------------'
	#print '                     MEASURING SHAPES                           '

	im2_entrada, im2_psf = psfex.compute_psf_4im2shape(psfex_salida, hdul_gx, corrida)
	
	im2_conf, im2_salida = im2_config_file(image, im2_entrada, im2_psf, corrida)	# defalut: chains=2,niter=200	
	call_im2 = './im2shape '+im2_conf#+' > im2_output'
	print clr.OKBLUE + call_im2 + clr.ENDC
	#os.system(call_im2)
	im2_txt = np.loadtxt(im2_salida, skiprows=1)

	#print '----------------------------------------------------------------'
	#print '                 COMPUTE SOME PARAMETERS                        '

	alfa 	= hdul_gx[2].data['ALPHA_J2000']
	delta 	= hdul_gx[2].data['DELTA_J2000']
	posx 	= hdul_gx[2].data['X_IMAGE']
	posy 	= hdul_gx[2].data['Y_IMAGE']
	e1 		= im2_txt[:,18] 
	e2 		= im2_txt[:,19]

	D_ang	 = cd.angular_diameter_distance(z_nobj, z0=0, **cosmo)
	kpcscale = np.deg2rad(D_ang/3600.) * 1000.
	r = np.rad2deg(np.arccos(np.sin(np.deg2rad(delta))*np.sin(np.deg2rad(DELTA0))+np.cos(np.deg2rad(delta))*np.cos(np.deg2rad(DELTA0))*np.cos(np.deg2rad(alfa-ALFA0))))
	r_kpc = r*3600.0*kpcscale
	peso = np.ones(len(alfa))
	et,ex = rot(posx,posy,e1,e2,xy0)
	nax = np.newaxis
	im2_txt = np.hstack((im2_txt, peso[:,nax],
						et[:,nax], ex[:,nax],
						r[:,nax], r_kpc[:,nax]))
	entrada = np.hstack((entrada, kpcscale))

	print '----------------------------------------------------------------'
	print '                     MERGE CATALOGS                             '

	hdul_gx = merge_gx_catalog(hdul_gx, im2_salida, entrada)

	N_gx	= len(hdul_gx.data)
	N_star  = len(hdul_stars[2].data)
	
	print clr.OKGREEN + 'Salida: ',iden, N_gx, N_star, clr.ENDC
	salida = [iden, hdul_gx, N_gx, N_star]

	return salida
