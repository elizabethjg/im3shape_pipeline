import sys
import numpy as np
from cosmolopy import *
from scipy.optimize import curve_fit
from scipy import integrate
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}

cvel=299792458;   # Speed of light (m.s-1)
G= 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc= 3.085678e16; # 1 pc (m)
Msun=1.989e30 # Solar mass (kg)

def chi_red(ajuste,data,err,gl):
	from scipy import integrate
	BIN=len(data)
	chi=((((ajuste-data)**2)/(err**2)).sum())/float(BIN-1-gl)
	return chi


def delta_fi(fi):
	return 1./np.sqrt(np.cos(fi)**2+(f**2)*np.sin(fi)**2)


def esis_profile_sigma(RDfi,f,disp=250.):
	R,fi = RDfi
	Rm=R*1.e6*pc
	R0 = ((disp*1.e3)**2)/G
	b = (Rm/R0)*np.sqrt((f**2)*np.cos(fi)**2+np.sin(fi)**2)
	D_Sigma = (np.sqrt(f)/(2.*b))*(pc**2/Msun)
	Rout = R*np.sqrt((f**2)*np.cos(fi)**2+np.sin(fi)**2)
	return Rout,D_Sigma
	
		
def e_SIS_stack_fit(R,fi,D_Sigma,err):	

	e_SIS_out=curve_fit(esis_profile_sigma2,(R,fi),D_Sigma,sigma=err,absolute_sigma=True)
	pcov=e_SIS_out[1]
	perr = np.sqrt(np.diag(pcov))
	e_f=perr[0]
	e_disp=perr[1]
	f=e_SIS_out[0][0]
	disp=e_SIS_out[0][1]		




def esis_profile_sigma_mod_per(R,fi,f,disp=250.):
	Rm = R*1.e6*pc
	R0 = ((disp*1.e3)**2)/G
	x2 = lambda fi: 1./np.sqrt((f**2)*np.cos(fi)**2+np.sin(fi)**2)
	integral1 = integrate.quad(x2, 0.25*np.pi, 0.75*np.pi)[0]
	integral2 = integrate.quad(x2, 1.25*np.pi, 1.75*np.pi)[0]
	D_Sigma = (R0/Rm)*(np.sqrt(f)/np.pi)*(integral1+integral2)
	Rout = R
	return Rout,D_Sigma	
	
	
def esis_profile_sigma_mod_par(R,fi,f,disp=250.):
	Rm = R*1.e6*pc
	R0 = ((disp*1.e3)**2)/G
	x2 = lambda fi: 1./np.sqrt((f**2)*np.cos(fi)**2+np.sin(fi)**2)
	integral1 = integrate.quad(x2, 0.*np.pi, 0.25*np.pi)[0]
	integral2 = integrate.quad(x2, 0.75*np.pi, 1.25*np.pi)[0]
	integral3 = integrate.quad(x2, 1.75*np.pi, 2.*np.pi)[0]
	D_Sigma = (R0/Rm)*(np.sqrt(f)/np.pi)*(integral1+integral2+integral3)
	Rout = R
	return Rout,D_Sigma	


def SIS_stack_fit(R,D_Sigma,err):
	
	
	# R en Mpc, D_Sigma M_Sun/pc2
	def sis_profile_sigma(R,sigma):
		Rm=R*1.e6*pc
		return (((sigma*1.e3)**2)/(2.*G*Rm))*(pc**2/Msun)
	

	sigma,err_sigma_cuad=curve_fit(sis_profile_sigma,R,D_Sigma,sigma=err,absolute_sigma=True)
	
	ajuste=sis_profile_sigma(R,sigma)
	
	chired=chi_red(ajuste,D_Sigma,err,1)	
	
	xplot=np.arange(0.001,R.max()+1.,0.001)
	yplot=sis_profile_sigma(xplot,sigma)
	
	return sigma[0],np.sqrt(err_sigma_cuad)[0][0],chired,xplot,yplot
	
def NFW_stack_fit(R,D_Sigma,err,z,roc):
	# R en Mpc, D_Sigma M_Sun/pc2
	#Ecuacion 15 (g(x)/2)
	roc_mpc=roc*((pc*1.0e6)**3.0)

	def NFW_profile(R,R200):
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*0.7
		c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
		
		####################################################
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=(R*c)/R200
		m1=x< 1.0
		atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
		jota=np.zeros(len(x))
		jota[m1]=(4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
			+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
			+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))
		m2=x> 1.0     
		atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
		jota[m2]=(4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
			+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
			+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
		m3=(x == 1.0)
		jota[m3]=2.0*np.log(0.5)+5.0/3.0
		rs_m=(R200*1.e6*pc)/c
		kapak=(2.*rs_m*deltac*roc)*(pc**2/Msun)
		return kapak*jota

	NFW_out=curve_fit(NFW_profile,R,D_Sigma,sigma=err,absolute_sigma=True)
	e_R200=np.sqrt(NFW_out[1][0][0])
	R200=NFW_out[0][0]
	
	ajuste=NFW_profile(R,R200)
	
	chired=chi_red(ajuste,D_Sigma,err,1)	

	xplot=np.arange(0.001,R.max()+1.,0.001)
	yplot=NFW_profile(xplot,R200)

	#calculo de c usando la relacion de Duffy et al 2008
	M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*0.7
	c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
	####################################################
	
	return R200,e_R200,chired,xplot,yplot,c

def SIS_fit(R,shear,err,beta,zlens):
	# R en Mpc
	D_ang=cd.angular_diameter_distance(zlens, z0=0, **cosmo)

	def sis_profile_sigma(R,sigma):
			#parameters
		return (4.*np.pi*(sigma**2)*beta*D_ang)/(((cvel/1000.)**2.)*2.*R)

	sigma,err_sigma_cuad=curve_fit(sis_profile_sigma,R,shear,sigma=err,absolute_sigma=True)

	ajuste=sis_profile_sigma(R,sigma)
	
	chired=chi_red(ajuste,shear,err,1)	
	
	xplot=np.arange(0.001,R.max()+1.,0.001)
	yplot=sis_profile_sigma(xplot,sigma)
	
	return sigma[0],np.sqrt(err_sigma_cuad)[0][0],chired,xplot,yplot
	


	
def NFW_fit_c(R,shear,err,roc,zlens,sigmac):
	#Ecuacion 15 (g(x)/2)
	
	H=cd.hubble_z(zlens,**cosmo) #H at z_cluster s-1
	roc=(3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_cluster (kg.m-3)

	#Ecuacion 15 (g(x)/2)
	# R en Mpc

	def NFW_profile(R,R200,c):
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=(R*c)/R200
		m1=x< 1.0
		atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
		jota=np.zeros(len(x))
		jota[m1]=(4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
			+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
			+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))
		m2=x> 1.0     
		atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
		jota[m2]=(4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
			+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
			+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
		m3=(x == 1.0)
		jota[m3]=2.0*np.log(0.5)+5.0/3.0
		rs_m=(R200*1.e6*pc)/c
		kapak=(2.*rs_m*deltac*roc)/sigmac
		return kapak*jota

	NFW_out=curve_fit(NFW_profile,R,shear,sigma=err,absolute_sigma=True)
	pcov=NFW_out[1]
	perr = np.sqrt(np.diag(pcov))
	e_R200=perr[0]
	e_c=perr[1]
	R200=NFW_out[0][0]
	c=NFW_out[0][1]
	
	
	ajuste=NFW_profile(R,R200,c)
	
	chired=chi_red(ajuste,shear,err,2)	

	xplot=np.arange(0.001,R.max()+1.,0.001)
	yplot=NFW_profile(xplot,R200,c)

	
	return R200,e_R200,c,e_c,chired,xplot,yplot


def NFW_fit(R,shear,err,zlens,sigmac):
	#Ecuacion 15 (g(x)/2)
	
	H=cd.hubble_z(zlens,**cosmo) #H at z_cluster s-1
	roc=(3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_cluster (kg.m-3)
	roc_mpc=roc*((pc*1.0e6)**3.0)

	def NFW_profile(R,R200):
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*0.7
		c=5.71*((M/2.e12)**-0.084)*((1.+zlens)**-0.47)
		####################################################
	
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=(R*c)/R200
		m1=x< 1.0
		atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
		jota=np.zeros(len(x))
		jota[m1]=(4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
			+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
			+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))
		m2=x> 1.0     
		atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
		jota[m2]=(4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
			+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
			+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
		m3=(x == 1.0)
		jota[m3]=2.0*np.log(0.5)+5.0/3.0
		rs_m=(R200*1.e6*pc)/c
		kapak=(2.*rs_m*deltac*roc)/sigmac
		return kapak*jota

	NFW_out=curve_fit(NFW_profile,R,shear,sigma=err,absolute_sigma=True)
	e_R200=np.sqrt(NFW_out[1][0][0])
	R200=NFW_out[0][0]
	ajuste=NFW_profile(R,R200)
	
	chired=chi_red(ajuste,shear,err,1)	

	xplot=np.arange(0.001,R.max()+1.,0.001)
	yplot=NFW_profile(xplot,R200)
	
	
	#calculo de c usando la relacion de Duffy et al 2008
	M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*0.7
	c=5.71*((M/2.e12)**-0.084)*((1.+zlens)**-0.47)
	####################################################
	
	
	return R200,e_R200,chired,xplot,yplot,c
	
def shear_map(x,y,e,theta,npix):
	stepx=(x.max()-x.min())/npix
	stepy=(y.max()-y.min())/npix
	xbin=np.zeros(npix**2,float)
	ybin=np.zeros(npix**2,float)
	ex=np.zeros(npix**2,float)
	ey=np.zeros(npix**2,float)
	ngx=np.zeros(npix**2,float)
	#~ plt.plot(x,y,'k.')
	inx=x.min()
	
	ind=0
	print len(ex)
	for j in range(npix):
		#~ plt.plot(x,y,'k.')
		maskx=(x>inx)*(x<(inx+stepx))
		#~ plt.plot(x[maskx],y[maskx],'r.')
		iny=y.min()
		for i in range(npix):
			masky=(y[maskx]>iny)*(y[maskx]<(iny+stepy))
			#~ plt.plot(x[maskx][masky],y[maskx][masky],'b.')
			#~ print ind,len(e[maskx][masky]),iny,iny+stepy
			ex[ind]=e[maskx][masky].mean()*np.cos(theta[maskx][masky].mean())
			ey[ind]=e[maskx][masky].mean()*np.sin(theta[maskx][masky].mean())
			xbin[ind]=x[maskx][masky].mean()
			ybin[ind]=y[maskx][masky].mean()
			ngx[ind]=len(y[maskx][masky])
			ind=ind+1
			iny=iny+stepy
		inx=inx+stepx	
		#~ plt.show()
	return xbin,ybin,ex,ey,ngx
	
def shear_map2(x,y,e1,e2,npix):
	stepx=(x.max()-x.min())/npix
	stepy=(y.max()-y.min())/npix
	xbin=np.zeros(npix**2,float)
	ybin=np.zeros(npix**2,float)
	ex=np.zeros(npix**2,float)
	ey=np.zeros(npix**2,float)
	ngx=np.zeros(npix**2,float)
	#~ plt.plot(x,y,'k.')
	inx=x.min()
	
	ind=0
	print len(ex)
	for j in range(npix):
		#~ plt.plot(x,y,'k.')
		maskx=(x>inx)*(x<(inx+stepx))
		#~ plt.plot(x[maskx],y[maskx],'r.')
		iny=y.min()
		for i in range(npix):
			masky=(y[maskx]>iny)*(y[maskx]<(iny+stepy))
			#~ plt.plot(x[maskx][masky],y[maskx][masky],'b.')
			#~ print ind,len(e[maskx][masky]),iny,iny+stepy
			ex[ind]=e1[maskx][masky].mean()
			ey[ind]=e2[maskx][masky].mean()
			xbin[ind]=x[maskx][masky].mean()
			ybin[ind]=y[maskx][masky].mean()
			ngx[ind]=len(y[maskx][masky])
			ind=ind+1
			iny=iny+stepy
		inx=inx+stepx	
		#~ plt.show()
	return xbin,ybin,ex,ey,ngx


def Sigma_NFW(R,R200,c,zlens):
	# Ecuacion 11 Wright y Brainerd 2000
	# x= R/Rs
	deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
	x=(R*c)/R200
	RS=(R200/c)*1.e6
	H=cd.hubble_z(zlens,**cosmo) #H at z_cluster s-1
	roc=(3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_cluster (kg.m-3)
	roc_pc=roc*(pc**3.0)/Msun
	if x < 1.:
		factor=(2.*RS*deltac*roc_pc)/(x**2-1.)
		sigma_nfw=factor*(1.-(2./np.sqrt(1.-x**2))*np.arctanh(np.sqrt((1.-x)/(1.+x))))
	elif x == 1.:
		factor=(2.*RS*deltac*roc_pc)/3.
		sigma_nfw=factor
	else:
		factor=(2.*RS*deltac*roc_pc)/(x**2-1.)
		sigma_nfw=factor*(1.-(2./np.sqrt(x**2-1.))*np.arctan(np.sqrt((x-1.)/(1.+x))))
	#~ print RS, deltac,roc_pc,factor,(2.*RS*deltac*roc_mpc)/3.
	return sigma_nfw
