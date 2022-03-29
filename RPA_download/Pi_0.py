
#Compute Pi_0(omega,q), have integrate omega analytically
#haoting xu, Feb, 2022
#--------------------------------------------------------

import csv
import numpy as np
import libtetrabz
#--------------------------------------------------------
# preset parameter, assuming lattice spacing a=1
I = 1j            #imaginary unit
mu =-1      #chemical potential/fermi energy
t=1               #hopping energy, in the unit of eV
k_up = np.pi  #boundary of BZ
T=40             #Temprature of the system, K
kb=8.617e-5       #Boltzman constant, in the unit of eV/K
u=1           #electron interaction energy, in the unit of eV
N = 100 #number of grid in each direction
N_omega = 123 #number of data point of omega
b1 = [1,0,0]
b2 = [0,1,0]
b3 = [0,0,1]
bvec = np.array([b1,b2,b3]) #lattice vectors
data_folder = "March_16/"
#--------------------------------------------------------

#from numba import jit
#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
from scipy.integrate import dblquad
import numpy as np 
#import warnings
#warnings.filterwarnings("ignore")
#--------------------------------------------------------

#some definitions of the model
def beta(T):
    return 1/(kb*T)

def f(z):  #fermi-dirac distribution
    if z>0:
        return np.exp(-beta(T)*z)/(1+np.exp(-beta(T)*z))
    if z<0:
        return 1/(1+np.exp(beta(T)*z))
    if z==0:
    	return 0.5
    
def inte(func,wght):   #multivariable integral done by simpson's method, func is an array
    inte = 0
    for i in range(0,N):
        for j in range(0,N):                
            inte = inte+wght[i,j]*func[i,j] #i,j are the label of kx, ky
    return inte

    
def E(kx,ky):  #energy of Hubbard model 
    return -t*(np.cos(kx)+np.cos(ky))-mu

def integrand(omega,Ep,Eqp,T): 
    integrand =np.empty([N,N],dtype = np.cfloat)
    for i in range(0,N):
        for j in range(0,N):
            E1 = Eqp[i,j]
            E2 = Ep[i,j]
            integrand[i,j] = (beta(T))*(u**2)*(f(E2)-f(E1))/(I*omega+E2-E1)
    return integrand

    
def Pi_02(qx,qy,omega,T):  #integrate over px,py to give the final Pi0
    om = np.array([omega])
    px = np.linspace(-k_up,k_up,N)
    py = np.linspace(-k_up,k_up,N)
    Px, Py = np.meshgrid(px,py)
    Ep = E(Px,Py)
    qpx = qx+px
    qpy = qy+py
    QPx,QPy = np.meshgrid(qpx,qpy)
    Eqp = E(QPx,QPy)
    
    Ep_3D = np.empty([N, N, 1, 1], dtype=np.float_) #numpy.empty(grid in x,y,z,number of band)
    Eqp_3D = np.empty([N, N, 1, 1], dtype=np.float_)
    for i in range(0,N):
        for j in range(0,N):
            Ep_3D[i,j,0,0] = Ep[i,j]
            Eqp_3D[i,j,0,0] = Eqp[i,j]
    
    wght_5D = libtetrabz.polcmplx(bvec, Ep_3D, Eqp_3D, I*om)
    wght = wght_5D[:,:,0,0,0]
    func = integrand(omega,Ep,Eqp,T)
    return inte(func,wght)
#--------------------------------------------------------


omega1=np.linspace(-10,10,N_omega)      #determine the omega range

def Pi(qx,qy,T):  #compute Pi as an array when input is an array
    Pi_r = np.zeros(np.size(omega1))
    Pi_i = np.zeros(np.size(omega1))
    for i in range(0,np.size(Pi_i)):
        Pi_r[i] = np.real(Pi_02(qx,qy,omega1[i],T))[0]
        Pi_i[i] = np.imag(Pi_02(qx,qy,omega1[i],T))[0]
    return Pi_r,Pi_i

def Pi_om_in_file(qqx):
	Pi_r,Pi_i = Pi(qqx,0,T)
	file = open(data_folder+"qx="+str(round(qqx,2))+",real.csv", "w")
	writer = csv.writer(file)

	for w in range(np.size(omega1)): #from 0-3
		writer.writerow([omega1[w], Pi_r[w]])
	file.close()

	file = open(data_folder+"qx="+str(round(qqx,2))+",imag.csv", "w")
	writer = csv.writer(file)
	for w in range(np.size(omega1)): #from 0-3
		writer.writerow([omega1[w], Pi_i[w]])

	file.close()

Pi_om_in_file(np.pi)
Pi_om_in_file(np.pi/4)
Pi_om_in_file(-np.pi/4)
Pi_om_in_file(np.pi/3)
