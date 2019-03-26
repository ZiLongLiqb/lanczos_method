import numpy as np 
from numba import jit

gamma=0.06
omega=np.sqrt(1-gamma**2)
l=len(np.arange(-6,6,0.01))
result_real=np.zeros((l,l))
result_imag=np.zeros((l,l))

def accu_cal():
	global result_real
	global result_imag
	i=0
	j=0
	for w_tt in np.arange(-6,6,0.01):
		j=0
		for w_t in np.arange(-6,6,0.01):
			part1=1/(3*gamma+(w_t-omega)*1j)*1/(gamma+(w_tt-omega)*1j)
			part2=-1/(3*gamma+(w_t+omega)*1j)*1/(gamma+(w_tt+omega)*1j)
			part3=-1/(3*gamma+(w_t-3*omega)*1j)*1/(gamma+(w_tt-omega)*1j)*1/2
			part4=1/(3*gamma+(w_t+3*omega)*1j)*1/(gamma+(w_tt+omega)*1j)*1/2
			part5=-1/(3*gamma+(w_t+omega)*1j)*1/(gamma+(w_tt-omega)*1j)*1/2
			part6=1/(3*gamma+(w_t-omega)*1j)*1/(gamma+(w_tt+omega)*1j)*1/2
			part7=1/(3*gamma+(w_t-omega)*1j)*1/(2*gamma+(w_tt)*1j)
			part8=-1/(3*gamma+(w_t+omega)*1j)*1/(2*gamma+(w_tt)*1j)
			part9=-1/(3*gamma+(w_t-3*omega)*1j)*1/(2*gamma+(w_tt-2*omega)*1j)*1/2
			part10=1/(3*gamma+(w_t+3*omega)*1j)*1/(2*gamma+(w_tt+2*omega)*1j)*1/2
			part11=-1/(3*gamma+(w_t+omega)*1j)*1/(2*gamma+(w_tt+2*omega)*1j)*1/2
			part12=1/(3*gamma+(w_t-omega)*1j)*1/(2*gamma+(w_tt-2*omega)*1j)*1/2
			part13=(1/(gamma-(omega-w_t)*1j)-1/(gamma+(omega+w_t)*1j))*1/(2*omega)
			full_part=(part1+part2+part3+part4+part5+part6+part7+part8+part9+part10+part11+part12)*part13/(8*omega**3)*0.08**3
			#full_part=(part1+part2+part3+part4+part5+part6)
			result_real[i,j]=abs(full_part)
			result_imag[i,j]=full_part.imag 
			j=j+1
		i=i+1

accu_cal()
import matplotlib.pyplot as plt
plt.figure()
im=plt.imshow(result_real,cmap="jet",origin="lower")
plt.colorbar(im)
plt.xlabel(r"$\omega_t$")
plt.ylabel(r"$\omega_\tau$")
#im.set_clim(-0.02,0.02)
plt.xticks(np.arange(0,1300,100),np.arange(-6,7,1))
plt.yticks(np.arange(0,1300,100),np.arange(-6,7,1))
plt.figure()
im=plt.imshow(result_imag,cmap="jet",origin="lower")
plt.colorbar(im)
plt.xlabel(r"$\omega_t$")
plt.ylabel(r"$\omega_\tau$")
im.set_clim(-0.03,0.03)
plt.xticks(np.arange(0,1300,100),np.arange(-6,7,1))
plt.yticks(np.arange(0,1300,100),np.arange(-6,7,1))
plt.show()