import numpy as np 

lam=0.03
delta=1

l=len(np.arange(-3,3,0.01))
spectra_real=np.zeros((l,l))
spectra_imag=np.zeros((l,l))
i=0
j=0
for w_t in np.arange(-3,3,0.01):
	i=0
	for w_tt in np.arange(-3,3,0.01):
		part1=2*(1/(lam+w_tt*1j))
		part2=(1/(lam+(w_tt-delta)*1j)+1/(lam+(w_tt+delta)*1j))
		part3=(1/(lam+(w_t-delta)*1j)-1/(lam+(w_t+delta)*1j))
		full_part=(part1+part2)*part3*1j
		spectra_real[i,j]=abs(full_part) 
		spectra_imag[i,j]=full_part.imag 
		i=i+1
	j=j+1

import matplotlib.pyplot as plt
plt.figure()
im=plt.imshow(spectra_real,cmap="jet",origin="lower")
plt.colorbar(im)
plt.xlabel(r"$\omega_t$")
plt.ylabel(r"$\omega_\tau$")
#im.set_clim(-0.02,0.02)
plt.xticks(np.arange(0,700,100),np.arange(-3,4,1))
plt.yticks(np.arange(0,700,100),np.arange(-3,4,1))
plt.figure()
im=plt.imshow(spectra_imag,cmap="jet",origin="lower")
plt.colorbar(im)
plt.xlabel(r"$\omega_t$")
plt.ylabel(r"$\omega_\tau$")
#im.set_clim(-0.03,0.03)
plt.xticks(np.arange(0,700,100),np.arange(-3,4,1))
plt.yticks(np.arange(0,700,100),np.arange(-3,4,1))
plt.show()