import numpy as np 
import time 

#problem solved
#the even case is not that good!
#the final problem is the normalization function!!!!
#Cause now we consider plural

#----------------Global parameter-------------------
N=14
Gamma=10
J=0.9
h=0.1
#----------------Global parameter end---------------

#----------------function area----------------------
def get_base_2(num,N):
	n=num
	result_base_2=[]
	while True:
		if n==1 or n==0:
			result_base_2.append(n)
			break
		m=n%2
		n=n//2
		result_base_2.append(m)
	if len(result_base_2)<N:
		result_base_2=result_base_2+[0]*(N-len(result_base_2))
	result_base_2.reverse()
	return result_base_2

#get the dictionary which can convert a decimal number to binary number
#decimal_to_binary_dict
decimal_to_binary_dict={}
for i in np.arange(2**N):
	decimal_to_binary_dict[i]=get_base_2(i,N)

#get the dictionary which can convert a binary number to decimal
#binary_to_decimal_dict
binary_to_decimal_dict={}
for item in decimal_to_binary_dict:
	binary_to_decimal_dict[tuple(decimal_to_binary_dict[item])]=item

#return a list
def flip(num_list,i):
	new_num_list=[]
	new_num_list=new_num_list+num_list
	new_num_list[-i-1]=int(not new_num_list[-i-1])
	return new_num_list

#return a list
def flip_all(num_list):
	num_list_new=[]
	num_list_new=num_list_new+num_list
	L=len(num_list_new)
	for i in range(L):
		num_list_new=flip(num_list_new,i)
	return num_list_new

#return a array
def normolize(phi):
	sq=np.dot(np.conjugate(phi),phi)
	return np.array([item/(np.sqrt(sq)) for item in phi]),np.sqrt(sq)

def find_min(l):
	a=l[0]
	j=0
	for i in range(len(l)):
		if l[i]<a:
			j=i 
			a=l[i]
	return j,a

def check_rep(num):
	num_group=[num]
	num_binary=decimal_to_binary_dict[num]
	num_binary_unchanged=[]+num_binary
	n=0
	while True:
		n=n+1
		num_binary_next=num_binary[1:]+[num_binary[0]]
		if num_binary_next==num_binary_unchanged:
			break
		else:
			num_group.append(binary_to_decimal_dict[tuple(num_binary_next)])
			num_binary=num_binary_next
	num_group_p=[]
	for m in num_group:
		num_group_p.append(binary_to_decimal_dict[tuple(flip_all(decimal_to_binary_dict[m]))])
	n_1,min_1=find_min(num_group)
	n_2,min_2=find_min(num_group_p)
	if np.min([min_1,min_2])==num:
		is_rep=1
		r=0
		p=0
	elif np.min([min_1,min_2])==binary_to_decimal_dict[tuple(flip_all(decimal_to_binary_dict[num]))]:
		is_rep=0
		r=0
		p=1
	else:
		is_rep=0
		if min_2<min_1:
			p=1
			r=n-n_2
		else:
			p=0
			r=n-n_1
	if min_1==min_2:
		n=n/2
		t=1
	else:
		t=0
	return [is_rep,n,r,p,np.min([min_1,min_2]),t]
 
#test
#get the check_rep_dict
check_rep_dict={}
for i in range(2**N):
	check_rep_dict[i]=check_rep(i)

def find_state(l,mi):
	i_min=0
	i_max=len(l)-1
	while True:
		b=i_min+int((i_max-i_min)/2)
		if mi<l[b]:
			i_max=b-1
		elif mi>l[b]:
			i_min=b+1
		else:
			return b 
		if i_min>i_max:
			return -1


def sub_hoperation(phi,N,J,h,yita,phi_base,k_): 
	phi_next=[]
	for i in np.arange(len(phi)):
		midle=0
		result=[] 
		phi_i=decimal_to_binary_dict[phi_base[i]]
		for k in np.arange(N):
			j=(k+1)%(N)
			if phi_i[k]==phi_i[j]: 
				midle+=-J
			else:
				midle+=J
		midle=midle*phi[i]
		another_item=0
		for j in range(N):
			l=binary_to_decimal_dict[tuple(flip(phi_i,j))]
			is_rep_2,n_2,r_2,p_2,min_2,t_2=check_rep_dict[l]
			is_rep_1,n_1,r_1,p_1,min_1,t_1=check_rep_dict[phi_base[i]]
			l_location=find_state(phi_base,min_2)
			if l_location==-1:
				continue
			another_item+=-h*phi[l_location]*np.sqrt(n_1/n_2)
			#*(np.exp(yita*np.pi*1j)**p_2)*np.exp(-2*np.pi*k_*r_2*1j/N)
		phi_next.append(midle+another_item)
	return np.array(phi_next)

def M_x(phi,phi_base,N=N):
	phi_next=[]
	for i in np.arange(len(phi)):
		phi_i=decimal_to_binary_dict[phi_base[i]]
		item_1=0
		for j in range(N):
			l=binary_to_decimal_dict[tuple(flip(phi_i,j))]
			is_rep_2,n_2,r_2,p_2,min_2,t_2=check_rep_dict[l]
			is_rep_1,n_1,r_1,p_1,min_1,t_1=check_rep_dict[phi_base[i]]
			l_location=find_state(phi_base,min_2)
			item_1+=phi[l_location]*np.sqrt(n_1/n_2)
		phi_next.append(item_1)
	return np.array(phi_next)


def lanczos_method(N,J,h,Gamma,yita,k):
	phi_base=[]
	for i in range(2**N):
		if check_rep_dict[i][0]==1 and (2*check_rep_dict[i][1]*k+check_rep_dict[i][5]*yita*N)%(2*N)==0:
			phi_base.append(i) 
	phi_all=np.zeros((len(phi_base),Gamma+1),dtype=complex)
	n=[]
	a=[]
	phi_0=np.random.random(len(phi_base)) 
	phi_0,n_0=normolize(phi_0)
	phi_all[:,0]=phi_0
	n.append(n_0)
	phi_1=sub_hoperation(phi_0,N,J,h,yita,phi_base,k)
	a_0=np.dot(np.conjugate(phi_0),np.transpose(phi_1))
	a.append(a_0)
	phi_1=phi_1-a_0*phi_0
	phi_1,n_1=normolize(phi_1)
	phi_all[:,1]=phi_1
	n.append(n_1)
	for i in range(Gamma): 
		phi_2=sub_hoperation(phi_1,N,J,h,yita,phi_base,k)
		a.append(np.dot(np.conjugate(phi_1),np.transpose(phi_2)))
		phi_2=phi_2-a[-1]*phi_1-n[-1]*phi_0
		phi_2,n_2=normolize(phi_2)
		if i+2<=Gamma:
			phi_all[:,i+2]=phi_2
		n.append(n_2)
		phi_0,phi_1=phi_1,phi_2	
	return n,a,phi_all,phi_base


def main_operation_2(n,a,Gamma=10):
	H=np.zeros((Gamma+1,Gamma+1),dtype=complex)
	for i in range(Gamma+1):
		for j in range(Gamma+1):
			if i==j:
				H[i,j]=a[i] 
			if j==i+1:
				H[i,j]=n[j]
			if j==i-1:
				H[i,j]=n[i]
	a,b=np.linalg.eig(H)
	return a,b


def main_operation(N,J,h,Gamma,yita=0,k=0):
	n,a,phi_all,phi_base=lanczos_method(N,J,h,Gamma,yita,k)
	H=np.zeros((Gamma+1,Gamma+1),dtype=complex)
	for i in range(Gamma+1):
		for j in range(Gamma+1):
			if i==j:
				H[i,j]=a[i] 
			if j==i+1:
				H[i,j]=n[j]
			if j==i-1:
				H[i,j]=n[i]
	a,b=np.linalg.eig(H)
	return a,b,phi_all,phi_base


def get_phi_0_positive(alpha,J,h,N=N,Gamma=Gamma):
	a,b,phi_all,phi_base=main_operation(N,J,h,Gamma)
	energy_min=np.min(a)
	b_min=np.zeros(Gamma+1,dtype=complex)
	for i in range(len(a)):
		if a[i]==energy_min:
			b_min=b[:,i]
	phi_g_s=np.dot(phi_all,b_min)
	M_x_phi=M_x(phi_g_s,phi_base,N)
	M_x_phi_2=M_x(M_x_phi,phi_base,N)
	phi_0_positive=phi_g_s-alpha*M_x_phi*1j-alpha**2/2*M_x_phi_2
	phi_0_positive,n=normolize(phi_0_positive)
	return phi_0_positive,phi_base

begin_time=time.time()
phi_0_positive,phi_base=get_phi_0_positive(0.1,J,h)
M_x_average=[]
phi_t_0=phi_0_positive
end_time=time.time()
time_inteval=end_time-begin_time
print(time_inteval)
#print(np.dot(np.transpose(np.conjugate(phi_t_0)),sub_hoperation(phi_t_0,N,J,h,0,phi_base,0)))
 
def lanczos_method_with_phi_0(phi_0,phi_base,J,h,N=N,Gamma=10,yita=0,k=0):
	n=[]
	a=[]
	phi_all=np.zeros((len(phi_base),Gamma+1),dtype=complex)
	phi_0=phi_0
	phi_all[:,0]=phi_0
	n_0=1
	n.append(n_0)
	phi_1=sub_hoperation(phi_0,N,J,h,yita,phi_base,k)
	a_0=np.dot(np.conjugate(phi_0),np.transpose(phi_1))
	a.append(a_0)
	phi_1=phi_1-a_0*phi_0
	phi_1,n_1=normolize(phi_1)
	phi_all[:,1]=phi_1
	n.append(n_1)
	for i in range(Gamma): 
		phi_2=sub_hoperation(phi_1,N,J,h,yita,phi_base,k)
		a.append(np.dot(np.conjugate(phi_1),np.transpose(phi_2)))
		phi_2=phi_2-a[-1]*phi_1-n[-1]*phi_0
		phi_2,n_2=normolize(phi_2)
		n.append(n_2)
		if i>0:
			for p in range(i-1):
				q=np.dot(np.transpose(np.conjugate(phi_all[:,p])),phi_2)
				phi_2=phi_2-q*phi_all[:,p]
			phi_2,n_2=normolize(phi_2)
			#phi_2,n_2=normolize(phi_2)
		if i+2<=Gamma:
			phi_all[:,i+2]=phi_2
		phi_0,phi_1=phi_1,phi_2	
	return n,a,phi_all	


def error(p,l,delta_t):
	if p==2:
		return l*delta_t
	else:
		return l*delta_t/(p-1)*error(p-1,l,delta_t)


def get_after_delta_t(phi_t_0,phi_base,delta_t=0.1,J=J,h=h,N=N):
	H_value=np.dot(np.conjugate(np.transpose(phi_t_0)),sub_hoperation(phi_t_0,N,J,h,0,phi_base,0))
	M_x_value=np.dot(np.conjugate(np.transpose(phi_t_0)),M_x(phi_t_0,phi_base))
	n,a,phi_all=lanczos_method_with_phi_0(phi_t_0,phi_base,J,h,Gamma=10)
	diag_before,u=main_operation_2(n,a,Gamma=10)
	diag=np.diag(diag_before)
	phi_new_appearance=np.zeros(11)
	phi_new_appearance[0]=1
	u_dag=(np.transpose(u))
	diag_len=len(diag)
	diag_new=np.zeros((diag_len,diag_len),dtype=complex) 
	for i in range(diag_len):
		diag_new[i,i]=np.exp(-diag[i,i]*delta_t*1j) 
	matrix_1=np.dot(u_dag,np.transpose(phi_new_appearance))  
	matrix_2=np.dot(diag_new,matrix_1)
	matrix_3=np.dot(u,matrix_2)
	phi_after_delta_t=np.dot(phi_all,matrix_3)
	#print(np.dot(np.transpose(np.conjugate(phi_after_delta_t)),sub_hoperation(phi_after_delta_t,12,J,h,0,phi_base,0)))
	e=error(10,max(abs(diag_before)),delta_t)
	#print(np.dot(np.transpose(np.conjugate(phi_all)),phi_all))
	return H_value,M_x_value,e,phi_after_delta_t

#get_after_delta_t(phi_t_0,phi_base)
'''H_diag_num=len(phi_base)
H_diag=np.zeros((H_diag_num,H_diag_num),dtype=complex)
for i in range(H_diag_num):
	i_th=np.zeros(H_diag_num)
	i_th[i]=1
	H_i_th=sub_hoperation(i_th,N,J,h,0,phi_base,0)
	for j in range(H_diag_num):
		j_th=np.zeros(H_diag_num)
		j_th[j]=1
		H_diag[i,j]=np.dot(np.transpose(np.conjugate(j_th)),H_i_th)
H_value,H_diag_vector=np.linalg.eig(H_diag)

M_x_average=[]
for t in np.arange(0,20,0.1):
	phi_now=np.zeros(H_diag_num,dtype=complex)
	for k in range(len(H_value)):
		phi_now=phi_now+(np.exp(-H_value[k]*t*1j)*np.dot(np.transpose(np.conjugate(H_diag_vector[:,k])),phi_t_0))*H_diag_vector[:,k]
	M_x_average.append(np.dot(np.transpose(np.conjugate(phi_now)),M_x(phi_now,phi_base)))

import pickle
f=open("H_x_average_1.txt","wb")
pickle.dump(M_x_average,f)
f.close()'''
'''import matplotlib.pyplot as plt 
plt.figure()
plt.scatter(np.arange(0,20,0.1),M_x_average,s=5)
plt.figure()
plt.plot(np.fft.fft(M_x_average).imag)
plt.show()'''

Energy=[]
Error=[]
for i in range(100):
	H_value,M_x_value,e,phi_t_0=get_after_delta_t(phi_t_0,phi_base)
	M_x_average.append(M_x_value)
	Energy.append(H_value)
	Error.append(e)

import pickle
f=open("Energy.txt","wb")
pickle.dump(Energy,f)
f.close()
f=open("Error.txt","wb")
pickle.dump(Error,f)
f.close()
f=open("H_x_average.txt","wb")
pickle.dump(M_x_average,f)
f.close()


'''import matplotlib.pyplot as plt 
plt.figure()
plt.scatter(np.arange(0,0.2,0.002),M_x_average,s=5)
plt.figure()
plt.scatter(np.arange(0,0.2,0.002),)'''


#---------------test code------------------------------
'''jj=[]
M_x_value_all=[]
for h in np.arange(0.1,1,0.05):
	a,b,phi_all,phi_base=main_operation(N,1-h,h,Gamma)
	energy_min=np.min(a)
	b_min=np.zeros(Gamma+1,dtype=complex)
	for i in range(len(a)):
		if a[i]==energy_min:
			b_min=b[:,i]
	phi_g_s=np.dot(phi_all,b_min)
	M_x_value=np.dot(np.conjugate(phi_g_s),M_x(phi_g_s,phi_base,N))
	M_x_value_all.append(M_x_value)
	jj.append(h)

print(M_x_value_all)

import matplotlib.pyplot as plt 
plt.figure()
plt.scatter(jj,M_x_value_all)
plt.show()'''

'''a,b,phi_all,phi_base=main_operation(N,0.1,0.9,Gamma)
energy_min=min(a)
for i in range(len(a)):
	if a[i]==energy_min:
		b_min=b[:,i]
phi_g_s=np.dot(phi_all,b_min)
print(np.dot(np.conjugate(phi_g_s),M_x(phi_g_s,phi_base,N)))'''

'''phi_base=[]
for i in range(2**N):
	if check_rep_dict[i][0]==1:
		phi_base.append(i)
nn=len(phi_base)
H=np.zeros((nn,nn))
for j in range(nn):
	a=np.zeros(nn)
	a[j]=1
	phi_next=M_x(a,phi_base,N)
	for k in range(nn):
		b=np.zeros(nn)
		b[k]=1
		H[k,j]=np.dot(np.conjugate(b),phi_next)
c,d=np.linalg.eig(H)
print(c)'''