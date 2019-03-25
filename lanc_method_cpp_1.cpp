#include <iostream>
#include <assert.h>
#include <math.h>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <Eigen/Dense>
#include <vector>
#define NUMBER 12
#define GAMMA 20
using namespace std;
using namespace Eigen;
//This program is to accomplish Lanczos method algorithm.
//All global parameters are defined firstly.

//The flip function is to flip the ith digit of a binary number,N is the digits we need.
int flip(int num,int i,int N){
	assert(i<=N&&i>=0);
	if(num>>(i)&1){
		return num-pow(2,i);
	}
	else{
		return num+pow(2,i);
	}
}

Matrix<int,1,2> find_min(vector<int> l){
	int min_item=l[0];
	int j=0;
	for(int i=0;i<l.size();++i){
		if(l[i]<min_item){
			j=i;
			min_item=l[i];
		}
	}
	Matrix<int,1,2> result;
	result<<j,min_item;
	return result; 
}


int rotate(int num){
	int final_num=num&1;
	int next_num=num>>1;
	int begin_num=num>>NUMBER&1;
	if(begin_num!=final_num){
		next_num=flip(next_num,NUMBER-1,NUMBER);
	}
	return next_num;
}

//The flip_all function is to flip all of the digits of a binary number.
int flip_all(int num,int N){
	int not_num=~num;
	int max_num=pow(2,N-1);
	return (not_num&max_num);
}

Matrix<int,1,5> check_rep(int num){
	int count=0;
	int num_begin=num;
	vector<int> num_group;
	num_group.push_back(num);
	while(1){
		count++;
		int num_next=rotate(num);
		if(num_next==num_begin)
			break;
		else
		{
			num_group.push_back(num_next);
			num=num_next;
		}
	}
	vector<int> num_group_p;
	for(int i=0;i<num_group.size();++i){
		int item=flip_all(num_group[i],NUMBER);
		num_group_p.push_back(item);
	}
	Matrix<int,1,2> find_min_1=find_min(num_group);
	Matrix<int,1,2> find_min_2=find_min(num_group_p);
	int min_all=min(find_min_1(0,1),find_min_2(0,1));
	int r=0;
	int p=0;
	int is_rep=0;
	if(min_all==num_begin)
		is_rep=1;
	else{
		int num_flip_all=flip_all(num_begin,NUMBER);
		if(min_all==num_flip_all)
			p=1;
		else{
			if(find_min_1(0,1)>find_min_2(0,1)){
				p=1;
				r=count-find_min_2(0,1);
			}
			else
				r=count-find_min_1(0,1);
		}
	}
	Matrix<int,1,5> result;
	result<<is_rep,count,r,p,min_all;
	return result;
}
//If we offer a  state which is expressed by a one-column matrix,the Normalize function 
//can help us normalize the state.
complex<double> Normalize(Matrix<complex<double>,1,Dynamic> phi_){
	complex<double> norm_para=phi_.dot(phi_.transpose().conjugate());
	norm_para=sqrt(norm_para);
	return norm_para;
}

Matrix<complex<double>,1,Dynamic> H(Matrix<complex<double>,1,Dynamic> phi,int N,double J,double h){
	int max_i=pow(2,N);
	Matrix<complex<double>,1,Dynamic> phi_next=phi;
	for(int i=0;i<max_i;++i){
		complex<double> H_z_value_item=0;
		for(int k=0;k<N;++k){
			int j=(k+1)%N;
			if((i>>k&1)==(i>>j&1)){
				H_z_value_item+=-J;
			}
			else{
				H_z_value_item+=J;
			}
		}
		H_z_value_item=H_z_value_item*phi(0,i);
		complex<double> H_x_value_item=0;
		for(int j=0;j<N;++j){
			int phi_after_flip=flip(i,j,N);
			H_x_value_item+=-h*phi(0,phi_after_flip);
		}
		phi_next(0,i)=H_x_value_item+H_z_value_item;
	}
	return phi_next;
}


Matrix<complex<double>,Dynamic,1> lanczos_method(Matrix<complex<double>,1,Dynamic> cin_phi_0)
{
	const int NUM=pow(2,NUMBER);
	Matrix<complex<double>,1,Dynamic> phi_0=cin_phi_0;
	Matrix<complex<double>,2,GAMMA+1> parameter;
	complex<double> norm=Normalize(phi_0);
	parameter(0,0)=norm;
	phi_0=phi_0/norm;
	Matrix<complex<double>,1,Dynamic> phi_1=H(phi_0,NUMBER,0.2,0.8);
	complex<double> d=phi_0.dot(phi_1.conjugate().transpose());
	parameter(1,0)=d;
	phi_1=phi_1-d*phi_0;
	norm=Normalize(phi_1);
	parameter(0,1)=norm;
	phi_1=phi_1/norm;
	for(int i=0;i<GAMMA-1;++i){
		cout<<"once"<<endl;
		Matrix<complex<double>,1,Dynamic> phi_2=H(phi_1,NUMBER,0.2,0.8);
		parameter(1,i+1)=phi_1.dot(phi_2.conjugate().transpose());
		phi_2=phi_2-parameter(1,i+1)*phi_1-parameter(0,i+1)*phi_0;
		parameter(0,i+2)=Normalize(phi_2);
		phi_2=phi_2/parameter(0,i+2);
		phi_0=phi_1;
		phi_1=phi_2;
	}
	Matrix<complex<double>,GAMMA,GAMMA> H;
	for(int i=0;i<GAMMA;++i){
		for(int j=0;j<GAMMA;++j){
			if(i==j)
				H(i,j)=parameter(1,i);
			if(j==i+1)
				H(i,j)=parameter(0,j);
			if(j==i-1)
				H(i,j)=parameter(0,i);
		}
	}
	ComplexEigenSolver<MatrixXcd>es(H);
	return es.eigenvalues();
}

Matrix<complex<double>,1,Dynamic> Judge_if_phi_0(Matrix<complex<double>,1,Dynamic> cin_phi_0)
{
	const int NUM=pow(2,NUMBER);
	Matrix<complex<double>,1,NUM> phi_0;
	if(cin_phi_0.size()==0){
		srand(unsigned(time(NULL)));
		complex<double> * p=new complex<double> [NUM];
		for(int i=0;i<NUM;++i){
			p[i]=double(rand()%1000)/1000;
			//p[i]=1;
		}
		Map<Matrix<complex<double>,1,NUM>> phi_0(p);
		return lanczos_method(phi_0);
	}
	else
		return lanczos_method(cin_phi_0);
}


int main(){
	/*Matrix<complex<double>,1,0> pi;
	Matrix<complex<double>,Dynamic,1> eigenvalue=Judge_if_phi_0(pi);
	cout<<eigenvalue<<endl;*/
	return 0;
}