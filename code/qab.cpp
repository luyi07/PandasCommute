
#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"3j6j9j.h"

double r2(int n,int l1,int l2)
{
	n=2*n+l1;
	double y;
	if(l1==l2)
		y=n+1.5;
	else if(l1==l2+2)
		y=-sqrt(n+l2+3)*sqrt(n-l2);  //warning
	else if(l1==l2-2)
		y=-sqrt(n+l2+1)*sqrt(n-l2+2); //warning
	else 
		y=0;
	return y;
}

double gamma(int x)//x is already doubled
{
	if(x<0||x==0){
		cout<<"error: invalid x for gamma(x)"<<endl;
		return 1E6;
	}
	else if(x%2==0)return factorial(x/2-1);
	else if(x==1) return sqrt(M_PI);
	else return (double)(x-2)/2*gamma(x-2);
}

double rr(int n1,int l1, int n2, int l2, int lambda) //   /pow(alpha,lambda)
{
	double y;
	double temp;
	temp=factorial(n1)*factorial(n2)/gamma(2*(l1+n1)+3)/gamma(2*(l2+n2)+3);      
	temp=theta(n1+n2)*sqrt(temp);
	y=temp*gamma(l1-l2+lambda+2)*gamma(l2-l1+lambda+2);
	int uplimit=n1;
	if(n2<n1)uplimit=n2;
	int lowlimit=2*n2-(l1-l2+lambda);
	if(lowlimit<2*n1-(l2-l1+lambda))
		lowlimit=2*n1-(l2-l1+lambda);
	if(lowlimit%2==0)lowlimit/=2;					  
	else lowlimit=(lowlimit+1)/2;
	if(lowlimit<0)lowlimit=0;

	temp=0;								   
	for(int k=lowlimit;k<=uplimit;k++)                                
	{
		temp+=gamma(l1+l2+lambda+3+2*k)/factorial(k)/factorial(n1-k)
			/factorial(n2-k)/gamma(l1-l2+lambda+2*(k+1-n2))
			/gamma(l2-l1+lambda+2*(k+1-n1));
	}                                                            
	y=y*temp;
	return y;
}      

/*
   rr() is symmetric in respect to the exchange of (n1l1j1, n2l2j2) orbits.
   in qab() only the phase factor theta(j1-1/2) is not.
   So qab(ab)=(-1)^{j1-j2}q(ba)=(-1)^{j1+j2+1}q(ba) (Conden-Shortly)
 */
double qab(int n1,int l1,int n2,int l2,int j1,int j2,int lambda)//spins j1 and j2 are all doubled
{
	double y;
	if(j1+j2>=2*lambda
	&&abs(j1-j2)<=2*lambda){
		y=theta(l1+l2+1)*(1+theta(l1+l2+lambda))/2*theta((j1-1)/2)
			*sqrt(j1+1)*sqrt(j2+1)/sqrt(4*M_PI)/sqrt(2*lambda+1);
		y=y*calCG(j1,1,j2,-1,lambda*2,0);
		y*=rr(n1,l1,n2,l2,lambda);
		//	if(lambda==1)y=n1*n2;
		return y;
	}
	else{
		return 0;
	}
}

/*
   M1 operator: \vec{L}--structure coefficients
   j1, j2, lambda are doubled
   n1, l1, n2, l2 are not
 */
double qM_L(int n1, int l1, int j1, int n2, int l2, int j2, int lambda){

	if(l1!=l2)return 0;

	double y1=2/(double)(lambda+1)*theta((j1+j2)/2+l1+l2+lambda+1);
	y1*=(1+theta(l1+l2+lambda-1))/2;
	y1*=sqrt(j2*(j2+2)/(16*M_PI));
	y1*=sqrt(j1+1)*sqrt(j2+1)*sqrt(2*(lambda-1)+1);
	y1*=calCG(j1,1,2*(lambda-1),0,j2,1);
	y1*=cal6j(j1,2*lambda,j2,2,j2,2*(lambda-1));

	double y2=-2/(double)(lambda+1)*sqrt(1/(16*M_PI*lambda));
	y2*=theta(lambda)*(1+theta(l1+l2+lambda-1))/2;
	y2*=sqrt(j1+1)/sqrt(2*lambda+1);
	y2*=lambda-(theta((j1+1)/2+l1)*(j1+1)+theta((j2+1)/2+l2)*(j2+1))/2;
	y2*=calCG(j1,1,2*lambda,0,j2,1);

	double y=y1+y2;
	y*=sqrt(lambda*(2*lambda+1))*rr(n1,l1,n2,l2,lambda-1);

	return y;
}

/*
   M1 operator: \vec{S}--structure coefficients
   j1, j2, are doubled
   n1, l1, n2, l2, lambda are not
 */
double qM_S(int n1, int l1, int j1, int n2, int l2, int j2, int lambda){

	if(l1!=l2)return 0;

	double y=sqrt(1/(16*M_PI*lambda))*theta(lambda)*(1+theta(l1+l2+lambda-1))/2;
	y*=sqrt(j1+1)/sqrt(2*lambda+1);
	y*=lambda-(theta((j1+1)/2+l1)*(j1+1)+theta((j2+1)/2+l2)*(j2+1))/2;
	y*=calCG(j1,1,2*lambda,0,j2,1);
	y*=sqrt(lambda*(2*lambda+1))*rr(n1,l1,n2,l2,lambda-1);
	return y;
}
