
#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"3j6j9j.h"

const double stepTab[15]={1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800.0, 479001600.0, 6227020800.0,87178291200.0};
#define max(x,z) ((x)>(z) ? (x):(z))
#define min(x,z) ((x)<(z) ? (x):(z))

double theta(int x)
{
	if (x%2==0)return 1;
	else return -1;    
}

double factorial(int n)
{   
	if(n<0)return 0;
	if(n<=14)return stepTab[n];
	else return n*factorial(n-1);

}

double delta(int a,int b,int c)
{
	int s=a+b+c;
	double z;
	z=sqrt(factorial(s/2-a)*factorial(s/2-b)*factorial(s/2-c)/factorial(s/2+1));
	return z;  
}

double j3symbol(int j1,int j2,int j3,int m1,int m2,int m3)
{
	double y;int t,t1=0,t2;double z=0;
	if(j3<fabs(j1-j2)||j3>(j1+j2)||(m1+m2)!=-m3)return 0;
	if(m1<-j1||m1>j1||m2<-j2||m2>j2||m3<-j3||m3>j3)return 0;
	if ((j1-m1)%2 || (j2-m2)%2 || (j3-m3)%2) return 0;
	if(t1<(-j3+j2-m1)/2)t1=(-j3+j2-m1)/2;
	if(t1<(-j3+j1+m2)/2)t1=(-j3+j1+m2)/2;

	t2=(j1+j2-j3)/2;
	if(t2>(j1-m1)/2)t2=(j1-m1)/2;
	if(t2>(j2+m2)/2)t2=(j2+m2)/2;

	for(t=t1;t<=t2;t++)z+=theta((j1-j2-m3)/2+t)/(factorial(t)*factorial((j1+j2-j3)/2-t)*factorial((j3-j2+m1)/2+t)
			*factorial((j3-j1-m2)/2+t)*factorial((j1-m1)/2-t)*factorial((j2+m2)/2-t));
	t=0;
	y=pow((factorial((j1+j2-j3)/2)*factorial((j2+j3-j1)/2)*factorial((j3+j1-j2)/2)/factorial((j1+j2+j3)/2+1)),0.5)
		*pow((factorial((j1+m1)/2)*factorial((j1-m1)/2)*factorial((j2+m2)/2)*factorial((j2-m2)/2)*factorial((j3-m3)/2)*factorial((j3+m3)/2)),0.5)
		*z;
	return y;
}

/*
   double calCG()----------Clebsch-Gordan coefficients

Input: doubled angular j1, m1, j2, m2, j3, m3

Output: double float
 */
double calCG(int j1, int m1, int j2, int m2, int j3, int m3)//spins are all doubled
{
	double y;
	y=theta((j1-j2+m3)/2)*sqrt(j3+1)*j3symbol(j1,j2,j3,m1,m2,-m3);
	return y;
}

int table_num;
double ******j66;

bool flag_j6table=false;

double cal6j(int j1,int j2,int j3,int l1,int l2,int l3)
{
if(flag_j6table){
cout<<"Oh no ----------------------------------Oh no-------------------------------------------Oh no---------------------"<<endl;
cout<<j1<<"\t"<<j2<<"\t"<<j3<<"\t"<<l1<<"\t"<<l2<<"\t"<<l3<<endl;
}
	if(j3<abs(j1-j2)||j3>j1+j2)return 0;//triangle j1,j2,j3
	if(l2<abs(j3-l1)||l2>j3+l1)return 0;//triangle l1,l2,j3
	if(l3<abs(j2-l1)||l3>j2+l1)return 0;//triangle l1,j2,l3
	if((j3-j1-j2)%2||(l2-j3-l1)%2||(l3-j2-l1)%2)return 0;

	//int m1,m2,m11;

	double z=0,co;
	int y[7];
	y[0]=j1+j2+j3;
	y[1]=j1+l2+l3;
	y[2]=l1+j2+l3;
	y[3]=l1+l2+j3;
	y[4]=j1+j2+l1+l2;
	y[5]=j2+j3+l2+l3;
	y[6]=j3+j1+l3+l1;
	if (y[0]%2 || y[1]%2 || y[2]%2 || y[3]%2 || y[4]%2 || y[5]%2 || y[6]%2) return 0;
	if (j1+j2<j3 || abs(j1-j2)>j3 || l1+j3<l2 || abs(l1-j3)>l2 || j2+l1<l3 || abs(j2-l1)>l3 || j1+l3<l2 || abs(j1-l3)>l2) return 0;
	int tmin,tmax,t;
	tmin=max(max(y[0],y[1]),max(y[2],y[3]));
	tmax=min(y[4],min(y[5],y[6]));
	//if(j1==12&&8==j2&&4==j3&&4==l1&&0==l2&&4==l3)cout<<"tmin="<<tmin<<" tmax="<<tmax<<endl;
	for(t=tmin;t<=tmax;t+=2)
	{
		//if(j1==12&&8==j2&&4==j3&&4==l1&&0==l2&&4==l3) 
		//cout<<(factorial((t-y[0])/2)*factorial((t-y[1])/2)*factorial((t-y[2])/2)
		//*factorial((t-y[3])/2)*factorial((y[4]-t)/2)*factorial((y[5]-t)/2)*factorial((y[6]-t)/2))<<endl;
		z+=theta(t/2)*factorial(t/2+1)/(factorial((t-y[0])/2)*factorial((t-y[1])/2)*factorial((t-y[2])/2)
				*factorial((t-y[3])/2)*factorial((y[4]-t)/2)*factorial((y[5]-t)/2)*factorial((y[6]-t)/2));
	}
	//cout<<setprecision(20)<<z<<endl;
	co=delta(j1,j2,j3)*delta(j1,l2,l3)*delta(l1,j2,l3)*delta(l1,l2,j3);
	//if(j1==12&&8==j2&&4==j3&&4==l1&&0==l2&&4==l3)cout<<delta(j1,j2,j3)<<" "<<delta(j1,l2,l3)<<" "<<delta(l1,j2,l3)<<" "<<delta(l1,l2,j3)<<endl;
	z*=co;
	//cout<<""
	/*for(m1=-j1;m1<=j1;m1+=2)
	  for(m2=-j2;m2<=j2;m2+=2)
	  for(m11=-l1;m11<=l1;m11+=2)
	  y+=theta((j1+j2+j3+l1+l2+l3+3*m11+m1+2*m2)/2)*j3symbol(j1,j2,j3,m1,m2,-m1-m2)
	 *j3symbol(j1,l2,l3,-m1,m11+m1+m2,-m11-m2)*j3symbol(l1,j2,l3,-m11,-m2,m11+m2)
	 *j3symbol(l1,l2,j3,m11,-m11-m1-m2,m1+m2);*/

	return z;

}

/*
   j6table()-------make a table of the 6j coefficients in the memory

Input: void
Output: void
Contents: open memory for j66[][][][][][], and store 6j coefficients
 */
/*
J	K	Jp
jd	je	jf
*/
void j6table(int K){

	j66=new double *****[2*table_num];
	for(int j1=0;j1<2*table_num;j1++){
		j66[j1]=new double ****[1];
		for(int j2=0;j2<1;j2++){
			j66[j1][j2]=new double ***[2*table_num];
			for(int j3=0;j3<2*table_num;j3++){
				j66[j1][j2][j3]=new double **[table_num];
				for(int l1=0;l1<table_num;l1++){
					j66[j1][j2][j3][l1]=new double *[table_num];
					for(int l2=0;l2<table_num;l2++){
						j66[j1][j2][j3][l1][l2]=new double [table_num];
					}
				}
			}
		}
	}

	for(int j1=0;j1<2*table_num;j1++){
		for(int j2=0;j2<1;j2++){
			for(int j3=0;j3<2*table_num;j3++){
				for(int l1=0;l1<table_num;l1++){
					for(int l2=0;l2<table_num;l2++){                  
						for(int l3=0;l3<table_num;l3++){       
							j66[j1][j2][j3][l1][l2][l3]=cal6j(2*j1,K,2*j3,1+2*l1,1+2*l2,1+2*l3);
						}
					}
				}
			}
		}
	}
}
