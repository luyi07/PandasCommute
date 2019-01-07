/*
In usdb.int V(1234;JT) has constraints like j1>=j2, j3>=j4 (j1, j2 can be just indices), because of symmetry. 
nn, pp parts is the same as T=1 part, except the indices
pn should be collected from both T=0 and T=1 interactions
*/
#include<iostream>
using namespace std;
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<cmath>

#include"input.cpp"

int theta(int a){
	if(a%2==0)return 1;
	else return -1;
}

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)

int main (){

	char file_sp[80]="../input/pn.sp";
//      read in single particle orbits
        int nj_p, nj;
        read_nj(file_sp, nj_p, nj);

        int *nn=new int[nj];
        int *ll=new int[nj];
        int *sh=new int[nj];

        read_sp(file_sp, nj, nn, ll, sh);

//      label sequences of orbits, according their j value, from large to small
        int *seq=new int[nj];
        for(int i=0;i<nj;i++){
                seq[i]=i;//0,1,2,...
        }
        cout<<"seq:";
        for(int i=0;i<nj;i++)cout<<seq[i]<<",";
        cout<<endl;
        for(int i=0;i<nj_p;i++){
                for(int j=i;j<nj_p-1;j++){
                        if(sh[seq[j]]<sh[seq[j+1]]){
                                int temp=seq[j];
                                seq[j]=seq[j+1];
                                seq[j+1]=temp;
                        }
                }
        }
        for(int i=nj_p;i<nj;i++){
                for(int j=i;j<nj-1;j++){
                        if(sh[seq[j]]<sh[seq[j+1]]){
                                int temp=seq[j];
                                seq[j]=seq[j+1];
                                seq[j+1]=temp;
                        }
                }
        }
        cout<<"seq:";
        for(int i=0;i<nj;i++)cout<<seq[i]<<",";
        cout<<endl;

//------------read general interactions------------
        char file_GME[80]="../input/GME.int";
        int num_GME;  //number of matrix elements
        read_GME_num(file_GME, num_GME);

	double *spe=new double [nj_p];
        int *GME_a=new int [num_GME];                   //sequence of ja, starts from 1
        int *GME_b=new int [num_GME];                   //sequence of jb
        int *GME_c=new int [num_GME];                   //sequence of jc
        int *GME_d=new int [num_GME];                   //sequence of jd
        int *GME_I=new int [num_GME];                   //I
        int *GME_T=new int [num_GME];                   //T
        double *GME_V=new double [num_GME];             //V(abcd;JT)
        read_GME(file_GME, nj_p, spe, num_GME, GME_a, GME_b, GME_c, GME_d, GME_I, GME_T, GME_V);

//      determine J_min, J_max for proton-neutron interactions
        int J_min=1E6, J_max=0;
        for(int i=0;i<nj_p;i++)
        for(int j=nj_p;j<nj;j++){
                if(J_min>abs(sh[i]-sh[j]))J_min=abs(sh[i]-sh[j]);
                if(J_max<sh[i]+sh[j])J_max=sh[i]+sh[j];
        }
        cout<<"J_min= "<<J_min<<" J_max="<<J_max<<endl;
	
//      count number of proton-neutron interactions
	int num_pn=0;
        for(int J=J_min;J<=J_max;J+=2){
                for(int a=0;a<nj_p;a++)
                for(int b=nj_p;b<nj;b++)
                for(int c=0;c<nj_p;c++)
                for(int d=nj_p;d<nj;d++)
                if(	min(a,b-nj_p) < min(c, d-nj_p)
			||
			( min(a,b-nj_p) == min(c,d-nj_p)
			&& max(a,b-nj_p) < max(c,d-nj_p))
			||
			( min(a,b-nj_p) == min(c,d-nj_p)
			&& max(a,b-nj_p) == max(c,d-nj_p)
			&& (	(a==c&&b==d)
				||
			    	(a<=c)
			   )
			)
		)
                if(J>=abs(sh[seq[a]]-sh[seq[b]])
                &&J<=sh[seq[a]]+sh[seq[b]]//angular triangles
                &&J>=abs(sh[seq[c]]-sh[seq[d]])
                &&J<=sh[seq[c]]+sh[seq[d]])
                if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2){//parity conservation
		num_pn++;
		}
	}

	cout<<"num_pn="<<num_pn<<endl;

//------get Vpns with only symmetry V(abcd)=V(cdab)--------------
        int *GME_a_pn=new int [num_pn];//sequence of ja
        int *GME_b_pn=new int [num_pn];//sequence of jb
        int *GME_c_pn=new int [num_pn];//sequence of jc
        int *GME_d_pn=new int [num_pn];//sequence of jd
        int *GME_I_pn=new int [num_pn];//I
        double *GME_V_pn=new double [num_pn];//V(abcd;JT)
	for(int i=0;i<num_pn;i++){
		GME_V_pn[i]=0;
	}
	
	int temp=0;
        for(int J=J_min;J<=J_max;J+=2){
                for(int a=0;a<nj_p;a++)
                for(int b=nj_p;b<nj;b++)
                for(int c=0;c<nj_p;c++)
                for(int d=nj_p;d<nj;d++)
                if(	min(a,b-nj_p) < min(c, d-nj_p)
			||
			( min(a,b-nj_p) == min(c,d-nj_p)
			&& max(a,b-nj_p) < max(c,d-nj_p))
			||
			( min(a,b-nj_p) == min(c,d-nj_p)
			&& max(a,b-nj_p) == max(c,d-nj_p)
			&& (	(a==c&&b==d)
				||
			    	(a<=c)
			   )
			)
		)
                if(J>=abs(sh[seq[a]]-sh[seq[b]])
                &&J<=sh[seq[a]]+sh[seq[b]]//angular triangles
                &&J>=abs(sh[seq[c]]-sh[seq[d]])
                &&J<=sh[seq[c]]+sh[seq[d]])
                if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2){//parity conservation
		GME_a_pn[temp]=seq[a]+1;
		GME_b_pn[temp]=seq[b]+1;
		GME_c_pn[temp]=seq[c]+1;
		GME_d_pn[temp]=seq[d]+1;
		GME_I_pn[temp]=J/2;
                temp++;
                }
        }

	for(int i=0;i<num_GME;i++){
                int a=GME_a[i];
                int b=GME_b[i];
                int c=GME_c[i];
                int d=GME_d[i];
                int I=2*GME_I[i];
                int T=2*GME_T[i];
                double V=GME_V[i];

		for(int k=0;k<num_pn;k++){
			if(GME_a_pn[k]==a	//V(abcd;I)
			&&GME_b_pn[k]==b+nj_p
			&&GME_c_pn[k]==c
			&&GME_d_pn[k]==d+nj_p
			&&GME_I_pn[k]==I/2){
				GME_V_pn[k]+=V;
			}
			
			if(a!=b){
				if(GME_a_pn[k]==b	//V(bacd;I) omitted in usdb.int is collected in GME_V_pn
				&&GME_b_pn[k]==a+nj_p
				&&GME_c_pn[k]==c
				&&GME_d_pn[k]==d+nj_p
				&&GME_I_pn[k]==I/2){
					GME_V_pn[k]+=V*theta((sh[a-1]+sh[b-1]+I+T)/2);
				}

				if(c!=d){
					if(GME_a_pn[k]==b	//V(badc;I) omitted in usdb.int is collected in GME_V_pn
					&&GME_b_pn[k]==a+nj_p
					&&GME_c_pn[k]==d
					&&GME_d_pn[k]==c+nj_p
					&&GME_I_pn[k]==I/2){
						GME_V_pn[k]+=V*theta((sh[a-1]+sh[b-1]+sh[c-1]+sh[d-1])/2);
					}
					if(GME_a_pn[k]==a       //V(abdc;I) omitted in usdb.int is collected in GME_V_pn
					&&GME_b_pn[k]==b+nj_p
					&&GME_c_pn[k]==d
					&&GME_d_pn[k]==c+nj_p
					&&GME_I_pn[k]==I/2){
						GME_V_pn[k]+=V*theta((sh[c-1]+sh[d-1]+I+T)/2);
					}
				}
			}
			else if(c!=d){
				if(GME_a_pn[k]==a	//V(aadc;I) omitted in usdb.int is collected in GME_V_pn
				&&GME_b_pn[k]==b+nj_p
				&&GME_c_pn[k]==d
				&&GME_d_pn[k]==c+nj_p
				&&GME_I_pn[k]==I/2){
					GME_V_pn[k]+=V*theta((sh[c-1]+sh[d-1]+I+T)/2);
				}
			}
		}
	}

	int num_GMEpn=0;
	for(int i=0;i<num_GME;i++){
		if(GME_T[i]==1){
			num_GMEpn+=2;
		}
	}
	num_GMEpn+=num_pn;
/*
	FILE *fp1;
	if((fp1=fopen("unnormalized.int","w"))==NULL)
	cout<<"error: failed to open unnormalized.int"<<endl;
	for(int i=0;i<num_GME;i++){
                int a=GME_a[i];
                int b=GME_b[i];
                int c=GME_c[i];
                int d=GME_d[i];
                int I=2*GME_I[i];
                int T=2*GME_T[i];
                double V=GME_V[i];
		if(GME_T[i]==1)fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
				a, b, c, d, I/2, T/2, V);
	}
	for(int i=0;i<num_GME;i++){
                int a=GME_a[i];
                int b=GME_b[i];
                int c=GME_c[i];
                int d=GME_d[i];
                int I=2*GME_I[i];
                int T=2*GME_T[i];
                double V=GME_V[i];
		if(GME_T[i]==1)fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
				a+nj, b+nj, c+nj, d+nj, I/2, T/2, V);
	}
	for(int i=0;i<num_GME;i++){
                int a=GME_a[i];
                int b=GME_b[i];
                int c=GME_c[i];
                int d=GME_d[i];
                int I=2*GME_I[i];
                int T=2*GME_T[i];
                double V=GME_V[i];
		fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
			a, b+nj, c, d+nj, I/2, T/2, V);
		if(a!=b)fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
			b, a+nj, c, d+nj, I/2, T/2, theta((sh[a-1]+sh[b-1]+I+T)/2)*V);
		if(c!=d)fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
			a, b+nj, d, c+nj, I/2, T/2, theta((sh[c-1]+sh[d-1]+I+T)/2)*V);
		if(a!=b&&c!=d)fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
			b, a+nj, d, c+nj, I/2, T/2, theta((sh[a-1]+sh[b-1]+sh[c-1]+sh[d-1])/2)*V);
	}
	fclose(fp1);
*/
	FILE *fp;
	char file_GMEpn[80]="../input/GMEpn.int";
	if((fp=fopen(file_GMEpn,"w"))==NULL)
	cout<<"error: failed to open "<<file_GMEpn<<endl;
	fprintf(fp,"! orbits 1-%d protons, %d-%d neutrons\n",nj_p,nj_p+1,nj);
	fprintf(fp,"!SP ORB  N  L  J    Tz\n");
	for(int i=0;i<nj_p;i++){
		fprintf(fp,"!SP   %d  %d  %d  %d/2    %d\n",i+1,nn[i],ll[i],sh[i],1);
	}
	for(int i=nj_p;i<nj;i++){
		fprintf(fp,"!SP   %d  %d  %d  %d/2    %d\n",i+1,nn[i],ll[i],sh[i],-1);
	}
	fprintf(fp,"%d   ",num_GMEpn);
	for(int i=0;i<nj_p;i++){
		fprintf(fp,"%lf   ",spe[i]);
	}
	for(int i=0;i<nj_p;i++){
		fprintf(fp,"%lf   ",spe[i]);
	}
	fprintf(fp,"\n");
	for(int i=0;i<num_GME;i++){
		if(GME_T[i]==1){
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t1\t%lf\n",GME_a[i],GME_b[i],GME_c[i],GME_d[i],GME_I[i],GME_V[i]);//Vpp=V(T=1)
		}
	}
	for(int i=0;i<num_GME;i++){
		if(GME_T[i]==1){//Vnn=V(T=1)
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t1\t%lf\n",GME_a[i]+nj_p,GME_b[i]+nj_p,GME_c[i]+nj_p,GME_d[i]+nj_p,GME_I[i],GME_V[i]);
		}
	}
	for(int i=0;i<num_pn;i++){
		if(GME_a_pn[i]==GME_b_pn[i]-nj_p)GME_V_pn[i]*=sqrt(2);//normalize proton-neutron interactions
		if(GME_c_pn[i]==GME_d_pn[i]-nj_p)GME_V_pn[i]*=sqrt(2);
		GME_V_pn[i]/=2;
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\t1\t%lf\n",GME_a_pn[i],GME_b_pn[i],GME_c_pn[i],GME_d_pn[i],GME_I_pn[i],GME_V_pn[i]);
	}
	fclose(fp);
	return 0;
}
