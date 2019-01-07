#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<fstream>
using namespace std;

#include"3j6j9j.h"
#include"input.h"
#include"qab.h"

double beta_coef(int a, int b, int sign, int L){

	char file_sp[80]="input/pn.sp";
//--------------read shell info--------------
        int nj_p, nj;
        read_nj(file_sp, nj_p, nj);
        int nj_n=nj-nj_p;

        int *nn=new int [nj];
        int *ll=new int [nj];
        int *sh=new int [nj];
        read_shell(file_sp, nj, nn, ll, sh);

        int *parity=new int [nj];
        for(int i=0;i<nj;i++){
                if(ll[i]%2==0)parity[i]=1;
                else parity[i]=-1;
        }
//-----------------------------------------------
	double y=0;
	if(sign==1){
		if(L==0){
			y=-sqrt(sh[a]+1);
			if(sh[a]!=sh[b]) y=0;
		}
		else if(L==2){
			if(nn[a]!=nn[b])return 0;
			if(ll[a]!=ll[b])return 0;
			y=-sqrt(2)*sqrt(sh[a]+1)*sqrt(sh[b]+1)*theta(ll[a]+(sh[a]+3)/2);
			y*=cal6j(1,1,2,sh[a],sh[b],2*ll[a]);
		}
	}
	else if(sign==-1){
		if(L==0){
			y=-sqrt(sh[a]+1);
			if(sh[a]!=sh[b]) y=0;
		}
		else if(L==2){
			if(nn[a]!=nn[b])return 0;
			if(ll[a]!=ll[b])return 0;
			y=-sqrt(2)*sqrt(sh[a]+1)*sqrt(sh[b]+1)*theta(ll[a]+(sh[a]+3)/2);
			y*=cal6j(1,1,2,sh[a],sh[b],2*ll[a]);
		}
	}
   
	delete[] nn;
	delete[] ll;
	delete[] sh;
	delete[] parity;

	return y;
}

void cal_Fcoef(char * file_sp, char * file_F){

//--------------read shell info--------------
        int nj_p, nj;
        read_nj(file_sp, nj_p, nj);
        int nj_n=nj-nj_p;

        int *nn=new int [nj];
        int *ll=new int [nj];
        int *sh=new int [nj];
        read_shell(file_sp, nj, nn, ll, sh);

        int *parity=new int [nj];
        for(int i=0;i<nj;i++){
                if(ll[i]%2==0)parity[i]=1;
                else parity[i]=-1;
        }
//	for test only
	/*
	for(int i=0;i<nj;i++){
		cout<<"n="<<nn[i]<<" l="<<ll[i]<<" j="<<sh[i]<<endl;
	}
	*/
//-------------produce Fcoef-------------------

	int Ftype=0;
	while(Ftype<1||Ftype>3){
		cout<<"\tWhat kind of one-body transition would you like to do?\n";
		cout<<"\t\t\t'1' for Electronic Multipole transitions."<<endl;
		cout<<"\t\t\t'2' for Magnetic Multipole transitions."<<endl;
		cout<<"\t\t\t'3' for Charge Exchange transitions."<<endl;
		cout<<"\t";
		cin>>Ftype;
		if(Ftype<1||Ftype>3)cout<<"\tI don't understand. Allow me to ask you again...\n";
	}

	int K=-1;
	while(K<0){
		cout<<"\tThe angular momentum of the one-body operator K=?\n\t";
		cin>>K;
		if(Ftype==3&&K>1)K=-1;
		if(K<0)cout<<"\tSorry, I can not do forbidden decay yet.\n";
	}

	int pF;
	if(Ftype==1){
		if(K%2==0)pF=1;
		else pF=-1;
	}
	else if(Ftype==2){
		if(K%2==0)pF=-1;
		else pF=1;
	}
	else if(Ftype==3){
		pF=1;
	}

	if(Ftype==1||Ftype==2){
		int switch_formalism=2;
		while(switch_formalism!=0
			&&switch_formalism!=1){

			cout<<"\tThe one-body operator has good isospin/or not?(1/0)\n\t";
			cin>>switch_formalism;

			if(switch_formalism!=0
			&&switch_formalism!=1){
				cout<<"\tinvalid input, sir/madam."<<endl;
			}
		}

		if(switch_formalism==1){
			int T_F=2;
			int Tz_F=2;
			while(T_F<0||T_F>1||Tz_F<-T_F||Tz_F>T_F){
				cout<<"\tThe isospin of the operator dT= ? (0/1)\n\t";
				cin>>T_F;
				if(T_F<0||T_F>1){
					cout<<"\tinvalid input..."<<endl;
				}
				if(T_F>0){
					cout<<"\tdTz= ? (-dT <= dTz <= dT)\n\t";
					cin>>Tz_F;
					if(Tz_F<-T_F||Tz_F>T_F){
						cout<<"\tinvalid input..."<<endl;
					}
				}
				else{
					Tz_F=0;
				}
			}

			double ee;
			if(Ftype==1){
				cout<<"\teffective charge = ";
				cin>>ee;
			}
			double gl,gs;
			if(Ftype==2){
				cout<<"\tgl = ";
				cin>>gl;
				cout<<"\tgs = ";
				cin>>gs;
			}

			FILE *fp;
			if((fp=fopen(file_F,"w"))==NULL)
			cerr<<"error: failed to open "<<file_F<<endl;

			fprintf(fp, "! orbits 1-%d protons, %d-%d neutrons\n", nj_p, nj_p+1, nj);
			fprintf(fp, "!SP ORB  N  L  J     Tz\n");
			for(int i=0;i<nj;i++){
				fprintf(fp, "!SP   %d  %d  %d  %d/2\n", i+1, nn[i], ll[i], sh[i]);
			}
			fprintf(fp,"K= %d\tpF= %d\n",K,pF);

			for(int u=0;u<4;u++){//pp pn np nn
				for(int i=0;i<nj_p;i++){
					for(int j=i;j<nj_p;j++){
						double y;
						if(Ftype==0){
							y=0;
						}
						else if(Ftype==1){
							y=ee*qab(nn[i],ll[i],nn[j],ll[j],sh[i],sh[j],K)*sqrt(2*K+1);
						}
						else if(Ftype==2){
							y=(gl*qM_L(nn[i],ll[i],sh[i],nn[j],ll[j],sh[j],K)+gs*qM_S(nn[i],ll[i],sh[i],nn[j],ll[j],sh[j],K))*sqrt(2*K+1);
						}

						if(T_F==0){
							if(u>0&&u<3)y=0;
						}
						else{
							if(Tz_F==0){
								if(u>0&&u<3)y=0;
								else if(u==0)y*=-1;
							}
							else if(Tz_F==-1){
								if(u!=1)y=0;
							}
							else if(Tz_F==1){
								if(u!=2)y=0;
								else y*=-1;
							}
						}

						fprintf(fp, "%d\t%d\t%lf\n", i+1 + nj_p*(u/2), j+1 + nj_p*(u%2), y);
					}
				}
			}
			fclose(fp);
		}//if(switch_formalism)
		else{	
			FILE *fp;
			if((fp=fopen(file_F,"w"))==NULL)
			cout<<"error: failed to open "<<file_F<<endl;

			fprintf(fp, "! orbits 1-%d protons, %d-%d neutrons\n", nj_p, nj_p+1, nj);
			fprintf(fp, "!SP ORB  N  L  J     Tz\n");
			for(int i=0;i<nj;i++){
				fprintf(fp, "!SP   %d  %d  %d  %d/2\n", i+1, nn[i], ll[i], sh[i]);
			}
			fprintf(fp,"K= %d\tpF= %d\n",K,pF);
		//-------------produce Fpp, Fpn, Fnp, Fnn---------------------
			double ee;
			double gl,gs;
			double y;
			if(Ftype==1){
				cout<<"\t proton effective charge = ";
				cin>>ee;
			}
			else{
				cout<<"\t proton orbital g factor gl_p = ";
				cin>>gl;
				cout<<"\t proton spin g factor gs_p = ";
				cin>>gs;
			}
			for(int i=0;i<nj_p;i++){
				for(int j=0;j<nj_p;j++){
					if(Ftype==1)y=ee*qab(nn[i],ll[i],nn[j],ll[j],sh[i],sh[j],K)*sqrt(2*K+1);
					else y=( gl*qM_L(nn[i],ll[i],sh[i],nn[j],ll[j],sh[j],K)
							+
						 gs*qM_S(nn[i],ll[i],sh[i],nn[j],ll[j],sh[j],K)
						)*sqrt(2*K+1);
					if(fabs(y)>1E-6)fprintf(fp, "%d\t%d\t%lf\n", i+1, j+1, y);
				}
			}
//----------------------output Fnn
			if(Ftype==1){
				cout<<"\t neutron effective charge = ";
				cin>>ee;
			}
			else{
				cout<<"\t neutron orbital g factor gl_n = ";
				cin>>gl;
				cout<<"\t neutron spin g factor gs_n = ";
				cin>>gs;
			}

			for(int i=0;i<nj_n;i++){
				for(int j=0;j<nj_n;j++){
					if(Ftype==1)y=ee*qab(nn[nj_p+i], ll[nj_p+i], nn[nj_p+j],
							ll[nj_p+j], sh[nj_p+i], sh[nj_p+j], K)*sqrt(2*K+1);
					else y=( gl*qM_L(nn[nj_p+i],ll[nj_p+i],sh[nj_p+i],nn[nj_p+j],ll[nj_p+j],sh[nj_p+j],K)
							+
						 gs*qM_S(nn[nj_p+i],ll[nj_p+i],sh[nj_p+i],nn[nj_p+j],ll[nj_p+j],sh[nj_p+j],K)
						)*sqrt(2*K+1);
					if(fabs(y)>1E-6)fprintf(fp, "%d\t%d\t%lf\n", i+1+nj_p, j+1+nj_p, y);
				}
			}

			fclose(fp);
		}//else(switch_formalism)
	}
	else{
		FILE *fp;
		if((fp=fopen(file_F,"w"))==NULL)
		cout<<"error: failed to open "<<file_F<<endl;

		fprintf(fp, "! orbits 1-%d protons, %d-%d neutrons\n", nj_p, nj_p+1, nj);
		fprintf(fp, "!SP ORB  N  L  J     Tz\n");
		for(int i=0;i<nj;i++){
			fprintf(fp, "!SP   %d  %d  %d  %d/2\n", i+1, nn[i], ll[i], sh[i]);
		}
		fprintf(fp,"K= %d\tpF= %d\n",K,pF);

		for(int i=0;i<nj_p;i++){
			for(int j=0;j<nj_n;j++){
				double y=beta_coef(i,j,1,2*K)*sqrt(2*K+1);
				if(fabs(y)>1E-6)fprintf(fp, "%d\t%d\t%lf\n", i+1, j+1 + nj_p, y);
			}
		}

		for(int i=0;i<nj_n;i++){
			for(int j=0;j<nj_p;j++){
				double y=beta_coef(i,j,1,2*K)*sqrt(2*K+1);
				if(fabs(y)>1E-6)fprintf(fp, "%d\t%d\t%lf\n", i+1 + nj_p, j+1, y);
			}
		}
		fclose(fp);
	}
}
