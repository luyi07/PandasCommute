

#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"struct.h"
#include"3j6j9j.h"
#include"input.h"

/*
nj:			number of single-j orbits
sh_input:		j values (doubled) of single-j orbits
pab_start:		start coordinate of "pab"
pab_end:		end coordinate of "pab"
num_pab:		total number of "pab"
pa_input/pb_input:	a/b value for "pab"
index_start:		start coordinate of index
K:			angularmomentum (doubled) of transition operator
Fcoef_input:		matrix elements of transition operator
TBME:			two-body matrix elements (2BME) of the hamiltonian
NME:			new 2BME of the EWSR double commutator
*/

void W1(int nj_input, int *sh_input, 
	int pab_start_input, int pab_end_input, 
	int num_pab_input, int *pa_input, int *pb_input, int index_start_input, 
	int K_input, double **Fcoef_input, 
	VJ **TBME, VJ **NME){
/*
	W1(abcd;J) = xxxx V(abef) xxxx, so pab, pcd, pef are all pp(nn/pn-np) types, so pab_start/pab_end define their ranges
	As a result, if Hpp is even, the double commutator can only end up to be H'pp-even by W1.
	In these two cases:
	Hpp-even, Fpp, Fpp -> H'pp-even
	Hpp-odd, Fpp, Fpp -> H'pp-odd
	I only need to set pab_start, pab_end, index_start accordingly	2018/1/14
*/
	#pragma parallel schedule(dynamic,1)
	#pragma omp parallel for
        for(int i=pab_start_input;i<pab_end_input;i++){// collect W1
//------make local variables, to be used only in one thread, avoiding clashes between different omp threads-------------
		int nj;
		int num_pab, pab_start, pab_end;
		int index_start;
		int K;
		#pragma omp criticle
		{
			nj=nj_input;
			num_pab=num_pab_input;
			pab_start=pab_start_input;
			pab_end=pab_end_input;
			index_start=index_start_input;
			K=K_input;
		}
		int *sh=new int[nj];
		double **Fcoef=new double *[nj];
		for(int j=0;j<nj;j++)Fcoef[j]=new double [nj];
		int *pa=new int[num_pab];
		int *pb=new int[num_pab];
		#pragma omp criticle
		{
			for(int j=0;j<nj;j++)sh[j]=sh_input[j];
			for(int j=0;j<nj;j++)
			for(int k=0;k<nj;k++){
				Fcoef[j][k]=Fcoef_input[j][k];
			}
			for(int j=0;j<num_pab;j++)pa[j]=pa_input[j];
			for(int j=0;j<num_pab;j++)pb[j]=pb_input[j];
		}
//----------------------------------------------------------------------------------------------------------------------
		int pab=i;
		if(i!=0&&i%100==0)cout<<"\t\ti="<<i<<"\tin\t["<<pab_start<<","<<pab_end<<"]\n";
                for(int j=i;j<pab_end;j++){//pef

			int case_ij_num=1;
			if(i!=j)case_ij_num=2;// in V(abef), pef < pab is also possible

			int index= index_start + (i - pab_start) + (j - pab_start +1)*(j - pab_start)/2;//every index is read only by one thread
			int nJ=(*TBME[index]).numJ;
			for(int k=0;k<nJ;k++){//J       here a specific V(abef;J) is located
				int J= (*TBME[index]).Jmin + k * (*TBME[index]).dJ;
				double V=(*TBME[index]).V[k];

				for(int case_ij=0;case_ij<case_ij_num;case_ij++){

					int a,b,e,f;
					int case_ef_num=1;
					if(case_ij==0){
						a=pa[i];
						b=pb[i];
						e=pa[j];
						f=pb[j];
						if(e!=f)case_ef_num=2;
					}
					else{
						e=pa[i];
						f=pb[i];
						a=pa[j];
						b=pb[j];
						if(e!=f)case_ef_num=2;
					}
					double y1=V;
					if(e==f)y1=y1*sqrt(2);//zeta(ef)
					if(J%4!=0)y1=-y1;
					
					for(int case_ef=0;case_ef<case_ef_num;case_ef++){//recover symmetry Pef V(abef)

						if(case_ef==1){
							int temp=e;
							e=f;
							f=temp;
							y1=y1*theta(1+(sh[e]+sh[f]+J)/2);
						}

						for(int l=pab_start;l<pab_end;l++){// pcd    a specific W(abcd;J) is located
							int c=pa[l];
							int d=pb[l];
							if(sh[c]+sh[d]>=J
							&&abs(sh[c]-sh[d])<=J
							&&(c!=d||J%4==0)){//triangle jc,jd,J

								double y2=y1;
								if(c==d)y2=y2*sqrt(2);// zeta(cd)^-1 * (1+PcdJ)
								if(case_ij==0&&l==i)y2=y2*2;//W1 + PacPbdW1, for W1(abab;J) it's 2W1(abab;J)
								if(case_ij==1&&l==j)y2=y2*2;

								double z=0;

								int case_cd_num=1;
								if(c!=d)case_cd_num=2;
								for(int case_cd=0;case_cd<case_cd_num;case_cd++){// (1+PcdJ) xxx

									if(case_cd==1){
										int temp=c;
										c=d;
										d=temp;
										y2=y2*theta(1+(sh[c]+sh[d]+J)/2);
									}

									if(abs(sh[d]-sh[e])<=J+K//two triangle rules: jd,je,Jp; J,K,Jp
									&&abs(J-K)<=sh[d] + sh[e]){//

										int Jp_min=abs(sh[d]-sh[e]);
										if(Jp_min<abs(J-K))Jp_min=abs(J-K);
										int Jp_max=sh[d]+sh[e];
										if(Jp_max>J+K)Jp_max=J+K;

										double y3=0;
										for(int Jp=Jp_min;Jp<=Jp_max;Jp+=2)//Jp
										if(d!=e||Jp%4==0){
											double temp;
											#pragma atomic
											temp=j66[J/2][0][Jp/2][sh[d]/2][sh[e]/2][sh[f]/2]*j66[ J/2 ][ 0 ][ Jp/2 ][ sh[e]/2 ][ sh[d]/2 ][ sh[c]/2 ];
											if(Jp%4==0)y3=y3+(Jp+1)*temp;
											else y3=y3-(Jp+1)*temp;
										}
										z=z-0.5*y2*y3*Fcoef[e][c]*Fcoef[f][d];
									}
								}//if(case_cd)

								int index_new=index_start;

								if(case_ij==0){
									if(l>=i)index_new = index_start + (i - pab_start) + (l - pab_start +1)*(l - pab_start)/2;
									else index_new = index_start + (l - pab_start) + (i - pab_start +1)*(i - pab_start)/2;
								}
								else{
									if(l>=j)index_new = index_start + (j - pab_start) + (l - pab_start +1)*(l - pab_start)/2;
									else index_new = index_start + (l - pab_start) + (j - pab_start +1)*(j - pab_start)/2;
								}

								//-----------------------------------------------
								//determine which V in (*NME[index_new]) to visit
								int JJmin=0, JJmax=1E6, dJJ=2;
								if(JJmin<abs(sh[a]-sh[b]))JJmin=abs(sh[a]-sh[b]);
								if(JJmin<abs(sh[c]-sh[d]))JJmin=abs(sh[c]-sh[d]);
								if(JJmax>sh[a]+sh[b])JJmax=sh[a]+sh[b];
								if(JJmax>sh[c]+sh[d])JJmax=sh[c]+sh[d];
								if(a==b||c==d){
									if(JJmin%4!=0)JJmin+=2;
									dJJ=4;
								}
								if( J>=JJmin && J<=JJmax && (J-JJmin)%dJJ ==0){
									int new_address= (J-JJmin)/dJJ;
									#pragma omp atomic
									(*NME[index_new]).V[new_address]+=z;
								}
								//------------------------------------------------
							}//if triangle jc,jd,J
						}//for(l)
					}//for(case_ef)
				}//for(case_ij)
			}//for(k)
		}//for(j)
		for(int i=0;i<nj;i++)delete [] Fcoef[i];
		delete [] Fcoef;
		delete [] sh;
		delete [] pa;
		delete [] pb;
	}//for(i)
}
