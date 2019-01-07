

#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"struct.h"
#include"3j6j9j.h"
#include"input.h"

/*
nj_start:               starting point of single-j orbits
nj_end:			ending point of single-j orbits
sh:            		j values (doubled) of single-j orbits
pab_start:              start coordinate of "pab"
pab_end:                end coordinate of "pab"
pa/pb:			a/b value for "pab"
index_start:            start coordinate of index
K:                      angularmomentum (doubled) of transition operator
Fcoef_input:            matrix elements of transition operator
TBME:                   two-body matrix elements (2BME) of the hamiltonian
NME:                    new 2BME of the EWSR double commutator
*/

void W2(int nj_input, int nj_start_input, int nj_end_input, int *sh_input, 
	int num_pab_input, int pab_start_input, int pab_end_input, 
	int *pa_input, int *pb_input, 
	int index_start_input, 
	int K_input, double **Fcoef_input, 
	VJ **TBME, VJ **NME){
/*
	W2(abcd;J) = xxxx V(abce) xxxx, so pab, pcd, pce are all pp(nn/pn-np) types, so pab_start/pab_end define their ranges
	index_start is 0 for Vpp-even, but can be num_Vpp_even for Vpp-odd
*/
	#pragma parallel schedule(dynamic,1)
	#pragma omp parallel for
        for(int i=pab_start_input;i<pab_end_input;i++){

//--------------make local variables, to be used only in one thread, avoiding clashes between different threads
		int nj, nj_start, nj_end;
		int num_pab, pab_start, pab_end;
		int index_start;
		int K;
		#pragma omp critical
		{
			nj=nj_input;
			nj_start=nj_start_input;
			nj_end=nj_end_input;

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
//----------------------------------------------------------------------------------------------------------

                if(i!=0&&i%100==0)cout<<"\t\ti="<<i<<"\tin\t["<<pab_start<<","<<pab_end<<"]\n";
                for(int j=i;j<pab_end;j++){
		//pce: V(abce) can be out of the recorded V(abcd;J) with ellipsises, so j can be both </>= i

			int index=index_start + (i - pab_start) + (j - pab_start +1)*(j - pab_start)/2;
			int nJ=(*TBME[index]).numJ;
			for(int k=0;k<nJ;k++){
				int J= (*TBME[index]).Jmin + k * (*TBME[index]).dJ;
                                double V=(*TBME[index]).V[k];

				int case_ij_num=1;
				if(i!=j)case_ij_num=2;
				for(int case_ij=0;case_ij<case_ij_num;case_ij++){

					int pab, pce;
					int a,b,c,e;
					if(case_ij==0){
						pab=i;
						a=pa[i];
						b=pb[i];
						pce=j;
						c=pa[j];
						e=pb[j];
					}
					else{
						pab=j;
						a=pa[j];
						b=pb[j];
						pce=i;
						c=pa[i];
						e=pb[i];
					}
					int case_ce_num=1;
					if(c!=e)case_ce_num=2;
					for(int case_ce=0;case_ce<case_ce_num;case_ce++){

						double phase_ce=1;
						if(case_ce==1){
							int temp=c;
							c=e;
							e=temp;
							phase_ce*=theta(1+(sh[c]+sh[e]+J)/2);
						}
			
						for(int l=pab_start;l<pab_end;l++)// when l>=pab, produce W2(abcd;J); when l<pab, produce W5(cdab;J)=W2(abcd;J)
						if(pa[l]==c || pb[l]==c){// W2=(1 + PcdJ)xxx; so when pa[l]==c, collect xxx; when pb[l]==c, collect PcdJ xxx.
							double phase_cd=1;
							int d;
							if(pa[l]==c){
								d=pb[l];
							}
							else{
								d=pa[l];
								phase_cd=theta(1+(sh[c]+sh[d]+J)/2);
							}

							int index_new=index_start;
							if(l>=pab)index_new += pab - pab_start + (l - pab_start +1)*(l - pab_start)/2;//index in the new matrix elements
							else index_new += l - pab_start + (pab - pab_start +1)*(pab - pab_start)/2;//W5(cdab;J)
							for(int f=nj_start;f<nj_end;f++)
							if(abs(sh[c]-sh[f])<=J+K//two triangle rules: jc,jf,Jp; J,K,Jp
							&&abs(J-K)<=sh[c] + sh[f]
							&&abs(sh[e]-sh[f])<=K
							&&sh[e]+sh[f]>=K//triangle je,jf,K
							&&abs(sh[d]-sh[f])<=K
							&&sh[d]+sh[f]>=K){//triangle jd,jf,K;
								int Jp_min=abs(sh[c]-sh[f]);
								if(Jp_min<abs(J-K))Jp_min=abs(J-K);
								int Jp_max=sh[c]+sh[f];
								if(Jp_max>J+K)Jp_max=J+K;//two triangles rules: jc,jf,Jp;J,K,Jp
								for(int Jp=Jp_min;Jp<=Jp_max;Jp+=2)
								if(c!=f||Jp%4==0){
									double y1=V*(Jp+1);
									if(c==e)y1*=sqrt(2);
									if(c==d)y1/=sqrt(2);
									#pragma omp atomic	//fetching j66[][] should be atomic
									y1*=-0.5*Fcoef[e][f]*Fcoef[d][f]*j66[J/2][0][Jp/2][sh[f]/2][sh[c]/2][sh[e]/2]*j66[J/2][0][Jp/2][sh[f]/2][sh[c]/2][sh[d]/2];
									y1*=phase_ce;
									y1*=phase_cd;
									if(c==d)y1*=2;// (1+PcdJ)xxx
									if(l==pab)y1*=2;//W1 + PacPbdW1, so for W1(abab;J) it's 2W1(abab;J)

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
										(*NME[index_new]).V[new_address]+=y1;
									}
									//------------------------------------------------
								}
							}//for(f)
						}//for(l)
					}//for(case_ce)
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
