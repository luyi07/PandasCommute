

#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"struct.h"
#include"3j6j9j.h"
#include"input.h"

void W3(int nj_input, int *sh_input, 
	int pbe_start_input, int pbe_end_input, int index_start_bedf_input, 
	int num_pab1_input, int *pa1_input, int *pb1_input,
	int pab_start_input, int pab_end_input, int index_start_abcd_input,
	int num_pab2_input, int *pa2_input, int *pb2_input,
	int K_input, double **Fcoef_input, VJ **TBME, VJ **NME){
/*
	W3(abcd;J) = xxxx V(bedf) xxxx, so pbe, pdf are pp(nn/pn-np) types, so pab_start/pab_end define their ranges
	index_start is 0 for Vpp-even, but can be num_Vpp_even for Vpp-odd
*/
        #pragma parallel schedule(dynamic,1)
        #pragma omp parallel for
        for(int i=pbe_start_input;i<pbe_end_input;i++){// i=pbe
//-------------make parameters variable, to avoid omp threads clashing when fetching the same memory unit
		int nj;
		int pbe_start, pbe_end, index_start_bedf;
		int num_pab1;
		int pab_start, pab_end, index_start_abcd;
		int num_pab2;
		int K;
		#pragma omp critical
		{
			nj=nj_input;

			pbe_start=pbe_start_input;
			pbe_end=pbe_end_input;
			index_start_bedf=index_start_bedf_input;
			num_pab1=num_pab1_input;
			
			pab_start=pab_start_input;
			pab_end=pab_end_input;
			index_start_abcd=index_start_abcd_input;
			num_pab2=num_pab2_input;

			K=K_input;
		}
                int *sh=new int[nj];
                double **Fcoef=new double *[nj];
                for(int j=0;j<nj;j++)Fcoef[j]=new double [nj];
                int *pa1=new int[num_pab1];
                int *pb1=new int[num_pab1];
                int *pa2=new int[num_pab2];
                int *pb2=new int[num_pab2];
                #pragma omp criticle
                {
                        for(int j=0;j<nj;j++)sh[j]=sh_input[j];
                        for(int j=0;j<nj;j++)
                        for(int k=0;k<nj;k++){
                                Fcoef[j][k]=Fcoef_input[j][k];
                        }
                        for(int j=0;j<num_pab1;j++)pa1[j]=pa1_input[j];
                        for(int j=0;j<num_pab1;j++)pb1[j]=pb1_input[j];
                        for(int j=0;j<num_pab2;j++)pa2[j]=pa2_input[j];
                        for(int j=0;j<num_pab2;j++)pb2[j]=pb2_input[j];
                }
//------------------------------------------------------------------------
		if(i!=0&&i%100==0)cout<<"\t\ti="<<i<<"\tin\t["<<pbe_start<<","<<pbe_end<<"]\n";
//----------------------------------------------------------------------------------
		for(int j=i;j<pbe_end;j++){
			int index= index_start_bedf + (i - pbe_start) + (j - pbe_start +1)*(j - pbe_start)/2;
			int nJ= (*TBME[index]).numJ;
			for(int k=0;k<nJ;k++){
				int Jp=(*TBME[index]).Jmin + k * (*TBME[index]).dJ;
				double V=(*TBME[index]).V[k];

				int pbe, pdf, b,e,d,f;
				
				int case_ij_num=1;
				if(i!=j)case_ij_num=2;

				for(int case_ij=0;case_ij<case_ij_num;case_ij++){
					if(case_ij==0){
						pbe=i;
						pdf=j;
					}
					else{
						pbe=j;
						pdf=i;
					}

					int case_be_num=1, case_df_num=1;
					if(pa1[pbe]!=pb1[pbe])case_be_num=2;
					if(pa1[pdf]!=pb1[pdf])case_df_num=2;

					double phase_be=1, phase_df=1;
					for(int case_be=0;case_be<case_be_num;case_be++){

						if(case_be==1){
                                                        b=pa1[pbe];
                                                        e=pb1[pbe];
                                                        phase_be=1;
						}
						else{
							b=pb1[pbe];
							e=pa1[pbe];
							phase_be=theta(1+(sh[b]+sh[e]+Jp)/2);
						}

						for(int case_df=0;case_df<case_df_num;case_df++){
							if(case_df==0){
								d=pa1[pdf];
								f=pb1[pdf];
								phase_df=1;
							}
							else{
								d=pb1[pdf];
								f=pa1[pdf];
								phase_df=theta(1+(sh[d]+sh[f]+Jp)/2);
							}

							for(int l=pab_start;l<pab_end;l++)
							if(pa2[l]==b||pb2[l]==b){
								
								int a;
								if(pb2[l]==b){
									a=pa2[l];
								}
								else{
									a=pb2[l];
								}

								if(abs(sh[a]-sh[b])<=Jp+K// triangles ja,jb,J;Jp,K,J
								&&abs(Jp-K)<=sh[a]+sh[b])
								for(int m=l;m<pab_end;m++)
								// m>=l because W3(abcd;J) is in the shrinked pool.
								if(pa2[m]==d||pb2[m]==d){
		
									int c;
									if(pb2[m]==d){
										c=pa2[m];
									}
									else{
										c=pb2[m];
									}

									if(abs(sh[c]-sh[d])<=Jp+K// triangles jc,jd,J;Jp,K,J
									&&abs(Jp-K)<=sh[c]+sh[d]){
										int index_new=index_start_abcd
											+(l-pab_start)
											+(m-pab_start+1)*(m-pab_start)/2;

										int J_min=abs(sh[a]-sh[b]);// triangle ja, jb, J; Jp, K, J
										if(J_min<abs(Jp-K))J_min=abs(Jp-K);
										if(J_min<abs(sh[c]-sh[d]))J_min=abs(sh[c]-sh[d]);

										int J_max=sh[a]+sh[b];
										if(J_max>Jp+K)J_max=Jp+K;
										if(J_max>sh[c]+sh[d])J_max=sh[c]+sh[d];

										for(int J=J_min;J<=J_max;J+=2)//J
										if(a!=b||J%4==0)
										if(c!=d||J%4==0){
											double y1=V*(Jp+1);
											y1*=phase_be;
											y1*=phase_df;
											if(b==e)y1*=sqrt(2);
											if(d==f)y1*=sqrt(2);
											if(a==b)y1/=sqrt(2);
											if(c==d)y1/=sqrt(2);
											#pragma omp atomic// protect the process of reading j66
											y1*=Fcoef[e][a]*Fcoef[f][c] * j66[J/2][0][Jp/2][sh[e]/2][sh[b]/2][sh[a]/2] * j66[J/2][0][Jp/2][sh[f]/2][sh[d]/2][sh[c]/2];
											if(a==b)y1*=2;// (1+PabJ)xxx
											if(a>b)y1*=theta(1+(sh[a]+sh[b]+J)/2);
											if(c==d)y1*=2;// (1+PcdJ)xxx
											if(c>d)y1*=theta(1+(sh[c]+sh[d]+J)/2);
				
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
										}//for(J)
									}//if(c,d,Jp,K)
								}//for(m) if(pa2[m]==d)
							}//for(l) if(pa2[l])
                                                }//for(case_df)
                                        }//for(case_be)
                                }//for(case_ij)
			}//for(k)
		}//for(j)
                for(int i=0;i<nj;i++)delete [] Fcoef[i];
                delete [] Fcoef;
                delete [] sh;
                delete [] pa1;
                delete [] pb1;
                delete [] pa2;
                delete [] pb2;
        }//for(i)
}
