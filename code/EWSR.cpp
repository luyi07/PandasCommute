
#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"struct.h"
#include"3j6j9j.h"
#include"input.h"
#include"W.h"

/*
cal_EWSR_operator(..)	calculate \hat{O}_{EWSR} in the paper
file_sp:	filename of single particle infomation
file_GMEpn:	filename of general matrix elements
file_F:		filename for Fcoef
file_NMEpn:	filename for new matrix elements of \hat{O}_{EWSR}-----the *OUTPUT*
*/
void cal_EWSR_operator(char *file_sp, char *file_GMEpn, char *file_F, char *file_NMEpn){
//--------------read shell info--------------
	int nj_p, nj;
	read_nj(file_sp, nj_p, nj);
	int nj_n=nj-nj_p;

	int *nn=new int [nj];
	int *ll=new int [nj];
	int *sh=new int [nj];
	read_shell(file_sp, nj, nn, ll, sh);

	int *tz=new int [nj];
	int *parity=new int [nj];
	for(int i=0;i<nj;i++){
		if(i<nj_p)tz[i]=-1;//doubled
		else tz[i]=1;

		if(ll[i]%2==0)parity[i]=1;
		else parity[i]=-1;
	}
//--------------read in F
        int pF;//parity of one-body operator F
        double **Fcoef=new double *[nj];
        for(int i=0;i<nj;i++){
                Fcoef[i]=new double [nj];
                for(int j=0;j<nj;j++){
                        Fcoef[i][j]=0;
                }
        }
        read_F_new(file_F, sh, K, pF, Fcoef);//therein K is doubled
	cout<<"hahaha pF="<<pF<<endl;
	// for test only
	/*
	cout<<"Fcoef:"<<endl;
	for(int i=0;i<nj;i++){
		for(int j=0;j<nj;j++){
			cout<<Fcoef[i][j]<<",";
		}
		cout<<endl;
	}
	*/

	bool flag_Fcoef_pp=false, flag_Fcoef_nn=false, flag_Fcoef_pn=false;
	for(int i=0;i<nj_p;i++)
	for(int j=0;j<nj_p;j++)
	if(fabs(Fcoef[i][j])>1E-9)flag_Fcoef_pp=true;

	for(int i=nj_p;i<nj;i++)
	for(int j=nj_p;j<nj;j++)
	if(fabs(Fcoef[i][j])>1E-9)flag_Fcoef_nn=true;

	for(int i=0;i<nj_p;i++)
	for(int j=nj_p;j<nj;j++)
	if(fabs(Fcoef[i][j])>1E-9)flag_Fcoef_pn=true;
//--------------read general interactions-----------------
        int num_GMEpn;            //number of matrix elements
        read_GMEpn_num(file_GMEpn, num_GMEpn);

	double *spe=new double [nj];
        int *GME_a=new int [num_GMEpn];                   //sequence of ja
        int *GME_b=new int [num_GMEpn];                   //sequence of jb
        int *GME_c=new int [num_GMEpn];                   //sequence of jc
        int *GME_d=new int [num_GMEpn];                   //sequence of jd
        int *GME_I=new int [num_GMEpn];                   //2*I
        int *GME_T=new int [num_GMEpn];                   //2*T
        double *GME_V=new double [num_GMEpn];             //V(abcd;JT)
        read_GMEpn(file_GMEpn, nj, num_GMEpn, spe, GME_a, GME_b, GME_c, GME_d, GME_I, GME_T, GME_V);
/*
	// for test only
	for(int i=0;i<num_GMEpn;i++){
		cout<<GME_a[i]<<"\t";
		cout<<GME_b[i]<<"\t";
		cout<<GME_c[i]<<"\t";
		cout<<GME_d[i]<<"\t";
		cout<<GME_I[i]<<"\t";
		cout<<GME_V[i]<<"\t";
		cout<<endl;
	}
*/
//	cout<<"num_GMEpn="<<num_GMEpn<<endl;
//	nushell uses unnormalized pn 2bme, while cwj uses normalized ones. see "help.pdf" p14 of nushellx@msu
//	so we address that with a bool flag_pn_norm, when it's true, normalized; when it's false, unnormalized
	if(!flag_pn_norm){
		for(int i=0;i<num_GMEpn;i++){//unnormalized proton-neutron interactions need to be reduced before used here
			if( (GME_a[i]-1)/nj_p==0 && (GME_b[i]-1)/nj_p>0 ){// p-n interaction
				if( GME_a[i] == GME_b[i]-nj_n ) GME_V[i]*=sqrt(2);
				if( GME_c[i] == GME_d[i]-nj_n ) GME_V[i]*=sqrt(2);
				GME_V[i]/=2;
				// pn(aa;cc) = upn(aa;cc)
				// pn(ab;cc) = 1/sqrt(2) upn(ab;cc)
				// pn(ab;cd) = 1/2 upn(ab;cd)
			}
		}
	}
	time_t t_start=time(NULL);
//**************calculate one-body "interactions"*************
	double **H1=new double *[nj];
	for(int i=0;i<nj;i++){
		H1[i]=new double [nj];	
		for(int j=0;j<nj;j++){
			if(i==j)H1[i][j]=spe[i];
			else H1[i][j]=0;
		}
	}

	//for test only
	/*
	cout<<"H1:"<<endl;
	for(int i=0;i<nj;i++){
		for(int j=0;j<nj;j++){
			cout<<H1[i][j]<<",";
		}
		cout<<endl;
	}
	*/

	double **g=new double *[nj];
	for(int i=0;i<nj;i++){
		g[i]=new double [nj];
		for(int j=0;j<nj;j++){
			g[i][j]=0;
		}
	}

	for(int a=0;a<nj;a++)
	for(int b=0;b<nj;b++)
	if(sh[a]==sh[b]){
		double y=0;
		for(int c=0;c<nj;c++)
		for(int d=0;d<nj;d++){
			y+=	- H1[a][c] * Fcoef[c][d] * Fcoef[b][d]
				+ Fcoef[a][c] * H1[c][d] * Fcoef[b][d]
				+ Fcoef[c][a] * H1[c][d] * Fcoef[d][b]
				- Fcoef[c][a] * Fcoef[c][d] * H1[d][b];
		}

		y *= 1.0/2/(sh[a]+1);
		g[a][b]=y;
	}
	cout<<"\tEWSR:\tg values finished\t";
	print_time();
	time_t t_gab=time(NULL);
//************** two-body parts of O_EWSR ****************
	int num_NME=0;
//*************************************************************
//	open memory for proton-proton matrix elements
//*************************************************************
//------	proton-proton: pab	----------------------
	int num_pab_pp = nj_p * (nj_p+1)/2;
	int *pa_pp=new int [num_pab_pp];
	int *pb_pp=new int [num_pab_pp];
	int temp=0;
	for(int i=0;i<nj_p;i++){//pab-even
		for(int j=i;j<nj_p;j++){
			if(parity[i]*parity[j]==1){
				pa_pp[temp]=i;
				pb_pp[temp]=j;
				temp++;
			}
		}
	}
	int num_pab_pp_even=temp;
	//cout<<"num_pab_pp_even="<<num_pab_pp_even<<endl;
	for(int i=0;i<nj_p;i++){//pab-odd
		for(int j=i;j<nj_p;j++){
			if(parity[i]*parity[j]==-1){
				pa_pp[temp]=i;
				pb_pp[temp]=j;
				temp++;
			}
		}
	}
	int num_pab_pp_odd = num_pab_pp - num_pab_pp_even;
//	for test only
/*
	cout<<endl;
	cout<<"proton-proton pab:"<<endl;
	for(int i=0;i<num_pab_pp;i++){
		cout<<i+1<<"\t"<<pa_pp[i]+1<<"\t"<<pb_pp[i]+1<<"\n";
	}
	cout<<endl;
*/
//------open memory for TBMEpp (Hpp)-----
	int num_Vpp_even = (num_pab_pp_even +1)*num_pab_pp_even/2;
	int num_Vpp =	(num_pab_pp_even +1)*num_pab_pp_even/2
			+ (num_pab_pp_odd + 1)*num_pab_pp_odd/2;

	VJ **TBMEpp=new VJ *[num_Vpp];// even-parity & odd-parity interactions
/*
	for(int i=0;i<num_Vpp;i++){
		TBMEpp[i].numJ=0;
	}
*/
	for(int odevity=1;odevity>=-1;odevity-=2){
		int pab_start,pab_end,index_start;
		if(odevity==1){
			pab_start=0;
			pab_end=num_pab_pp_even;
			index_start=0;
		}
		else{
			pab_start=num_pab_pp_even;
			pab_end=num_pab_pp;
			index_start=num_Vpp_even;
		}
		for(int i=pab_start;i<pab_end;i++){
		// index = i + (j+1)*j/2
		// (j+1)*j/2 - j*(j-1)/2 = j, so i+(j+1)*j/2 (i=0,...j) just fill up that gap
			int a=pa_pp[i], b=pb_pp[i];
			int Jabmin=abs(sh[a]-sh[b]), Jabmax=sh[a]+sh[b];
			
			for(int j=i;j<pab_end;j++){
				
				int index= index_start+ (i-pab_start) + (j - pab_start +1)*(j - pab_start)/2;

				int c=pa_pp[j], d=pb_pp[j];
				int Jcdmin=abs(sh[c]-sh[d]), Jcdmax=sh[c]+sh[d];
				
				int JJmin=0,JJmax=1E6;
				if(JJmin<Jabmin)JJmin=Jabmin;
				if(JJmax>Jabmax)JJmax=Jabmax;
				if(JJmin<Jcdmin)JJmin=Jcdmin;
				if(JJmax>Jcdmax)JJmax=Jcdmax;

				if((a==b||c==d)){
					if(JJmin%4!=0)JJmin+=2;//Pauli's principle
					if(JJmax%4!=0)JJmax-=2;
				}

				TBMEpp[index]=new VJ;//allocate dynamic memory: a new VJ

				if(JJmax>=JJmin){

					int numJJ=0;
					if((a==b||c==d)){
						numJJ=(JJmax-JJmin)/4+1;
					}
					else{
						numJJ=(JJmax-JJmin)/2+1;
					}

					(*TBMEpp[index]).numJ=numJJ;
					(*TBMEpp[index]).Jmin=JJmin;
					if(numJJ==1)(*TBMEpp[index]).dJ=0;
					else{
						(*TBMEpp[index]).dJ=(JJmax-JJmin)/(numJJ-1);
					}
					(*TBMEpp[index]).V=new float [numJJ];
					for(int k=0;k<numJJ;k++){
						(*TBMEpp[index]).V[k]=0;
					}
					num_NME+=numJJ;
				}
				else{
					(*TBMEpp[index]).numJ=0;
				}
			}
		}
	}
//	for test only
/*
	cout<<"memory addresses!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	cout<<&(*TBMEpp[0]).numJ<<"\t"<<&(*TBMEpp[0]).Jmin<<"\t"<<&(*TBMEpp[0]).dJ<<"\t"<<&(*TBMEpp[0]).V[0]<<endl;
*/
//------open memory for new matrix elements (proton-proton)-------
        //cout<<"num_Vpp="<<num_Vpp<<endl;
        VJ **NMEpp=new VJ * [num_Vpp];
        for(int i=0;i<num_Vpp;i++){
		NMEpp[i]=new VJ;
                int nJ=(*TBMEpp[i]).numJ;
                (*NMEpp[i]).numJ=nJ;
                (*NMEpp[i]).Jmin=(*TBMEpp[i]).Jmin;
                (*NMEpp[i]).dJ=(*TBMEpp[i]).dJ;
                if(nJ>0){
                        (*NMEpp[i]).V=new float [nJ];
                        for(int j=0;j<nJ;j++){
                                (*NMEpp[i]).V[j]=0;
                        }
                }//if(nJ)
        }
//*************************************************************
//****** open memory for proton-neutron matrix elements *******
//*************************************************************
//------        proton-neutron: pab      ----------------------
        int num_pab_pn = nj_p * nj_n;
        int *pa_pn=new int [num_pab_pn];
        int *pb_pn=new int [num_pab_pn];
        temp=0;
        for(int i=0;i<nj_p;i++){//pab-even
                for(int j=nj_p;j<nj;j++){
                        if(parity[i]*parity[j]==1){
                                pa_pn[temp]=i;
                                pb_pn[temp]=j;
                                temp++;
                        }
                }
        }
        int num_pab_pn_even=temp;
        //cout<<"\tnum_pab_pn_even="<<num_pab_pn_even<<endl;
        for(int i=0;i<nj_p;i++){//pab-odd
                for(int j=nj_p;j<nj;j++){
                        if(parity[i]*parity[j]==-1){
                                pa_pn[temp]=i;
                                pb_pn[temp]=j;
                                temp++;
                        }
                }
        }
        int num_pab_pn_odd = num_pab_pn - num_pab_pn_even;
//      for test only
	/*
        cout<<endl;
        cout<<"proton-neutron pab:"<<endl;
        for(int i=0;i<num_pab_pn;i++){
                cout<<i+1<<"\t"<<pa_pn[i]+1<<"\t"<<pb_pn[i]+1<<"\n";
        }
        cout<<endl;
	*/
//------open memory for TBMEpn (Hpn)-----
	int num_Vpn_even = (num_pab_pn_even +1)*num_pab_pn_even/2;
	int num_Vpn =	(num_pab_pn_even +1)*num_pab_pn_even/2
			+ (num_pab_pn_odd + 1)*num_pab_pn_odd/2;

	VJ **TBMEpn=new VJ * [num_Vpn];// even-parity & odd-parity interactions
	for(int odevity=1;odevity>=-1;odevity-=2){
		int pab_start,pab_end,index_start;// to be set as [0,num_pab_pn_even],0 or [num_pab_pn_even,num_pab_pn],num_Vpn_even
		if(odevity==1){
			pab_start=0;
			pab_end=num_pab_pn_even;
			index_start=0;
		}
		else{
			pab_start=num_pab_pn_even;
			pab_end=num_pab_pn;
			index_start=num_Vpn_even;
		}
		for(int i=pab_start;i<pab_end;i++){
		// index = i + (j+1)*j/2
		// (j+1)*j/2 - j*(j-1)/2 = j, so i+(j+1)*j/2 (i=0,...j) just fill up that gap
			int a=pa_pn[i], b=pb_pn[i];
			int Jabmin=abs(sh[a]-sh[b]), Jabmax=sh[a]+sh[b];
			
			for(int j=i;j<pab_end;j++){
				
				int index= index_start+ (i-pab_start) + (j - pab_start +1)*(j - pab_start)/2;

				int c=pa_pn[j], d=pb_pn[j];
				int Jcdmin=abs(sh[c]-sh[d]), Jcdmax=sh[c]+sh[d];
				
				int JJmin=0,JJmax=1E6;
				if(JJmin<Jabmin)JJmin=Jabmin;
				if(JJmax>Jabmax)JJmax=Jabmax;
				if(JJmin<Jcdmin)JJmin=Jcdmin;
				if(JJmax>Jcdmax)JJmax=Jcdmax;

				if((a==b||c==d)){
					if(JJmin%4!=0)JJmin+=2;//Pauli's principle
					if(JJmax%4!=0)JJmax-=2;
				}

				TBMEpn[index]=new VJ;
				if(JJmax>=JJmin){

					int numJJ=0;
					if((a==b||c==d)){
						numJJ=(JJmax-JJmin)/4+1;
					}
					else{
						numJJ=(JJmax-JJmin)/2+1;
					}

					(*TBMEpn[index]).numJ=numJJ;
					(*TBMEpn[index]).Jmin=JJmin;
					if(numJJ==1)(*TBMEpn[index]).dJ=0;
					else{
						(*TBMEpn[index]).dJ=(JJmax-JJmin)/(numJJ-1);
					}
					(*TBMEpn[index]).V=new float [numJJ];
					for(int k=0;k<numJJ;k++){
						(*TBMEpn[index]).V[k]=0;
				
						if(flag_pn_norm)num_NME++;
						else{
							int multiple=1;
							if(pa_pn[i]!=pb_pn[i]-nj_p && pa_pn[j]!=pb_pn[j]-nj_p){//(a,b) (c,d) are on different orbits
								multiple*=2;//T=0,1
								if(i!=j && ( pa_pn[i]==pb_pn[j]-nj_p && pa_pn[j]==pb_pn[i]-nj_p) )//V(abba), V(baab)
								multiple+=2;
							}
							num_NME+=multiple;
						}
					}
				}
				else{
					(*TBMEpn[index]).numJ=0;
				}
			}
		}
	}
	//cout<<"\tmemory for TBMEpn is opened"<<endl;
//------open memory for new matrix elements (neutron-neutron)-----
        //cout<<"\tnum_Vpn="<<num_Vpn<<endl;
        VJ **NMEpn=new VJ * [num_Vpn];
        for(int i=0;i<num_Vpn;i++){
		NMEpn[i]=new VJ;
                int nJ=(*TBMEpn[i]).numJ;
                (*NMEpn[i]).numJ=nJ;
                (*NMEpn[i]).Jmin=(*TBMEpn[i]).Jmin;
                (*NMEpn[i]).dJ=(*TBMEpn[i]).dJ;
                if(nJ>0){
                        (*NMEpn[i]).V=new float [nJ];
                        for(int j=0;j<nJ;j++){
                                (*NMEpn[i]).V[j]=0;
                        }
                }//if(nJ)
        }
//*************************************************************
//****** open memory for neutron-neutron matrix elements *******
//*************************************************************
//------        neutron-neutron: pab      ----------------------
        int num_pab_nn = nj_n * (nj_n+1)/2;
        int *pa_nn=new int [num_pab_nn];
        int *pb_nn=new int [num_pab_nn];
        temp=0;
        for(int i=nj_p;i<nj;i++){//pab-even
                for(int j=i;j<nj;j++){
                        if(parity[i]*parity[j]==1){
                                pa_nn[temp]=i;
                                pb_nn[temp]=j;
                                temp++;
                        }
                }
        }
        int num_pab_nn_even=temp;
        //cout<<"\tnum_pab_nn_even="<<num_pab_nn_even<<endl;
        for(int i=nj_p;i<nj;i++){//pab-odd
                for(int j=i;j<nj;j++){
                        if(parity[i]*parity[j]==-1){
                                pa_nn[temp]=i;
                                pb_nn[temp]=j;
                                temp++;
                        }
                }
        }
        int num_pab_nn_odd = num_pab_nn - num_pab_nn_even;
//      for test only
/*
        cout<<endl;
        cout<<"neutron-neutron pab:"<<endl;
        for(int i=0;i<num_pab_nn;i++){
                cout<<i+1<<"\t"<<pa_nn[i]+1<<"\t"<<pb_nn[i]+1<<"\n";
        }
        cout<<endl;
*/
//------open memory for TBMEnn (Hnn)-----
	int num_Vnn_even = (num_pab_nn_even +1)*num_pab_nn_even/2;
	int num_Vnn =	(num_pab_nn_even +1)*num_pab_nn_even/2
			+ (num_pab_nn_odd + 1)*num_pab_nn_odd/2;

	VJ **TBMEnn=new VJ * [num_Vnn];// even-parity & odd-parity interactions
	for(int odevity=1;odevity>=-1;odevity-=2){
		int pab_start,pab_end,index_start;// to be set as [0,num_pab_nn_even],0 or [num_pab_nn_even,num_pab_nn],num_Vnn_even
		if(odevity==1){
			pab_start=0;
			pab_end=num_pab_nn_even;
			index_start=0;
		}
		else{
			pab_start=num_pab_nn_even;
			pab_end=num_pab_nn;
			index_start=num_Vnn_even;
		}
		for(int i=pab_start;i<pab_end;i++){
		// index = i + (j+1)*j/2
		// (j+1)*j/2 - j*(j-1)/2 = j, so i+(j+1)*j/2 (i=0,...j) just fill up that gap
			int a=pa_nn[i], b=pb_nn[i];
			int Jabmin=abs(sh[a]-sh[b]), Jabmax=sh[a]+sh[b];
			
			for(int j=i;j<pab_end;j++){
				
				int index= index_start+ (i-pab_start) + (j - pab_start +1)*(j - pab_start)/2;

				int c=pa_nn[j], d=pb_nn[j];
				int Jcdmin=abs(sh[c]-sh[d]), Jcdmax=sh[c]+sh[d];
				
				int JJmin=0,JJmax=1E6;
				if(JJmin<Jabmin)JJmin=Jabmin;
				if(JJmax>Jabmax)JJmax=Jabmax;
				if(JJmin<Jcdmin)JJmin=Jcdmin;
				if(JJmax>Jcdmax)JJmax=Jcdmax;

				if((a==b||c==d)){
					if(JJmin%4!=0)JJmin+=2;//Pauli's principle
					if(JJmax%4!=0)JJmax-=2;
				}

				TBMEnn[index]=new VJ;
				if(JJmax>=JJmin){

					int numJJ=0;
					if((a==b||c==d)){
						numJJ=(JJmax-JJmin)/4+1;
					}
					else{
						numJJ=(JJmax-JJmin)/2+1;
					}

					(*TBMEnn[index]).numJ=numJJ;
					(*TBMEnn[index]).Jmin=JJmin;
					if(numJJ==1)(*TBMEnn[index]).dJ=0;
					else{
						(*TBMEnn[index]).dJ=(JJmax-JJmin)/(numJJ-1);
					}
					(*TBMEnn[index]).V=new float [numJJ];
					for(int k=0;k<numJJ;k++){
						(*TBMEnn[index]).V[k]=0;
					}
					num_NME+=numJJ;
				}
				else{
					(*TBMEnn[index]).numJ=0;
				}
			}
		}
	}
	//cout<<"\tmemory for TBMEnn is opened"<<endl;
//------open memory for new matrix elements (neutron-neutron)-------
        //cout<<"\tnum_Vnn="<<num_Vnn<<endl;
        VJ **NMEnn=new VJ *[num_Vnn];
        for(int i=0;i<num_Vnn;i++){
		NMEnn[i]=new VJ;
                int nJ=(*TBMEnn[i]).numJ;
                (*NMEnn[i]).numJ=nJ;
                (*NMEnn[i]).Jmin=(*TBMEnn[i]).Jmin;
                (*NMEnn[i]).dJ=(*TBMEnn[i]).dJ;
                if(nJ>0){
                        (*NMEnn[i]).V=new float [nJ];
                        for(int j=0;j<nJ;j++){
                                (*NMEnn[i]).V[j]=0;
                        }
                }//if(nJ)
        }
//********************************************************************
//	assign interactions read from file_GMEpn into TBMEpp, TBMEnn, TBMEpn
//********************************************************************
	cout<<"\tEWSR:\tInput Matrix elements are all collected\t";
	print_time();
	time_t t_GME=time(NULL);
//	******** assign values for TBMEpp ***********
	for(int i=0;i<num_GMEpn;i++){// complexity = num_GMEpn * num_pab
		int index_pab, index_pcd, index;
		int a=GME_a[i]-1;
		int b=GME_b[i]-1;
		int c=GME_c[i]-1;
		int d=GME_d[i]-1;
		int J=2*GME_I[i];
		double V=GME_V[i];
		int parity_V=parity[a]*parity[b];
		if(a<nj_p){
			if(b<nj_p){
				double phase=1;
				for(int j=0;j<num_pab_pp;j++){// find pab
					if(a==pa_pp[j]&&b==pb_pp[j]){
						index_pab=j;
					}
					else if(a==pb_pp[j]&&b==pa_pp[j]){//if b<a, I exchange them and induce an additional phase factor
						index_pab=j;
						phase=theta(1+(sh[a]+sh[b]+J)/2);
					}
				}
				for(int j=0;j<num_pab_pp;j++){// find pcd
					if(c==pa_pp[j]&&d==pb_pp[j]){
						index_pcd=j;
					}
					else if(c==pb_pp[j]&&d==pa_pp[j]){
						index_pcd=j;
						phase*=theta(1+(sh[c]+sh[d]+J)/2);
					}
				}
				if(index_pab>index_pcd){// in the input file GMEpn.int index_pab is not necessarily less than or equal to index_pcd, so we have to exchange them
					int temp=index_pab;
					index_pab=index_pcd;
					index_pcd=temp;
				}
				if(parity_V==1)index = index_pab + (index_pcd +1)*index_pcd/2;// even-parity interactions
				else index= num_Vpp_even + 
					(index_pab - num_pab_pp_even) 
					+ (index_pcd - num_pab_pp_even + 1)*(index_pcd - num_pab_pp_even)/2;// odd-parity interactions
				int numJJ=(*TBMEpp[index]).numJ;
				int JJmin=(*TBMEpp[index]).Jmin;
				int dJJ=(*TBMEpp[index]).dJ;
				if(J==JJmin) (*TBMEpp[index]).V[0]=V*phase;
				else if( (J-JJmin)%dJJ !=0){ 
					cout<<"Ah!!!!!! FIRE!!!! THE input GMEpn.int does not respect Pauli's Principle!"<<endl;
				}
				else if(J>JJmin && J< JJmin + dJJ * numJJ){
					int j=(J-JJmin)/dJJ;
					(*TBMEpp[index]).V[j]=V*phase;
				}
			}
			else{
				for(int j=0;j<num_pab_pn;j++){//find pab
					if(a==pa_pn[j]&&b==pb_pn[j]){
						index_pab=j;
					}
				}
				for(int j=0;j<num_pab_pn;j++){//find pcd
					if(c==pa_pn[j]&&d==pb_pn[j]){
						index_pcd=j;
					}
				}
				bool flag_abba=false;
                                if(index_pab>index_pcd){// in the input file GMEpn.int index_pab is not necessarily less than or equal to index_pcd, so we have to exchange them
                                        int temp=index_pab;
                                        index_pab=index_pcd;
                                        index_pcd=temp;
                                        if(!flag_pn_norm && a==d-nj_p && c==b-nj_p)flag_abba=true;
					//in NushellX format, p-n matrix elements V(abba;J) and V(baab;J) both appear in the interaction file.
					//we should collect only one of those, so I set flag_abba to block the other of those.
                                }
				if(!flag_abba){//if V(abba) is collected, V(baab) is blocked outside
					if(parity_V==1)index = index_pab + (index_pcd +1)*index_pcd/2;// find index
					else index= num_Vpn_even + 
						(index_pab - num_pab_pn_even) 
						+ (index_pcd - num_pab_pn_even + 1)*(index_pcd - num_pab_pn_even)/2;
					//odd-parity interaction
					int numJJ=(*TBMEpn[index]).numJ;
					int JJmin=(*TBMEpn[index]).Jmin;
					int dJJ=(*TBMEpn[index]).dJ;
					if(J==JJmin) (*TBMEpn[index]).V[0]+=V;//in NushellX format, one (TBMEpn[index]).V needs to collect several matrix elements, so here it is a "+=" rather than a "=". 2018/2/7
					else if( (J-JJmin)%dJJ !=0){
						cout<<"Ah!!!!!! FIRE!!!! THE input GMEpn.int does not respect Pauli's Principle!"<<endl;
					}
					else if(J>JJmin && J< JJmin + dJJ * numJJ){
						int j=(J-JJmin)/dJJ;
						(*TBMEpn[index]).V[j]+=V;//similar to the upper annotation, it is a "+=" here, rather than "=".
					}
					else{
						cout<<"Ah!!!!!! FIRE!!!! THE input GMEpn.int does not respect triangle rule!"<<endl;
					}
				}
			}
		}
		else{
			double phase=1;
			for(int j=0;j<num_pab_nn;j++){// find pab
				if(a==pa_nn[j]&&b==pb_nn[j]){
					index_pab=j;
				}
				else if(a==pb_nn[j]&&b==pa_nn[j]){//if b<a, I exchange them and induce an additional phase factor
					index_pab=j;
					phase=theta(1+(sh[a]+sh[b]+J)/2);
				}
			}
			for(int j=0;j<num_pab_nn;j++){// find pcd
				if(c==pa_nn[j]&&d==pb_nn[j]){
					index_pcd=j;
				}
				else if(c==pb_nn[j]&&d==pa_nn[j]){
					index_pcd=j;
					phase*=theta(1+(sh[c]+sh[d]+J)/2);
				}
			}
			if(index_pab>index_pcd){// in the input file GMEpn.int index_pab is not necessarily less than or equal to index_pcd, so we have to exchange them
				int temp=index_pab;
				index_pab=index_pcd;
				index_pcd=temp;
			}
			if(parity_V==1)index = index_pab + (index_pcd +1)*index_pcd/2;// even-parity interactions
			else index= num_Vnn_even +
				(index_pab - num_pab_nn_even)
				+ (index_pcd - num_pab_nn_even + 1)*(index_pcd - num_pab_nn_even)/2;// odd-parity interactions
			int numJJ=(*TBMEnn[index]).numJ;
			int JJmin=(*TBMEnn[index]).Jmin;
			int dJJ=(*TBMEnn[index]).dJ;
			if(J==JJmin) (*TBMEnn[index]).V[0]=V*phase;
			else if( (J-JJmin)%dJJ !=0){
				cout<<"Ah!!!!!! FIRE!!!! THE input GMEpn.int does not respect Pauli's Principle!"<<endl;
			}
			else if(J>JJmin && J< JJmin + dJJ * numJJ){
				int j=(J-JJmin)/dJJ;
				(*TBMEnn[index]).V[j]=V*phase;
			}
		}
	}
/*
	cout<<"check Vpp"<<endl;
	for(int i=0;i<num_Vpp;i++){
		for(int j=0;j<TBMEpp[i].numJ;j++){
			cout<<"\t"<<TBMEpp[i].V[j]<<"\n";
		}
	}
	cout<<"check Vnn"<<endl;
	for(int i=0;i<num_Vnn;i++){
		for(int j=0;j<(*TBMEnn[i]).numJ;j++){
			cout<<"\t"<<(*TBMEnn[i]).V[j]<<"\n";
		}
	}
	cout<<"check Vpn"<<endl;
	for(int i=0;i<num_Vpn;i++){
		for(int j=0;j<(*TBMEpn[i]).numJ;j++){
			cout<<"\t"<<(*TBMEpn[i]).V[j]<<"\n";
		}
	}
*/
//	assign value for table_num = jmax + 1/2
	table_num=1;
	for(int i=0;i<nj;i++)
	if(table_num<(sh[i]+1)/2)table_num=(sh[i]+1)/2;
	j6table(K);
	flag_j6table=true;
	cout<<"\tEWSR:\t6j table is done, start to collect W1&W4\t"; print_time();
	time_t t_6j=time(NULL);
//	collect W1	----------------------------------------
//------ calculate TBMEpp -> NMEpp -----------------------------
	if(flag_Fcoef_pp){
		W1(nj, sh, 0, num_pab_pp_even, num_pab_pp, pa_pp, pb_pp, 0, K, Fcoef, TBMEpp, NMEpp);
		//even-parity interactions:	pab_start=0
		//				pab_end=num_pab_pp_even
		//				index_start=0
		W1(nj, sh, num_pab_pp_even, num_pab_pp, num_pab_pp, pa_pp, pb_pp, num_Vpp_even, K, Fcoef, TBMEpp, NMEpp);
		//odd-parity interactions:	pab_start=num_pab_pp_even
		//				pab_end=num_pab_pp
		//				index_start=num_Vpp_even
	}
	cout<<"\t\tW1: Hpp -> H'pp finished\t"; print_time();
	if(flag_Fcoef_nn){
		W1(nj, sh, 0, num_pab_nn_even, num_pab_nn, pa_nn, pb_nn, 0, K, Fcoef, TBMEnn, NMEnn);
		W1(nj, sh, num_pab_nn_even, num_pab_nn, num_pab_nn, pa_nn, pb_nn, num_Vnn_even, K, Fcoef, TBMEnn, NMEnn);
	}
	cout<<"\t\tW1: Hnn -> H'nn finished\t"; print_time();
	W1(nj, sh, 0, num_pab_pn_even, num_pab_pn, pa_pn, pb_pn, 0, K, Fcoef, TBMEpn, NMEpn);
	W1(nj, sh, num_pab_pn_even, num_pab_pn, num_pab_pn, pa_pn, pb_pn, num_Vpn_even, K, Fcoef, TBMEpn, NMEpn);
	cout<<"\t\tW1: Hpn -> H'pn finished\t"; print_time();
	cout<<"\tEWSR:\tW1&W4 is collected,start to collect W2&W5\t"; print_time();
	time_t t_W1=time(NULL);
//	collect W2	----------------------------------------
	//	Hpp -> H'pp
	if(flag_Fcoef_pp){
		W2(nj, 0, nj_p, sh, num_pab_pp, 0, num_pab_pp_even, pa_pp, pb_pp, 0, K, Fcoef, TBMEpp, NMEpp);
		W2(nj, 0, nj_p, sh, num_pab_pp, num_pab_pp_even, num_pab_pp, pa_pp, pb_pp, num_Vpp_even, K, Fcoef, TBMEpp, NMEpp);
	}
	if(flag_Fcoef_pn){
		W2(nj, nj_p, nj, sh, num_pab_pp, 0, num_pab_pp_even, pa_pp, pb_pp, 0, K, Fcoef, TBMEpp, NMEpp);
		W2(nj, nj_p, nj, sh, num_pab_pp, num_pab_pp_even, num_pab_pp, pa_pp, pb_pp, num_Vpp_even, K, Fcoef, TBMEpp, NMEpp);
	}
	cout<<"\t\tW2: Hpp -> H'pp finished\t"; print_time();
	//	Hnn -> H'nn
	if(flag_Fcoef_nn){
		W2(nj, nj_p, nj, sh, num_pab_nn, 0, num_pab_nn_even, pa_nn, pb_nn, 0, K, Fcoef, TBMEnn, NMEnn);
		W2(nj, nj_p, nj, sh, num_pab_nn, num_pab_nn_even, num_pab_nn, pa_nn, pb_nn, num_Vnn_even, K, Fcoef, TBMEnn, NMEnn);
	}
	if(flag_Fcoef_pn){
		W2(nj, 0, nj_p, sh, num_pab_nn, 0, num_pab_nn_even, pa_nn, pb_nn, 0, K, Fcoef, TBMEnn, NMEnn);
		W2(nj, 0, nj_p, sh, num_pab_nn, num_pab_nn_even, num_pab_nn, pa_nn, pb_nn, num_Vnn_even, K, Fcoef, TBMEnn, NMEnn);
	}
	cout<<"\t\tW2: Hnn -> H'nn finished\t"; print_time();
	//	Hpn -> H'pn
	//cout<<"\tnum_pab_pn_even="<<num_pab_pn_even<<endl;
	//cout<<"\tnum_pab_pn="<<num_pab_pn<<endl;
	W2(nj, 0, nj, sh, num_pab_pn, 0, num_pab_pn_even, pa_pn, pb_pn, 0, K, Fcoef, TBMEpn, NMEpn);//even-parity
	W2(nj, 0, nj, sh, num_pab_pn, num_pab_pn_even, num_pab_pn, pa_pn, pb_pn, num_Vpn_even, K, Fcoef, TBMEpn, NMEpn);//odd-parity
	cout<<"\t\tW2: Hpn -> H'pn finished\t"; print_time();
	cout<<"\tEWSR:\tW2&W5 is collected, start to collect W3\t"; print_time();
	time_t t_W2=time(NULL);

//-----------	collect W3	--------------------------------
	if(pF==1){
//------------	Hpp -> H'pp/H'pn	------------------------
//		Hpp -> H'pp
		if(flag_Fcoef_pp){
			W3(nj, sh, 
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp, 
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpp, NMEpp);

			W3(nj, sh,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpp, NMEpp);
		}
		//cout<<"\tW3:	Hpp -> H'pp is finished"<<endl;
//		Hpp -> H'pn
		if(flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEpp, NMEpn);

			W3(nj, sh,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEpp, NMEpn);
		}
		cout<<"\t\tW3: Hpp -> H'pp/H'pn finished\t"; print_time();
		//cout<<"\tW3:	Hpp -> H'pn is finished"<<endl;
//------------	Hnn -> H'nn/H'pn	------------------------
//		Hnn -> H'nn
		if(flag_Fcoef_nn){
			W3(nj, sh, 
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEnn, NMEnn);

			W3(nj, sh,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEnn, NMEnn);
		}
		//cout<<"\tW3:	Hnn -> H'nn is finished"<<endl;
//		Hnn -> H'pn
		if(flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEnn, NMEpn);

			W3(nj, sh,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEnn, NMEpn);
		}
		cout<<"\t\tW3: Hnn -> H'nn/H'pn finished\t"; print_time();
		//cout<<"\tW3:	Hnn -> H'pn is finished"<<endl;
//------------	Hpn -> H'pp/H'nn/H'pn	------------------------
//		Hpn -> H'pn
		W3(nj, sh, 
		0, num_pab_pn_even, 0, 
		num_pab_pn, pa_pn, pb_pn,
		0, num_pab_pn_even, 0, 
		num_pab_pn, pa_pn, pb_pn,
		K, Fcoef, TBMEpn, NMEpn);

		W3(nj, sh,
		num_pab_pn_even, num_pab_pn, num_Vpn_even, 
		num_pab_pn, pa_pn, pb_pn,
		num_pab_pn_even, num_pab_pn, num_Vpn_even, 
		num_pab_pn, pa_pn, pb_pn,
		K, Fcoef, TBMEpn, NMEpn);
		//cout<<"\tW3:	Hpn -> H'pn is finished"<<endl;

//		Hpn -> H'pp
		if(flag_Fcoef_pp||flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpn, NMEpp);

			W3(nj, sh,
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpn, NMEpp);
		}
		//cout<<"\tW3:	Hpn -> H'pp is finished"<<endl;
//		Hpn -> H'nn
		if(flag_Fcoef_nn||flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEpn, NMEnn);

			W3(nj, sh,
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEpn, NMEnn);
		}
		cout<<"\t\tW3: Hpn -> H'pp/H'nn/H'pn finished\t"; print_time();
	}
	else{
//		Hpp -> H'pp
		if(flag_Fcoef_pp){
			W3(nj, sh, 
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpp, NMEpp);

			W3(nj, sh,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpp, NMEpp);
		}
		//cout<<"\tW3:	Hpp -> H'pp is finished"<<endl;
//		Hpp -> H'pn
		if(flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp, 
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEpp, NMEpn);

			W3(nj, sh,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn, 
			K, Fcoef, TBMEpp, NMEpn);
		}
		cout<<"\t\tW3: Hpp -> H'pp/H'pn finished\t"; print_time();
//		Hnn -> H'nn
		if(flag_Fcoef_nn){
			W3(nj, sh, 
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEnn, NMEnn);

			W3(nj, sh,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEnn, NMEnn);
		}
		//cout<<"\tW3:	Hnn -> H'nn is finished"<<endl;
//		Hnn -> H'pn
		if(flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn, 
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEnn, NMEpn);

			W3(nj, sh,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			K, Fcoef, TBMEnn, NMEpn);
		}
		cout<<"\t\tW3: Hnn -> H'nn/H'pn finished\t"; print_time();
//		Hpn -> H'pn
		W3(nj, sh, 
		0, num_pab_pn_even, 0, 
		num_pab_pn, pa_pn, pb_pn,
		num_pab_pn_even, num_pab_pn, num_Vpn_even, 
		num_pab_pn, pa_pn, pb_pn,
		K, Fcoef, TBMEpn, NMEpn);

		W3(nj, sh, 
		num_pab_pn_even, num_pab_pn, num_Vpn_even, 
		num_pab_pn, pa_pn, pb_pn,
		0, num_pab_pn_even, 0, 
		num_pab_pn, pa_pn, pb_pn,
		K, Fcoef, TBMEpn, NMEpn);
		//cout<<"\tW3:	Hpn -> H'pn is finished"<<endl;
//		Hpn -> H'pp
		if(flag_Fcoef_pp||flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			num_pab_pp_even, num_pab_pp, num_Vpp_even, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpn, NMEpp);

			W3(nj, sh, 
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			0, num_pab_pp_even, 0, 
			num_pab_pp, pa_pp, pb_pp,
			K, Fcoef, TBMEpn, NMEpp);
		}
		//cout<<"\tW3:	Hpn -> H'pp is finished"<<endl;
//		Hpn -> H'nn
		if(flag_Fcoef_nn||flag_Fcoef_pn){
			W3(nj, sh, 
			0, num_pab_pn_even, 0, 
			num_pab_pn, pa_pn, pb_pn,
			num_pab_nn_even, num_pab_nn, num_Vnn_even, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEpn, NMEnn);

			W3(nj, sh, 
			num_pab_pn_even, num_pab_pn, num_Vpn_even, 
			num_pab_pn, pa_pn, pb_pn,
			0, num_pab_nn_even, 0, 
			num_pab_nn, pa_nn, pb_nn,
			K, Fcoef, TBMEpn, NMEnn);
		}
		cout<<"\t\tW3: Hpn -> H'pp/H'nn/H'pn finished\t"; print_time();
	}
	cout<<"\tEWSR:\tW3 is collected\t";
	print_time();
	time_t t_W3=time(NULL);
	//for test only
/*
	cout<<"NMEpp:"<<endl;
	for(int i=0;i<num_Vpp;i++){
		cout<<NMEpp[i].numJ<<"\t";
		for(int j=0;j<NMEpp[i].numJ;j++){
			cout<<NMEpp[i].J[j]<<",";
		}
		cout<<"\t";
		for(int j=0;j<NMEpp[i].numJ;j++){
			cout<<NMEpp[i].V[j]<<",";
		}
		cout<<endl;
	}
*/
//----------------------------------------------------------------
/*
	ME2B *NME=new ME2B [num_2body];
	int count=0;
	tempint=0;
        for(int a=0;a<nj;a++)//this block consumes several seconds
        for(int b=0;b<nj;b++)
        for(int J=abs(sh[a]-sh[b]);J<=sh[a]+sh[b];J+=2)
        if(a!=b||J%4==0)
        for(int c=0;c<nj;c++)
        for(int d=0;d<nj;d++)
        if(J>=abs(sh[c]-sh[d])&&J<=sh[c]+sh[d])
        if(c!=d||J%4==0)
	if( (ll[a]+ll[b])%2==(ll[c]+ll[d])%2 )
	if(tz[a]+tz[b]==tz[c]+tz[d]){
		count++;
		NME[tempint].a=a;
		NME[tempint].b=b;
		NME[tempint].c=c;
		NME[tempint].d=d;
		NME[tempint].J=J;
		NME[tempint].V2=TBME[tempint].V2;
		if(!flag_pn_norm){
			if( (a/nj_p==0 && b/nj_p>0)){
				if( a== b-nj_n ) NME[tempint].V2/=sqrt(2);
				if( c== d-nj_n ) NME[tempint].V2/=sqrt(2);
				NME[tempint].V2*=2;//put new "interactions" into unnormalized form again, for possible use in NushellX
			}
		}
		tempint++;
        }
*/
	//------------------Output NMEpn-----
	FILE *fp;
        if((fp=fopen(file_NMEpn,"w"))==NULL)
        cout<<"error: failed to open NMEpn.int"<<endl;
	fprintf(fp,"!Operator for energy-weighted sum rules\n");
	fprintf(fp,"! orbits 1-%d protons, %d-%d neutrons\n",
		nj_p, nj_p+1, 2*nj_p);
	fprintf(fp,"!SP ORB  N  L  J    Tz\n");
	for(int i=0;i<nj_p;i++){
		fprintf(fp,"!SP   %d  %d  %d  %d/2    %d\n",
			i+1,nn[i],ll[i],sh[i],-1);
	}
	for(int i=0;i<nj_p;i++){
		fprintf(fp,"!SP   %d  %d  %d  %d/2    %d\n",
			i+1+nj_p,nn[i],ll[i],sh[i],1);
	}
	fprintf(fp,"!One-body parts\n");
        int num_g=0;
        for(int i=0;i<nj;i++){
                for(int j=i;j<nj;j++){
                        if(fabs(g[i][j])>1E-6){
                                num_g++;
                        }
                }
        }
        fprintf(fp,"%d\n",num_g);
        for(int i=0;i<nj;i++){
                for(int j=i;j<nj;j++){
                        if(fabs(g[i][j])>1E-6){
                                fprintf(fp,"%d\t%d\t%lf\n",i+1,j+1,g[i][j]);
                        }
                }
        }
	fprintf(fp,"!Two-body parts\n");
        fprintf(fp,"%d\n", num_NME);
	int pab_start, pab_end, index_start;
//------	print out NMEpp	----------------------
	for(int odevity=1;odevity>=-1;odevity-=2){
		if(odevity==1){
			pab_start=0;
			pab_end=num_pab_pp_even;
			index_start=0;
		}
		else{
			pab_start=num_pab_pp_even;
			pab_end=num_pab_pp;
			index_start=num_Vpp_even;
		}
		for(int i=pab_start;i<pab_end;i++){
			for(int j=i;j<pab_end;j++){
				int index=index_start
					+ i - pab_start
					+ (j - pab_start + 1)*(j - pab_start)/2;
				if((*NMEpp[index]).numJ>0)
				for(int k=0;k<(*NMEpp[index]).numJ;k++){
					int JJ = (*NMEpp[index]).Jmin + (*NMEpp[index]).dJ * k;
					fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
					pa_pp[i]+1, pb_pp[i]+1, pa_pp[j]+1, pb_pp[j]+1, JJ/2, 1, (*NMEpp[index]).V[k]);
				}
			}
		}
	}
//------	print out NMEnn	----------------------
	for(int odevity=1;odevity>=-1;odevity-=2){
		if(odevity==1){
			pab_start=0;
			pab_end=num_pab_nn_even;
			index_start=0;
		}
		else{
			pab_start=num_pab_nn_even;
			pab_end=num_pab_nn;
			index_start=num_Vnn_even;
		}
		for(int i=pab_start;i<pab_end;i++){
			for(int j=i;j<pab_end;j++){
				int index=index_start
					+ i - pab_start
					+ (j - pab_start + 1)*(j - pab_start)/2;
				if((*NMEnn[index]).numJ>0)
				for(int k=0;k<(*NMEnn[index]).numJ;k++){
					int JJ = (*NMEnn[index]).Jmin + (*NMEnn[index]).dJ * k;
					fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
					pa_nn[i]+1, pb_nn[i]+1, pa_nn[j]+1, pb_nn[j]+1, JJ/2, 1, (*NMEnn[index]).V[k]);
				}
			}
		}
	}
//------        print out NMEpn ----------------------
        for(int odevity=1;odevity>=-1;odevity-=2){
                if(odevity==1){
                        pab_start=0;
                        pab_end=num_pab_pn_even;
                        index_start=0;
                }
                else{
                        pab_start=num_pab_pn_even;
                        pab_end=num_pab_pn;
                        index_start=num_Vpn_even;
                }
                for(int i=pab_start;i<pab_end;i++){
                        for(int j=i;j<pab_end;j++){
                                int index=index_start
                                        + i - pab_start
                                        + (j - pab_start + 1)*(j - pab_start)/2;

                                int case_i_num=1;
                                if(pa_pn[i]!=pb_pn[i]-nj_p)case_i_num=2;
                                int case_j_num=1;
                                if(pa_pn[j]!=pb_pn[j]-nj_p)case_j_num=2;

                                if((*NMEpn[index]).numJ>0)
                                for(int k=0;k<(*NMEpn[index]).numJ;k++){
                                        int JJ = (*NMEpn[index]).Jmin + (*NMEpn[index]).dJ * k;
                                        double V = (*NMEpn[index]).V[k];

                                        if(flag_pn_norm){
                                                fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
                                                pa_pn[i]+1, pb_pn[i]+1, pa_pn[j]+1, pb_pn[j]+1, JJ/2, 1, (*NMEpn[index]).V[k]);
                                        }
                                        else{
                                                if(case_i_num * case_j_num ==1){
                                                        int T=0;
                                                        if(JJ%4==0)T=1;
                                                        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
                                                        pa_pn[i]+1, pb_pn[i]+1, pa_pn[j]+1, pb_pn[j]+1, JJ/2, T, (*NMEpn[index]).V[k]);
                                                }
                                                if(case_i_num * case_j_num ==2){
                                                        int T=0;
                                                        if(JJ%4==0)T=1;
                                                        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
                                                        pa_pn[i]+1, pb_pn[i]+1, pa_pn[j]+1, pb_pn[j]+1, JJ/2, T, sqrt(2)*(*NMEpn[index]).V[k]);
                                                }//sqrt(2)*V(aa;cd), or sqrt(2)*V(ab;cc), because NushellX needs unnormalized matrix elements
                                                else{
                                                        double V_sym;
                                                        for(int l=pab_start;l<pab_end;l++)
                                                        if(pa_pn[l]==pb_pn[i]-nj_p
                                                        &&pb_pn[l]==pa_pn[i]+nj_p){

                                                                int index_sym = index_start;
                                                                if(l<=j)index_sym += l - pab_start
                                                                                + (j - pab_start +1)*(j - pab_start)/2;
                                                                else index_sym += j - pab_start
                                                                                + (l - pab_start +1)*(l - pab_start)/2;
                                                                V_sym = (*NMEpn[index_sym]).V[k];

                                                                for(int T=0;T<2;T++){
                                                                        double V_final = (V + V_sym * theta((sh[pa_pn[i]]+sh[pb_pn[i]]+JJ)/2+T))/2;
                                                                        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
                                                                        pa_pn[i]+1, pb_pn[i]+1, pa_pn[j]+1, pb_pn[j]+1, JJ/2, T, 2*V_final);// 2*V(ab;cd), in NushellX format
                                                                        if(i!=j&& pa_pn[i]==pb_pn[j]-nj_p && pa_pn[j]==pb_pn[i]-nj_p)//NushellX prints out both V(abba;J) and V(baab;J), but NMEpn only stores one of them, so these 3 lines prints out the other.
                                                                        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",
                                                                        pa_pn[j]+1, pb_pn[j]+1, pa_pn[i]+1, pb_pn[i]+1, JJ/2, T, 2*V_final);// 2*V(ab;cd) in NushellX format
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
//---------------------------------------------------
	fclose(fp);
//	delete dynamic memories for TBME and NME
//------	pp
        for(int i=0;i<num_Vpp;i++){
                if( (*TBMEpp[i]).numJ >0 ) delete [] (*TBMEpp[i]).V;
                delete TBMEpp[i];
        }
        delete [] TBMEpp;

        for(int i=0;i<num_Vpp;i++){
                if( (*NMEpp[i]).numJ >0 ) delete [] (*NMEpp[i]).V;
                delete NMEpp[i];
        }
        delete [] NMEpp;
//------	nn
        for(int i=0;i<num_Vnn;i++){
                if( (*TBMEnn[i]).numJ >0 ) delete [] (*TBMEnn[i]).V;
                delete TBMEnn[i];
        }
        delete [] TBMEnn;

        for(int i=0;i<num_Vnn;i++){
                if( (*NMEnn[i]).numJ >0 ) delete [] (*NMEnn[i]).V;
                delete NMEnn[i];
        }
        delete [] NMEnn;
//------	pn
        for(int i=0;i<num_Vpn;i++){
                if( (*TBMEpn[i]).numJ >0 ) delete [] (*TBMEpn[i]).V;
                delete TBMEpn[i];
        }
        delete [] TBMEpn;

        for(int i=0;i<num_Vpn;i++){
                if( (*NMEpn[i]).numJ >0 ) delete [] (*NMEpn[i]).V;
                delete NMEpn[i];
        }
        delete [] NMEpn;
	time_t t_EWSR=time(NULL);
	cout<<"time to calculate gab: "<<t_gab-t_start<<"s"<<endl;
	cout<<"time to input general matrix elements: "<<t_GME-t_gab<<"s"<<endl;
	cout<<"time to make 6j-coefficients tables: "<<t_6j-t_GME<<"s"<<endl;
	cout<<"time to collect W1&W4: "<<t_W1-t_6j<<"s"<<endl;
	cout<<"time to collect W2&W5: "<<t_W2-t_W1<<"s"<<endl;
	cout<<"time to collect W3: "<<t_W3-t_W2<<"s"<<endl;
	cout<<"total time to calculate EWSR: "<<t_EWSR-t_start<<"s"<<endl;
}
