
#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"global.h"
#include"struct.h"
#include"3j6j9j.h"
#include"input.h"

/*
This file contains functions to calculate non-energy-weighted sum rules of one-body operator transitions in the shell model spaces.
*/

/*
Goal: calculate matrix elements of O_NEWSR

input: 
char *file_sp	filename of single particle info
char *file_F	filename of Fcoef, i.e. the one-body operator
char *file_GMEpn	filename of input matrix elements
char *file_NEWSR	filename of output: matrix elements of \hat{O}_{NEWSR}
*/
void cal_NEWSR_operator(char *file_sp, char *file_F, char *file_GMEpn, char *file_NEWSR){
//--------------read shell info--------------
        int nj_p, nj;	// nj_p: number of single-j orbits for protons
			// nj: number of single-j orbits
        read_nj(file_sp, nj_p, nj);
        int nj_n=nj-nj_p;// nj_n: number of single-j orbits for neutrons
        int *nn=new int [nj];// n value, 0 for 0d5/2 in sd shell
        int *ll=new int [nj];// l value, 2 for 0d5/2 in sd shell
        int *sh=new int [nj];// doubled j-value of orbits
        read_shell(file_sp, nj, nn, ll, sh);
        int *parity=new int [nj];// parity of orbits, 1 for 0d5/2, -1 for 0p3/2
        for(int i=0;i<nj;i++){
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
	//	for test only
/*
	cout<<"Fcoef:"<<endl;
	for(int i=0;i<nj;i++){
		for(int j=0;j<nj;j++){
			cout<<Fcoef[i][j]<<",";
		}
		cout<<endl;
	}
*/
//--------------calculate gab
	double **g=new double *[nj];
	for(int i=0;i<nj;i++){
		g[i]=new double [nj];
		for(int j=0;j<nj;j++){
			g[i][j]=0;
		}
	}

	for(int a=0;a<nj;a++){
		for(int b=0;b<nj;b++)
		if(sh[a]==sh[b]&&parity[a]==parity[b]){
			double y=0;
			for(int c=0;c<nj;c++){
				y+=Fcoef[c][a]*Fcoef[c][b];
			}
			y/=sh[a]+1;
			g[a][b]=y;
		}
	}
	cout<<"\tNEWSR:\tg calculated\t";
	print_time();
//------read general interactions, cause I need the sequences of the 2body interactions therein. Pandawarrior relies on that-------
        int num_GMEpn;            //number of matrix elements
        read_GMEpn_num(file_GMEpn, num_GMEpn);
//	cout<<"num_GMEpn="<<num_GMEpn<<endl;

        double *spe=new double [nj];
        int *GME_a=new int [num_GMEpn];                   //sequence of ja
        int *GME_b=new int [num_GMEpn];                   //sequence of jb
        int *GME_c=new int [num_GMEpn];                   //sequence of jc
        int *GME_d=new int [num_GMEpn];                   //sequence of jd
        int *GME_I=new int [num_GMEpn];                   //2*I
        int *GME_T=new int [num_GMEpn];                   //2*T
        double *GME_V=new double [num_GMEpn];             //V(abcd;JT)
//	cout<<"dynamic memories opened"<<endl;
        read_GMEpn(file_GMEpn, nj, num_GMEpn, spe, GME_a, GME_b, GME_c, GME_d, GME_I, GME_T, GME_V);
//************** two-body parts of O_NEWSR ****************
	int num_NME=0;
//*************************************************************
//      proton-proton matrix elements
//*************************************************************
//------        proton-proton: pab      ----------------------
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
        //cout<<"\tnum_pab_pp_even="<<num_pab_pp_even<<endl;
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
//------open memory and assign values for NMEpp (Hpp)-----
        int num_Vpp_even = (num_pab_pp_even +1)*num_pab_pp_even/2;
        int num_Vpp =   (num_pab_pp_even +1)*num_pab_pp_even/2
                        + (num_pab_pp_odd + 1)*num_pab_pp_odd/2;

        VJ **NMEpp=new VJ *[num_Vpp];// even-parity & odd-parity interactions
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

                                NMEpp[index]=new VJ;//allocate dynamic memory: a new VJ

                                if(JJmax>=JJmin){

                                        int numJJ=0;
                                        if((a==b||c==d)){
                                                if(JJmin%4!=0)JJmin+=2;//Pauli's principle
                                                if(JJmax%4!=0)JJmax-=2;
                                                numJJ=(JJmax-JJmin)/4+1;
                                        }
                                        else{
                                                numJJ=(JJmax-JJmin)/2+1;
                                        }

                                        (*NMEpp[index]).numJ=numJJ;
                                        (*NMEpp[index]).Jmin=JJmin;
                                        if(numJJ==1)(*NMEpp[index]).dJ=0;
                                        else{
                                                (*NMEpp[index]).dJ=(JJmax-JJmin)/(numJJ-1);
                                        }
                                        (*NMEpp[index]).V=new float [numJJ];
                                        for(int k=0;k<numJJ;k++){
						int J=(*NMEpp[index]).Jmin + k * (*NMEpp[index]).dJ;

						double y1=-2*Fcoef[c][b]*Fcoef[a][d]*cal6j(sh[a], sh[d], K, sh[c], sh[b], J);
						if(a==b)y1/=sqrt(2);
						if(c==d)y1/=sqrt(2);

						if(a==b)y1=y1*2;
						else{
							double y2=-2*theta(1+(sh[a]+sh[b]+J)/2)*Fcoef[c][a]*Fcoef[b][d]*cal6j(sh[b], sh[d], K, sh[c], sh[a], J);
							if(a==b)y2=y2/sqrt(2);
							if(c==d)y2=y2/sqrt(2);
							
							y1=y1+y2;
						}

                                                (*NMEpp[index]).V[k]=y1;
                                        }
					num_NME+=numJJ;
                                }
                                else{
                                        (*NMEpp[index]).numJ=0;
                                }
                        }
                }
        }
        cout<<"\tNEWSR:\tNMEpp is finished\t";
	print_time();
//*************************************************************
//****** neutron-neutron matrix elements *******
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
//------open memory and assign values for NMEnn (Hnn)-----
        int num_Vnn_even = (num_pab_nn_even +1)*num_pab_nn_even/2;
        int num_Vnn =   (num_pab_nn_even +1)*num_pab_nn_even/2
                        + (num_pab_nn_odd + 1)*num_pab_nn_odd/2;

        VJ **NMEnn=new VJ * [num_Vnn];// even-parity & odd-parity interactions
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

                                NMEnn[index]=new VJ;
                                if(JJmax>=JJmin){

                                        int numJJ=0;
                                        if((a==b||c==d)){
                                                if(JJmin%4!=0)JJmin+=2;//Pauli's principle
                                                if(JJmax%4!=0)JJmax-=2;
                                                numJJ=(JJmax-JJmin)/4+1;
                                        }
                                        else{
                                                numJJ=(JJmax-JJmin)/2+1;
                                        }

                                        (*NMEnn[index]).numJ=numJJ;
                                        (*NMEnn[index]).Jmin=JJmin;
                                        if(numJJ==1)(*NMEnn[index]).dJ=0;
                                        else{
                                                (*NMEnn[index]).dJ=(JJmax-JJmin)/(numJJ-1);
                                        }
                                        (*NMEnn[index]).V=new float [numJJ];
                                        for(int k=0;k<numJJ;k++){


                                                int J=(*NMEnn[index]).Jmin + k * (*NMEnn[index]).dJ;

                                                double y1=-2*Fcoef[c][b]*Fcoef[a][d]*cal6j(sh[a], sh[d], K, sh[c], sh[b], J);
                                                if(a==b)y1/=sqrt(2);
                                                if(c==d)y1/=sqrt(2);

                                                if(a==b)y1=y1*2;
                                                else{
                                                        double y2=-2*theta(1+(sh[a]+sh[b]+J)/2)*Fcoef[c][a]*Fcoef[b][d]*cal6j(sh[b], sh[d], K, sh[c], sh[a], J);
                                                        if(a==b)y2=y2/sqrt(2);
                                                        if(c==d)y2=y2/sqrt(2);

                                                        y1=y1+y2;
                                                }

                                                (*NMEnn[index]).V[k]=y1;
                                        }
					num_NME+=numJJ;
                                }
                                else{
                                        (*NMEnn[index]).numJ=0;
                                }
                        }
                }
        }
        cout<<"\tNEWSR:\tNMEnn is finished\t";
	print_time();
//*************************************************************
//****** proton-neutron matrix elements *******
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
//------open memory and assign values for NMEpn (Hpn)-----
        int num_Vpn_even = (num_pab_pn_even +1)*num_pab_pn_even/2;
        int num_Vpn =   (num_pab_pn_even +1)*num_pab_pn_even/2
                        + (num_pab_pn_odd + 1)*num_pab_pn_odd/2;

        VJ **NMEpn=new VJ * [num_Vpn];// even-parity & odd-parity interactions
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

                                NMEpn[index]=new VJ;
                                if(JJmax>=JJmin){

                                        int numJJ=0;
                                        if((a==b||c==d)){
                                                if(JJmin%4!=0)JJmin+=2;//Pauli's principle
                                                if(JJmax%4!=0)JJmax-=2;
                                                numJJ=(JJmax-JJmin)/4+1;
                                        }
                                        else{
                                                numJJ=(JJmax-JJmin)/2+1;
                                        }

                                        (*NMEpn[index]).numJ=numJJ;
                                        (*NMEpn[index]).Jmin=JJmin;
                                        if(numJJ==1)(*NMEpn[index]).dJ=0;
                                        else{
                                                (*NMEpn[index]).dJ=(JJmax-JJmin)/(numJJ-1);
                                        }
                                        (*NMEpn[index]).V=new float [numJJ];
                                        for(int k=0;k<numJJ;k++){

                                                int J=(*NMEpn[index]).Jmin + k * (*NMEpn[index]).dJ;

                                                double y1=-2*Fcoef[c][b]*Fcoef[a][d]*cal6j(sh[a], sh[d], K, sh[c], sh[b], J);
                                                if(a==b)y1/=sqrt(2);
                                                if(c==d)y1/=sqrt(2);

                                                if(a==b)y1=y1*2;
                                                else{
                                                        double y2=-2*theta(1+(sh[a]+sh[b]+J)/2)*Fcoef[c][a]*Fcoef[b][d]*cal6j(sh[b], sh[d], K, sh[c], sh[a], J);
                                                        if(a==b)y2=y2/sqrt(2);
                                                        if(c==d)y2=y2/sqrt(2);

                                                        y1=y1+y2;
                                                }

                                                (*NMEpn[index]).V[k]=y1;
//                                               (*NMEpn[index]).V[k]=0;
						if(flag_pn_norm)num_NME++;
						else{
							int multiple=1;
							if(pa_pn[i]!=pb_pn[i]-nj_p && pa_pn[j]!=pb_pn[j]-nj_p){
								multiple*=2;
								if(i!=j && ( pa_pn[i]==pb_pn[j]-nj_p && pa_pn[j]==pb_pn[i]-nj_p) )
								multiple+=2;
							}
							num_NME+=multiple;
						}
                                        }
                                }
                                else{
                                        (*NMEpn[index]).numJ=0;
                                }
                        }
                }
        }

        cout<<"\tNEWSR:\tNMEpn is finished\t";
	print_time();

//***************************************************
//***************************************************
//****** Print gab and W(abcd;I) into file_NEWSR
//------ print gab into file_NEWSR ----------
        FILE *fp;
        if((fp=fopen(file_NEWSR,"w"))==NULL)
        cout<<"error: failed to open ../output/g.int"<<endl;
        fprintf(fp, "!Operator for non-energy weighted sum rules\n");
	fprintf(fp, "!SP ORB  N  L  J\n");
	for(int i=0;i<nj;i++){
		fprintf(fp, "!SP   %d  %d  %d  %d/2\n",i+1,nn[i],ll[i],sh[i]);
	}
	fprintf(fp, "!one-body parts of the NEWSR operator\n");
	int num_g=0;
	for(int a=0;a<nj;a++){
		for(int b=a;b<nj;b++){
			if(fabs(g[a][b])>1E-6)num_g++;
		}
	}
        fprintf(fp,"%d\n",num_g);
	for(int a=0;a<nj;a++){
		for(int b=a;b<nj;b++){
			if(fabs(g[a][b])>1E-6)fprintf(fp,"%d\t%d\t%lf\n",a+1,b+1,g[a][b]);
		}
	}
//	print W(abcd;I) into file_NEWSR
	fprintf(fp, "!two-body parts of the NEWSR operator\n");
        fprintf(fp,"%d\n", num_NME);
        int pab_start, pab_end, index_start;
//------        print out NMEpp ----------------------
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
//------        print out NMEnn ----------------------
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
                                                        &&pb_pn[l]==pa_pn[i]+nj_p){//given V(abcd;J), find V(bacd;J)

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
//******delete dynamic memories, before they turn into zombies
//------delete TBME 
//      pp
        for(int i=0;i<num_Vpp;i++){
		if( (*NMEpp[i]).numJ >0 ) delete [] (*NMEpp[i]).V;
		delete NMEpp[i];
	}
        delete [] NMEpp;
//      nn
        for(int i=0;i<num_Vnn;i++){
		if( (*NMEnn[i]).numJ >0 ) delete [] (*NMEnn[i]).V;
		delete NMEnn[i];
	}
        delete [] NMEnn;
//      pn
        for(int i=0;i<num_Vpn;i++){
		if( (*NMEpn[i]).numJ >0 ) delete [] (*NMEpn[i]).V;
		delete NMEpn[i];
	}
        delete [] NMEpn;
//------delete GMEpn
	delete [] spe;
	delete [] GME_a;
	delete [] GME_b;
	delete [] GME_c;
	delete [] GME_d;
	delete [] GME_I;
	delete [] GME_T;
	delete [] GME_V;
//------delete gab
	for(int i=0;i<nj;i++)delete [] g[i];
	delete [] g;
//------delete Fcoef
	for(int i=0;i<nj;i++)delete [] Fcoef[i];
	delete [] Fcoef;
//------delete sp info
	delete [] nn;
	delete [] ll;
	delete [] sh;
	delete [] parity;
}

