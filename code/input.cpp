
#include<iostream>
using namespace std;

#include<fstream>
#include<cmath>

#include"3j6j9j.h"

#define temps_scale 3072
/*
   read_nj()--------------read number of orbits

Input: 
char *filename---------filename
int &nj_p--------------number of proton single-j orbits
int &nj----------------total number of proton, neutron single-j orbits

thus nj is read in
*/
void read_nj(char* filename, int &nj_p, int &nj){
	char temps[temps_scale];
	char tempc;
	FILE *fp;
	if((fp=fopen(filename,"r"))==NULL)
		cout<<"error: failed to open "<<filename<<endl;
	temps[0]='!';
	while(temps[0]=='!'||temps[0]=='#'){
		for(int i=0;i<temps_scale;i++)temps[i]=0;
		fscanf(fp,"%[^\n]%c",temps,&tempc);
	}
	int A_core, Z_core;
	fscanf(fp, "%d %d", &A_core, &Z_core);
	fscanf(fp, "%d", &nj);
	fscanf(fp, "%s %d", temps, &nj_p);
	fclose(fp);
}

/*
   read_sp()--------------------read single particle orbits

Input:
char *filename------------------filename
int &nj-------------------------number of orbits
int *sh-------------------------2*j
int *parity---------------------parity of orbits
double *spe---------------------single particle energy
int *nn-------------------------n value of orbit
int *ll-------------------------l value of orbit
 */
void read_shell(char* filename, int nj, int *nn, int *ll, int *sh){

	char temps[temps_scale];
	char tempc;
	int tempint;
	FILE *fp;
	if((fp=fopen(filename,"r"))==NULL)
		cout<<"error: failed to open "<<filename<<endl;
	temps[0]='!';
	while(temps[0]=='!'||temps[0]=='#'){
		for(int i=0;i<temps_scale;i++)temps[i]=0;
		fscanf(fp,"%[^\n]%c",temps,&tempc);
	}
	int A_core, Z_core;
	fscanf(fp, "%d %d", &A_core, &Z_core);
	fscanf(fp, "%d", &nj);
	for(int i=0;i<3;i++)fscanf(fp,"%d",&tempint);
	
	for(int i=0;i<nj;i++){
		fscanf(fp,"%d %d %d %d",&tempint,&nn[i],&ll[i],&sh[i]);
	}
	bool flag_n=false;
	for(int i=0;i<nj;i++){
		if(nn[i]==0)flag_n=true;
	}
	if(!flag_n){
		for(int i=0;i<nj;i++){
			nn[i]--;
		}
	}
	fclose(fp);
}

/*
read_GMEpn_num()-----------read number of general two-body interactions in pn formalism

Input:
char *filename-------------filename
int & num_GMEpn------------number of general two-body interactions
*/
void read_GMEpn_num(char *filename, int &num_GMEpn){

        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
        char temps[temps_scale];//temps_scale is the upper limit
//-------- skip comment lines headed with '!' --------------------
	temps[0]='!';
	char c;
        while(temps[0]=='!'||temps[0]=='#'){
		for(int i=0;i<temps_scale;i++)temps[i]=0;
		fscanf(fp,"%[^\n]%c",temps,&c);
/*
		bool flag_blankline=true;
		for(int i=0;i<temps_scale;i++)
		if(temps[i]!=0)flag_blankline=false;
		if(flag_blankline)temps[0]='!';
*/
        }
//-------- find the first integer, as num_GMEPn-------------------
	int length=0;
	for(int i=0;i<temps_scale;i++){
		if(temps[i]>=48&&temps[i]<=57){
			length++;//meet a figure, count its digits
		}
		if(length>0&&(temps[i]<48||temps[i]>57)){
			int y=0;
			for(int j=1;j<=length;j++){
				y+=pow(10,j-1)*(temps[i-j]-48);
			}
			num_GMEpn=y;//picks up the 1st positive integer in the 1st line after comments
			break;
		}
	}
        if(length==0){
		cout<<"Error:\tfailed to find number of two-body matrix elements in GMEpn.int\nAn integer is expected in the first line unheaded with ! or # \n";
	}
        fclose(fp);
}

/*
   read_GMEpn--------------read general interactions

Input:
char *filename-------------filename

int &num_GMEpn-------------number of general interactions in the pn form
double *spe----------------single particle energies
int *GME_a-----------------orbit sequence in orbit shop ordered by a
int *GME_b-----------------orbit sequence in orbit shop ordered by b
int *GME_c-----------------orbit sequence in orbit shop ordered by c
int *GME_d-----------------orbit sequence in orbit shop ordered by d
int *GME_I-----------------2*I of this interaction
int *GME_T-----------------isospin of this interacction
double *GME_V--------------strength of this interaction

It skips the first lines headed with "!" or "#", read the first integer as number of two-body matrix elements (2BME), and floats as single particle energies, then start reading in 2BMEs. So no matter the pn matrix elements are normalized or not, this function is applicable as well.
*/
void read_GMEpn(char *filename, int nj, int num_GMEpn, double *spe,
int *GME_a, int *GME_b, int *GME_c, int *GME_d,
int *GME_I, int *GME_T, double *GME_V){

        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
        char temps[temps_scale];//temps_scale is the upper limit
//-------- skip comment lines headed with '!' --------------------
	temps[0]='!';
	char c;
        while(temps[0]=='!'||temps[0]=='#'){
		for(int i=0;i<temps_scale;i++)temps[i]=0;
		fscanf(fp,"%[^\n]%c",temps,&c);
        }
//-------- find the 1st integer and assign it to num_GMEpn--------
//-------- find the following floats and assign them to *spe-----
	for(int i=0;i<nj;i++)spe[i]=0;
	int realnum=0;
	while(realnum<nj+1){
		bool flaginteger=false;
		bool flagfloat=false;
		double y=0;
		int floatdigit=0;
		char sign=1;
		for(int i=0;i<temps_scale;i++){
			if(temps[i]=='-'){
				flaginteger=true;
				sign=-1;
			}
			else if(temps[i]>=48&&temps[i]<=57){
				if(flagfloat==false){
					if(!flaginteger)flaginteger=true;//flip flaginteger to true
					y=y*10+temps[i]-48;//picks up this digit
				}
				else{//OK, this digit is after a pointer, its the float part of a number
					y=y+(temps[i]-48)*pow(0.1,floatdigit+1);
					floatdigit++;
				}
			}
			else if(temps[i]=='.'){
				flaginteger=false;//integer part ends
				flagfloat=true;//float part begins
			}
			else if(flaginteger||flagfloat){//ends picking up digits
				if(realnum==0){//the 1st real number is assianged to num_GMEpn
					num_GMEpn=y;
					realnum++;
					flaginteger=false;
					flagfloat=false;
					y=0;
					sign=1;
				}
				else if(realnum<=nj){//s.p.e. real numbers have not ended yet
					spe[realnum-1]=sign*y;
					realnum++;
					flaginteger=false;
					flagfloat=false;
					floatdigit=0;
					y=0;
					sign=1;
				}
				else{
					break;
				}
			}
                        else{
                                flaginteger=false;
                                flagfloat=false;
                        }
		}
		if(realnum<nj+1){//to tolerate a line break among single-particle-energies
			for(int i=0;i<temps_scale;i++)temps[i]=0;
			fscanf(fp,"%[^\n]%c",temps,&c);
		}
	}
//-------- read in GME ------------
	double tempdouble;
        for(int i=0;i<num_GMEpn;i++){
                fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf", &GME_a[i], &GME_b[i], &GME_c[i], &GME_d[i], &GME_I[i], &GME_T[i], &GME_V[i]);
        }
        fclose(fp);
}

/*
read_F_new

What'sNew: this version does not preassume proton or neutron single particle states, but treat them as same particles on different single particle states.

Goal: read in K, the angular momentum of one-body operator F
pF, the parity of this operator,
and the structure coefficients of F.

Input:
char * file_F: file name
int nj_p: number of proton j-orbits
int *sh
int & K
int & pF
double **F[nj][nj]
*/
void read_F_new(char *file_F, int *sh, int &K, int &pF, double **Fcoef){

        FILE *fp;
        char temps[temps_scale];
	char tempchar;
	int tempint;
        if((fp=fopen(file_F,"r"))==NULL)
        cout<<"error: failed to open "<<file_F<<endl;
	
	temps[0]='!';
	while(temps[0]=='!'||temps[0]=='#'){
		for(int i=0;i<temps_scale;i++)temps[i]=0;
		fscanf(fp,"%[^\n]%c",temps,&tempchar);
	}
//----- find the 1st integer, as K and the 2nd as pF------------
        int num_integer=0;
        while(num_integer<2){
                bool flaginteger=false;
                int y=0;
                int sign=1;
                for(int i=0;i<200;i++){
                        if(temps[i]=='-'){
                                flaginteger=true;
                                sign=-1;
                        }
                        else if(temps[i]>=48&&temps[i]<=57){
                                if(!flaginteger)flaginteger=true;
                                y=y*10+temps[i]-48;
                        }
                        else if(flaginteger){//ends an integer
                                if(fabs(y)>0){
                                        if(num_integer==0){
                                                K=2*y;//K must be positive, so K does not pick up sign. K is doubled.
                                                num_integer++;
                                        }
                                        else{
                                                pF=y*sign;
                                                num_integer++;
                                                break;
                                        }
                                        flaginteger=false;
                                        y=0;
                                        sign=1;
                                }
                                else{
                                        flaginteger=false;
                                        sign=1;
                                }
                        }
                }
        }//while
cout<<"K="<<K<<"\tpF="<<pF<<endl;
	int a,b; 
	double y;
	while( fscanf(fp,"%d\t%d\t%lf",&a,&b,&y)!=EOF ){
			Fcoef[a-1][b-1]=y;
			if(a!=b)Fcoef[b-1][a-1]=theta(1+(sh[a-1]+sh[b-1])/2)*y;
        }
        fclose(fp);
}
