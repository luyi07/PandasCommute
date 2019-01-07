#define temps_scale 3072
/*
   read_nj()--------------read number of orbits

Input: 
char *filename---------filename
int &nj_p--------------number of proton single-j orbits
int &nj----------------total number of proton, neutron single-j orbits

it expects in file "filename":
	comment lines headed with "!"
	a line of anything
	integer integer
	nj
	string nj_p

and nj, nj_p are read and recorded
 */
void read_nj(char* filename, int &nj_p, int &nj){
        char temps[temps_scale];
        char tempc;
        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
        temps[0]='!';
        while(temps[0]=='!')fscanf(fp,"%[^\n]%c",temps,&tempc);
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
int *jj-------------------------2*j
int *parity---------------------parity of orbits
double *spe---------------------single particle energy
int *nn-------------------------n value of orbit
int *ll-------------------------l value of orbit

it expects in the file "filename":
	comment lines headed with "!"
	a line of anything
	integer integer
	nj
	integer	integer	integer
	integer	nn[0]	ll[0]	jj[0]
	...
	integer	nn[nj-1] ll[nj-1] jj[nj-1]
 */
void read_sp(char* filename, int nj, int *nn, int *ll, int *jj){

        char temps[temps_scale];
        char tempc;
        int tempint;
        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
        temps[0]='!';
        while(temps[0]=='!')fscanf(fp,"%[^\n]%c",temps,&tempc);
        int A_core, Z_core;
        fscanf(fp, "%d %d", &A_core, &Z_core);
        fscanf(fp, "%d", &nj);
        for(int i=0;i<3;i++)fscanf(fp,"%d",&tempint);

        for(int i=0;i<nj;i++){
                fscanf(fp,"%d %d %d %d",&tempint,&nn[i],&ll[i],&jj[i]);
                nn[i]--;
        }
        fclose(fp);
}


/*
read_GME_num(): read number of 2BMEs

Input:
char *filename-------------filename
int & num_GME------------number of 2BMEs
*/
void read_GME_num(char *filename, int &num_GME){

        char temps[temps_scale];
        int tempint;
        char tempchar;
        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
//-------- skip comment lines headed with '!' -----------------
        temps[0]='!';
        while(temps[0]=='!'){
                fscanf(fp,"%[^\n]%c",temps,&tempchar);
        }
//-------- find the first integer, as num_GME----------------
        //-------------------
        int length=0;
        for(int i=0;i<temps_scale;i++){
                if(temps[i]>=48&&temps[i]<=57){
                        length++;
                }
                if(length>0&&(temps[i]<48||temps[i]>57)){
                        int y=0;
                        for(int j=1;j<=length;j++){
                                y+=pow(10,j-1)*(temps[i-j]-48);
                        }
                        num_GME=y;
                        break;
                }
        }
        cout<<"num_GME="<<num_GME<<endl;
	fclose(fp);
}

/*
   read_GME--------------read general interactions

Input:
char *filename-------------filename
int nj_p-------------------number of orbits
double *spe----------------single particle energies
int &num_GME---------------number of general interactions
int *GME_a-----------------orbit sequence in orbit shop ordered by a
int *GME_b-----------------orbit sequence in orbit shop ordered by b
int *GME_c-----------------orbit sequence in orbit shop ordered by c
int *GME_d-----------------orbit sequence in orbit shop ordered by d
int *GME_I-----------------2*I of this interaction
int *GME_T-----------------isospin of this interacction
double *GME_V--------------strength of this interaction
*/
void read_GME(char *filename, 
int nj_p, double *spe,
int num_GME, 
int *GME_a, int *GME_b, int *GME_c, int *GME_d,
int *GME_I, int *GME_T, double *GME_V){

        FILE *fp;
        if((fp=fopen(filename,"r"))==NULL)
                cout<<"error: failed to open "<<filename<<endl;
        char temps[temps_scale];
        temps[0]=0;
//-------- skip comment lines headed with '!' -------------------
        temps[0]='!';
        char c;
        while(temps[0]=='!'){
                fscanf(fp,"%[^\n]%c",temps,&c);
        }
//-------- find the 1st integer and assign it to num_GMEpn-------
//-------- find the following floats and assign them to *spe-----
        for(int i=0;i<nj_p;i++)spe[i]=0;
        int realnum=0;
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
                                num_GME=y;
                                realnum++;
                                flaginteger=false;
                                flagfloat=false;
                                y=0;
                                sign=1;
                        }
                        else if(realnum<=nj_p){//s.p.e. real numbers have not ended yet
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
        }
//-------- read in GME ------------
        double tempdouble;
        for(int i=0;i<num_GME;i++){
                fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf", &GME_a[i], &GME_b[i], &GME_c[i], &GME_d[i], &GME_I[i], &GME_T[i], &GME_V[i]);
        }
        fclose(fp);
}
