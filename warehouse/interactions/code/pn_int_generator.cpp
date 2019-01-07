/*
This file contains: 
void pn_int_norm_generator()

still needs to add comment lines and this line:
num_GME, spe[0], spe[1], spe[2], ...

pn_int_norm_generator reads orbits from file_sp, and generate pn-formalism interactions in the convention of th BigStick.
*/
void pn_int_generator(char *file_sp, char *file_pn_int){
//	read in single particle orbits
	int nj_p, nj;
	read_nj(file_sp, nj_p, nj);
	
	int *nn=new int[nj];
	int *ll=new int[nj];
	int *jj=new int[nj];

	read_sp(file_sp, nj, nn, ll, jj);
//	label sequences of orbits, according their j value, from large to small
	int *seq=new int[nj];
	for(int i=0;i<nj;i++){
		seq[i]=i;//0,1,2,...
	}
	cout<<"seq:";
	for(int i=0;i<nj;i++)cout<<seq[i]<<",";
	cout<<endl;
	for(int i=0;i<nj_p;i++){
		for(int j=i;j<nj_p-1;j++){
			if(jj[seq[j]]<jj[seq[j+1]]){//bubble algorithm to sort proton orbits
				int temp=seq[j];
				seq[j]=seq[j+1];
				seq[j+1]=temp;
			}
		}
	}
	for(int i=nj_p;i<nj;i++){
		for(int j=i;j<nj-1;j++){
			if(jj[seq[j]]<jj[seq[j+1]]){//bubble algorithm to sort neutron orbits
				int temp=seq[j];
				seq[j]=seq[j+1];
				seq[j+1]=temp;
			}
		}
	}
	cout<<"seq:";
	for(int i=0;i<nj;i++)cout<<seq[i]<<",";
	cout<<endl;
//***************************************************************
//	calculate	J_pp_min, J_pp_max
//			J_nn_min, J_nn_max
//			J_pn_min, J_pn_max
//***************************************************************
//	determint J_pp_min, J_pp_max for proton-proton interactions
	int J_pp_min=1E6, J_pp_max=0;
	for(int i=0;i<nj_p;i++)
	for(int j=0;j<nj_p;j++){
		if(J_pp_min>abs(jj[i]-jj[j]))J_pp_min=abs(jj[i]-jj[j]);
		if(J_pp_max<jj[i]+jj[j])J_pp_max=jj[i]+jj[j];
	}
	cout<<"J_pp_min= "<<J_pp_min<<" J_pp_max="<<J_pp_max<<endl;
//	determine J_nn_min, J_nn_max for proton-proton interactions
	int J_nn_min=1E6, J_nn_max=0;
	for(int i=nj_p;i<nj;i++)
	for(int j=nj_p;j<nj;j++){
		if(J_nn_min>abs(jj[i]-jj[j]))J_nn_min=abs(jj[i]-jj[j]);
		if(J_nn_max<jj[i]+jj[j])J_nn_max=jj[i]+jj[j];
	}
	cout<<"J_nn_min= "<<J_nn_min<<" J_nn_max="<<J_nn_max<<endl;
//	determine J_pn_min, J_pn_max for proton-neutron interactions
	int J_pn_min=1E6, J_pn_max=0;
	for(int i=0;i<nj_p;i++)
	for(int j=nj_p;j<nj;j++){
		if(J_pn_min>abs(jj[i]-jj[j]))J_pn_min=abs(jj[i]-jj[j]);
		if(J_pn_max<jj[i]+jj[j])J_pn_max=jj[i]+jj[j];
	}
	cout<<"J_pn_min= "<<J_pn_min<<" J_pn_max="<<J_pn_max<<endl;
//***************************************************************
//	calculate	num_pp, num_nn, num_pn
//***************************************************************
//	num_pp
	int T=1;
	int num_pp=0;
	for(int J=J_pp_min;J<=J_pp_max;J+=2){
                for(int a=0;a<nj_p;a++)
                for(int b=a;b<nj_p;b++)//kick out V(bacd;JT)
                for(int c=0;c<nj_p;c++)
                for(int d=c;d<nj_p;d++)//kick out V(abdc;JT)
                if(a<c||(a==c&&b<=d))//kick out V(cdab;JT)
                if(J>=abs(jj[seq[a]]-jj[seq[b]])
                &&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
                &&J>=abs(jj[seq[c]]-jj[seq[d]])
                &&J<=jj[seq[c]]+jj[seq[d]])
                if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2)//parity conservation
                if(a!=b||((jj[seq[a]]+jj[seq[b]]+J)/2+T)%2==0)//in j-j coupling, J=even/odd when T=1/0
                if(c!=d||((jj[seq[c]]+jj[seq[d]]+J)/2+T)%2==0){
			num_pp++;
                }
        }
//	num_nn
	T=1;
	int num_nn=0;
        for(int J=J_nn_min;J<=J_nn_max;J+=2){
                for(int a=nj_p;a<nj;a++)
                for(int b=a;b<nj;b++)//kick out V(bacd;JT)
                for(int c=nj_p;c<nj;c++)
                for(int d=c;d<nj;d++)//kick out V(abdc;JT)
                if(a<c||(a==c&&b<=d))//kick out V(cdab;JT)
                if(J>=abs(jj[seq[a]]-jj[seq[b]])
                &&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
                &&J>=abs(jj[seq[c]]-jj[seq[d]])
                &&J<=jj[seq[c]]+jj[seq[d]])
                if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2)//parity conservation
                if(a!=b||((jj[seq[a]]+jj[seq[b]]+J)/2+T)%2==0)//in j-j coupling, J=even/odd when T=1/0
                if(c!=d||((jj[seq[c]]+jj[seq[d]]+J)/2+T)%2==0){
			num_nn++;
                }
        }
//	num_pn
	int num_pn=0;
        for(int J=J_pn_min;J<=J_pn_max;J+=2){
                for(int a=0;a<nj_p;a++)
                for(int b=nj_p;b<nj;b++)
                for(int c=0;c<nj_p;c++)
                for(int d=nj_p;d<nj;d++)
                if(a<c||(a==c&&b<=d))
                if(J>=abs(jj[seq[a]]-jj[seq[b]])
                &&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
                &&J>=abs(jj[seq[c]]-jj[seq[d]])
                &&J<=jj[seq[c]]+jj[seq[d]])
                if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2){//parity conservation
			num_pn++;
                }
        }
//***************************************************************
//	print out interactions
//***************************************************************
	FILE *fp;
	if((fp=fopen(file_pn_int,"w"))==NULL)
	cout<<"error: failed to open "<<file_pn_int<<endl;

	int num_int =num_pp + num_nn + num_pn;
	fprintf(fp,"!pn-formalism interactions generated by pn_int_norm_generator()\n");
	fprintf(fp,"%d  ",num_int);
	double *spe=new double [nj];
	for(int i=0;i<nj/2;i++)spe[i]=-3*(double)rand()/RAND_MAX;
	for(int i=nj/2;i<nj;i++)spe[i]=spe[i-nj/2];
	for(int i=0;i<nj;i++){
		fprintf(fp,"%lf  ",spe[i]);
	}
	fprintf(fp,"\n");
//	print proton-proton interactions
	for(int J=J_pp_min;J<=J_pp_max;J+=2){
		for(int a=0;a<nj_p;a++)
		for(int b=a;b<nj_p;b++)//kick out V(bacd;JT)
		for(int c=0;c<nj_p;c++)
		for(int d=c;d<nj_p;d++)//kick out V(abdc;JT)
		if(a<c||(a==c&&b<=d))//kick out V(cdab;JT)
		if(J>=abs(jj[seq[a]]-jj[seq[b]])
		&&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
		&&J>=abs(jj[seq[c]]-jj[seq[d]])
		&&J<=jj[seq[c]]+jj[seq[d]])
		if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2)//parity conservation
		if(a!=b||((jj[seq[a]]+jj[seq[b]]+J)/2+T)%2==0)//in j-j coupling, J=even/odd when T=1/0
		if(c!=d||((jj[seq[c]]+jj[seq[d]]+J)/2+T)%2==0){
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", seq[a]+1, seq[b]+1, seq[c]+1, seq[d]+1, J/2, T, -2*(double)rand()/RAND_MAX);
		}
	}
//	print neutron-neutron interactions
	for(int J=J_nn_min;J<=J_nn_max;J+=2){
		for(int a=nj_p;a<nj;a++)
		for(int b=a;b<nj;b++)//kick out V(bacd;JT)
		for(int c=nj_p;c<nj;c++)
		for(int d=c;d<nj;d++)//kick out V(abdc;JT)
		if(a<c||(a==c&&b<=d))//kick out V(cdab;JT)
		if(J>=abs(jj[seq[a]]-jj[seq[b]])
		&&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
		&&J>=abs(jj[seq[c]]-jj[seq[d]])
		&&J<=jj[seq[c]]+jj[seq[d]])
		if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2)//parity conservation
		if(a!=b||((jj[seq[a]]+jj[seq[b]]+J)/2+T)%2==0)//in j-j coupling, J=even/odd when T=1/0
		if(c!=d||((jj[seq[c]]+jj[seq[d]]+J)/2+T)%2==0){
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", seq[a]+1, seq[b]+1, seq[c]+1, seq[d]+1, J/2, T, -2*(double)rand()/RAND_MAX);
		}
	}
//	print proton-neutron interactions
	for(int J=J_pn_min;J<=J_pn_max;J+=2){
		for(int a=0;a<nj_p;a++)
		for(int b=nj_p;b<nj;b++)
		for(int c=0;c<nj_p;c++)
		for(int d=nj_p;d<nj;d++)
		if(a<c||(a==c&&b<=d))//kick out V(cdab;J)
		if(J>=abs(jj[seq[a]]-jj[seq[b]])
		&&J<=jj[seq[a]]+jj[seq[b]]//angular triangles
		&&J>=abs(jj[seq[c]]-jj[seq[d]])
		&&J<=jj[seq[c]]+jj[seq[d]])
		if( (ll[seq[a]]+ll[seq[b]])%2==(ll[seq[c]]+ll[seq[d]])%2){//parity conservation
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", seq[a]+1, seq[b]+1, seq[c]+1, seq[d]+1, J/2, T, -2*(double)rand()/RAND_MAX);
		}
	}
	fclose(fp);
}
