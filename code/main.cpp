#include<iostream>
#include<cstdlib>
#include<fstream>
using namespace std;

#include"SR.h"
#include"Fcoef_producer.h"

bool flag_pn_norm=true;
void print_time();
int K;

int main(){

	cout<<"\n\t\t/-----------------------------------\\ \n";
	cout<<"\t        |  PandasCommute 3.1       by Yi Lu |\n";
	cout<<"\t\t\\-----------------------------------/ \n";
	cout<<"\nGreetings. PandasCommute is at your service.\n";
	cout<<"\nI'm a small code to calculate scalar operators, which helps evaluate non-energy-weighted and energy-weighted sum rules of E&M, beta transition strengths in the shell model.\n";
	cout<<"\nI need three files as input, one for single particle orbits i.e. to define the valence space, one for general interactions, and one file for the one-body transition operator.\n";

	char file_NEWSR[80]="output/O_NEWSR.int";
	char file_EWSR[80]="output/O_EWSR.int";

        string sp;
        cout<<"\nPlease input the filename for s.p. info, as examplified in example/sd-shell/input/pn.sp \n\t";
        cin>> sp;
	char file_sp[80];
	for(int i=0;i<80;i++)file_sp[i]=sp[i];

        string GMEpn;
        cout<<"\nPlease input the filename for general interactions, as examplified in example/sd-shell/input/GMEpn.int \n\t";
        cin>> GMEpn;
	char file_GMEpn[80];
	for(int i=0;i<80;i++)file_GMEpn[i]=GMEpn[i];

	char switch_pn_norm='h';
	while(switch_pn_norm!='y'
	&&switch_pn_norm!='n'){
		cout<<"Are the pn matrix elements normalized or not? (y/n)\n";
		cout<<"\ty if you're using interaction files adapted to BigStick\n";
		cout<<"\tn if you're using interaction files adapted to NushellX\n\t";
		cin>>switch_pn_norm;
	}

	char switch_F_book='h';
	while(switch_F_book!='y'
	&&switch_F_book!='n'){
		cout<<"\nI can generate a one-body transition operator automatically, would you like to book one? y/n\nInput 'n' if you have already prepared one, and I'll ask for its name later.\n\t";
		cin>>switch_F_book;
		if(switch_F_book != 'y' && switch_F_book != 'n')
			cout<<"Sorry I don't understand, allow me to ask again\n";
	}
	
	char file_F[80];
	if(switch_F_book == 'n'){
		string Fname;
		cout<<"\nPlease input the filename for one-body transition operator, as examplified in example/sd-shell/input/F.coef \n\t";
		cin>> Fname;
		for(int i=0;i<80;i++)file_F[i]=Fname[i];
	}
	else{
		cout<<"\nPlease input the filename for one-body transition operator, I'll write info into it.\n\t";
		string Fname;
		cin>> Fname;
		for(int i=0;i<80;i++)file_F[i]=Fname[i];
		cal_Fcoef(file_sp, file_F);
	}

	if(switch_pn_norm=='y')flag_pn_norm=true;
	else flag_pn_norm=false;
	cout<<"Start to calculate the NEWSR operator...\t";
	print_time();
	cal_NEWSR_operator(file_sp, file_F, file_GMEpn, file_NEWSR);
	cout<<"Finished NEWSR, start to calculate EWSR ...\t";
	print_time();
	cal_EWSR_operator(file_sp, file_GMEpn, file_F, file_EWSR);
	cout<<"Finished EWSR. \t";
	print_time();
	cout<<"-----------------------------------------------------------------\n";

	return 0;
}

/*
print_time()----------------print out the time at a moment
*/
void print_time(void)
{
    time_t ttt = time( 0 );
    char tmp[64];
    strftime( tmp, sizeof(tmp), "%Y/%m/%d %X",localtime(&ttt) );
    puts( tmp ); //output the year-date-time
}
