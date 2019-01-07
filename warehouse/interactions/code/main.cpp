#include<iostream>
using namespace std;
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<time.h>

#include"input.cpp"
#include"iso_int_generator.cpp"
#include"pn_int_generator.cpp"

int main(){

	srand((unsigned)time(0));
	char file_sp[80]="../input/pn.sp";
	char file_iso_int[80]="../output/iso.int";
	char file_pn_int[80]="../output/pn.int";
//	iso_int_generator(file_sp, file_iso_int);
	pn_int_generator(file_sp, file_pn_int);
	return 0;
}
