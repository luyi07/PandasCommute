
void read_nj(char * filename, int & nj_p, int & nj);
// read number of single j orbits from filename
// nj_p: number of proton single-j orbits
// nj: number of single-j orbits, including neutron ones

void read_shell(char * filename, int nj, int * nn, int * ll, int * sh);
// read single-j orbits from filename, including the quantum numbers n, l, j of the orbits

void read_GMEpn_num(char * filename, int & num_GMEpn);
// read number of two-body interactions in filename

void read_GMEpn(char *filename, int nj, int num_GMEpn, double *spe,
int *GME_a, int *GME_b, int *GME_c, int *GME_d,
int *GME_I, int *GME_T, double *GME_V); 
// read general interactions

void read_F_new(char *file_F, int *sh, int &K, int &pF, double **Fcoef);
// read info of \hat{F}, the one-body transition operator
