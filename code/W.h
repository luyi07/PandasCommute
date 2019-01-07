
void W1(int nj_input, int *sh_input,
        int pab_start_input, int pab_end_input,
        int num_pab_input, int *pa_input, int *pb_input, int index_start_input,
        int K_input, double **Fcoef_input,
        VJ **TBME, VJ **NME);
// W1(abcd;J) in the paper

void W2(int nj_input, int nj_start_input, int nj_end_input, int *sh_input,
        int num_pab_input, int pab_start_input, int pab_end_input,
        int *pa_input, int *pb_input,
        int index_start_input,
        int K_input, double **Fcoef_input,
        VJ **TBME, VJ **NME);
// W2(abcd;J) in the paper

void W3(int nj_input, int *sh_input,
        int pbe_start_input, int pbe_end_input, int index_start_bedf_input,
        int num_pab1_input, int *pa1_input, int *pb1_input,
        int pab_start_input, int pab_end_input, int index_start_abcd_input,
        int num_pab2_input, int *pa2_input, int *pb2_input,
        int K_input, double **Fcoef_input, VJ **TBME, VJ **NME);
// W3(abcd;J) in the paper
