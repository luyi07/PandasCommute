double beta_coef(int a, int b, int sign, int L);
// beta transition coefficients, F[a][b], sign can be 1 or -1, for beta plus or minus, L can be 0 or 2, for Fermi or Gamow-Teller (angular momentum is doubled)

void cal_Fcoef(char *file_sp, char *file_F);
// calculate coefficients of transition operators
//		\hat{F} = \sum_{ab} [K]^{-1} Fcoef[a][b] (a^\dagger \otimes \tilde{b})_K
