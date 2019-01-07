
double r2(int n,int l1,int l2);// < || r2 || >
double gamma(int x);
double rr(int n1,int l1, int n2, int l2, int lambda);
double qab(int n1,int l1,int n2,int l2,int j1,int j2,int lambda);
double qM_L(int n1, int l1, int j1, int n2, int l2, int j2, int lambda);
// orbital M1 operator structure coefficients

/*
   M1 operator: \vec{S}--structure coefficients
   j1, j2, are doubled
   n1, l1, n2, l2, lambda are not
 */
double qM_S(int n1, int l1, int j1, int n2, int l2, int j2, int lambda);
