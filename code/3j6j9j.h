double theta(int x);// returns 1 if x is even, -1 if odd.

double factorial(int n);
// literally factorial n!

double delta(int a, int b, int c);
// a quantity to be used in j3symbol(...)

double j3symbol(int j1, int j2, int j3, int m1, int m2, int m3);
// calculate 3j-coefficient, j1,...,m3 are all doubled

double calCG(int j1, int m1, int j2, int m2, int j3, int m3);
// CG coefficient, j1, ..., m3 are all doubled

double cal6j(int j1, int j2, int j3, int l1, int l2, int l3);
// calculate a 6j-coefficient, j1,...,l3 are all doubled

void j6table(int K);
// make a table of 6j-coefficients in the memory
