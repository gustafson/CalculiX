#ifdef __cplusplus
extern "C"
#endif

int suitesparsecholmod(double *ad, double *au, double *adb, double *aub, double *sigma,
		       double *b, int *icol, int *irow, 
		       int *neq, int *nzs);

int suitesparseqr(double *ad, double *au, double *adb, double *aub, double *sigma,
		       double *b, int *icol, int *irow, 
		       int *neq, int *nzs);
