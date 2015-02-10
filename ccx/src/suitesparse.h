#ifdef __cplusplus
extern "C"
#endif

int suitesparsecholmod(double *ad, double *au, double *adb, double *aub, double *sigma,
		       double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);

int suitesparseqr(double *ad, double *au, double *adb, double *aub, double *sigma,
		       double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);
