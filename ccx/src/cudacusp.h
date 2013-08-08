#ifdef __cplusplus
extern "C"
#endif

int cudacusp(double *ad, double *au, double *adb, double *aub, double *sigma, 
	     double *b, int *icol, int *irow, int *neq, int *nzs, 
	     int *symmetryflag, int *inputformat, int *jq, int *nzs3);
