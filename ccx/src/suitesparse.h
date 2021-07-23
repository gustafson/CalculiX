#ifdef __cplusplus
extern "C"
#endif

#ifdef LONGLONG

#ifndef SuiteSparse_long

#ifdef _WIN64
#define SuiteSparse_long __int64
#define SuiteSparse_long_max _I64_MAX
#define SuiteSparse_long_idd "I64d"

#else

#define SuiteSparse_long long long
#define SuiteSparse_long_max LONG_LONG_MAX
#define SuiteSparse_long_idd "lld"

#endif
#endif
#endif


int suitesparsecholmod(double *ad, double *au, double *adb, double *aub, double *sigma,
		       double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);

// int suitesparseqr(double *ad, double *au, double *adb, double *aub, double *sigma,
// 		       double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);
