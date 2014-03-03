void exo(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne0,
	 double *v,double *stn,int *inum,int *nmethod,int *kode,
	 char *filab,double *een,double *t1,double *fn,double *time,
	 double *epn,int *ielmat,char *matname,double *enern,
	 double *xstaten,int *nstate_,int *istep,int *iinc,
	 int *ithermal,double *qfn,int *mode,int *noddiam,
	 double *trab,int *inotr,int *ntrans,double *orab,
	 int *ielorien,int *norien,char *description,int *ipneigh,
	 int *neigh,int *mi,double *stx,double *vr,double *vi,
	 double *stnr,double *stni,double *vmax,double *stnmax,
	 int *ngraph,double *veold,double *ener,int *ne,double *cs,
	 char *set,int *nset,int *istartset,int *iendset,int *ialset,
	 double *eenmax,double *fnr,double *fni,double *emn,
	 double *thicke,char *jobnamec,char *output,double *qfx);

void exosetfind(char *set, int *nset, int *ialset, int *istartset, int *iendset, 
		int *num_ns, int *num_ss, int *num_es, int *num_fs, int *node_map_inv,
		int exoid, 
		int store);

void exoset_warn(int val);

int exoset_check(int n, int *node_map_inv);

void exovector(double *v,int *iset,int *ntrans,char * filabl,int *nkcoords,
               int *inum,int *inotr,double *trab,double *co,
               int *istartset,int *iendset,int *ialset,int *mi,int *ngraph,
               int exoid, int time_step, int countvar, int nout);

int exoset_check(int n, int *node_map_inv);

void exoselect(double *field1,double *field2,int *iset,int *nkcoords,int *inum,
	       int *istartset,int *iendset,int *ialset,int *ngraph,int *ncomp,
	       int *ifield,int *icomp,int *nfield,int *iselect,int exoid,
	       int time_step, int countvar, int nout);
