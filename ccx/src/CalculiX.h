/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#define Linux 1
#define IRIX 2
#define IRIX64 3
#define HP 4

#if ARCH == Linux
#define FORTRAN(A,B) A##_  B
#elif ARCH == IRIX || ARCH == IRIX64
#define FORTRAN(A,B) A##_##B
#elif ARCH == HP
#define FORTRAN(A,B) A##B
#endif

#define NNEW(a,b) (a *)u_calloc((b),sizeof(a))
#define RENEW(a,b,c) a=(b *) realloc((b *)(a),(c)*sizeof(b))

void FORTRAN(allocation,(int *nload_,int *nforc_,int *nboun_,
             int *nk_,int *ne_,int *nmpc_,int *nset_,int *nalset_,
	     int *nmat_,int *ntmat_,int *npmat_,int *norien_,int *nam_,
             int *nprint_,int *mint_,int *ntrans_,
             char *set,int *meminset,
             int *rmeminset, int *ncs_, int *namtot_, int *ncmat_,
             int *memmpc_, int *ne1d, int *ne2d, int *nflow,
             char *jobnamec, int *irstrt, int *ithermal, int *nener,
             int *nstate_, int *istep, char *inpc,
             int *ipoinp, int *inp, int *ntie_, int *nbody_,
             int *nprop_, int *ipoinpc));

void FORTRAN(allocont,(int *ncont,int *ntie,char *tieset,int *nset,
             char *set,int *istartset,int *iendset,int *ialset,
	     char *lakon, int *ncone, double *tietol, int *ismallsliding));

void FORTRAN(applyboun,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon));

void arpack(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int *icol, int *jq, int *irow, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun, 
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *shcon, int *nshcon, double *cocon, int *ncocon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old,
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,   
	     int *kode, double *adb, double *aub,int *mei, double *fei,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mint_,
             int *ncmat_, int *nstate_, double *ener, char *jobnamec,
             char *output, char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     char *cbody, int *ibody,double *xbody, int *nbody);

void arpackbu(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b,int *nactdof, 
	     int *icol, int *jq, int *irow, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun, 
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs, 
	     int *kode, double *adb, double *aub,int *mei, double *fei,
             char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mint_,
             int *ncmat_, int *nstate_, double *ener, char *output, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     char *cbody, int *ibody, double *xbody, int *nbody);

void arpackcs(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int *icol, int *jq, int *irow, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun, 
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old,
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,int *mei, double *fei,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mint_,
             int *ics, double *cs, int *mpcend, int *ncmat_, int *nstate_,
             int *mcs, int *nkon, double *ener, char *jobnamec,
             char *output, char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     char *cbody, int *ibody,double *xbody, int *nbody, int *nevtot);

void FORTRAN(bodyforce,(char *cbody,int *ibody,int *ipobody,int *nbody,
             char *set,int *istartset,int *iendset,int *ialset,
             int *inewton,int *nset, int *ifreebody, int *k));

void calcresidual(int *nmethod, int *neq, double *b, double *fext, double *f,
        int *iexpl, int *nactdof, double *aux1, double *aux2, double *vold,
        double *vini, double *dtime, double *accold, int *nk, double *adb,
        double *aub, int *icol, int *irow, int *nzl, double *alpha,
	double *fextini, double *fini);

void FORTRAN(calinput,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	       int *nkon, int *ne,
	       int *nodeboun, int *ndirboun, double *xboun, int *nboun,
	       int *ipompc, int *nodempc, double *coefmpc, int *nmpc,
	       int *nmpc_,int *nodeforc, int *ndirforc, double *xforc,
	       int *nforc, int *nforc_, int *nelemload, char *sideload,
	       double *xload, int *nload, int *nload_,
	       int *nprint, char *prlab, char *prset,int *mpcfree, int *nboun_,
	       int *mei, char *set, int *istartset, 
	       int *iendset, int *ialset, int *nset, int *nalset,
	       double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	       double *alcon, int *nalcon, double *alzero, double *t0,
	       double *t1, char *matname,
	       int *ielmat, char *orname, double *orab, int *ielorien,
	       char *amname, double *amta, int *namta, int *nam,
	       int *nmethod, int *iamforc, int *iamload,
	       int *iamt1,int *ithermal, int *iperturb, 
	       int *istat, int *istep, int *nmat,
	       int *ntmat_, int *norien, double *prestr, int *iprestr,
	       int *isolver, double *fei,double *veold, double *tinc,
	       double *tper, double *xmodal,
               char *filab, int *jout, int *nlabel,
	       int *idrct, int *jmax, double *tmin, double *tmax,
	       int *iexpl, double *alpha, int *iamboun,
	       double *plicon, int *nplicon, double *plkcon, int *nplkcon,
	       int *iplas, int *npmat_, int *mint_, int *nk_,
	       double *trab, int *inotr, int *ntrans, int *ikboun,
               int *ilboun, int *ikmpc, int *ilmpc, int *ics,
	       double *dcs, int *ncs_, int *namtot_, double *cs,
               int *nstate_, int *ncmat_, int *iumat, int *mcs,
               char *labmpc, int *iponor,double *xnor,int *knor,
	       double *thickn,double *thicke,int *ikforc,int *ilforc,
               double *offset,int *iponoel,int *inoel,int *rig,
               int *infree, int *nshcon, double *shcon, double *cocon,
               int *ncocon,double *physcon, int *nflow, double *ctrl,
               int *memmpc_, int *maxlenmpc, int *ne1d, int *ne2d, int *nener, 
               double *vold, int *nodebounold,
               int *ndirbounold, double *xbounold, double *xforcold,
               double *xloadold, double *t1old, double *eme,
               double *sti, double *ener,
               double *xstate, char *jobnamec, int *nnn, int *irstrt,
               double *ttime,double *qaold,
               char *output, char *typeboun, char *inpc, int *nline,
               int *ipoinp,int *inp, char *tieset, double *tietol, 
               int *ntie, double *fmpc, char *cbody, int *ibody, double *xbody,
               int *nbody, int *nbody_, double *xbodyold, int *nam_,
               int *ielprop, int *nprop, int *nprop_, double *prop,
               int *itpamp, int *iviewfile, int *ipoinpc));    

void cascade(int *ipompc, double **coefmpcp, int **nodempcp, int *nmpc,
   int *mpcfree, int *nodeboun, int *ndirboun, int*nboun, int*ikmpc,
   int *ilmpc, int *ikboun, int *ilboun, int *mpcend, int *mpcmult,
   char *labmpc, int *nk, int *memmpc_, int *icascade, int *maxlenmpc,
   int *callfrommain, int *iperturb);

int cgsolver(double *A, double *x, double *b, int neq, int len, int *ia, int *iz, 
				double *eps, int *niter, int precFlg);

void checkconvergence(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	  int *ne, double *stn, int *nmethod, 
	  int *kode, char *filab, double *een, double *t1act,
          double *time, double *epn,int *ielmat,char *matname,
          double *enern, double *xstaten, int *nstate_, int *istep,
          int *iinc, int *iperturb, double *ener, int *mint_, char *output,
          int *ithermal, double *qfn, int *mode, int *noddiam, double *trab,
          int *inotr, int *ntrans, double *orab, int *ielorien, int *norien,
          char *description, double *sti,
	  int *icutb, int *iit, double *dtime, double *qa, double *vold,
          double *qam, double *ram1, double *ram2, double *ram,
          double *cam, double *uam, int *ntg, double *ttime,
          int *icntrl, double *theta, double *dtheta, double *veold,
          double *vini, int *idrct, double *tper,int *istab, double *tmax, 
	  int *nactdof, double *b, double *tmin, double *ctrl, double *amta,
          int *namta, int *itpamp, int *inext, double *dthetaref,int *itp,
          int *jprint, int *jout, int *uncoupled, double *t1, int *iitterm,
          int *nelemload, int *nload, int *nodeboun, int *nboun, int *itg,
          int *ndirboun, double *deltmx);

void checkconvgas(int *icutb, int *iin,
		  double *qamt, double *qamf, double *qamp, 
		  double *ram1t, double *ram1f, double *ram1p,
		  double *ram2t, double *ram2f, double *ram2p,
		  double *ramt, double *ramf, double *ramp,
		  int *icntrl, double *dtheta, double *ctrl);

void checkinclength(double *time,double *ttime,double *theta, double *dtheta,
          int *idrct, double *tper,double *tmax, double *tmin, double *ctrl, 
          double *amta,int *namta, int *itpamp, int *inext, double *dthetaref, 
	  int *itp,int *jprint, int *jout);

void FORTRAN(checktime,(int *itpamp,int *namta,double *tinc,double *ttime,
             double *amta,double *tmin, int *inext, int *itp));

void FORTRAN(closefile,());

  
void FORTRAN(compdt,(int *nk,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,
       int *ntmat_,int *iponoel,int *inoel,double *dtimef,int *iexplicit,
       int *ielmat,double *physcon,
       double *dh));

void compfluid(double *co, int *nk, int *ipkon, int *kon, char *lakon,
    int *ne, int *ipoface, char *sideface, int *ifreestream, 
    int *nfreestream, int *isolidsurf, int *neighsolidsurf,
    int *nsolidsurf, int *iponoel, int *inoel, int *nshcon, double *shcon,
    int *nrhcon, double *rhcon, double *vold, int *ntmat_,int *nodeboun, 
    int *ndirboun, int *nboun, int *ipompc,int *nodempc, int *nmpc,
    int *ikmpc, int *ilmpc, int *ithermal, int *ikboun, int *ilboun,
    int *turbulent, int *isolver, int *iexpl, double *voldtu, double *ttime,
    double *time, double *dtime, int *nodeforc,int *ndirforc,double *xforc,
    int *nforc, int *nelemload, char *sideload, double *xload,int *nload,
    double *xbody,int *ipobody,int *nbody, int *ielmat, char *matname,
    int *mint_, int *ncmat_, double *physcon, int *istep, int *iinc,
    int *ibody, double *xloadold, double *xboun,
    double *coefmpc, int *nmethod, double *xforcold, double *xforcact,
    int *iamforc,int *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, int *iamt1, double *amta,
    int *namta, int *nam, double *ampli, double *xbounold, double *xbounact,
    int *iamboun, int *itg, int *ntg, char *amname, double *t0, int *nelemface,
    int *nface, double *cocon, int *ncocon, double *xloadact, double *tper);

void contact(int *ncont, int *ntie, char *tieset, int *nset, char *set,
	     int *istartset, int *iendset, int *ialset, int *itietri,
	     char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
	     double *cg, double *straight, int *ifree, double *co,
	     double *vold, int *ielmat, double *cs, double *elcon,
             int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
             int *ifcont1, int *ifcont2, int *ne0, double *vini,
             int *nmethod);

void FORTRAN(cprint,(char *text, int *before, int *after));

void FORTRAN(datri,(double *au, double *aa, double *ad, int *icol, int *neq, 
	    int *flg));

void FORTRAN(dasol,(double *au, double *aa, double *ad, double *b, int *icol, 
	    int *neq, double *energy));
      
void FORTRAN(dgesv, (int *nteq, int *nhrs,double *ac,int *lda,int *ipiv,
                     double *bc,int *ldb,int *info)); 

void FORTRAN(drfftf,(int *ndata, double *r, double *wsave, int *isave));

void FORTRAN(drffti,(int *ndata, double *wsave, int *isave));

void FORTRAN(dsaupd,(int *ido, char *bmat, int *n, char *which, int *nev,
	     double *tol, double *resid, int *ncv, double *z, int *ldz,
	     int *iparam, int *ipntr, double *workd, double *workl,
	     int *lworkl, int *info));

void FORTRAN(dseupd,(int *, char *, int *, double *, double *,
	     int *, double *, char *, int *, char *, 
	     int *, double *, double *, int *, double *,
	     int *, int *, int *, double *,
	     double *, int *, int *));

void FORTRAN(dsort,(double *dx, int *iy, int *n, int *kflag));

void FORTRAN(dspgv,(int *itype, char *jobz, char *uplo, int *n, double *ap,
                    double *bp, double *w, double *z, int *ldz, double *work,
                    int *info));

void dyna(double **cop, int *nk, int **konp, int **ipkonp, char **lakonp, int *ne, 
	       int **nodebounp, int **ndirbounp, double **xbounp, int *nboun,
	       int **ipompcp, int **nodempcp, double **coefmpcp, char **labmpcp,
               int *nmpc, int *nodeforc,int *ndirforc,double *xforc, 
               int *nforc,int *nelemload, char *sideload,double *xload,
	       int *nload, 
	       int **nactdofp,int *neq, int *nzl,int *icol, int *irow, 
	       int *nmethod, int **ikmpcp, int **ilmpcp, int **ikbounp, 
	       int **ilbounp,double *elcon, int *nelcon, double *rhcon, 
	       int *nrhcon,double *cocon, int *ncocon,
               double *alcon, int *nalcon, double *alzero, 
               int **ielmatp,int **ielorienp, int *norien, double *orab, 
               int *ntmat_,double **t0, 
	       double **t1,int *ithermal,double *prestr, int *iprestr, 
	       double **voldp,int *iperturb, double *sti, int *nzs, 
	       double *tinc, double *tper, double *xmodal,
	       double **veoldp, char *amname, double *amta,
	       int *namta, int *nam, int *iamforc, int *iamload,
	       int **iamt1,int *jout,
	       int *kode, char *filab,double **emep, double *xforcold, 
	       double *xloadold,
               double **t1old, int **iambounp, double **xbounoldp, int *iexpl,
               double *plicon, int *nplicon, double *plkcon,int *nplkcon,
               double *xstate, int *npmat_, char *matname, int *mint_,
               int *ncmat_, int *nstate_, double **enerp, char *jobnamec,
               double *ttime, char *set, int *nset, int *istartset,
               int *iendset, int *ialset, int *nprint, char *prlab,
               char *prset, int *nener, double *trab, 
               int **inotrp, int *ntrans, double **fmpcp, char *cbody, int *ibody,
               double *xbody, int *nbody, double *xbodyold, int *istep,
	       int *isolver,int *jq, char *output, int *mcs, int *nkon,
               int *mpcend, int *ics, double *cs, int *ntie, char *tieset,
               int *idrct, int *jmax, double *tmin, double *tmax,
	       double *ctrl, int *itpamp, double *tietol);
 
void dynboun(double *amta,int *namta,int *nam,double *ampli, double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,int *iamboun,int *nboun,int *nodeboun,
             int *ndirboun, double *ad, double *au, double *adb,
             double *aub, int *icol, int *irow, int *neq, int *nzs,
             double *sigma, double *b, int *isolver,
             double *alpham, double *betam, int *nzl,
             int *init, double *bact, double *bmin, int *jq, char *amname);

void FORTRAN(dynresults,(int *nk,double *v,int *ithermal,int *nactdof,
             double *vold,int *nodeboun,int *ndirboun,double *xboun,
             int *nboun,int *ipompc,int *nodempc,double *coefmpc,
	     char *labmpc,int *nmpc,double *b));

void FORTRAN(dynsolv,(double *b, double *z, double *d, double *zeta, int *nev,
	      int *neq, double *tinc, int *jinc, int *jout, double *vold,
	      double *xforcold, int *nodeforc, int *ndirforc,
	      double *xforc, int *iamforc, int *nforc, double *xloadold,
	      double *xload, int *iamload, int *nelemload, char *sideload,
	      int *nload,
	      double *t1old,
	      double *t1, int *iamt1, int *nk, double *amta, int *namta, 
	      int *nam, double *ampli, double *aa, double *bb, double *bj, 
	      double *v, int *nodeboun, int *ndirboun, double *xboun, 
	      int *nboun, int *ipompc, int *nodempc, double *coefmpc, 
              char *labmpc,int *nmpc, int *nactdof, 
	      int *iperturb, int *nmethod, double *co, int *kon,
	      int *ipkon, char *lakon,int *ne, double *stn,double *stx,
	      double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	      double *alcon, int *nalcon, double *alzero, int *ielmat,
	      int *ielorien, int *norien, double *orab, int *ntmat_,
	      double *t0, int *ithermal, int *kode, double *cv,
	      double *cd, int *inum, double *prestr, int *iprestr,
	      int *ikmpc, int *ilmpc, int *ikboun, int *ilboun,
	      char *filab, double *eme, double *een,
	      double *sti, double *f, double *fn, double *xforcact,
	      double *xloadact,
              double *t1act, double *xbounold, double *xbounact,
              int *iamboun, int *iexpl, double *plicon, int *nplicon,
              double *plkcon,
	      int *nplkcon, double *xstateini,double *xstiff,
              double *xstate, int *npmat_, double *epn, char *matname,
              int *mint_, int *ncmat_, int *nstate_, double *stiini,
              double *vini, double *ener, double *enern, double *xstaten,
              double *ttime, double *eei, double *enerini,double *cocon,
              int *ncocon, char *set, int *nset, int *istartset,
              int *iendset, int *ialset, int *nprint, char *prlab,
              char *prset, double *qfx, double *qfn, double *trab,
	      int *inotr, int *ntrans, double *fmpc, double *veold,
	      char *cbody, int *ibody, double *xbody, int *nbody, 
              double *xbodyold, double *xbodyact, int *ipobody,
              double *cgr, double *xmodal,double *au,
              double *aub,double *vbounact, double *abounact, int *nzs));

void FORTRAN(envtemp,(int *itg,int *ieg,int *ntg,int *ntr,char *sideload,
                      int *nelemload,int *ipkon,int *kon,char *lakon,
                      int *ielmat,int *ne, int *nload, int *iptri,
                      int *kontri, int *ntri, int *nloadtr,
                      int *nflow, int *ndirboun, int *nactdog,
                      int *nodeboun, int *nacteq,
                      int *nboun, int *ielprop, double *prop, int *nteq,
                      double *v,int *network,double *physcon,
		      double *shcon,int *ntmat_,double *co, int *ipogn,
                      int *ign, double *vold,char *set,int *nshcon,
		      double *rhcon, int *nrhcon));

void FORTRAN(equationcheck,(double *ac,int *ntm,int *nteq,int *nactdog,
                            int *itg,int *ntg, int *nacteq, int *network));

void expand(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, int *nodeforc, int *ndirforc,double *xforc, 
             int *nforc, int *nelemload, char *sideload, double *xload,
             int *nload, int *nactdof, int *neq, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat,
	     double *t0,int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     double *adb, double *aub,char *filab, double *eme,
             double *plicon, int *nplicon, double *plkcon,int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mint_,
	     int *ics, double *cs, int *mpcend, int *ncmat_,
             int *nstate_, int *mcs, int *nkon, double *ener,
             char *jobnamec, char *output, char *set, int *nset,int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     int *nev, double *z, int *iamboun, double *xbounold,
	     int *nsegments, int *nm,int *icol,int *irow,int *nzl,int *nam,
	     int *ipompcold, int *nodempcold, double *coefmpcold,
             char *labmpcold, int *nmpcold, double *xloadold, int *iamload,
             double *t1old,double *t1,int *iamt1, double *xstiff);

double FORTRAN(fcrit,(double *time,double *tend,double *aai,double *bbi,
		      double *zetaj,double *dj,double *ddj,
		      double *h1,double *h2,double *h3,double *h4));

void FORTRAN(updatecfd,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
           int *nactdoh));

void FORTRAN (flowoutput,(int *itg,int *ieg,int *ntg,int *ntm,
			  double *bc,char *lakon,
			  int *ntmat_,double *v,double *shcon,int *nshcon,
			  int *ipkon,int *kon,double *co,int *nflow,
			  double *dtime,double *ttime,double *time,
			  int *ielmat,double *prop,
			  int *ielprop,int *nactdog,int *nacteq,int *iin,
			  double *physcon, double *camt, double *camf,double *camp,
			  double *uamt, double *uamf,double *uamp,
			  double *rhcon, int *nrhcon,
			  double *vold,char *jobnamef, char *set));

void FORTRAN(flowresult,(int *ntg, int *itg,double *cam,double *vold,
              double *v,
              int *nload, char *sideload, int *nelemload,
              double *xloadact, int *nactdog, int *network));

void FORTRAN(frddummy,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *vold,int *kode,double *time,
             int *ielmat,char *matname,int *nnstep,double *vtu,
             double *voldtu,double *voldaux));

void frdcyc(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne,double *v,
	    double *stn,int *inum,int *nmethod,int *kode,char *filab,
	    double *een,double *t1,double *fn,double *time,double *epn,
	    int *ielmat,char *matname, double *cs, int *mcs, int *nkon,
	    double *enern, double *xstaten, int *nstate_, int *istep,
            int *iinc, int *iperturb, double *ener, int *mint_, char *output,
            int *ithermal, double *qfn, int *ialset, int *istartset,
            int *iendset, double *trab, int *inotr, int *ntrans,double *orab,
	    int *ielorien, int *norien, double *sti);

void FORTRAN(frdphase,(int *kode,double *time,int *nk,int *inum,
             double *vr, double *vi, double *stnr, double *stni,
             char *filab, int *mode, int *noddiam, int *nmethod,
	     double *vmax, double *stnmax, int *nkcoords));

double FORTRAN(fsub,(double *time,double *tend,double *aai,double *bbi,
		     double *ddj,double *h1,double *h2,double *h3,double *h4));

double FORTRAN(fsuper,(double *time,double *tend,double *aai,double *bbi,
		       double *h1,double *h2,double *h3,double *h4,
		       double *h5,double *h6));

void FORTRAN(gasmechbc,(double *vold,int *nload,char *sideload,
			int *nelemload,double *xload));

void FORTRAN(gencontelem,(char *tieset,int *ntie,int *itietri,int *ne,
     int *ipkon,int *kon,char *lakon,char *set,int *istartset,int *iendset,
     int *ialset,double *cg,double *straight,int *ifree,int *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,int *nx,int *ny,int *nz, int *nset,
     int *ielmat,double *cs, double *elcon,int *istep, int *iinc, int *iit,
     int *ncmat_,int *ntmat_, int *ifcont1, int *ifcont2, int *ne0,
     double *vini, int *nmethod));

void FORTRAN(getfneig,(char *fneig));

void FORTRAN(identamta,(double *amta,double *reftime,int *istart,int *iend,
               int *id));

void FORTRAN(includefilename,(char *buff,char *includefn,int *lincludefn));

void inicont(int *ncont, int *ntie, char *tieset, int *nset, char *set,
             int *istartset, int *iendset, int *ialset, int **itietrip,
	     char *lakon, int *ipkon, int *kon, int **koncontp,
             int *ncone, double *tietol, int *ismallsliding);
  
void FORTRAN(initialcfd,(double *yy,int *nk,double *co,int *ne,int *ipkon,
       int *kon,char *lakon,double *x,double *y,double *z,double *x0,
       double *y0,double *z0,int *nx,int *ny,int *nz,int *isolidsurf,
       int *neighsolidsurf,double *xsolidsurf,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,double *voldaux,
       int *ntmat_,int *iponoel,int *inoel,int *iexplicit,
       int *ielmat, int *nsolidsurf, int *turbulent, double *physcon));

void FORTRAN(initialgas,(int *itg,int *ieg,int *ntm,int *ntg,double *ac,
                         double *bc,
                         char *lakon,double *v,int * ipkon,int *kon,
                         int *nflow,int *ikboun,int *nboun,double *prop,
                         int *ielprop,int *nactdog,int *nacteq,int *ndirboun,
                         int *nodeboun,double *xbounact, int *ielmat,
                         int *ntmat_,double *shcon,int *nshcon,
                         double *physcon, int *ipiv, int *nteq,
                         double *rhcon, int *nrhcon, int *ipobody, int *ibody,
                         double *xbody, double *co, int *nbody, int *network,
                         int *iin_abs, double *vold, char *set));

void insert(int *ipointer, int **mast1p, int **mast2p, int *i1,
	    int *i2, int *ifree, int *nzs_);

void FORTRAN(isortii,(int *ix, int *iy, int *n, int *kflag));

void FORTRAN(isortiid,(int *ix, int *iy, double *dy1, int *n, int *kflag));

void FORTRAN(isortiddc1,(int *ix, double *dy1, double *dy2, char *cy, int *n, 
                         int *kflag));

void FORTRAN(isortiddc2,(int *ix1, int *ix2,double *dy1, double *dy2, 
                         char *cy, int *n, int *kflag));

void FORTRAN(iter,(double *coef, int *jcoef, int *ndim, int *n, int *p,
	   int *ip, double *u, double *ubar, double *rhs, double *wksp,
	   int *iwksp, int *nw, int *inw, int *iparm, double *rparm));

void FORTRAN(keystart,(int *ifreeinp,int *ipoinp,int *inp,char *name,
           int *iline,int *ikey));

void FORTRAN(mafilldm,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	       int *ne, int *nodeboun, int *ndirboun, double *xboun, 
	       int *nboun,int *ipompc, int *nodempc, double *coefmpc, 
	       int *nmpc, int *nodeforc,int *ndirforc,
	       double *xforc, int *nforc, int *nelemload, char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr,
	       double *ad, double *au, int *nactdof, 
	       int *icol, int *jq, int *irow, int *neq, int *nzl, 
	       int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	       int *ilboun,
	       double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	       double *alcon, int *nalcon, double *alzero, int *ielmat,
	       int *ielorien, int *norien, double *orab, int *ntmat_,
	       double *t0, double *t1, int *ithermal,
	       double *prestr, int *iprestr, double *vold,
	       int *iperturb, double *sti, int *nzs, double *stx,
	       double *adb, double *aub, int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstiff, 
	       int *npmat_, double *dtime, char *matname, int *mint_,
               int *ncmat_,double *ttime, double *time,
               int *istep, int *kinc, int *ibody));

void FORTRAN(mafillgas,(int *itg,int *ieg,int *ntg,int *ntm,
			double *ac,int *nload,char *sideload,
			int *nelemload,double *xloadact,char *lakon,
			int *ntmat_,double *v,double *shcon,int *nshcon,
			int *ipkon,int *kon,double *co,int *nflow,
			int *iinc,int *istep,
			double *dtime,double *ttime,double *time,
			int *ielmat,int *nteq,double *prop,
			int *ielprop,int *nactdog,int *nacteq,
			double *physcon, double *rhcon, int *nrhcon,
			int *ipobody, int *ibody, double *xbody, int *nbody,
			double *vold, double *xloadold, double *reltime,
			int *nmethod, char *set));

void FORTRAN(mafillklhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
         int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
         int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdok,
         int *icolk,int *jqk,int *irowk,int *neqk,int *nzlk,int *ikmpc,
         int *ilmpc,int *ikboun,int *ilboun,int *nzsk,double *adbk,
         double *aubk));

void FORTRAN(mafillkrhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
        int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
        int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nelemface,
        char *sideface,int *nface,int *nactdok,int *neqk,int *nmethod,
        int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,double *rhcon,
        int *nrhcon,int *ielmat,int *ntmat_,double *vold,double *voldaux,
        int *nzsk,double *dtimef,char *matname,int *mint_,int *ncmat_,
        double *shcon,int *nshcon,double *v,double *theta1,double *bk,
        double *bt,double *voldtu,int *isolidsurf,int *nsolidsurf,
        int *ifreestream,int *nfreestream,double *xsolidsurf));

void FORTRAN(mafillplhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
          int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
          int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
          int *icolp,int *jqp,int *irowp,int *neqp,int *nzlp,int *ikmpc,
          int *ilmpc,int *ikboun,int *ilboun,int *nzsp,double *adbp,
	  double *aubp,int *nmethod,int *iexplicit));

void FORTRAN(mafillprhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nelemface,
       char *sideface,int *nface,double *b,int *nactdoh,int *icolp,int *jqp,
       int *irowp,int *neqp,int *nzlp,int *nmethod,int *ikmpc,int *ilmpc,
       int *ikboun,int *ilboun,double *rhcon,int *nrhcon,int *ielmat,
       int *ntmat_,double *vold,double *voldaux,int *nzsp,double *dtimef,
       char *matname,int *mint_,int *ncmat_,double *shcon,int *nshcon,
       double *v,double *theta1,int *iexplicit,
       double *physcon));

void FORTRAN(mafillsm,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	       int *ne, int *nodeboun, int *ndirboun, double *xboun, 
	       int *nboun,int *ipompc, int *nodempc, double *coefmpc, 
	       int *nmpc, int *nodeforc,int *ndirforc,
	       double *xforc, int *nforc, int *nelemload, char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr,
	       double *ad, double *au, double *bb, int *nactdof, 
	       int *icol, int *jq, int *irow, int *neq, int *nzl, 
	       int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	       int *ilboun,
	       double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	       double *alcon, int *nalcon, double *alzero, int *ielmat,
	       int *ielorien, int *norien, double *orab, int *ntmat_,
	       double *t0, double *t1, int *ithermal,
	       double *prestr, int *iprestr, double *vold,
	       int *iperturb, double *sti, int *nzs, double *stx,
	       double *adb, double *aub, int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstiff, 
	       int *npmat_, double *dtime, char *matname, int *mint_,
               int *ncmat_,int *mass,int *stiffness,int *buckling,int *rhs,
               int *intscheme, double *physcon, double *shcon, int *nshcon,
               double *cocon, int *ncocon, double *ttime, double *time,
               int *istep, int *kinc, int *coriolis, int *ibody,
               double *xloadold, double *reltime));

void FORTRAN(mafillsmcs,(double *co, int *nk, int *kon, int *ipkon, 
               char *lakon,
	       int *ne, int *nodeboun, int *ndirboun, double *xboun, 
	       int *nboun,int *ipompc, int *nodempc, double *coefmpc, 
	       int *nmpc, int *nodeforc,int *ndirforc,
	       double *xforc, int *nforc, int *nelemload, char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr, 
	       double *ad, double *au, double *bb, int *nactdof, 
	       int *icol, int *jq, int *irow, int *neq, int *nzl, 
	       int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	       int *ilboun,
	       double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	       double *alcon, int *nalcon, double *alzero, int *ielmat,
	       int *ielorien, int *norien, double *orab, int *ntmat_,
	       double *t0, double *t1, int *ithermal,
	       double *prestr, int *iprestr, double *vold,
	       int *iperturb, double *sti, int *nzs, double *stx,
	       double *adb, double *aub, int *iexpl,double *plicon,
               int *nplicon,double *plkcon,int *nplkcon,double *xstiff, 
	       int *npmat_, double *dtime, char *matname, int *mint_,
               int *ics, double *cs, int *nm, int *ncmat_,char *labmpc,
               int *mass,int *stiffness,int *buckling,int *rhs,
               int *intscheme, int *mcs, int *coriolis, int *ibody,
               double *xloadold, double *reltime, int *ielcs));

void FORTRAN(mafilltlhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
       int *icolt,int *jqt,int *irowt,int *neqt,int *nzlt,int *ikmpc,
       int *ilmpc,int *ikboun,int *ilboun,int *nzst,double *adbt,
       double *aubt));
	  
void FORTRAN(mafilltrhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
             int *ipompc,int *nodempc,double *coefmpc,int *nmpc,
             int *nodeforc,int *ndirforc,double *xforc,int *nforc,
             int *nelemload,char *sideload,double *xload,int *nload,
             double *xbody,int *ipobody,int *nbody,double *b,int *nactdoh,
             int *neqt,int *nmethod,int *ikmpc,int *ilmpc,int *ikboun,
             int *ilboun,double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,
             double *t0,int *ithermal,double *vold,double *voldaux,int *nzst,
             double *dtimef,char *matname,int *mint_,int *ncmat_,
             double *physcon,double *shcon,int *nshcon,double *ttime,
             double *timef,int *istep,int *iinc,int *ibody,double *xloadold,
             double *reltime,double *cocon,int *ncocon,int *nelemface,
             char *sideface,int *nface));

void FORTRAN(mafillvlhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
       int *icolv,int *jqv,int *irowv,int *neqv,int *nzlv,int *ikmpc,
       int *ilmpc,int *ikboun,int *ilboun,int *nzsv,double *adbv,
       double *aubv));

void FORTRAN(mafillv1rhs,(double *co,int *nk,int *kon,int *ipkon,
         char *lakon,int *ne,int *nodeboun,int *ndirboun,
	 double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
         int *nmpc,int *nodeforc,int *ndirforc,double *xforc,
	 int *nforc,int *nelemload,char *sideload,double *xload,
         int *nload,double *xbody,int *ipobody,int *nbody,
         double *b,int *nactdoh,int *icolv,int *jqv,int *irowv,
         int *neqv,int *nzlv,int *nmethod,int *ikmpc,int *ilmpc,
         int *ikboun,int *ilboun,double *rhcon,int *nrhcon,int *ielmat,
         int *ntmat_,double *t0,int *ithermal,double *vold,
         double *voldaux,int *nzsv,double *dtimef,char *matname,
         int *mint_,int *ncmat_,double *physcon,double *shcon,int *nshcon,
         double *ttime,double *timef,int *istep,int *iinc,int *ibody,
         double *xloadold,int *turbulent,double *voldtu,
         double *yy,int *nelemface,char *sideface,int *nface));

void FORTRAN(mafillv2rhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,double *
       b,int *nactdoh,int *icolv,int *jqv,int *irowv,int *neqv,int *nzlv,
       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
       double *vold,int *nzsv,double *dtimef,double *v,double *theta2,
       int *iexplicit));

void mastruct(int *nk, int *kon, int *ipkon, char*lakon, int *ne,
	      int *nodeboun, int *ndirboun, int *nboun, int *ipompc,
	      int *nodempc, int *nmpc, int *nactdof, int *icol,
	      int *jq, int **mast1p, int **irowp, int *isolver, int *neq,
	      int *nnn, int *ikmpc, int *ilmpc,
	      int *ipointer, int *nzs, int *nmethod,
              int *ithermal, int *ikboun, int *ilboun, int *iperturb);

void mastructcs(int *nk, int *kon, int *ipkon, char *lakon,
	       int *ne, int *nodeboun,
	       int *ndirboun, int *nboun, int *ipompc, int *nodempc,
	       int *nmpc, int *nactdof, int *icol, int *jq, int **mast1p,
	       int **irowp, int *isolver, int *neq, int *nnn, 
	       int *ikmpc, int *ilmpc, int *ipointer,
	       int *nzs, int *nmethod, int *ics, double *cs,
               char *labmpc, int *mcs);

void mastructf(int *nk, int *kon, int *ipkon, char *lakon, int *ne,
	      int *nodeboun, int *ndirboun, int *nboun, int *ipompc,
	      int *nodempc, int *nmpc, int *nactdoh, int *icolt,
	      int *icolv, int *icolp, int *icolk,int *jqt, int *jqv, int *jqp,
	      int *jqk,int **mast1p, int **irowtp, int **irowvp, int **irowpp, 
	      int **irowkp,int *isolver, int *neqt, int *neqv, int *neqp,
	      int *neqk,int *ikmpc, int *ilmpc,int *ipointer, 
	      int *nzst, int *nzsv, int *nzsp, int *nzsk, 
	      int *ithermal, int *ikboun, int *ilboun, int *turbulent,
              int *nactdok, int *ifreestream, int *nfreestream,
              int *isolidface, int *nsolidface, int *nzs);

void FORTRAN(mult,(double *matrix, double *trans, int *n));

void FORTRAN(nident,(int *x, int *px, int *n, int *id));

void nonlingeo(double **co, int *nk, int **konp, int **ipkonp, char **lakonp,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int **nodempcp, double **coefmpcp, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int **icolp, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int **ielmatp,
	     int **ielorienp, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double **vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub, 
	     char *filab, int *idrct,
	     int *jmax, int *jout, double *tinc, double *tper,
	     double *tmin, double *tmax, double *eme, double *xbounold,
	     double *xforcold, double *xloadold,
             double *veold, double *accold,
             char *amname, double *amta, int *namta, int *nam,
             int *iamforc, int *iamload,
             int *iamt1, double *alpha, int *iexpl,
	     int *iamboun, double *plicon, int *nplicon, double *plkcon,
	     int *nplkcon,
             double *xstate, int *npmat_, int *istep, double *ttime,
	     char *matname, double *qaold, int *mint_,
             int *isolver, int *ncmat_, int *nstate_, int *iumat,
             double *cs, int *mcs, int *nkon, double *ener, int *mpcinfo,
             int *nnn, char *output,
             double *shcon, int *nshcon, double *cocon, int *ncocon,
             double *physcon, int *nflow, double *ctrl, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener,int *ikforc, int *ilforc, double *trab, 
             int *inotr, int *ntrans, double *fmpc, char *cbody,
             int *ibody, double *xbody, int *nbody, double *xbodyold,
             int *ielprop, double *prop, int *ntie, char *tieset,
	     int *itpamp, int *iviewfile, char *jobnamec, double *tietol);

void FORTRAN(nonlinmpc,(double *co,double *vold,int *ipompc,int *nodempc,
		   double *coefmpc,char *labmpc,int *nmpc,int *ikboun,
		   int *ilboun,int *nboun,double *xbounact,double *aux,
		   int *iaux, int *maxlenmpc, int *ikmpc, int *ilmpc,
                   int *icascade, int *kon, int *ipkon, char *lakon,
		   int *ne, double *reltime, int *newstep, double *xboun,
		   double *fmpc, int *newinc, int *idiscon, int *ncont,
		   double *trab, int *ntrans));

void FORTRAN(op,(int *, double *, double *, double *, double *, double *, int *,
	 int *, int *));

void FORTRAN(opcs,(int *, double *, double *, double *, double *, double *, int *,
	 int *, int *));

void FORTRAN(openfile,(char *jobname, char *output));

void FORTRAN(out,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	  int *ne, double *v, 
	  double *stn, int *inum, int *nmethod, 
	  int *kode, char *filab, double *een, double *t1,
          double *fn, double *time, double *epl,int *ielmat,char *matname,
          double *enern, double *xstaten, int *nstate_, int *istep,
          int *iinc, int *iperturb, double *ener, int *mint_, char *output,
          int *ithermal, double *qfn, int *mode, int *noddiam, double *trab,
          int *inotr, int *ntrans, double *orab, int *ielorien, int *norien,
	  char *description, int *ipneigh,int *neigh,double *sti,
          double *vr,double *vi,double *stnr,double *stni,
	  double *vmax,double *stnmax,int *ngraph));

void FORTRAN(precfd,(int *nelemface,char *sideface,int *nface,int *ipoface,
        int *nodface,int *cfd,int *ne,int *ipkon,int *kon,char *lakon,
        int *ikboun,int *ilboun,double *xboun,int *nboun,int *nk,
        int *isolidsurf,int *nsolidsurf,int *ifreestream,int *nfreestream,
        int *neighsolidsurf,int *iponoel,int *inoel,int *inoelfree,
        int *nef, double *co));

void prediction(double *uam, int *nmethod, double *bet, double *gam, double *dtime,
               int *ithermal, int *nk, double *veold, double *accold, double *v,
	       int *iinc, int *idiscon, double *vold, int *nactdof);

void preiter(double *ad, double **aup, double *b, int **icolp, int **irowp, 
	     int *neq, int *nzs, int *isolver, int *iperturb);
  
void prespooles(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int **icolp, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs, 
	     int *kode, double *adb, double *aub, 
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *isolver,
	     int *mint_, int *ncmat_, int *nstate_, double *cs, int *mcs,
             int *nkon, double *ener, double *xbounold,
	     double *xforcold, double *xloadold,
             char *amname, double *amta, int *namta,
             int *nam, int *iamforc, int *iamload,
             int *iamt1, int *iamboun, double *ttime, char *output, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, double *trab, 
             int *inotr, int *ntrans, double *fmpc, char *cbody, int *ibody,
	     double *xbody, int *nbody, double *xbodyold);

void radcyc(int *nk,int *kon,int *ipkon,char *lakon,int *ne,
	    double *cs, int *mcs, int *nkon,int *ialset, int *istartset,
            int *iendset,int **kontrip,int *ntri,
            double **cop, double **voldp,int *ntrit, int *inocs);

void radflowload(int *itg, int *ieg, int *ntg, int *ntr,int *ntm,
       double *ac,double *bc,int *nload,
       char *sideload,int *nelemload,double *xloadact,char *lakon,int *ipiv,
       int *ntmat_,double *vold,double *shcon,int *nshcon,int *ipkon,
       int *kon,double *co,double *pmid,double *e1,double *e2,double *pnor,
       int *iptr,int *kontri,int *ntri,
       int *nloadtr,double *tarea,double *tenv,double *physcon,double *erad,
       double *fij,double *dist,int *idist,double *area,
       int *nflow,int *ikboun,double *xboun,int *nboun,int *ithermal,
       int *iinc,int *iit,double *cs, int *mcs, int *inocs, int *ntrit,
       int *nk, double *fenv,int *istep,double *dtime,double *ttime,
       double *time,int *ilboun,int *ikforc,int *ilforc,double *xforc,
       int *nforc,double *cam, int *ielmat, int *nteq, double *prop,
       int *ielprop, int *nactdog, int *nacteq, int *nodeboun, int *ndirboun,
       int *network, double *rhcon, int *nrhcon,
       int *ipobody, int *ibody, double *xbody, int *nbody, int *iviewfile,
       char *jobnamef, double *ctrl, double *xloadold, double *reltime,
       int *nmethod,char *set);

void FORTRAN (radmatrix,(int *ntr,int *ntm,double *ac,double *bc,
       char *sideload,int *nelemload,double *xloadact,char *lakon,
       double *vold,int *ipkon,int *kon,double *co,double *pmid,double *e1,
       double *e2,double *e3,int *iptri,int *kontri,int *ntri,int *nloadtr,
       double *tarea,double *tenv,double *physcon,double *erad,
       double *f,double *dist,int *idist,double *area,int *ithermal,int *iinc,
       int *iit,double *cs,int *mcs,int *inocs,int *ntrit,int *nk,
       double *fenv,int *istep,
       double *dtime,double *ttime,double *time,int *iviewfile,
       char *jobnamef, double *xloadold, double *reltime, int *nmethod));

void FORTRAN (radresult,(int *ntr,double *xloadact,int *ntm,double *bc,
       int *nloadtr,double *tarea,double * tenv,double *physcon,double *erad,
       double *f,double *fenv));

void readinput(char *jobnamec,char **inpcp,int *nline,int *nset, int *ipoinp,
        int **inpp, int **ipoinpcp, int *ithermal); 

void FORTRAN(rearrange,(double *au,int *irow,int *icol,int *ndim,int *neq));


void FORTRAN(rectcyl,(double *co,double *v,double *fn,double *stn,
		      double *qfn, double *een,double *cs,int *nk, 
                      int *icntrl,double *t, char *filab, int *imag));

void remastruct(int *ipompc, double **coefmpcp, int **nodempcp, int *nmpc,
              int *mpcfree, int *nodeboun, int *ndirboun, int *nboun,
              int *ikmpc, int *ilmpc, int *ikboun, int *ilboun,
              char *labmpc, int *nk,
              int *memmpc_, int *icascade, int *maxlenmpc,
              int *kon, int *ipkon, char *lakon, int *ne, int *nnn,
              int *nactdof, int *icol, int *jq, int **irowp, int *isolver,
              int *neq, int *nzs,int *nmethod, double **fp,
              double **fextp, double **bp, double **aux2p, double **finip,
              double **fextinip,double **adbp, double **aubp, int *ithermal,
              int *iperturb);

void FORTRAN(renumber,(int *nk, int *kon, int *ipkon, char *lakon, int *ne, 
	       int *ipompc, int *nodempc, int *nmpc, int *nnn, int *npn, 
	       int *adj,int *xadj, int *iw, int *mmm, int *xnpn, int *inum1, 
               int *inum2));

void FORTRAN(restartshort,(int *nset,int *nload,int *nbody, int *nforc,
    int *nboun,
    int *nk,int *ne,int *nmpc,int *nalset,int *nmat,int *ntmat,int *npmat,
    int *norien,int *nam,int *nprint,int *mint,int *ntrans,int *ncs,
    int *namtot,int *ncmat,int *memmpc,int *ne1d,int *ne2d,int *nflow,
    char *set,int *meminset,int *rmeminset,char *jobnamec,int *irestartstep,
    int *icntrl,int *ithermal,int *nener,int *nstate_,int *ntie));

void FORTRAN(restartwrite,(int *istep, int *nset, int*nload, int *nforc, 
  int * nboun,int *nk, int *ne, int *nmpc, int *nalset, int *nmat, int *ntmat_, 
  int *npmat_, int *norien, int *nam, int *nprint, int *mint_, 
  int *ntrans, int *ncs_, int *namtot_, int *ncmat_, int *mpcend, 
  int *maxlenmpc, int *ne1d, 
  int *ne2d, int *nflow, int *nlabel, int *iplas, int *nkon, int *ithermal, 
  int *nmethod,int *iperturb, int *nstate_,int *nener,char *set, 
  int *istartset, int *iendset, int *ialset, double *co, int *kon, int *ipkon, 
  char *lakon, int *nodeboun, int *ndirboun, int *iamboun, double *xboun, 
  int *ikboun, int *ilboun, int *ipompc, int *nodempc, double *coefmpc, 
  char *labmpc, int *ikmpc, int *ilmpc, int *nodeforc, int *ndirforc, 
  int *iamforc, double *xforc, int *ikforc, int *ilforc, int *nelemload, 
  int *iamload, char *sideload, double *xload,  
  double *elcon, int *nelcon, double *rhcon, int *nrhcon, double *alcon, 
  int *nalcon, double *alzero, double *plicon, int *nplicon, double *plkcon, 
  int *nplkcon, char *orname, double *orab, int *ielorien, double *trab, 
  int *inotr, char *amname, double *amta, int *namta, double *t0, double *t1, 
  int *iamt1, double *veold, int *ielmat,char *matname, 
  char *prlab, char *prset, char *filab, double *vold, 
  int *nodebounold, int *ndirbounold, double *xbounold, double *xforcold, 
  double *xloadold, double *t1old, double *eme, int *iponor, 
  double *xnor, int *knor, double *thickn, double *thicke, double *offset, 
  int *iponoel, int *inoel, int *rig, 
  double *shcon, int *nshcon, double *cocon, int *ncocon, 
  int *ics, double *sti, double *ener, double *xstate, 
  char *jobnamec, int *infree, int *nnn, double *prestr,int *iprestr,
  char *cbody, int *ibody, double *xbody, int *nbody, double *xbodyold,
  double *ttime,double *qaold,double *cs,
  int *mcs, char *output,double *physcon, double *ctrl, char *typeboun,
  double *fmpc,char *tieset, int *ntie));

void FORTRAN(resultgas,(int *itg,int *ieg,int *ntg,int *ntm,
                        double *bc,int *nload,char *sideload,
                        int *nelemload,double *xloadact,char *lakon,
                        int *ntmat_,double *v,double *shcon,int *nshcon,
                        int *ipkon,int *kon,double *co,int *nflow,
			int *iinc,int *istep,
                        double *dtime,double *ttime,double *time,
			int *ikforc,int *ilforc,
                        double *xforcact,int *nforc,
                        int *ielmat,int *nteq,double *prop,
                        int *ielprop,int *nactdog,int *nacteq,int *iin,
                        double *physcon, double *camt, double *camf,
                        double *camp, double *rhcon, int *nrhcon,
			int *ipobody, int *ibody, double *xbody, int *nbody,
                        double *dtheta, double *vold, double *xloadold,
                        double *reltime, int *nmethod,char *set));

void FORTRAN(results,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, double *v, double *stn, int *inum, 
	     double *stx,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1,int *ithermal,double *prestr, 
             int *iprestr, char *filab, double *eme, 
             double *een, int *iperturb,double *f, double *fn, int *nactdof,
             int *iout, double *qa,
	     double *vold, double *b, int *nodeboun, int *ndirboun,
	     double *xboun, int *nboun, int *ipompc, int *nodempc,
	     double *coefmpc, char *labmpc, int *nmpc, int *nmethod, 
             double *vmax,int *neq, double *veold, double *accold,
	     double *beta, double *gamma, double *dtime, double *time,
             double *ttime, double *plicon,
             int *nplicon,double *plkcon,int *nplkcon,
             double *xstateini,double *xstiff,double *xstate, int *npmat_,
	     double *epl, char *matname, int *mint_, int *ielas,
	     int *icmd, int *ncmat_, int *nstate_, double *stiini,
	     double *vini, int *ikboun, int *ilboun, double *ener,
	     double *enern, double *sti, double *xstaten, double *eei,
             double *enerini, double *cocon, int *ncocon, char *set, 
             int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, double *qfx, double *qfn, double *trab,
             int *inotr, int *ntrans, double *fmpc, int *nelemload,
             int *nload));

void FORTRAN(resultsk,(int *nk,int *nactdoh,double *vtu,double *solk,
      double *solt,int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultsp,(int *nk,int *nactdoh,double *v,double *sol,
      int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultst,(int *nk,int *nactdoh,double *v,double *sol,
      int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultsv1,(int *nk,int *nactdoh,double *v,double *sol,
      int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultsv2,(int *nk,int *nactdoh,double *v,double *sol,
      int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(rhs,(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	       int *ne,int *ipompc, int *nodempc, double *coefmpc, 
	       int *nmpc, int *nodeforc,int *ndirforc,
	       double *xforc, int *nforc, int *nelemload, char *sideload,
	       double *xload,int *nload, double *xbody, int *ipobody,
               int *nbody, double *cgr, double *bb, int *nactdof, int *neq, 
	       int *nmethod, int *ikmpc, int *ilmpc,
	       double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	       double *alcon, int *nalcon, double *alzero, int *ielmat,
	       int *ielorien, int *norien, double *orab, int *ntmat_,
	       double *t0, double *t1, int *ithermal, 
               int *iprestr, double *vold,int *iperturb,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               int *npmat_, double *ttime, double *time, int *istep,
               int *iinc, double *dtime, double *physcon, int *ibody,
               double *xbodyold, double *reltime));

void FORTRAN(solveeq,(double *adbv,double *aubv,double *adl,double *addiv,
            double *b,double *sol,double *aux,int *icolv,int *irowv,
            int *jqv,int *neqv,int *nzsv,int *nzlv));

void FORTRAN(spcmatch,(double *xboun, int *nodeboun, int *ndirboun, int *nboun,
	       double *xbounold, int *nodebounold, int *ndirbounold,
	       int *nbounold, int *ikboun, int *ilboun, double *vold,
	       double *reorder, int *nreorder));

void FORTRAN(splitline,(char *text,char *textpart,int *n));

void spooles(double *ad, double *au, double *adb, double *aub,
             double *sigma, double *b,
	     int *icol, int *irow, int *neq, int *nzs, int *symmtryflag,
             int *inputformat);

void FORTRAN(springforc,(double *xl,int *konl,double *vl,int *imat,
             double *elcon,int *nelcon,double *elas,double *fnl,int *ncmat_,
             int *ntmat_,int *nope,char *lakonl,double *t0l,double *t1l,
             int *kodem,double *elconloc,double *plicon,int *nplicon,
	     int *npmat_,double *veloldl));

void steadystate(double **co, int *nk, int **kon, int **ipkon, char **lakon, int *ne, 
	  int **nodeboun, int **ndirboun, double **xboun, int *nboun,
	  int **ipompcp, int **nodempcp, double **coefmpcp, char **labmpcp, int *nmpc, 
	  int *nodeforc,int *ndirforc,double *xforc, int *nforc, 
	  int *nelemload, char *sideload,double *xload,
	  int *nload, 
	  int **nactdof, int *neq, int *nzl,int *icol, int *irow, 
	  int *nmethod, int **ikmpcp, int **ilmpcp, int **ikboun, 
	  int **ilboun,
	  double *elcon, int *nelcon, double *rhcon, int *nrhcon,
          double *cocon, int *ncocon,
	  double *alcon, int *nalcon, double *alzero, int **ielmat,
	  int **ielorien, int *norien, double *orab, int *ntmat_,
	  double **t0, 
	  double **t1,int *ithermal,double *prestr, int *iprestr, 
	  double **voldp,int *iperturb, double *sti, int *nzs, 
	  double *tinc, double *tper, double *xmodal,
	  double *veold, char *amname, double *amta,
	  int *namta, int *nam, int *iamforc, int *iamload,
	  int **iamt1, int *jout,int *kode, char *filab,
	  double **emep, double *xforcold, double *xloadold,
          double **t1old, int **iamboun,
          double **xbounold, int *iexpl, double *plicon, int *nplicon,
          double *plkcon,int *nplkcon,
          double *xstate, int *npmat_, char *matname, int *mint_,
          int *ncmat_, int *nstate_, double **enerp, char *jobnamec,
          double *ttime, char *set, int *nset, int *istartset,
          int *iendset, int *ialset, int *nprint, char *prlab,
          char *prset, int *nener, double *trab, 
          int **inotr, int *ntrans, double **fmpcp, char *cbody, int *ibody,
          double *xbody, int *nbody, double *xbodyold,int *istep,
          int *isolver, int *jq, char *output, int *mcs, int *nkon,
          int *ics, double *cs, int *mpcend);

void FORTRAN(stop,());

void FORTRAN(storeresidual,(int *nactdof, double *b, double *fn, char *filab,
             int *ithermal, int *nk, double *sti, double *stn,
             int *ipkon, int *inum, int *kon, char *lakon,
             int *ne, int *mint_, double *orab, int *ielorien,
             double *co, int *nelemload, int *nload, int *nodeboun,
             int *nboun, int *itg, int *ntg, double *vold, int *ndirboun));

int strcmp1(const char *s1, const char *s2);

int strcpy1(char *s1, const char *s2, int length);

void FORTRAN(subspace,(double *d,double *aa,double *bb,double *cc,
             double *alpham,double *betam,int *nev,
             double *xini,double *cd,double *cv,double *time,
             double *rwork,int *lrw, int *k, int *jout, double *rpar,
             double *bj, int *iwork, int *liw, int *iddebdf));

void FORTRAN(tempload,(double *xforcold, double *xforc, double *xforcact,
               int *iamforc, int *nforc, double *xloadold, double *xload,
               double *xloadact, int *iamload, int *nload, int *ibody,
               double *xbody, int *nbody,double *xbodyold,double *xbodyact, 
               double *t1old, double *t1, double *t1act, int *iamt1,
               int *nk, double *amta, int *namta, int *nam, double *ampli,
               double *time, double *reltime, double *ttime, double *dtime,
               int *ithermal, int *nmethod,
	       double *xbounold, double *xboun, double *xbounact,
	       int *iamboun, int *nboun,int *nodeboun,
               int *ndirboun,int *nodeforc,int *ndirforc,int *istep,
               int *iint,double *co,double *vold,int *itg,int *ntg,
               char *amname,int *ikboun,int *ilboun,int *nelemload,
               char *sideload));

void FORTRAN(temploadmodal,(double *amta,int *namta,int *nam,double *ampli,
         double *timemin,double *ttimemin,double *dtime,double *xbounold,
         double *xboun,double *xbounmin,int *iamboun,int *nboun,
         int *nodeboun,int *ndirboun, char *amname));

void FORTRAN(transformatrix,(double *xab, double *p, double *a));

void FORTRAN(triangucont,(int *ncont,int *ntie,char *tieset,int *nset,
          char *set,int *istartset,int *iendset,int *ialset,int *itietri,
          char *lakon,int *ipkon,int *kon,int *koncont));

void FORTRAN(uout,(double *v));

void FORTRAN(updatecont,(int *koncont,int *ncont,double *co,double *vold,
			double *cg,double *straight));

void *u_calloc(size_t num,size_t size);

void FORTRAN(writeboun,(int *nodeboun,int *ndirboun,double *xboun,
      char *typeboun,int *nboun));

void FORTRAN(writebv,(double *, int *));

void writec(char *ifield, int nfield, FILE *f1);

void FORTRAN(writeev,(double *, int *, double *, double *));

void FORTRAN(writehe,(int *));

void FORTRAN(writeevcs,(double *, int *, int *, double *, double *));

void writef(double *ifield, int nfield, FILE *f1);

void writei(int *ifield, int nfield, FILE *f1);

void FORTRAN(writeim,());

void FORTRAN(writeinput,(char *inpc,int *ipoinp,int *inp,int *nline,int *ninp));

void FORTRAN(writempc,(int *, int *, double *, char *,int *));

void FORTRAN(writepf,(double *d, double *bjr, double *bji, double *freq , 
       int *nev));

void FORTRAN(writere,());

void FORTRAN(writesummary,(int *istep, int *j, int *icutb, int *l, double *ttime,
		   double *time, double *dtime));
