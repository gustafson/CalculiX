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

#define DMEMSET(a,b,c,d) for(im=b;im<c;im++)a[im]=d

void FORTRAN(addimdnodecload,(int *nodeforc,int *i,int *imdnode, 
             int *nmdnode,double *xforc,int *ikmpc,int *ilmpc,
             int *ipompc,int *nodempc,int *nmpc,int *imddof,int *nmddof,
             int *nactdof,int *mi,int *imdmpc,int *nmdmpc,int *imdboun,
	     int *nmdboun,int *ikboun,int *nboun,int *ilboun,int *ithermal,
             int *usercload));

void FORTRAN(addimdnodedload,(int *nelemload,char *sideload,int *ipkon,
             int *kon,char *lakon,int *i,int *imdnode,int *nmdnode,
             int *ikmpc,int *ilmpc,
             int *ipompc,int *nodempc,int *nmpc,int *imddof,int *nmddof,
             int *nactdof,int *mi,int *imdmpc,int *nmdmpc,int *imdboun,
	     int *nmdboun,int *ikboun,int *nboun,int *ilboun,int *ithermal));

void FORTRAN(addizdofcload,(int *nodeforc,int *ndirforc,int *nactdof,
	     int *mi,int *izdof,int *nzdof,int *i,int *iznode,int *nznode,
	     int *nk,int *imdnode,int *nmdnode,double *xforc));

void FORTRAN(addizdofdload,(int *nelemload,char *sideload,int *ipkon,
             int *kon,char *lakon,int *nactdof,int *izdof,int *nzdof,
	     int *mi,int *i,int *iznode,int *nznode,int *nk,
             int *imdnode,int *nmdnode));

void add_rect(double *au_1,int * irow_1,int * jq_1,int n_1,int m_1,
	       double *au_2,int * irow_2,int * jq_2,int n_2,int m_2,
	       double **au_rp,int **irow_rp,int * jq_r);

void FORTRAN(allocation,(int *nload_,int *nforc_,int *nboun_,
             int *nk_,int *ne_,int *nmpc_,int *nset_,int *nalset_,
	     int *nmat_,int *ntmat_,int *npmat_,int *norien_,int *nam_,
             int *nprint_,int *mi,int *ntrans_,
             char *set,int *meminset,
             int *rmeminset,int *ncs_,int *namtot_,int *ncmat_,
             int *memmpc_,int *ne1d,int *ne2d,int *nflow,
             char *jobnamec,int *irstrt,int *ithermal,int *nener,
             int *nstate_,int *istep,char *inpc,
             int *ipoinp,int *inp,int *ntie_,int *nbody_,
	     int *nprop_,int *ipoinpc,int *nevdamp,int *npt_,
	     int *nslavsm,int *nkon_));

void FORTRAN(allocont,(int *ncont,int *ntie,char *tieset,int *nset,
             char *set,int *istartset,int *iendset,int *ialset,
	     char *lakon,int *ncone,double *tietol,int *ismallsliding,
	     char *kind1,char *kind2,int *mortar));

void FORTRAN(applyboun,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon));

void FORTRAN(applyboundum,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *compressible));

void FORTRAN(applybounk,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *iponoel,double *vold,int *ipompc,int *nodempc,
       double *coefmpc,int *nmpc,int *nfreestream,int *ifreestream,
       int *nsolidsurf,int *isolidsurf,double *xsolidsurf,int *inoel,
       double *physcon,int *compressible,int *ielmat,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *voldtu,int *ntmat_,
       char *lampc,int *inomat,int *mi,int *ithermal));

void FORTRAN(applybounmpc,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *compressible,int *ismooth,int *nmpc,int *nodempc,
       int *ipompc,double *coefmpc,int *inomat,int *mi,char *matname));

void FORTRAN(applybounp,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *ipompc,int *nodempc,double *coefmpc,int *nmpc,
       int *inomat,int *mi));

void FORTRAN(applybounpgas,(int *nodeboun,int *ndirboun,int *nboun,
          double *xbounact,int *iponoel,double *vold,int *ipompc,
	  int *nodempc,double *coefmpc,int *nmpc,int *inomat,
          char *matname,int *nshcon,double *shcon,int *nrhcon,
	  double *rhcon,double *physcon,int *ntmat_,double *voldaux,
          int *mi));

void FORTRAN(applybount,(int *nodeboun,int *ndirboun,int *nboun,
          double *xbounact,int *iponoel,double *vold,int *ipompc,
	  int *nodempc,double *coefmpc,int *nmpc,int *inomat,
          char *matname,int *nshcon,double *shcon,int *nrhcon,double *rhcon,
	  double *physcon,int *compressible,int *ntmat_,double *voldaux,
	  int *mi,int *ithermal));

void FORTRAN(applybounv,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *voldtu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *compressible,int *ismooth,int *nmpc,int *nodempc,
       int *ipompc,double *coefmpc,int *inomat,int *mi));

void arpack(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     double *ad,double *au,double *b,int *nactdof, 
	     int *icol,int *jq,int *irow,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *shcon,int *nshcon,double *cocon,int *ncocon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old,
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs,   
	     int *kode,double *adb,double *aub,int *mei,double *fei,
	     char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *mi,
             int *ncmat_,int *nstate_,double *ener,char *jobnamec,
             char *output,char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	    char *cbody,int *ibody,double *xbody,int *nbody,double *thicke);

void arpackbu(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     double *ad,double *au,double *b,int *nactdof, 
	     int *icol,int *jq,int *irow,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs, 
	     int *kode,double *adb,double *aub,int *mei,double *fei,
             char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *mi,
             int *ncmat_,int *nstate_,double *ener,char *output, 
             char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     char *cbody,int *ibody,double *xbody,int *nbody,double *thicke);

void arpackcs(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     double *ad,double *au,double *b,int *nactdof, 
	     int *icol,int *jq,int *irow,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old,
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs,  
	     int *kode,double *adb,double *aub,int *mei,double *fei,
	     char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *mi,
             int *ics,double *cs,int *mpcend,int *ncmat_,int *nstate_,
             int *mcs,int *nkon,double *ener,char *jobnamec,
             char *output,char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     char *cbody,int *ibody,double *xbody,int *nbody, 
             int *nevtot,double *thicke);

void bdfill(int **irowbdp,int *jqbd,
        double **bdup,double *bdd,int *nzs_,int *ntie,int *ipkon,int *kon, 
        char *lakon,int *nslavnode,int *nmastnode,int *imastnode,
        int *islavnode,int *islavsurf,int *imastsurf,double *pmastsurf, 
	int *itiefac,int *neq,int *nactdof,double *co,double *vold,
	int *iponoels,int *inoels,int *mi,double *gapmints,double *gap,
        double *pslavsurf,double* pslavdual);

void FORTRAN(bodyforce,(char *cbody,int *ibody,int *ipobody,int *nbody,
             char *set,int *istartset,int *iendset,int *ialset,
             int *inewton,int *nset,int *ifreebody,int *k));

void FORTRAN(calcmach,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,double *factor));

void calcresidual(int *nmethod,int *neq,double *b,double *fext,double *f,
        int *iexpl,int *nactdof,double *aux1,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,int *nk,double *adb,
        double *aub,int *icol,int *irow,int *nzl,double *alpha,
	double *fextini,double *fini,int *islavnode,int *nslavnode,
        int *imastnode,int *nmastnode,int *mortar,int *ntie,
	double *f_cm,double *f_cs,int *mi);

void calcshapef(int *nvar_,int *ipvar,double **var,int *ne,
	     char *lakon,double *co,int *ipkon,int *kon,
             int *nelemface,char *sideface,int *nface,
             int *nvarf_,int *ipvarf,double **varfp);

void FORTRAN(calcview,(char *sideload,double *vold,double *co,
             double *pmid,double *e1,double *e2,double *e3,
	     int *kontri,int *nloadtr,double *adview,double *auview,
             double *dist,int *idist,double *area,int *ntrit,int *mi,int *jqrad,
	     int *irowrad,int *nzsrad,double *sidemean,int *ntria, 
             int *ntrib));

void *calcviewmt(void *i);

void FORTRAN(calinput,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	       int *nkon,int *ne,
	       int *nodeboun,int *ndirboun,double *xboun,int *nboun,
	       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,
	       int *nmpc_,int *nodeforc,int *ndirforc,double *xforc,
	       int *nforc,int *nforc_,int *nelemload,char *sideload,
	       double *xload,int *nload,int *nload_,
	       int *nprint,char *prlab,char *prset,int *mpcfree,int *nboun_,
	       int *mei,char *set,int *istartset, 
	       int *iendset,int *ialset,int *nset,int *nalset,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,double *t0,
	       double *t1,char *matname,
	       int *ielmat,char *orname,double *orab,int *ielorien,
	       char *amname,double *amta,int *namta,int *nam,
	       int *nmethod,int *iamforc,int *iamload,
	       int *iamt1,int *ithermal,int *iperturb, 
	       int *istat,int *istep,int *nmat,
	       int *ntmat_,int *norien,double *prestr,int *iprestr,
	       int *isolver,double *fei,double *veold,double *tinc,
	       double *tper,double *xmodal,
               char *filab,int *jout,int *nlabel,
	       int *idrct,int *jmax,double *tmin,double *tmax,
	       int *iexpl,double *alpha,int *iamboun,
	       double *plicon,int *nplicon,double *plkcon,int *nplkcon,
	       int *iplas,int *npmat_,int *mi,int *nk_,
	       double *trab,int *inotr,int *ntrans,int *ikboun,
               int *ilboun,int *ikmpc,int *ilmpc,int *ics,
	       double *dcs,int *ncs_,int *namtot_,double *cs,
               int *nstate_,int *ncmat_,int *iumat,int *mcs,
               char *labmpc,int *iponor,double *xnor,int *knor,
	       double *thickn,double *thicke,int *ikforc,int *ilforc,
               double *offset,int *iponoel,int *inoel,int *rig,
               int *infree,int *nshcon,double *shcon,double *cocon,
               int *ncocon,double *physcon,int *nflow,double *ctrl,
               int *maxlenmpc,int *ne1d,int *ne2d,int *nener, 
               double *vold,int *nodebounold,
               int *ndirbounold,double *xbounold,double *xforcold,
               double *xloadold,double *t1old,double *eme,
               double *sti,double *ener,
               double *xstate,char *jobnamec,int *nnn,int *irstrt,
               double *ttime,double *qaold,
               char *output,char *typeboun,char *inpc,int *nline,
               int *ipoinp,int *inp,char *tieset,double *tietol, 
               int *ntie,double *fmpc,char *cbody,int *ibody,double *xbody,
               int *nbody,int *nbody_,double *xbodyold,int *nam_,
               int *ielprop,int *nprop,int *nprop_,double *prop,
	       int *itpamp,int *iviewfile,int *ipoinpc,int *cfd,
	       int *nslavs,double *t0g,double *t1g));    

void cascade(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
   int *mpcfree,int *nodeboun,int *ndirboun,int*nboun,int*ikmpc,
   int *ilmpc,int *ikboun,int *ilboun,int *mpcend,int *mpcmult,
   char *labmpc,int *nk,int *memmpc_,int *icascade,int *maxlenmpc,
   int *callfrommain,int *iperturb,int *ithermal);

void FORTRAN(cfdconv,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co, 
	   double *factor,double *voldconini,double *dtimef));

int cgsolver(double *A,double *x,double *b,int neq,int len,int *ia,int *iz, 
				double *eps,int *niter,int precFlg);

void checkconvergence(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	  int *ne,double *stn,int *nmethod, 
	  int *kode,char *filab,double *een,double *t1act,
          double *time,double *epn,int *ielmat,char *matname,
          double *enern,double *xstaten,int *nstate_,int *istep,
          int *iinc,int *iperturb,double *ener,int *mi,char *output,
          int *ithermal,double *qfn,int *mode,int *noddiam,double *trab,
          int *inotr,int *ntrans,double *orab,int *ielorien,int *norien,
          char *description,double *sti,
	  int *icutb,int *iit,double *dtime,double *qa,double *vold,
          double *qam,double *ram1,double *ram2,double *ram,
          double *cam,double *uam,int *ntg,double *ttime,
          int *icntrl,double *theta,double *dtheta,double *veold,
          double *vini,int *idrct,double *tper,int *istab,double *tmax, 
	  int *nactdof,double *b,double *tmin,double *ctrl,double *amta,
          int *namta,int *itpamp,int *inext,double *dthetaref,int *itp,
          int *jprint,int *jout,int *uncoupled,double *t1,int *iitterm,
          int *nelemload,int *nload,int *nodeboun,int *nboun,int *itg,
	  int *ndirboun,double *deltmx,int *iflagact,char *set,int *nset,
	  int *istartset,int *iendset,int *ialset,double *emn,double *thicke);

void checkconvnet(int *icutb,int *iin,
		  double *qamt,double *qamf,double *qamp, 
		  double *ram1t,double *ram1f,double *ram1p,
		  double *ram2t,double *ram2f,double *ram2p,
		  double *ramt,double *ramf,double *ramp,
		  int *icntrl,double *dtheta,double *ctrl,
                  double *uama,double *ram1a,double *ram2a,double *rama,
                  double *vamt,double *vamf,double *vamp,double *vama);

void checkinclength(double *time,double *ttime,double *theta,double *dtheta,
          int *idrct,double *tper,double *tmax,double *tmin,double *ctrl, 
          double *amta,int *namta,int *itpamp,int *inext,double *dthetaref, 
	  int *itp,int *jprint,int *jout);

void FORTRAN(checktime,(int *itpamp,int *namta,double *tinc,double *ttime,
             double *amta,double *tmin,int *inext,int *itp));

void FORTRAN(closefile,());

  
void FORTRAN(compdt,(int *nk,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,
       int *ntmat_,int *iponoel,int *inoel,double *dtimef,int *iexplicit,
       int *ielmat,double *physcon,double *dh,double *cocon,int *ncocon,
       int *ithermal,int *mi,int *ipkon,int *kon,char *lakon,
       double *dtl,int *ne,double *v,double *co, int *turbulent, 
       double *voldtu));

void compfluid(double *co,int *nk,int *ipkon,int *kon,char *lakon,
    int *ne,int *ipoface,char *sideface,int *ifreestream, 
    int *nfreestream,int *isolidsurf,int *neighsolidsurf,
    int *nsolidsurf,int *iponoel,int *inoel,int *nshcon,double *shcon,
    int *nrhcon,double *rhcon,double *vold,int *ntmat_,int *nodeboun, 
    int *ndirboun,int *nboun,int *ipompc,int *nodempc,int *nmpc,
    int *ikmpc,int *ilmpc,int *ithermal,int *ikboun,int *ilboun,
    int *turbulent,int *isolver,int *iexpl,double *voldtu,double *ttime,
    double *time,double *dtime,int *nodeforc,int *ndirforc,double *xforc,
    int *nforc,int *nelemload,char *sideload,double *xload,int *nload,
    double *xbody,int *ipobody,int *nbody,int *ielmat,char *matname,
    int *mi,int *ncmat_,double *physcon,int *istep,int *iinc,
    int *ibody,double *xloadold,double *xboun,
    double *coefmpc,int *nmethod,double *xforcold,double *xforcact,
    int *iamforc,int *iamload,double *xbodyold,double *xbodyact,
    double *t1old,double *t1,double *t1act,int *iamt1,double *amta,
    int *namta,int *nam,double *ampli,double *xbounold,double *xbounact,
    int *iamboun,int *itg,int *ntg,char *amname,double *t0,int *nelemface,
    int *nface,double *cocon,int *ncocon,double *xloadact,double *tper,
    int *jmax,int *jout,char *set,int *nset,int *istartset,
    int *iendset,int *ialset,char *prset,char *prlab,int *nprint,
    double *trab,int *inotr,int *ntrans,char *filab,char *labmpc,
    double *sti,int *norien,double *orab);

void complexfreq(double **cop,int *nk,int **konp,int **ipkonp,char **lakonp,int *ne, 
	       int **nodebounp,int **ndirbounp,double **xbounp,int *nboun,
	       int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,
               int *nmpc,int *nodeforc,int *ndirforc,double *xforc, 
               int *nforc,int *nelemload,char *sideload,double *xload,
	       int *nload, 
	       int **nactdofp,int *neq,int *nzl,int *icol,int *irow, 
	       int *nmethod,int **ikmpcp,int **ilmpcp,int **ikbounp, 
	       int **ilbounp,double *elcon,int *nelcon,double *rhcon, 
	       int *nrhcon,double *cocon,int *ncocon,
               double *alcon,int *nalcon,double *alzero, 
               int **ielmatp,int **ielorienp,int *norien,double *orab, 
               int *ntmat_,double **t0p, 
	       double **t1p,int *ithermal,double *prestr,int *iprestr, 
	       double **voldp,int *iperturb,double **stip,int *nzs, 
	       double *tinc,double *tper,double *xmodal,
	       double **veoldp,char *amname,double *amta,
	       int *namta,int *nam,int *iamforc,int *iamload,
	       int **iamt1p,int *jout,
	       int *kode,char *filab,double **emep,double *xforcold, 
	       double *xloadold,
               double **t1oldp,int **iambounp,double **xbounoldp,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstate,int *npmat_,char *matname,int *mi,
               int *ncmat_,int *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,int *nset,int *istartset,
               int *iendset,int **ialsetp,int *nprint,char *prlab,
               char *prset,int *nener,double *trab, 
               int **inotrp,int *ntrans,double **fmpcp,char *cbody,int *ibody,
               double *xbody,int *nbody,double *xbodyold,int *istep,
               int *isolver,int *jq,char *output,int *mcs,int *nkon,
               int *mpcend,int *ics,double *cs,int *ntie,char *tieset,
               int *idrct,int *jmax,double *tmin,double *tmax,
	       double *ctrl,int *itpamp,double *tietol,int *nalset,
	       int **nnnp,int *ikforc,int *ilforc,double *thicke);

void FORTRAN(con2phys,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void contact(int *ncont,int *ntie,char *tieset,int *nset,char *set,
	     int *istartset,int *iendset,int *ialset,int *itietri,
	     char *lakon,int *ipkon,int *kon,int *koncont,int *ne,
	     double *cg,double *straight,int *ifree,double *co,
	     double *vold,int *ielmat,double *cs,double *elcon,
             int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
             int *ne0,double *vini,
             int *nmethod,int *nmpc,int *mpcfree,int *memmpc_,
             int **ipompcp,char **labmpcp,int **ikmpcp,int **ilmpcp,
             double **fmpcp,int **nodempcp,double **coefmpcp,
             int *iperturb,int *ikboun,int *nboun,int *mi,int *imastop,
             int *nslavnode,int *islavnode,int *islavsurf,int *itiefac,
             double *areaslav,int *iponoels,int *inoels,double *springarea,
             double *tietol,double *reltime,int *imastnode,int *nmastnode,
             double *xmastnor,double *xnormastface,char*filab);

void contactmortar(int *ncont,int *ntie,char *tieset,int *nset,char *set,
        int *istartset,int *iendset,int *ialset,int *itietri,
        char *lakon,int *ipkon,int *kon,int *koncont,int *ne,
        double *cg,double *straight,double *co,
        double *vold,int *ielmat,double *cs,double *elcon,
        int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
        int *ne0,double *vini,
        int *nmethod,int *neq,int *nzs,int *nactdof,int *itiefac,
        int *islavsurf,int *islavnode,int *imastnode,
        int *nslavnode,int *nmastnode,int *ncone,double *ad,
        double **aup,double *b,int **irowp,int *icol,int *jq,int *imastop,
        int *iponoels,int *inoels,int *nzsc,double **aucp,
        double *adc,int **irowcp,int *jqc,int *islavact,
        double *gap,double *bdd,double **auqdtp,int **irowqdtp,
	int *jqqdt,int *nzsqdt,int *nzl,double *slavnor,double *bhat,
	int *icolc,double **aubdp,int **irowbdp,int *jqbd,int *mi,
	int *ipe,int *ime,double *tietol);

void contactstress(double *bhat,double *adc,double *auc,int *jqc, 
    int *irowc,int *neq,double *gap,double *bdd,double *b,int *islavact,
    double *auqdt,int *irowqdt,int *jqqdt,int *ntie,int *nslavnode,
    int *islavnode,double *slavnor,int *icol,int *nzl,
    int *nactdof,int* flagact,double* cstress,int *mi,double *cdisp,
    double *alambda,int *iit);

void FORTRAN(coriolissolve,(double *cc,int *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
             int *iter,double *d,double *temp));

void FORTRAN(cprint,(char *text,int *before,int *after));

void FORTRAN(createbdentry,(int *itie,int *ipkon,int *kon,char *lakon, 
     int *nodem,int *nodes,int *islavsurf,int *imastsurf,
     double *pmastsurf,int *itiefac,double *contribution,double *co,
     double *vold,int *iponoels,int *inoels,int *mi,double* pslavsurf,
     double* pslavdual));

void FORTRAN(createddentry,(int *itie,int *ipkon,int *kon,int *nodes,
       char *lakon,int *islavsurf,int *itiefac,double *contribution,
       double *co,double *vold,int *iponoels,int *inoels,int *mi, 
       double* pslavdual));

void FORTRAN(creategap,(int *itie,int *ipkon,int *kon,char *lakon,int *nodes,
       int *islavsurf,int *itiefac,
       double *co,double *vold,int *iponoels,int *inoels,int *mi, 
       double *pslavsurf,double* pslavdual,double *gapmints,double *gapcont));

void FORTRAN(createinum,(int *ipkon,int *inum,int *kon,char *lakon,int *nk,
             int *ne,char *cflag,int *nelemload,int *nload,int *nodeboun,
             int *nboun,int *ndirboun,int *ithermal));

void FORTRAN(createmddof,(int *imddof,int *nmddof,int *istartset,
       int *iendset,int *ialset,int *nactdof,int *ithermal,int *mi,
       int *imdnode,int *nmdnode,int *ikmpc,int *ilmpc,int *ipompc,
       int *nodempc,int *nmpc,int *imdmpc,
       int *nmdmpc,int *imdboun,int *nmdboun,int *ikboun,int *nboun,
       int *nset,int *ntie,char *tieset,char *set,char *lakon,int *kon,
       int *ipkon,char *labmpc,int *ilboun,char *filab,char *prlab,
       char *prset,int *nprint,int *ne));

void FORTRAN(createmdelem,(int *imdnode,int *nmdnode,double *xforc,
             int *ikmpc,int *ilmpc,int *ipompc,int *nodempc,int *nmpc,
             int *imddof,int *nmddof,int *nactdof,int *mi,int *imdmpc,
             int *nmdmpc,int *imdboun,int *nmdboun,int *ikboun,int *nboun,
             int *ilboun,int *ithermal,int *imdelem,int *nmdelem,
             int *iponoel,int *inoel,char *prlab,char *prset,int *nprint,
             char *lakon,char *set,int *nset,int *ialset,int *ipkon,
             int *kon,int *istartset,int *iendset,int *nforc,
             int *ikforc,int *ilforc));

void FORTRAN(datri,(double *au,double *aa,double *ad,int *icol,int *neq, 
	    int *flg));

void FORTRAN(dasol,(double *au,double *aa,double *ad,double *b,int *icol, 
	    int *neq,double *energy));

void dfdbj(double *bcont,double **dbcontp,int *neq,int *nope,
	   int *konl,int *nactdof,double *s,double *z,int *ikmpc,
	   int *ilmpc,int *ipompc,int *nodempc,int *nmpc,
	   double *coefmpc,double *fnl,int *nev,
	   int **ikactcontp,int **ilactcontp,int *nactcont,int *nactcont_,
           int *mi,int *cyclicsymmetry,int *izdof,int *nzdof);
      
void FORTRAN(dgesv, (int *nteq,int *nhrs,double *ac,int *lda,int *ipiv,
                     double *bc,int *ldb,int *info)); 

void FORTRAN(dgetrs, (char *trans,int *nteq,int *nrhs,double *ac,int *lda,
		      int *ipiv,double *bc,int *ldb,int *info));

void FORTRAN(drfftf,(int *ndata,double *r,double *wsave,int *isave));

void FORTRAN(drffti,(int *ndata,double *wsave,int *isave));

void FORTRAN(dsaupd,(int *ido,char *bmat,int *n,char *which,int *nev,
	     double *tol,double *resid,int *ncv,double *z,int *ldz,
	     int *iparam,int *ipntr,double *workd,double *workl,
	     int *lworkl,int *info));

void FORTRAN(dseupd,(int *,char *,int *,double *,double *,
	     int *,double *,char *,int *,char *, 
	     int *,double *,double *,int *,double *,
	     int *,int *,int *,double *,
	     double *,int *,int *));

void FORTRAN(dsort,(double *dx,int *iy,int *n,int *kflag));

void FORTRAN(dspgv,(int *itype,char *jobz,char *uplo,int *n,double *ap,
                    double *bp,double *w,double *z,int *ldz,double *work,
                    int *info));

void dyna(double **cop,int *nk,int **konp,int **ipkonp,char **lakonp,int *ne, 
	       int **nodebounp,int **ndirbounp,double **xbounp,int *nboun,
	       int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,
               int *nmpc,int *nodeforc,int *ndirforc,double *xforc, 
               int *nforc,int *nelemload,char *sideload,double *xload,
	       int *nload, 
	       int **nactdofp,int *neq,int *nzl,int *icol,int *irow, 
	       int *nmethod,int **ikmpcp,int **ilmpcp,int **ikbounp, 
	       int **ilbounp,double *elcon,int *nelcon,double *rhcon, 
	       int *nrhcon,double *cocon,int *ncocon,
               double *alcon,int *nalcon,double *alzero, 
               int **ielmatp,int **ielorienp,int *norien,double *orab, 
               int *ntmat_,double **t0p, 
	       double **t1p,int *ithermal,double *prestr,int *iprestr, 
	       double **voldp,int *iperturb,double **stip,int *nzs, 
	       double *tinc,double *tper,double *xmodal,
	       double **veoldp,char *amname,double *amta,
	       int *namta,int *nam,int *iamforc,int *iamload,
	       int **iamt1p,int *jout,
	       int *kode,char *filab,double **emep,double *xforcold, 
	       double *xloadold,
               double **t1oldp,int **iambounp,double **xbounoldp,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstate,int *npmat_,char *matname,int *mi,
               int *ncmat_,int *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,int *nset,int *istartset,
               int *iendset,int **ialsetp,int *nprint,char *prlab,
               char *prset,int *nener,double *trab, 
               int **inotrp,int *ntrans,double **fmpcp,char *cbody,int *ibody,
               double *xbody,int *nbody,double *xbodyold,int *istep,
               int *isolver,int *jq,char *output,int *mcs,int *nkon,
               int *mpcend,int *ics,double *cs,int *ntie,char *tieset,
               int *idrct,int *jmax,double *tmin,double *tmax,
	       double *ctrl,int *itpamp,double *tietol,int *nalset,
       	       int **nnnp,int *ikforc,int *ilforc,double *thicke);

void dynacont(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne, 
	      int *nodeboun,int *ndirboun,double *xboun,int *nboun,
	      int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
	      int *nmpc,int *nodeforc,int *ndirforc,double *xforc, 
	      int *nforc,int *nelemload,char *sideload,double *xload,
	      int *nload, 
	      int *nactdof,int *neq,int *nzl,int *icol,int *irow, 
	      int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	      int *ilboun,double *elcon,int *nelcon,double *rhcon, 
	      int *nrhcon,double *cocon,int *ncocon,
	      double *alcon,int *nalcon,double *alzero, 
	      int *ielmat,int *ielorien,int *norien,double *orab, 
	      int *ntmat_,double *t0, 
	      double *t1,int *ithermal,double *prestr,int *iprestr, 
	      double *vold,int *iperturb,double *sti,int *nzs, 
	      double *tinc,double *tper,double *xmodal,
	      double *veold,char *amname,double *amta,
	      int *namta,int *nam,int *iamforc,int *iamload,
	      int *iamt1,int *jout,char *filab,double *eme,double *xforcold, 
	      double *xloadold,
	      double *t1old,int *iamboun,double *xbounold,int *iexpl,
	      double *plicon,int *nplicon,double *plkcon,int *nplkcon,
	      double *xstate,int *npmat_,char *matname,int *mi,
	      int *ncmat_,int *nstate_,double *ener,char *jobnamec,
	      double *ttime,char *set,int *nset,int *istartset,
	      int *iendset,int *ialset,int *nprint,char *prlab,
	      char *prset,int *nener,double *trab, 
	      int *inotr,int *ntrans,double *fmpc,char *cbody,int *ibody,
	      double *xbody,int *nbody,double *xbodyold,int *istep,
	      int *isolver,int *jq,char *output,int *mcs,int *nkon,
	      int *mpcend,int *ics,double *cs,int *ntie,char *tieset,
	      int *idrct,int *jmax,double *tmin,double *tmax,
	      double *ctrl,int *itpamp,double *tietol,int *iit,
	      int *ncont,int *ne0,double *reltime,double *dtime,
	      double *bcontini,double *bj,double *aux,int *iaux,
	      double *bcont,int *nev,double *v,
              int *nkon0,double *deltmx,double *dtheta,double *theta,
              int *iprescribedboundary,int *mpcfree,int *memmpc_,
              int *itietri,int *koncont,double *cg,double *straight,
              int *iinc,double *vini,
              double *aa,double *bb,double *aanew,double *d,double *z, 
	      double *zeta,double *b,double *time0,double *time1, 
	      int *ipobody,
              double *xforcact,double *xloadact,double *t1act, 
              double *xbounact,double *xbodyact,double *cd,double *cv,
              double *ampli,double *dthetaref,double *bjp,double *bp,
              double *cstr,int *imddof,
              int *nmddof,int **ikactcontp,int *nactcont,int *nactcont_,
              double *aamech,double *bprev,int *iprev,int *inonlinmpc,
              int **ikactmechp,int *nactmech,int *imdnode,int *nmdnode,
              int *imdboun,int *nmdboun,int *imdmpc,int *nmdmpc,
              int *itp,int *inext,int *ifricdamp,double *aafric,
              double *bfric,int *imastop,int *nslavnode,int *islavnode,
              int *islavsurf,
              int *itiefac,double *areaslav,int *iponoels,int *inoels,
              double *springarea,int *izdof,int *nzdof,double *fn,
	      int *imastnode,int *nmastnode,double *xmastnor,
              double *xnormastface);
 
void dynboun(double *amta,int *namta,int *nam,double *ampli,double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,int *iamboun,int *nboun,int *nodeboun,
             int *ndirboun,double *ad,double *au,double *adb,
             double *aub,int *icol,int *irow,int *neq,int *nzs,
             double *sigma,double *b,int *isolver,
             double *alpham,double *betam,int *nzl,
             int *init,double *bact,double *bmin,int *jq,char *amname,
             double *bv,double *bprev,double *bdiff,
             int *nactmech,int *icorrect,int *iprev);

void FORTRAN(dynresults,(int *nk,double *v,int *ithermal,int *nactdof,
             double *vold,int *nodeboun,int *ndirboun,double *xboun,
             int *nboun,int *ipompc,int *nodempc,double *coefmpc,
	     char *labmpc,int *nmpc,double *b,double *bp,double *veold,
	     double *dtime,int *mi,int *imdnode,int *nmdnode,int *imdboun,
	     int *nmdboun,int *imdmpc,int *nmdmpc,int *nmethod));

void FORTRAN(dynsolv,(double *b,double *z,double *d,double *zeta,int *nev,
	      int *neq,double *tinc,int *jinc,int *jout,double *vold,
	      double *xforcold,int *nodeforc,int *ndirforc,
	      double *xforc,int *iamforc,int *nforc,double *xloadold,
	      double *xload,int *iamload,int *nelemload,char *sideload,
	      int *nload,
	      double *t1old,
	      double *t1,int *iamt1,int *nk,double *amta,int *namta, 
	      int *nam,double *ampli,double *aa,double *bb,double *bj, 
	      double *v,int *nodeboun,int *ndirboun,double *xboun, 
	      int *nboun,int *ipompc,int *nodempc,double *coefmpc, 
              char *labmpc,int *nmpc,int *nactdof, 
	      int *iperturb,int *nmethod,double *co,int *kon,
	      int *ipkon,char *lakon,int *ne,double *stn,double *stx,
	      double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	      double *alcon,int *nalcon,double *alzero,int *ielmat,
	      int *ielorien,int *norien,double *orab,int *ntmat_,
	      double *t0,int *ithermal,int *kode,double *cv,
	      double *cd,int *inum,double *prestr,int *iprestr,
	      int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
	      char *filab,double *eme,double *een,
	      double *sti,double *f,double *fn,double *xforcact,
	      double *xloadact,
              double *t1act,double *xbounold,double *xbounact,
              int *iamboun,int *iexpl,double *plicon,int *nplicon,
              double *plkcon,
	      int *nplkcon,double *xstateini,double *xstiff,
              double *xstate,int *npmat_,double *epn,char *matname,
              int *mi,int *ncmat_,int *nstate_,double *stiini,
              double *vini,double *ener,double *enern,double *xstaten,
              double *ttime,double *eei,double *enerini,double *cocon,
              int *ncocon,char *set,int *nset,int *istartset,
              int *iendset,int *ialset,int *nprint,char *prlab,
              char *prset,double *qfx,double *qfn,double *trab,
	      int *inotr,int *ntrans,double *fmpc,double *veold,
	      char *cbody,int *ibody,double *xbody,int *nbody, 
              double *xbodyold,double *xbodyact,int *ipobody,
              double *cgr,double *xmodal,double *au,
              double *aub,double *vbounact,double *abounact,int *nzs));

void FORTRAN(elementpernode,(int *iponoel,int *inoel,char *lakon,int *ipkon,
              int *kon,int *ne));

void FORTRAN(envtemp,(int *itg,int *ieg,int *ntg,int *ntr,char *sideload,
                      int *nelemload,int *ipkon,int *kon,char *lakon,
                      int *ielmat,int *ne,int *nload,
                      int *kontri,int *ntri,int *nloadtr,
                      int *nflow,int *ndirboun,int *nactdog,
                      int *nodeboun,int *nacteq,
                      int *nboun,int *ielprop,double *prop,int *nteq,
                      double *v,int *network,double *physcon,
		      double *shcon,int *ntmat_,double *co,
                      double *vold,char *set,int *nshcon,
		      double *rhcon,int *nrhcon,int *mi,int *nmpc,
                      int *nodempc,int *ipompc,char *labmpc,int *ikboun));

void FORTRAN(equationcheck,(double *ac,int *nteq,int *nactdog,
                            int *itg,int *ntg,int *nacteq,int *network));

void expand(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc,int *nodeforc,int *ndirforc,double *xforc, 
             int *nforc,int *nelemload,char *sideload,double *xload,
             int *nload,int *nactdof,int *neq, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs,  
	     double *adb,double *aub,char *filab,double *eme,
             double *plicon,int *nplicon,double *plkcon,int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *mi,
	     int *ics,double *cs,int *mpcend,int *ncmat_,
             int *nstate_,int *mcs,int *nkon,double *ener,
             char *jobnamec,char *output,char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     int *nev,double **z,int *iamboun,double *xbounold,
             int *nsectors,int *nm,int *icol,int *irow,int *nzl,int *nam,
             int *ipompcold,int *nodempcold,double *coefmpcold,
             char *labmpcold,int *nmpcold,double *xloadold,int *iamload,
             double *t1old,double *t1,int *iamt1,double *xstiff,int **icolep,
	     int **jqep,int **irowep,int *isolver,
	     int *nzse,double **adbep,double **aubep,int *iexpl,int *ibody,
	     double *xbody,int *nbody,double *cocon,int *ncocon,
	     char* tieset,int* ntie,int **nnnp,int *imddof,int *nmddof,
	     int *imdnode,int *nmdnode,int *imdboun,int *nmdboun,
             int *imdmpc,int *nmdmpc,int **izdofp,int *nzdof);

void FORTRAN(extrapolate,(double *sti,double *stn,int *ipkon,int *inum,
             int *kon,char *lakon,int *nfield,int *nk,int *ne,int *mi,
             int *ndim,double *orab,int *ielorien,double *co,int *iorienglob,
	     char *cflag,int *nelemload,int *nload,int *nodeboun,int *nboun,
             int *ndirboun,double *vold,int *ithermal,int *force,
	     int *cfd,int *ielmat,double *thicke));

void FORTRAN(fcrit,(double *time,double *tend,double *aai,double *bbi,
		      double *zetaj,double *dj,double *ddj,
		      double *h1,double *h2,double *h3,double *h4,
                      double *func,double *funcp));

void FORTRAN(findsurface,(int *ipoface,int *nodface,int *ne,int *ipkon,int *kon,
                     char *lakon,int *ntie,char *tieset));

void FORTRAN (flowoutput,(int *itg,int *ieg,int *ntg,int *nteq,
			  double *bc,char *lakon,
			  int *ntmat_,double *v,double *shcon,int *nshcon,
			  int *ipkon,int *kon,double *co,int *nflow,
			  double *dtime,double *ttime,double *time,
			  int *ielmat,double *prop,
			  int *ielprop,int *nactdog,int *nacteq,int *iin,
			  double *physcon,double *camt,double *camf,double *camp,
			  double *uamt,double *uamf,double *uamp,
			  double *rhcon,int *nrhcon,
			  double *vold,char *jobnamef,char *set,int *istartset,
                          int *iendset,int *ialset,int *nset,int *mi));

void FORTRAN(flowresult,(int *ntg,int *itg,double *cam,double *vold,
              double *v,
              int *nload,char *sideload,int *nelemload,
	      double *xloadact,int *nactdog,int *network,int *mi));

void FORTRAN(forcesolve,(double *zc,int *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
	     int *iter,double *d,int *neq,double *z,int *istartnmd,
	     int *iendnmd,int *nmd,int *cyclicsymmetry,int *neqact));

void frdcyc(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne,double *v,
	    double *stn,int *inum,int *nmethod,int *kode,char *filab,
	    double *een,double *t1,double *fn,double *time,double *epn,
	    int *ielmat,char *matname,double *cs,int *mcs,int *nkon,
	    double *enern,double *xstaten,int *nstate_,int *istep,
            int *iinc,int *iperturb,double *ener,int *mi,char *output,
            int *ithermal,double *qfn,int *ialset,int *istartset,
            int *iendset,double *trab,int *inotr,int *ntrans,double *orab,
	    int *ielorien,int *norien,double *sti,double *veold,int *noddiam,
            char *set,int *nset,double *emn, double *thicke);

void FORTRAN(frdfluid,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *vold,int *kode,double *time,
             int *ielmat,char *matname,int *nnstep,double *vtu,
	     double *voldtu,double *voldaux,double *physcon,char *filab,
	     int *inomat,int *ntrans,int *inotr,double *trab,int *mi,
             double *stn));

void FORTRAN(frditeration,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *time,int *ielmat,char *matname,
             int *mi,int *istep,int *iinc));

void FORTRAN(frdphase,(int *kode,double *time,int *nk,int *inum,
             double *vr,double *vi,double *stnr,double *stni,
             char *filab,int *mode,int *noddiam,int *nmethod,
	     double *vmax,double *stnmax,int *nkcoords));

void FORTRAN(fridaforc,(double *xl,int *konl,double *vl,int *imat,
       double *elcon,int *nelcon,double *elas,double *fnl,int *ncmat_,
       int *ntmat_,int *nope,char *lakonl,double *t0l,double *t1l,
       int *kode,double *elconloc,double *plicon,int *nplicon,int *npmat_,
       double *veoldl,double *senergy,int *iener,double *cstr,int *mi,
       double *springarea));

void FORTRAN(fsub,(double *time,double *tend,double *aai,double *bbi,
		   double *ddj,double *h1,double *h2,double *h3,double *h4,
                   double *func,double *funcp));

void FORTRAN(fsuper,(double *time,double *tend,double *aai,double *bbi,
		       double *h1,double *h2,double *h3,double *h4,
		       double *h5,double *h6,double *func,double *funcp));

void FORTRAN(gasmechbc,(double *vold,int *nload,char *sideload,
			int *nelemload,double *xload,int *mi));

void FORTRAN(gencontelem,(char *tieset,int *ntie,int *itietri,int *ne,
     int *ipkon,int *kon,char *lakon,
     double *cg,double *straight,int *ifree,int *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,int *nx,int *ny,int *nz,
     int *ielmat,double *cs,double *elcon,int *istep,int *iinc,int *iit,
     int *ncmat_,int *ntmat_,int *ne0,
     double *vini,int *nmethod,int *mi,int *imastop,int *nslavnode,
     int *islavnode,int *islavsurf,int *itiefac,double *areaslav,
     int *iponoels,int *inoels,double *springarea,int *ikmpc,
     int *ilmpc,int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
     char *set,int *nset,int *istartset,int *iendset,int *ialset,
     double *tietol,double *reltime,double *xmastnor,double *xnormastface,
     int *imastnode,int *nmastnode,char* filab));

void gencontmpc(int *ne,int *ne0,char *lakon,int *ipkon,int *kon,
		int *nmpc,int **ikmpc,int **ilmpc,int **ipompc, 
                int *mpcfree,
                double **fmpc,char **labmpc,int **nodempc,int *memmpc_,
                double **coefmpc,int *nmpc_,int *ikboun,int *nboun);
       
void FORTRAN(gencontrel,(char *tieset,int *ntie,int *itietri,int *ipkon,
        int *kon,char *lakon,char *set,double *cg,double *straight,
        int *koncont,double *co,double *vold,double *xo,
        double *yo,double *zo,double *x,double *y,double *z,int *nx,
        int *ny,int *nz,int *nset,
        int *iinc,int *iit,
        int *islavsurf,int *imastsurf,double *pmastsurf,
        int *itiefac,int *islavnode,int *nslavnode,double *slavnor,
	double *slavtan,int *imastop,
	int *mi,int *ncont,int *ipe,int *ime,double *pslavsurf,
        double* pslavdual));

void FORTRAN(genfirstactif,(char *tieset,int *ntie,int *itietri,int *ne,
     int *ipkon,int *kon,char *lakon,
     double *cg,double *straight,int *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,int *nx,int *ny,int *nz,
     int *ielmat,double *cs,double *elcon,int *istep,int *iinc,int *iit,
     int *ncmat_,int *ntmat_,int *ne0,
     double *vini,int *nmethod,int *mi,int *imastop,int *nslavnode,
     int *islavnode,int *islavsurf,int *itiefac,double *areaslav,
     int *iponoels,int *inoels,
     char *set,int *nset,int *istartset,int *iendset,int *ialset,
     int *islavact,int *ifree,double *tietol));

void FORTRAN(genislavactdof,(int *ntie,int *neq,int *nactdof,
             int *nslavnode,int *islavact,int *islavactdof,
			     int *islavnode,int *mi));

void  FORTRAN(gennactdofinv,(int *nactdof,int *nactdofinv,int *nk,
       int *mi,int *nodorig,int *ipkon,char *lakon,int *kon,int *ne));

void FORTRAN(gentiedmpc,(char *tieset,int *ntie,int *itietri,
          int *ipkon,int *kon,char *lakon,char *set,int *istartset,
          int *iendset,int *ialset,double *cg,double *straight,
	  int *koncont,double *co,double *xo,double *yo,double *zo,
          double*x,double *y,double *z,int *nx,int *ny,int *nz,int *nset,
          int *ifaceslave,int *istartfield,int *iendfield,int *ifield,
          int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nmpc_,
          int *mpcfree,int *ikmpc,int *ilmpc,char *labmpc,int *ithermal,
	  double *tietol,int *icfd,int *ncont,int *imastop));

void FORTRAN(geomview,(double *vold,double *co,double *pmid,double *e1,
             double *e2,double *e3,int *kontri,double *area,double *cs,
             int *mcs,int *inocs,int *ntrit,int *nk,int *mi,double *sidemean));

void FORTRAN(getfneig,(char *fneig));

void FORTRAN(identamta,(double *amta,double *reftime,int *istart,int *iend,
               int *id));

void FORTRAN(identifytiedface,(char *tieset,int *ntie,char *set,int *nset,
           int *faceslave));

void FORTRAN(includefilename,(char *buff,char *includefn,int *lincludefn));

void inicont(int* nk,int *ncont,int *ntie,char *tieset,int *nset,char *set,
             int *istartset,int *iendset,int *ialset,int **itietrip,
	     char *lakon,int *ipkon,int *kon,int **koncontp,
             int *ncone,double *tietol,int *ismallsliding,int **itiefacp,
	     int **islavsurfp,int **islavnodep,int **imastnodep,
	     int **nslavnodep,int **nmastnodep,int *mortar,
	     int **imastopp,int *nkon,int **iponoels,int **inoelsp,
             int **ipep,int **imep,int *ne,int *ifacecount,
             int *nmpc,int *mpcfree,int *memmpc_,
             int **ipompcp,char **labmpcp,int **ikmpcp,int **ilmpcp,
             double **fmpcp,int **nodempcp,double **coefmpcp,
             int *iperturb,int *ikboun,int *nboun);
  
void FORTRAN(initialcfd,(double *yy,int *nk,double *co,int *ne,int *ipkon,
       int *kon,char *lakon,double *x,double *y,double *z,double *x0,
       double *y0,double *z0,int *nx,int *ny,int *nz,int *isolidsurf,
       int *neighsolidsurf,double *xsolidsurf,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,double *voldaux,
       int *ntmat_,int *iponoel,int *inoel,int *iexplicit,
       int *ielmat,int *nsolidsurf,int *turbulent,double *physcon,
       int *compressible,char *matname,int *inomat,double *voldtu,
       int *mi,int *euler,int *ithermal));

void FORTRAN(initialchannel,(int *itg,int *ieg,int *ntg,double *ac,double *bc,
                         char *lakon,double *v,int * ipkon,int *kon,
                         int *nflow,int *ikboun,int *nboun,double *prop,
                         int *ielprop,int *nactdog,int *ndirboun,
                         int *nodeboun,double *xbounact,int *ielmat,
                         int *ntmat_,double *shcon,int *nshcon,
                         double *physcon,int *ipiv,int *nteq,
                         double *rhcon,int *nrhcon,int *ipobody,int *ibody,
                         double *xbody,double *co,int *nbody,int *network,
                         int *iin_abs,double *vold,char *set,int *istep,
                         int *iit,int *mi,int *ineighe,int *ilboun));

void FORTRAN(initialnet,(int *itg,int *ieg,int *ntg,double *ac,double *bc,
                         char *lakon,double *v,int * ipkon,int *kon,
                         int *nflow,int *ikboun,int *nboun,double *prop,
                         int *ielprop,int *nactdog,int *ndirboun,
                         int *nodeboun,double *xbounact,int *ielmat,
                         int *ntmat_,double *shcon,int *nshcon,
                         double *physcon,int *ipiv,int *nteq,
                         double *rhcon,int *nrhcon,int *ipobody,int *ibody,
                         double *xbody,double *co,int *nbody,int *network,
                         int *iin_abs,double *vold,char *set,int *istep,
                         int *iit,int *mi,int *ineighe,int *ilboun,
                         int *channel));

void insert(int *ipointer,int **mast1p,int **mast2p,int *i1,
	    int *i2,int *ifree,int *nzs_);

void insertas(int **irowp,int **mast1p,int *i1,
      int *i2,int *ifree,int *nzs_,double *contribution,double **bdp);

void insertas_ws(int **irowp,int *i1,
      int *i2,int *ifree,int *nzs_,double *contribution,double **bdp);

void FORTRAN(isortid,(int *ix,double *dy,int *n,int *kflag));

void FORTRAN(isortii,(int *ix,int *iy,int *n,int *kflag));

void FORTRAN(isortiid,(int *ix,int *iy,double *dy1,int *n,int *kflag));

void FORTRAN(isortiddc1,(int *ix,double *dy1,double *dy2,char *cy,int *n, 
                         int *kflag));

void FORTRAN(isortiddc2,(int *ix1,int *ix2,double *dy1,double *dy2, 
                         char *cy,int *n,int *kflag));

void FORTRAN(iter,(double *coef,int *jcoef,int *ndim,int *n,int *p,
	   int *ip,double *u,double *ubar,double *rhs,double *wksp,
	   int *iwksp,int *nw,int *inw,int *iparm,double *rparm));

void FORTRAN(keystart,(int *ifreeinp,int *ipoinp,int *inp,char *name,
           int *iline,int *ikey));

void FORTRAN(lump,(double *adb,double *aub,double *adl,int *irow,int *jq,
                   int *neq));

void FORTRAN(mafillcorio,(double *co,int *nk,int *kon,int *ipkon, 
               char *lakon,
	       int *ne,int *nodeboun,int *ndirboun,double *xboun, 
	       int *nboun,int *ipompc,int *nodempc,double *coefmpc, 
	       int *nmpc,int *nodeforc,int *ndirforc,
	       double *xforc,int *nforc,int *nelemload,char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr,
	       double *ad,double *au,int *nactdof, 
	       int *icol,int *jq,int *irow,int *neq,int *nzl, 
	       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	       int *ilboun,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,int *ielmat,
	       int *ielorien,int *norien,double *orab,int *ntmat_,
	       double *t0,double *t1,int *ithermal,
	       double *prestr,int *iprestr,double *vold,
	       int *iperturb,double *sti,int *nzs,double *stx,
	       double *adb,double *aub,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstiff, 
	       int *npmat_,double *dtime,char *matname,int *mi,
               int *ncmat_,double *ttime,double *time,
               int *istep,int *kinc,int *ibody));

void FORTRAN(mafilldm,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	       int *ne,int *nodeboun,int *ndirboun,double *xboun, 
	       int *nboun,int *ipompc,int *nodempc,double *coefmpc, 
	       int *nmpc,int *nodeforc,int *ndirforc,
	       double *xforc,int *nforc,int *nelemload,char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr,
	       double *ad,double *au,int *nactdof, 
	       int *icol,int *jq,int *irow,int *neq,int *nzl, 
	       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	       int *ilboun,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,int *ielmat,
	       int *ielorien,int *norien,double *orab,int *ntmat_,
	       double *t0,double *t1,int *ithermal,
	       double *prestr,int *iprestr,double *vold,
	       int *iperturb,double *sti,int *nzs,double *stx,
	       double *adb,double *aub,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstiff, 
	       int *npmat_,double *dtime,char *matname,int *mi,
               int *ncmat_,double *ttime,double *time,
               int *istep,int *kinc,int *ibody));

void FORTRAN(mafillnet,(int *itg,int *ieg,int *ntg,
			double *ac,int *nload,char *sideload,
			int *nelemload,double *xloadact,char *lakon,
			int *ntmat_,double *v,double *shcon,int *nshcon,
			int *ipkon,int *kon,double *co,int *nflow,
			int *iinc,int *istep,
			double *dtime,double *ttime,double *time,
			int *ielmat,int *nteq,double *prop,
			int *ielprop,int *nactdog,int *nacteq,
			double *physcon,double *rhcon,int *nrhcon,
			int *ipobody,int *ibody,double *xbody,int *nbody,
			double *vold,double *xloadold,double *reltime,
			int *nmethod,char *set,int *mi,int *nmpc,
                        int *nodempc,int *ipompc,double *coefmpc,char *labmpc));

void FORTRAN(mafillklhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
         int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
         int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdok,
         int *icolk,int *jqk,int *irowk,int *neqk,int *nzlk,int *ikmpc,
         int *ilmpc,int *ikboun,int *ilboun,int *nzsk,double *adbk,
	 double *aubk,int *ipvar,double *var));

void FORTRAN(mafillkrhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
        int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
        int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nelemface,
        char *sideface,int *nface,int *nactdok,int *neqk,int *nmethod,
        int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,double *rhcon,
        int *nrhcon,int *ielmat,int *ntmat_,double *vold,double *voldaux,
        int *nzsk,double *dtimef,char *matname,int *mi,int *ncmat_,
        double *shcon,int *nshcon,double *v,double *theta1,double *bk,
        double *bt,double *voldtu,int *isolidsurf,int *nsolidsurf,
        int *ifreestream,int *nfreestream,double *xsolidsurf,double *yy,
	int *compressible,int *turbulent,int *ithermal,int *ipvar,
	double *var,int *ipvarf,double *varf, int *nea, int *neb));

void *mafillkrhsmt(void *i);

void FORTRAN(mafillplhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
          int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
          int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
          int *icolp,int *jqp,int *irowp,int *neqp,int *nzlp,int *ikmpc,
          int *ilmpc,int *ikboun,int *ilboun,int *nzsp,double *adbp,
	  double *aubp,int *nmethod,int *iexplicit,int *ipvar,double *var));

void FORTRAN(mafillprhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nelemface,
       char *sideface,int *nface,double *b,int *nactdoh,int *icolp,int *jqp,
       int *irowp,int *neqp,int *nzlp,int *nmethod,int *ikmpc,int *ilmpc,
       int *ikboun,int *ilboun,double *rhcon,int *nrhcon,int *ielmat,
       int *ntmat_,double *vold,double *voldaux,int *nzsp,double *dtl,
       char *matname,int *mi,int *ncmat_,double *shcon,int *nshcon,
       double *v,double *theta1,int *iexplicit,
       double *physcon,int *nea,int *neb,double *dtimef,int *ipvar, 
       double *var,int *ipvarf,double *varf));

void *mafillprhsmt(void *i);

void FORTRAN(mafillsm,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	       int *ne,int *nodeboun,int *ndirboun,double *xboun, 
	       int *nboun,int *ipompc,int *nodempc,double *coefmpc, 
	       int *nmpc,int *nodeforc,int *ndirforc,
	       double *xforc,int *nforc,int *nelemload,char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr,
	       double *ad,double *au,double *bb,int *nactdof, 
	       int *icol,int *jq,int *irow,int *neq,int *nzl, 
	       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	       int *ilboun,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,int *ielmat,
	       int *ielorien,int *norien,double *orab,int *ntmat_,
	       double *t0,double *t1,int *ithermal,
	       double *prestr,int *iprestr,double *vold,
	       int *iperturb,double *sti,int *nzs,double *stx,
	       double *adb,double *aub,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               double *xstiff, 
	       int *npmat_,double *dtime,char *matname,int *mi,
               int *ncmat_,int *mass,int *stiffness,int *buckling,int *rhs,
               int *intscheme,double *physcon,double *shcon,int *nshcon,
               double *cocon,int *ncocon,double *ttime,double *time,
               int *istep,int *kinc,int *coriolis,int *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,int *nstate_,double *xstateini,
	       double *xstate,double *thicke,double *xnormastface));

void FORTRAN(mafillsmcs,(double *co,int *nk,int *kon,int *ipkon, 
               char *lakon,
	       int *ne,int *nodeboun,int *ndirboun,double *xboun, 
	       int *nboun,int *ipompc,int *nodempc,double *coefmpc, 
	       int *nmpc,int *nodeforc,int *ndirforc,
	       double *xforc,int *nforc,int *nelemload,char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
	       int *nbody,double *cgr, 
	       double *ad,double *au,double *bb,int *nactdof, 
	       int *icol,int *jq,int *irow,int *neq,int *nzl, 
	       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	       int *ilboun,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,int *ielmat,
	       int *ielorien,int *norien,double *orab,int *ntmat_,
	       double *t0,double *t1,int *ithermal,
	       double *prestr,int *iprestr,double *vold,
	       int *iperturb,double *sti,int *nzs,double *stx,
	       double *adb,double *aub,int *iexpl,double *plicon,
               int *nplicon,double *plkcon,int *nplkcon,double *xstiff, 
	       int *npmat_,double *dtime,char *matname,int *mi,
               int *ics,double *cs,int *nm,int *ncmat_,char *labmpc,
               int *mass,int *stiffness,int *buckling,int *rhs,
               int *intscheme,int *mcs,int *coriolis,int *ibody,
	       double *xloadold,double *reltime,int *ielcs,double *veold,
	       double *springarea,double *thicke,double *xnormastface));

void FORTRAN(mafilltlhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
       int *icolt,int *jqt,int *irowt,int *neqt,int *nzlt,int *ikmpc,
       int *ilmpc,int *ikboun,int *ilboun,int *nzst,double *adbt,
       double *aubt,int *ipvar,double *var));
	  
void FORTRAN(mafilltrhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
             int *ipompc,int *nodempc,double *coefmpc,int *nmpc,
             int *nodeforc,int *ndirforc,double *xforc,int *nforc,
             int *nelemload,char *sideload,double *xload,int *nload,
             double *xbody,int *ipobody,int *nbody,double *b,int *nactdoh,
             int *neqt,int *nmethod,int *ikmpc,int *ilmpc,int *ikboun,
             int *ilboun,double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,
             double *t0,int *ithermal,double *vold,double *voldaux,int *nzst,
             double *dtl,char *matname,int *mi,int *ncmat_,
             double *physcon,double *shcon,int *nshcon,double *ttime,
             double *timef,int *istep,int *iinc,int *ibody,double *xloadold,
             double *reltime,double *cocon,int *ncocon,int *nelemface,
	     char *sideface,int *nface,int *compressible,double *v,
	     double *voldtu,double *yy,int *turbulent,int *nea,int *neb,
             double *dtimef,int *ipvar,double *var,int *ipvarf,double *varf));

void *mafilltrhsmt(void *i);

void FORTRAN(mafillvlhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nactdoh,
       int *icolv,int *jqv,int *irowv,int *neqv,int *nzlv,int *ikmpc,
       int *ilmpc,int *ikboun,int *ilboun,int *nzsv,double *adbv,
       double *aubv,int *ipvar,double *var));

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
         double *voldaux,int *nzsv,double *dtl,char *matname,
         int *mi,int *ncmat_,double *physcon,double *shcon,int *nshcon,
         double *ttime,double *timef,int *istep,int *iinc,int *ibody,
         double *xloadold,int *turbulent,double *voldtu,
	 double *yy,int *nelemface,char *sideface,int *nface,int *compressible,
	 int *nea,int *neb,double *dtimef,int *ipvar,double *var,
	 int *ipvarf,double *varf,double *sti));

void *mafillv1rhsmt(void *i);

void FORTRAN(mafillv2rhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,int *nodeboun,int *ndirboun,double *xboun,int *nboun,
       int *ipompc,int *nodempc,double *coefmpc,int *nmpc,double *
       b,int *nactdoh,int *icolv,int *jqv,int *irowv,int *neqv,int *nzlv,
       int *nmethod,int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
       double *vold,int *nzsv,double *dtl,double *v,double *theta2,
       int *iexplicit,int *nea,int *neb,int *mi,double *dtimef, 
       int *ipvar,double *var,int *ipvarf,double *varf));

void *mafillv2rhsmt(void *i);

void mastruct(int *nk,int *kon,int *ipkon,char*lakon,int *ne,
	      int *nodeboun,int *ndirboun,int *nboun,int *ipompc,
	      int *nodempc,int *nmpc,int *nactdof,int *icol,
	      int *jq,int **mast1p,int **irowp,int *isolver,int *neq,
	      int *nnn,int *ikmpc,int *ilmpc,
	      int *ipointer,int *nzs,int *nmethod,
              int *ithermal,int *ikboun,int *ilboun,int *iperturb,
              int *mi);

void mastructcs(int *nk,int *kon,int *ipkon,char *lakon,
	       int *ne,int *nodeboun,
	       int *ndirboun,int *nboun,int *ipompc,int *nodempc,
	       int *nmpc,int *nactdof,int *icol,int *jq,int **mast1p,
	       int **irowp,int *isolver,int *neq,int *nnn, 
	       int *ikmpc,int *ilmpc,int *ipointer,
	       int *nzs,int *nmethod,int *ics,double *cs,
	       char *labmpc,int *mcs,int *mi);

void mastructf(int *nk,int *kon,int *ipkon,char *lakon,int *ne,
	      int *nodeboun,int *ndirboun,int *nboun,int *ipompc,
	      int *nodempc,int *nmpc,int *nactdoh,int *icolt,
	      int *icolv,int *icolp,int *icolk,int *jqt,int *jqv,int *jqp,
	      int *jqk,int **mast1p,int **irowtp,int **irowvp,int **irowpp, 
	      int **irowkp,int *isolver,int *neqt,int *neqv,int *neqp,
	      int *neqk,int *ikmpc,int *ilmpc,int *ipointer, 
	      int *nzst,int *nzsv,int *nzsp,int *nzsk, 
	      int *ithermal,int *ikboun,int *ilboun,int *turbulent,
              int *nactdok,int *ifreestream,int *nfreestream,
	      int *isolidface,int *nsolidface,int *nzs,int *iexplicit,
	     int *ielmat,int *inomat,char *labmpc);

void mastructrad(int *ntr,int *nloadtr,char *sideload,int *ipointerrad,
              int **mast1radp,int **irowradp,int *nzsrad,
	      int *jqrad,int *icolrad);

void FORTRAN(mpcrem,(int *i,int *mpcfree,int *nodempc,int *nmpc,int *ikmpc,
                     int *ilmpc,char *labmpc,double *coefmpc,int *ipompc));

void FORTRAN(mult,(double *matrix,double *trans,int *n));
    
void multimortar(double *au,double *ad,int *irow,int *jq,int *nzs,
	   double *aubd,double *bdd,int *irowbd,int *jqbd,int *nzsbd,
	   double **aucp,double *adc,int **irowcp,int *jqc,int *nzsc,
           double *auqdt,int *irowqdt,int *jqqdt,int *nzsqdt,
	   int *neq,double *b,double *bhat,int* islavnode,int*imastnode,
           int*nactdof,int nslavs,int nmasts,int * mi);

void multi_rect(double *au_1,int * irow_1,int * jq_1,int n_1,int m_1,
	       double *au_2,int * irow_2,int * jq_2,int n_2,int m_2,
		double **au_rp,int **irow_rp,int * jq_r);

void multi_rectv(double *au_1,int * irow_1,int * jq_1,int n_1,int m_1,
		 double * b,double ** v_rp);

void multi_scal(double *au_1,int * irow_1,int * jq_1,
	       double *au_2,int * irow_2,int * jq_2,
	       int m,int n,double*value,int *flag);

void FORTRAN(nident,(int *x,int *px,int *n,int *id));

void FORTRAN(nodestiedface,(char *tieset,int *ntie,int *ipkon,int *kon,
       char *lakon,char *set,int *istartset,int *iendset,int *ialset,
       int *nset,int *faceslave,int *istartfield,int *iendfield,
       int *ifield,int *nconf,int *ncone));

void nonlingeo(double **co,int *nk,int **konp,int **ipkonp,char **lakonp,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     double *ad,double *au,double *b,int *nactdof, 
	     int **icolp,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int **ikmpcp,int **ilmpcp,int *ikboun, 
	     int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int **ielmatp,
	     int **ielorienp,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double **vold,int *iperturb,double *sti,int *nzs,  
	     int *kode,double *adb,double *aub, 
	     char *filab,int *idrct,
	     int *jmax,int *jout,double *tinc,double *tper,
	     double *tmin,double *tmax,double *eme,double *xbounold,
	     double *xforcold,double *xloadold,
             double *veold,double *accold,
             char *amname,double *amta,int *namta,int *nam,
             int *iamforc,int *iamload,
             int *iamt1,double *alpha,int *iexpl,
	     int *iamboun,double *plicon,int *nplicon,double *plkcon,
	     int *nplkcon,
             double **xstatep,int *npmat_,int *istep,double *ttime,
	     char *matname,double *qaold,int *mi,
             int *isolver,int *ncmat_,int *nstate_,int *iumat,
             double *cs,int *mcs,int *nkon,double **ener,int *mpcinfo,
             int *nnn,char *output,
             double *shcon,int *nshcon,double *cocon,int *ncocon,
             double *physcon,int *nflow,double *ctrl, 
             char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *ikforc,int *ilforc,double *trab, 
             int *inotr,int *ntrans,double **fmpcp,char *cbody,
             int *ibody,double *xbody,int *nbody,double *xbodyold,
             int *ielprop,double *prop,int *ntie,char *tieset,
	     int *itpamp,int *iviewfile,char *jobnamec,double *tietol,
             int *nslavs,double *thicke);

void FORTRAN(nonlinmpc,(double *co,double *vold,int *ipompc,int *nodempc,
		   double *coefmpc,char *labmpc,int *nmpc,int *ikboun,
		   int *ilboun,int *nboun,double *xbounact,double *aux,
		   int *iaux,int *maxlenmpc,int *ikmpc,int *ilmpc,
                   int *icascade,int *kon,int *ipkon,char *lakon,
		   int *ne,double *reltime,int *newstep,double *xboun,
		   double *fmpc,int *newinc,int *idiscon,int *ncont,
		   double *trab,int *ntrans,int *ithermal,int *mi));

void FORTRAN(op,(int *,double *,double *,double *,double *,double *,int *,
	 int *,int *));

void FORTRAN(op_corio,(int *,double *,double *,double *,double *,double *,
              int *,int *,int *));

void FORTRAN(openfile,(char *jobname,char *output));

void FORTRAN(opnonsym, (int *neq,double *aux,double *b,double *bhat, 
          double *bdd,double*bdu,int *jqbd,int *irowbd));

void FORTRAN(opnonsymt, (int *neq,double *aux,double *b,double *bhat, 
          double *bdd,double*bdu,int *jqbd,int *irowbd));

void FORTRAN(out,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	  int *ne0,double *v, 
	  double *stn,int *inum,int *nmethod, 
	  int *kode,char *filab,double *een,double *t1,
          double *fn,double *time,double *epl,int *ielmat,char *matname,
          double *enern,double *xstaten,int *nstate_,int *istep,
          int *iinc,int *iperturb,double *ener,int *mi,char *output,
          int *ithermal,double *qfn,int *mode,int *noddiam,double *trab,
          int *inotr,int *ntrans,double *orab,int *ielorien,int *norien,
	  char *description,int *ipneigh,int *neigh,double *stx,
          double *vr,double *vi,double *stnr,double *stni,
	  double *vmax,double *stnmax,int *ngraph,double *veold,int *ne,
          double *cs,char *set,int *nset,int *istartset,int *iendset,
	  int *ialset,double *eenmax,double *fnr,double *fni,
	  double *emn, double *thicke));

void FORTRAN(postview,(int *ntr,char *sideload,int *nelemload,int *kontri,
             int *ntri,int *nloadtr,double *tenv,double *adview,double *auview,
             double *area,double *fenv,int *jqrad,int *irowrad,int *nzsrad));

void FORTRAN(precfd,(int *nelemface,char *sideface,int *nface,int *ipoface,
        int *nodface,int *ne,int *ipkon,int *kon,char *lakon,
        int *ikboun,int *ilboun,double *xboun,int *nboun,int *nk,
        int *isolidsurf,int *nsolidsurf,int *ifreestream,int *nfreestream,
        int *neighsolidsurf,int *iponoel,int *inoel,int *inoelfree,
	int *nef,double *co,int *ipompc,int *nodempc,int *ikmpc,
        int *ilmpc,int *nmpc,char *set,int *istartset,int *iendset,
        int *ialset,int *nset, int *iturbulent));

void prediction(double *uam,int *nmethod,double *bet,double *gam,double *dtime,
               int *ithermal,int *nk,double *veold,double *accold,double *v,
	       int *iinc,int *idiscon,double *vold,int *nactdof,int *mi);

void preiter(double *ad,double **aup,double *b,int **icolp,int **irowp, 
	     int *neq,int *nzs,int *isolver,int *iperturb);

void FORTRAN(presgradient,(int *iponoel,int *inoel,double *sa,double *sav,
            int *nk,double *dt,double *shockcoef,double *dtimef,int *ipkon,
	    int *kon,char *lakon,double *vold,int *mi,int *compressible,
	    int *nmethod,double *dtl,int *isolidsurf,int *nsolidsurf,
	    double *co,int *euler));
  
void linstatic(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     double *ad,double *au,double *b,int *nactdof, 
	     int **icolp,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs, 
	     int *kode,double *adb,double *aub, 
	     char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *isolver,
	     int *mi,int *ncmat_,int *nstate_,double *cs,int *mcs,
             int *nkon,double *ener,double *xbounold,
	     double *xforcold,double *xloadold,
             char *amname,double *amta,int *namta,
             int *nam,int *iamforc,int *iamload,
             int *iamt1,int *iamboun,double *ttime,char *output, 
             char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,double *trab, 
             int *inotr,int *ntrans,double *fmpc,char *cbody,int *ibody,
	     double *xbody,int *nbody,double *xbodyold,double *tper, 
             double *thicke);

void FORTRAN(printout,(char *set,int *nset,int *istartset,int *iendset,
             int *ialset,int *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,int *ipkon,
             char *lakon,double *stx,double *eme,double *xstate,
             double *ener,int *mi,int *nstate_,int *ithermal,
             double *co,int *kon,double *qfx,double *ttime,double *trab,
             int *inotr,int *ntrans,double *orab,int *ielorien,
	     int *norien,int *nk,int *ne,int *inum,char *filab,double *vold,
             int *ikin));

void FORTRAN(printoutface,(double *co,double *rhcon,int *nrhcon,int *ntmat_,
            double *vold,double *shcon,int *nshcon,double *cocon,
            int *ncocon,int *compressible,int *istartset,int *iendset,
            int *ipkon,char *lakon,int *kon,int *ialset,char *prset,
	    double *timef,int *nset,char *set,int *nprint,char *prlab,
	    int *ielmat,int *mi));

void radcyc(int *nk,int *kon,int *ipkon,char *lakon,int *ne,
	    double *cs,int *mcs,int *nkon,int *ialset,int *istartset,
            int *iendset,int **kontrip,int *ntri,
            double **cop,double **voldp,int *ntrit,int *inocs,int *mi);

void radflowload(int *itg,int *ieg,int *ntg,int *ntr,double *adrad,
       double *aurad,
       double *bcr,int *ipivr,double *ac,double *bc,int *nload,
       char *sideload,int *nelemload,double *xloadact,char *lakon,int *ipiv,
       int *ntmat_,double *vold,double *shcon,int *nshcon,int *ipkon,
       int *kon,double *co,int *kontri,int *ntri,
       int *nloadtr,double *tarea,double *tenv,double *physcon,double *erad,
       double **adviewp,double **auviewp,
       int *nflow,int *ikboun,double *xboun,int *nboun,int *ithermal,
       int *iinc,int *iit,double *cs,int *mcs,int *inocs,int *ntrit,
       int *nk,double *fenv,int *istep,double *dtime,double *ttime,
       double *time,int *ilboun,int *ikforc,int *ilforc,double *xforc,
       int *nforc,double *cam,int *ielmat,int *nteq,double *prop,
       int *ielprop,int *nactdog,int *nacteq,int *nodeboun,int *ndirboun,
       int *network,double *rhcon,int *nrhcon,
       int *ipobody,int *ibody,double *xbody,int *nbody,int *iviewfile,
       char *jobnamef,double *ctrl,double *xloadold,double *reltime,
       int *nmethod,char *set,int *mi,int * istartset,int* iendset,
       int *ialset,int *nset,int *ineighe,int *nmpc,int *nodempc,
       int *ipompc,double *coefmpc,char *labmpc,int *iemchange,int *nam,
       int *iamload,int *jqrad,int *irowrad,int *nzsrad,int *icolrad);

void FORTRAN (radmatrix,(int *ntr,double *adrad,double *aurad,double *bcr,
       char *sideload,int *nelemload,double *xloadact,char *lakon,
       double *vold,int *ipkon,int *kon,double *co,int *nloadtr,
       double *tarea,double *tenv,double *physcon,double *erad,
       double *adview,double *auview,int *ithermal,int *iinc,
       int *iit,double *fenv,int *istep,
       double *dtime,double *ttime,double *time,int *iviewfile,
       double *xloadold,double *reltime,int *nmethod,
       int *mi,int *iemchange,int *nam,int *iamload,int *jqrad,
       int *irowrad,int *nzsrad));

void FORTRAN (radresult,(int *ntr,double *xloadact,double *bcr,
       int *nloadtr,double *tarea,double * tenv,double *physcon,double *erad,
       double *auview,double *fenv,int *irowrad,int *jqrad, 
       int *nzsrad,double *q));

void FORTRAN(readforce,(double *zc,int *neq,int *nk,int *nev,int *nactdof,
	     int *ikmpc,int *nmpc,int *ipompc,int *nodempc,int *mi, 
             double *coefmpc));

void readinput(char *jobnamec,char **inpcp,int *nline,int *nset,int *ipoinp,
        int **inpp,int **ipoinpcp,int *ithermal); 

void FORTRAN(readview,(int *ntr,double *adview,double *auview,double *fenv,
             int *nzsrad,int *ithermal,char *jobnamef));

void FORTRAN(rearrange,(double *au,int *irow,int *icol,int *ndim,int *neq));

void FORTRAN(rectcyl,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,int *nk, 
                      int *icntrl,double *t,char *filab,int *imag, 
                      int *mi,double *emn));

void FORTRAN(rectcylexp,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,int *nkt, 
		      int *icntrl,double *t,char *filab,int *imag,int *mi,
		      int *iznode,int *nznode,int *nsectors,int *nk,
                      double *emn));

void FORTRAN(rectcylvi,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,int *nk, 
		      int *icntrl,double *t,char *filab,int *imag,int *mi));

void reloadcontact(char *lakon,int *ipkon,int *kon,
	       int **nelemloadp,char **sideloadp,int **iamloadp, 
	       double **xloadp,int *nload,int *ne,double *t1,int *iamt1,
	       int *nam,int *ithermal,double *vold,int *mi);

void remastruct(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
              int *mpcfree,int *nodeboun,int *ndirboun,int *nboun,
              int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
              char *labmpc,int *nk,
              int *memmpc_,int *icascade,int *maxlenmpc,
              int *kon,int *ipkon,char *lakon,int *ne,int *nnn,
              int *nactdof,int *icol,int *jq,int **irowp,int *isolver,
              int *neq,int *nzs,int *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,int *ithermal,
	      int *iperturb,int *mass,int *mi);

void remcontmpc(int *nmpc,char *labmpc,int *mpcfree,int *nodempc,
		int *ikmpc,int *ilmpc,double *coefmpc,int *ipompc);

void remeshcontact(int *ntie,char *tieset,int *nset,char *set,
               int *istartset,int *iendset,int **ialsetp,
               char **lakonp,int **ipkonp,int **konp,
	       int *nalset,int *nmpc,int *mpcfree,int *memmpc_,
               int **ipompcp,char **labmpcp,int **ikmpcp,int **ilmpcp,
               double **fmpcp,int **nodempcp,double **coefmpcp,
	       double **cop,int *nmpc_,int *mi,int *nk,int *nkon, 
	       int *ne,int *nk_,int *ithermal,int **ielmat,
	       int **ielorien,double **t0p,double **voldp,double **veoldp,
               int *ncont,double **xstate,int *nstate_,double **prestr,
               int *iprestr);

void FORTRAN(remeshcontactel,(char *tieset,int *ntie,char *set,int *nset,
       int *istartset,int *iendset,int *ialset,int *ipkon,int *kon,int *nkon,
       char *lakon,int *nodface,int *ipoface,int *nk,int *ipompc,int *nodempc,
       int *ikmpc,int *ilmpc,int *nmpc,int *nmpc_,char *labmpc,double *coefmpc,
       int *mpcfree,int *nalset,double *co,int *ithermal,int *nk0,int *ne,
       int *ielmat,int *ielorien,int *mi,double *t0,double *vold,
       double *veold,int *iponoel,int *inoel,double *xstate,int *nstate_,
       double *prestr,int *iprestr));

void FORTRAN(remeshload,(int *ipkon,int *kon,char *lakon,
	     int *nelemload,char *sideload,int *iamload,double *xload,
	     int *nload,int *ne,double *t1,int *iamt1,int *nam,int *ithermal,
             double *vold,int *mi));

void FORTRAN(remeshsurf,(char *tieset,int *ntie,char *set,int *nset,
       int *istartset,int *iendset,int *ialset,int *ipkon,int *kon,
       char *lakon,int *nodface,int *ipoface,int *nface,int *nquadface,
       int *ninterface,int *ntotface,int *nk,int *ne,int *iponoel,int *inoel,
       int *ntets2remesh));

void FORTRAN(renumber,(int *nk,int *kon,int *ipkon,char *lakon,int *ne, 
	       int *ipompc,int *nodempc,int *nmpc,int *nnn,int *npn, 
	       int *adj,int *xadj,int *iw,int *mmm,int *xnpn,int *inum1, 
               int *inum2));

void FORTRAN(restartshort,(int *nset,int *nload,int *nbody,int *nforc,
    int *nboun,
    int *nk,int *ne,int *nmpc,int *nalset,int *nmat,int *ntmat,int *npmat,
    int *norien,int *nam,int *nprint,int *mint,int *ntrans,int *ncs,
    int *namtot,int *ncmat,int *memmpc,int *ne1d,int *ne2d,int *nflow,
    char *set,int *meminset,int *rmeminset,char *jobnamec,int *irestartstep,
    int *icntrl,int *ithermal,int *nener,int *nstate_,int *ntie,int *nslavs,
    int *nkon));

void FORTRAN(restartwrite,(int *istep,int *nset,int*nload,int *nforc, 
  int * nboun,int *nk,int *ne,int *nmpc,int *nalset,int *nmat,int *ntmat_, 
  int *npmat_,int *norien,int *nam,int *nprint,int *mi, 
  int *ntrans,int *ncs_,int *namtot_,int *ncmat_,int *mpcend, 
  int *maxlenmpc,int *ne1d, 
  int *ne2d,int *nflow,int *nlabel,int *iplas,int *nkon,int *ithermal, 
  int *nmethod,int *iperturb,int *nstate_,int *nener,char *set, 
  int *istartset,int *iendset,int *ialset,double *co,int *kon,int *ipkon, 
  char *lakon,int *nodeboun,int *ndirboun,int *iamboun,double *xboun, 
  int *ikboun,int *ilboun,int *ipompc,int *nodempc,double *coefmpc, 
  char *labmpc,int *ikmpc,int *ilmpc,int *nodeforc,int *ndirforc, 
  int *iamforc,double *xforc,int *ikforc,int *ilforc,int *nelemload, 
  int *iamload,char *sideload,double *xload,  
  double *elcon,int *nelcon,double *rhcon,int *nrhcon,double *alcon, 
  int *nalcon,double *alzero,double *plicon,int *nplicon,double *plkcon, 
  int *nplkcon,char *orname,double *orab,int *ielorien,double *trab, 
  int *inotr,char *amname,double *amta,int *namta,double *t0,double *t1, 
  int *iamt1,double *veold,int *ielmat,char *matname, 
  char *prlab,char *prset,char *filab,double *vold, 
  int *nodebounold,int *ndirbounold,double *xbounold,double *xforcold, 
  double *xloadold,double *t1old,double *eme,int *iponor, 
  double *xnor,int *knor,double *thickn,double *thicke,double *offset, 
  int *iponoel,int *inoel,int *rig, 
  double *shcon,int *nshcon,double *cocon,int *ncocon, 
  int *ics,double *sti,double *ener,double *xstate, 
  char *jobnamec,int *infree,int *nnn,double *prestr,int *iprestr,
  char *cbody,int *ibody,double *xbody,int *nbody,double *xbodyold,
  double *ttime,double *qaold,double *cs,
  int *mcs,char *output,double *physcon,double *ctrl,char *typeboun,
  double *fmpc,char *tieset,int *ntie,double *tietol,int *nslavs,
  double *t0g,double *t1g));

void FORTRAN(resultnet,(int *itg,int *ieg,int *ntg,
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
                        double *physcon,double *camt,double *camf,
                        double *camp,double *rhcon,int *nrhcon,
			int *ipobody,int *ibody,double *xbody,int *nbody,
                        double *dtheta,double *vold,double *xloadold,
                        double *reltime,int *nmethod,char *set,int *mi,
                        int *ineighe,double *cama,double *vamt,
                        double *vamf,double *vamp,double *vama,
                        int *nmpc,int *nodempc,int *ipompc,double *coefmpc,
                        char *labmpc));

void results(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne,double *v,double *stn,int *inum, 
	     double *stx,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,int *ithermal,double *prestr, 
             int *iprestr,char *filab,double *eme,double *emn,
             double *een,int *iperturb,double *f,double *fn,int *nactdof,
             int *iout,double *qa,
	     double *vold,double *b,int *nodeboun,int *ndirboun,
	     double *xboun,int *nboun,int *ipompc,int *nodempc,
	     double *coefmpc,char *labmpc,int *nmpc,int *nmethod, 
             double *vmax,int *neq,double *veold,double *accold,
	     double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             int *nplicon,double *plkcon,int *nplkcon,
             double *xstateini,double *xstiff,double *xstate,int *npmat_,
	     double *epl,char *matname,int *mi,int *ielas,
	     int *icmd,int *ncmat_,int *nstate_,double *stiini,
	     double *vini,int *ikboun,int *ilboun,double *ener,
	     double *enern,double *sti,double *xstaten,double *eei,
             double *enerini,double *cocon,int *ncocon,char *set, 
             int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             int *inotr,int *ntrans,double *fmpc,int *nelemload,
	     int *nload,int *ikmpc,int *ilmpc,int *istep,int *iinc,
	     double *springarea,double *reltime,int *ne0,double *xforc,
             int *nforc,double *thicke,double *xnormastface);

void  FORTRAN(resultsforc,(int *nk,double *f,double *fn,int *nactdof,
       int *ipompc,int *nodempc,double *coefmpc,char *labmpc,int *nmpc,
       int *mi,double *fmpc,int *calcul_fn,
       int *calcul_f));

void FORTRAN(resultsini,(int *nk,double *v,int *ithermal,char *filab,
       int *iperturb,double *f,double *fn,int *nactdof,int *iout,
       double *qa,double *vold,double *b,int *nodeboun,int *ndirboun,
       double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
       char *labmpc,int *nmpc,int *nmethod,double *cam,int *neq,
       double *veold,double *accold,double *bet,double *gam,double *dtime,
       int *mi,double *vini,int *nprint,char *prlab,int *intpointvar,
       int *calcul_fn,int *calcul_f,int *calcul_qa,int *calcul_cauchy,
       int *iener,int *ikin,int *intpointvart,double *xforc,int *nforc));

void FORTRAN(resultsk,(int *nk,int *nactdoh,double *vtu,double *solk,
      double *solt,int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultsmech,(double *co,int *kon,int *ipkon,char *lakon,int *ne,
          double *v,double *stx,double *elcon,int *nelcon,double *rhcon,
          int *nrhcon,double *alcon,int *nalcon,double *alzero,int *ielmat,
          int *ielorien,int *norien,double *orab,int *ntmat_,double *t0,
          double *t1,int *ithermal,double *prestr,int *iprestr,double *eme,
          int *iperturb,double *fn,int *iout,double *qa,double *vold,
          int *nmethod,double *veold,double *dtime,double *time,
          double *ttime,double *plicon,int *nplicon,double *plkcon,
          int *nplkcon,double *xstateini,double *xstiff,double *xstate,
          int *npmat_,char *matname,int *mi,int *ielas,int *icmd,int *ncmat_,
          int *nstate_,double *stiini,double *vini,double *ener,double *eei,
          double *enerini,int *istep,int *iinc,double *springarea,
          double *reltime,int *calcul_fn,int *calcul_qa,int *calcul_cauchy,
	  int *iener,int *ikin,int *nal,int *ne0,double *thicke,
	  double *xnormastface,int *nea,int *neb));

void *resultsmechmt(void *i);

void  FORTRAN(resultsprint,(double *co,int *nk,int *kon,int *ipkon,
       char *lakon,int *ne,double *v,double *stn,int *inum,double *stx,
       int *ielorien,int *norien,double *orab,double *t1,int *ithermal,
       char *filab,double *een,int *iperturb,double *fn,int *nactdof,
       int *iout,double *vold,int *nodeboun,int *ndirboun,int *nboun,
       int *nmethod,double *ttime,double *xstate,double *epn,int *mi,
       int *nstate_,double *ener,double *enern,double *xstaten,double *eei,
       char *set,int *nset,int *istartset,int *iendset,int *ialset,int *nprint,
       char *prlab,char *prset,double *qfx,double *qfn,double *trab,int *inotr,
       int *ntrans,int *nelemload,int *nload,int *ikin,int *ielmat,
       double *thicke,double *eme,double *emn));

void FORTRAN(resultsp,(int *nk,int *nactdoh,double *v,double *sol,
	     int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(resultst,(int *nk,int *nactdoh,double *v,double *sol,
	     int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(resultstherm,(double *co,int *kon,int *ipkon,
       char *lakon,int *ne,double *v,
       double *elcon,int *nelcon,double *rhcon,int *nrhcon,int *ielmat,
       int *ielorien,int *norien,double *orab,int *ntmat_,double *t0,
       int *iperturb,double *fn,int *iout,
       double *qa,double *vold,int *ipompc,int *nodempc,
       double *coefmpc,int *nmpc,double *dtime,
       double *time,double *ttime,double *plicon,int *nplicon,double *xstateini,
       double *xstiff,double *xstate,int *npmat_,char *matname,
       int *mi,int *ncmat_,int *nstate_,double *cocon,int *ncocon,
       double *qfx,int *ikmpc,int *ilmpc,int *istep,
       int *iinc,double *springarea,int *calcul_fn,int *calcul_qa,int *nal,
       int *nea,int *neb));

void *resultsthermmt(void *i);

void FORTRAN(resultsv1,(int *nk,int *nactdoh,double *v,double *sol,
	    int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(resultsv2,(int *nk,int *nactdoh,double *v,double *sol,
	     int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(shape3tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,int *iflag));

void FORTRAN(shape4q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,int *iflag));

void FORTRAN(shape4tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(shape6tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,int *iflag));

void FORTRAN(shape6w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(shape8h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(shape8q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,int *iflag));

void FORTRAN(shape10tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(shape15w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(shape20h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,int *iflag));

void FORTRAN(rhs,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	       int *ne,int *ipompc,int *nodempc,double *coefmpc, 
	       int *nmpc,int *nodeforc,int *ndirforc,
	       double *xforc,int *nforc,int *nelemload,char *sideload,
	       double *xload,int *nload,double *xbody,int *ipobody,
               int *nbody,double *cgr,double *bb,int *nactdof,int *neq, 
	       int *nmethod,int *ikmpc,int *ilmpc,
	       double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	       double *alcon,int *nalcon,double *alzero,int *ielmat,
	       int *ielorien,int *norien,double *orab,int *ntmat_,
	       double *t0,double *t1,int *ithermal, 
               int *iprestr,double *vold,int *iperturb,int *iexpl,
               double *plicon,int *nplicon,double *plkcon,int *nplkcon,
               int *npmat_,double *ttime,double *time,int *istep,
               int *iinc,double *dtime,double *physcon,int *ibody,
	       double *xbodyold,double *reltime,double *veold,
	       char *matname,int *mi,int *ikactmech,int *nactmech));
       
void FORTRAN(slavintmortar,(char *tieset,int *ntie,int *itietri,int *ipkon,
        int *kon,char *lakon,char *set,double *cg,double *straight,
        int *nintpoint,int *koncont,double *co,double *vold,double *xo,
        double *yo,double *zo,double *x,double *y,double *z,int *nx,
        int *ny,int *nz,int *nset,
        int *iinc,int *iit,
        int *islavsurf,int *imastsurf,double *pmastsurf,
        int *itiefac,int *islavnode,int *nslavnode,double *slavnor,
	double *slavtan,int *imastop,double *gap,int *islavact,
	int *mi,int *ncont,int *ipe,int *ime,double *pslavsurf,
        double* pslavdual,int *i,int *l,int *ntri));

void FORTRAN(smooth,(double *adbv,double *aubv,double *adl,
            double *sol,double *aux,int *icolv,int *irowv,
	    int *jqv,int *neqv,int *nzlv,double *csmooth));

void FORTRAN(smoothshock,(double *adbv,double *aubv,double *adl,
	    double *addiv,double *sol,double *aux,int *icolv,int *irowv,
	    int *jqv,int *neqv,int *nzlv,double *sa));

void FORTRAN(solveeq,(double *adbv,double *aubv,double *adl,double *addiv,
            double *b,double *sol,double *aux,int *icolv,int *irowv,
            int *jqv,int *neqv,int *nzsv,int *nzlv));

void FORTRAN(spcmatch,(double *xboun,int *nodeboun,int *ndirboun,int *nboun,
	       double *xbounold,int *nodebounold,int *ndirbounold,
	       int *nbounold,int *ikboun,int *ilboun,double *vold,
	       double *reorder,int *nreorder,int *mi));

void FORTRAN(splitline,(char *text,char *textpart,int *n));

void spooles(double *ad,double *au,double *adb,double *aub,
             double *sigma,double *b,
	     int *icol,int *irow,int *neq,int *nzs,int *symmtryflag,
             int *inputformat);

void FORTRAN(springforc,(double *xl,int *konl,double *vl,int *imat,
             double *elcon,int *nelcon,double *elas,double *fnl,int *ncmat_,
             int *ntmat_,int *nope,char *lakonl,double *t0l,double *t1l,
             int *kodem,double *elconloc,double *plicon,int *nplicon,
	     int *npmat_,double *veloldl,double *senergy,int *iener,
	     double *cstr,int *mi,double *springarea,int *nmethod,
             int *ne0,int *iperturb,int *nstate_,double *xstateini,
	     double *xstate,double *reltime,double *xnormastface));

void FORTRAN(springstiff,(double *xl,double *elas,int *konl,double *voldl,
             double *s,int *imat,double *elcon,int *nelcon,int *ncmat_,
             int *ntmat_,int *nope,char *lakonl,double *t0l,double *t1l,
             int *kode,double *elconloc,double *plicon,int *nplicon,
	     int *npmat_,int *iperturb,double *springarea,int *nmethod,
             int *mi,int *ne0,int *nstate_,double *xstateini,
	     double *xstate,double *reltime,double *xnormastface));

void steadystate(double **co,int *nk,int **kon,int **ipkon,char **lakon,int *ne, 
	  int **nodeboun,int **ndirboun,double **xboun,int *nboun,
	  int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,int *nmpc, 
	  int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	  int *nelemload,char *sideload,double *xload,
	  int *nload, 
	  int **nactdof,int *neq,int *nzl,int *icol,int *irow, 
	  int *nmethod,int **ikmpcp,int **ilmpcp,int **ikboun, 
	  int **ilboun,
	  double *elcon,int *nelcon,double *rhcon,int *nrhcon,
          double *cocon,int *ncocon,
	  double *alcon,int *nalcon,double *alzero,int **ielmat,
	  int **ielorien,int *norien,double *orab,int *ntmat_,
	  double **t0, 
	  double **t1,int *ithermal,double *prestr,int *iprestr, 
	  double **voldp,int *iperturb,double *sti,int *nzs, 
	  double *tinc,double *tper,double *xmodal,
	  double *veold,char *amname,double *amta,
	  int *namta,int *nam,int *iamforc,int *iamload,
	  int **iamt1,int *jout,int *kode,char *filab,
	  double **emep,double *xforcold,double *xloadold,
          double **t1old,int **iamboun,
          double **xbounold,int *iexpl,double *plicon,int *nplicon,
          double *plkcon,int *nplkcon,
          double *xstate,int *npmat_,char *matname,int *mi,
          int *ncmat_,int *nstate_,double **enerp,char *jobnamec,
          double *ttime,char *set,int *nset,int *istartset,
          int *iendset,int *ialset,int *nprint,char *prlab,
          char *prset,int *nener,double *trab, 
          int **inotr,int *ntrans,double **fmpcp,char *cbody,int *ibody,
          double *xbody,int *nbody,double *xbodyold,int *istep,
          int *isolver,int *jq,char *output,int *mcs,int *nkon,
	  int *ics,double *cs,int *mpcend,int **nnnp,double *ctrl,
	  int *ikforc,int *ilforc,double *thicke);

void FORTRAN(stop,());

void storecontactdof(int *nope,int *nactdof,int *mt,int *konl, 
          int **ikactcontp, 
          int *nactcont,int *nactcont_,double *bcont,double *fnl, 
          int *ikmpc,int *nmpc,int *ilmpc,int *ipompc,int *nodempc, 
	  double *coefmpc);

void FORTRAN(storeresidual,(int *nactdof,double *b,double *fn,char *filab,
             int *ithermal,int *nk,double *sti,double *stn,
             int *ipkon,int *inum,int *kon,char *lakon,
             int *ne,int *mi,double *orab,int *ielorien,
             double *co,int *nelemload,int *nload,int *nodeboun,
	     int *nboun,int *itg,int *ntg,double *vold,int *ndirboun,
             int *ielmat));

int strcmp1(const char *s1, const char *s2);

int strcpy1(char *s1, const char *s2,int length);

void FORTRAN(subspace,(double *d,double *aa,double *bb,double *cc,
             double *alpham,double *betam,int *nev,
             double *xini,double *cd,double *cv,double *time,
             double *rwork,int *lrw,int *k,int *jout,double *rpar,
	     double *bj,int *iwork,int *liw,int *iddebdf,double *bjp));

void FORTRAN(tempload,(double *xforcold,double *xforc,double *xforcact,
               int *iamforc,int *nforc,double *xloadold,double *xload,
               double *xloadact,int *iamload,int *nload,int *ibody,
               double *xbody,int *nbody,double *xbodyold,double *xbodyact, 
               double *t1old,double *t1,double *t1act,int *iamt1,
               int *nk,double *amta,int *namta,int *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               int *ithermal,int *nmethod,
	       double *xbounold,double *xboun,double *xbounact,
	       int *iamboun,int *nboun,int *nodeboun,
               int *ndirboun,int *nodeforc,int *ndirforc,int *istep,
               int *iint,double *co,double *vold,int *itg,int *ntg,
               char *amname,int *ikboun,int *ilboun,int *nelemload,
	       char *sideload,int *mi,int *ntrans,double *trab,
               int *inotr,double *veold));

void FORTRAN(temploaddiff,(double *xforcold,double *xforc,double *xforcact,
               int *iamforc,int *nforc,double *xloadold,double *xload,
               double *xloadact,int *iamload,int *nload,int *ibody,
               double *xbody,int *nbody,double *xbodyold,double *xbodyact, 
               double *t1old,double *t1,double *t1act,int *iamt1,
               int *nk,double *amta,int *namta,int *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               int *ithermal,int *nmethod,
	       double *xbounold,double *xboun,double *xbounact,
	       int *iamboun,int *nboun,int *nodeboun,
               int *ndirboun,int *nodeforc,int *ndirforc,int *istep,
               int *iint,double *co,double *vold,int *itg,int *ntg,
               char *amname,int *ikboun,int *ilboun,int *nelemload,
	       char *sideload,int *mi,double *xforcdiff,double *xloaddiff,
	       double *xbodydiff,double *t1diff,double *xboundiff,
	       int *icorrect,int *iprescribedboundary,int *ntrans,
               double *trab,int *inotr,double *veold,int *nactdof,
	       double *bcont,double *fn));

void FORTRAN(temploadmodal,(double *amta,int *namta,int *nam,double *ampli,
         double *timemin,double *ttimemin,double *dtime,double *xbounold,
         double *xboun,double *xbounmin,int *iamboun,int *nboun,
         int *nodeboun,int *ndirboun,char *amname));
    
void FORTRAN(tiefaccont,(char *lakon,int *ipkon,int *kon,int *ntie,
       char *tieset,int *nset,char *set,int *istartset,int *iendset,
       int *ialset,int *itiefac,int *islavsurf,int *islavnode,
       int *imastnode,int *nslavnode,int *nmastnode,int *nslavs,
       int *nmasts,int *ifacecount,int *iponoels,int *inoels,int *ifreenoels,
       int *mortar,int *ipoface,int *nodface,int *nk));   

void tiedcontact(int *ntie,char *tieset,int *nset,char *set,
               int *istartset,int *iendset,int *ialset,
               char *lakon,int *ipkon,int *kon,double *tietol,
               int *nmpc,int *mpcfree,int *memmpc_,
               int **ipompcp,char **labmpcp,int **ikmpcp,int **ilmpcp,
               double **fmpcp,int **nodempcp,double **coefmpcp,
	       int *ithermal,double *co,double *vold,int *cfd,
	       int *nmpc_,int *mi,int *nk);

void FORTRAN(transformatrix,(double *xab,double *p,double *a));

void FORTRAN(trianeighbor,(int *ipe,int *ime,int *imastop,int *ncont,
               int *koncont,int *ifreeme));

void FORTRAN(triangucont,(int *ncont,int *ntie,char *tieset,int *nset,
          char *set,int *istartset,int *iendset,int *ialset,int *itietri,
	  char *lakon,int *ipkon,int *kon,int *koncont,char *kind1,
          char *kind2));

#ifdef BAM
void FORTRAN(uexternaldb,(int *lop,int *lrestart,double *time,double *dtime,
                          int *kstep,int *kinc));
#endif

void FORTRAN(ufaceload,(double *co,int *ipkon,int *kon,char *lakon,
             int *nelemload,char *sideload,int *nload));

void FORTRAN(uout,(double *v,int *mi));

void FORTRAN(updatecfd,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void FORTRAN(updatecomp,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void FORTRAN(updatecon,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
	   double *physcon,int *iponoel,int *inoel,int *ithermal,
	   int *nactdoh,int *iit,int *compressible,int *ismooth,
	   double *voldtu,double *vtu,int *turbulent,int *inomat,
	   int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void FORTRAN(updatecont,(int *koncont,int *ncont,double *co,double *vold,
			 double *cg,double *straight,int *mi));

void FORTRAN(updatecontpen,(int *koncont,int *ncont,double *co,double *vold,
			 double *cg,double *straight,int *mi,int *imastnode,
                         int *nmastnode,double *xmastnor,int *ntie,
                         char *tieset,int *nset,char *set,int *istartset,
                         int *iendset,int *ialset,int *ipkon,char *lakon,
                         int *kon));

void *u_calloc(size_t num,size_t size);

void FORTRAN(writeboun,(int *nodeboun,int *ndirboun,double *xboun,
      char *typeboun,int *nboun));

void FORTRAN(writematrix,(double *au,double *ad,int *irow,int *jq,int *neq,
         int *number));

void FORTRAN(writebv,(double *,int *));

void writec(char *ifield,int nfield, FILE *f1);

void FORTRAN(writeev,(double *,int *,double *,double *));

void FORTRAN(writeevcomplex,(double *eigxx,int *nev));

void FORTRAN(writeevcs,(double *,int *,int *,double *,double *));

void FORTRAN(writeevcscomplex,(double *eigxx,int *nev,int *nm));

void FORTRAN(writehe,(int *));

void writef(double *ifield,int nfield, FILE *f1);

void writei(int *ifield,int nfield, FILE *f1);

void FORTRAN(writeim,());

void FORTRAN(writeinput,(char *inpc,int *ipoinp,int *inp,int *nline,int *ninp,
                         int *ipoinpc));

void FORTRAN(writempc,(int *,int *,double *,char *,int *));

void FORTRAN(writepf,(double *d,double *bjr,double *bji,double *freq , 
       int *nev));

void FORTRAN(writere,());

void FORTRAN(writesummary,(int *istep,int *j,int *icutb,int *l,double *ttime,
		   double *time,double *dtime));

void FORTRAN(writeview,(int *ntr,double *adview,double *auview,double *fenv,
            int *nzsrad,char *jobnamef));
