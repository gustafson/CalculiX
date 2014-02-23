/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                     */

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
	     int *nmdboun,int *ikboun,int *nboun,int *ilboun,int *ithermal));

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

void FORTRAN(adjustcontactnodes,(char *tieset,int *ntie,int *itietri,double *cg,
             double *straight,double *co,double *vold,double *xo,double *yo,
             double *zo,double *x,double *y,double *z,int *nx,int *ny,
             int *nz,int *istep,int *iinc,int *iit,int *mi,int *imastop,
             int *nslavnode,int *islavnode,char *set,int *nset,int *istartset,
	     int *iendset,int *ialset,double *tietol,double *clearini,
	     double *clearslavnode,int *itiefac,int *ipkon,int *kon,
             char *lakon,int *islavsurf));

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
	     int *nslavsm,int *nkon_,int *mcs,int *mortar,int *ifacecount,
             int *nintpoint));

void FORTRAN(allocont,(int *ncont,int *ntie,char *tieset,int *nset,
             char *set,int *istartset,int *iendset,int *ialset,
	     char *lakon,int *ncone,double *tietol,int *ismallsliding,
	     char *kind1,char *kind2,int *mortar,int *istep));

void FORTRAN(applyboun,(int *ifaext,int *nfaext,int *ielfa,int *ikboun,
             int *ilboun,int *nboun,char *typeboun,int *nelemload,
             int *nload,char *sideload,int *isolidsurf,int *nsolidsurf,
             int *ifabou,int *nfabou));

void arpack(double *co,int *nk,int **konp,int **ipkonp,char **lakonp,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     int *nactdof, 
	     int *icol,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *shcon,int *nshcon,double *cocon,int *ncocon,
	     double *alcon,int *nalcon,double *alzero,int **ielmatp,
	     int **ielorienp,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old,
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs,   
	     int *kode,int *mei,double *fei,
	     char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double **xstatep,int *npmat_,char *matname,int *mi,
             int *ncmat_,int *nstate_,double **enerp,char *jobnamec,
             char *output,char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     char *cbody,int *ibody,double *xbody,int *nbody,double *thicke,
	     int *nslavs,double *tietol,int *nkon,int *mpcinfo,int *ntie,
	     int *istep,int *mcs,int *ics,char *tieset,
             double *cs,
	     int *nintpoint,int *mortar,int *ifacecount,
	     int **islavsurfp,double **pslavsurfp,double **clearinip);

void arpackbu(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, 
	     int *nactdof, 
	     int *icol,int *jq,int *irow,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs, 
	     int *kode,int *mei,double *fei,
             char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double *xstate,int *npmat_,char *matname,int *mi,
             int *ncmat_,int *nstate_,double *ener,char *output, 
             char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     char *cbody,int *ibody,double *xbody,int *nbody,
	     double *thicke,char *jobnamec);

void arpackcs(double *co,int *nk,int **konp,int **ipkonp,char **lakonp,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, int *nactdof, 
	     int *icol,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun, 
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int **ielmatp,
	     int **ielorienp,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old,
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs,  
	     int *kode,int *mei,double *fei,
	     char *filab,double *eme,
             int *iexpl,double *plicon,int *nplicon,double *plkcon,
             int *nplkcon,
             double **xstatep,int *npmat_,char *matname,int *mi,
             int *ics,double *cs,int *mpcend,int *ncmat_,int *nstate_,
             int *mcs,int *nkon,double **enerp,char *jobnamec,
             char *output,char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *isolver,double *trab, 
             int *inotr,int *ntrans,double *ttime,double *fmpc,
	     char *cbody,int *ibody,double *xbody,int *nbody, 
             int *nevtot,double *thicke, int *nslavs, double *tietol, 
	     int *mpcinfo,int *ntie,int *istep,
	     char *tieset,
	     int *nintpoint,int *mortar,int *ifacecount,
	     int **islavsurfp,double **pslavsurfp,double **clearinip);

void FORTRAN(assigndomtonodes,(int *ne,char *lakon,int *ipkon,int *kon,
             int *ielmat,int *inomat,double *elcon,int *ncmat_,int *ntmat_,
             int *mi));

void FORTRAN(basis,(double *x,double *y,double *z,double *xo,double *yo,
                    double *zo,int *nx,int *ny,int *nz,double *planfa,
                    int *ifatet,int *nktet,int *netet,double *field,
                    int *nfield,double *cotet,int *kontyp,int *ipkon,
                    int *kon,int *iparent,double *xp,double *yp,double *zp,
                    double *value,double *ratio,int *iselect,int *nselect,
                    int *istartset,int *iendset,int *ialset,int *imastset,
                    int *ielemnr,int *nterms,int *konl));

void biosav(int *ipkon,int *kon,char *lakon,int *ne,double *co,
	    double *qfx,double *h0,int *mi,int *inomat,int *nk);

void FORTRAN(biotsavart,(int *ipkon,int *kon,char *lakon,int *ne,double *co,
                         double *qfx,double *h0,int *mi,int *nka,int *nkb));

void *biotsavartmt(int *i);

void FORTRAN(bodyforce,(char *cbody,int *ibody,int *ipobody,int *nbody,
             char *set,int *istartset,int *iendset,int *ialset,
             int *inewton,int *nset,int *ifreebody,int *k));

void FORTRAN(calcgradv,(int *ne,char *lakon,int *ipnei,double *vfa,
			double *area,double *xxn,double *gradv,int *neifa));
		      
void FORTRAN(calcmac,(int *neq,double *z,double *zz,int *nev,double *mac,
		      double* maccpx,int *istartnmd,int *iendnmd,int *nmd,
		      int *cyclicsymmetry,int *neqact,double *bett,
		      double *betm));

void FORTRAN(calcmass,(int *ipkon,char *lakon,int *kon,double *co,int *mi,
             int *nelem,int *ne,double *thicke,int *ielmat,
             int *nope,double *t0,double *t1,double *rhcon,
             int *nrhcon,int *ntmat_,int *ithermal,double *csmass));

void calcresidual(int *nmethod,int *neq,double *b,double *fext,double *f,
        int *iexpl,int *nactdof,double *aux1,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,int *nk,double *adb,
        double *aub,int *icol,int *irow,int *nzl,double *alpha,
	double *fextini,double *fini,int *islavnode,int *nslavnode,
        int *mortar,int *ntie,
        double *f_cm,double *f_cs,int *mi,int *nzs,int *nasym);

void calcresidual_em(int *nmethod,int *neq,double *b,double *fext,double *f,
        int *iexpl,int *nactdof,double *aux1,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,int *nk,double *adb,
        double *aub,int *icol,int *irow,int *nzl,double *alpha,
	double *fextini,double *fini,int *islavnode,int *nslavnode,
        int *mortar,int *ntie,
	double *f_cm,double *f_cs,int *mi,int *nzs,int *nasym,int *ithermal);

void FORTRAN(calcdvifa,(int *nface,double *vfa,double *shcon,
			int *nshcon,int *ielmat,int *ntmat_,int *ithermal,
                        int *mi,int *ielfa,double *umfa));

void FORTRAN(calcrhoel,(int *ne,char *lakon,double *vel,double *rhcon,
			int *nrhcon,int *ielmat,int *ntmat_,int *ithermal,
                        int *mi));

void FORTRAN(calcrhofa,(int *nface,double *vfa,double *rhcon,
			int *nrhcon,int *ielmat,int *ntmat_,int *ithermal,
                        int *mi,int *ielfa));

void FORTRAN(calcview,(char *sideload,double *vold,double *co,
             double *pmid,double *e1,double *e2,double *e3,
	     int *kontri,int *nloadtr,double *adview,double *auview,
             double *dist,int *idist,double *area,int *ntrit,int *mi,int *jqrad,
	     int *irowrad,int *nzsrad,double *sidemean,int *ntria, 
             int *ntrib));

void *calcviewmt(int *i);

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
               double *xstate,char *jobnamec,int *irstrt,
               double *ttime,double *qaold,
               char *output,char *typeboun,char *inpc,
               int *ipoinp,int *inp,char *tieset,double *tietol, 
               int *ntie,double *fmpc,char *cbody,int *ibody,double *xbody,
               int *nbody,int *nbody_,double *xbodyold,int *nam_,
               int *ielprop,int *nprop,int *nprop_,double *prop,
	       int *itpamp,int *iviewfile,int *ipoinpc,int *cfd,
	       int *nslavs,double *t0g,double *t1g,int *network,
	       int *cyclicsymmetry,int *idefforc,int *idefload,
               int *idefbody,int *mortar,int *ifacecount,int *islavsurf,
	       double *pslavsurf,double *clearini));    

void cascade(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
   int *mpcfree,int *nodeboun,int *ndirboun,int*nboun,int*ikmpc,
   int *ilmpc,int *ikboun,int *ilboun,int *mpcend,int *mpcmult,
   char *labmpc,int *nk,int *memmpc_,int *icascade,int *maxlenmpc,
   int *callfrommain,int *iperturb,int *ithermal);

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
	  int *istartset,int *iendset,int *ialset,double *emn,double *thicke,
	  char *jobnamec,int *mortar);

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

void FORTRAN(closefilefluid,());

void compfluid(double **cop,int *nk,int **ipkonp,int **konp,char **lakonp,
    int *ne,char **sideface,int *ifreestream, 
    int *nfreestream,int *isolidsurf,int *neighsolidsurf,
    int *nsolidsurf,int **iponoel,int **inoel,int *nshcon,double *shcon,
    int *nrhcon,double *rhcon,double **voldp,int *ntmat_,int *nodeboun, 
    int *ndirboun,int *nboun,int **ipompcp,int **nodempcp,int *nmpc,
    int **ikmpcp,int **ilmpcp,int *ithermal,int *ikboun,int *ilboun,
    int *turbulent,int *isolver,int *iexpl,double *vcontu,double *ttime,
    double *time,double *dtime,int *nodeforc,int *ndirforc,double *xforc,
    int *nforc,int *nelemload,char *sideload,double *xload,int *nload,
    double *xbody,int *ipobody,int *nbody,int **ielmatp,char *matname,
    int *mi,int *ncmat_,double *physcon,int *istep,int *iinc,
    int *ibody,double *xloadold,double *xboun,
    double **coefmpcp,int *nmethod,double *xforcold,double *xforcact,
    int *iamforc,int *iamload,double *xbodyold,double *xbodyact,
    double *t1old,double *t1,double *t1act,int *iamt1,double *amta,
    int *namta,int *nam,double *ampli,double *xbounold,double *xbounact,
    int *iamboun,int *itg,int *ntg,char *amname,double *t0,int **nelemface,
    int *nface,double *cocon,int *ncocon,double *xloadact,double *tper,
    int *jmax,int *jout,char *set,int *nset,int *istartset,
    int *iendset,int *ialset,char *prset,char *prlab,int *nprint,
    double *trab,int *inotr,int *ntrans,char *filab,char **labmpcp,
    double *sti,int *norien,double *orab,char *jobnamef,char *tieset,
    int *ntie,int *mcs,int *ics,double *cs,int *nkon,int *mpcfree,
    int *memmpc_,double **fmpcp,int *nef,int **inomat,double *qfx,
    int *neifa,int *neiel,int *ielfa,int *ifaext,double *vfa,double *vel,
    int *ipnei,int *nflnei,int *nfaext,char *typeboun);

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
	       int *ikforc,int *ilforc,double *thicke,
	       char *jobnamef,int *mei);

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
             double *xmastnor,char *filab,int *mcs,
             int *ics,int *nasym,double *xnoels,int *mortar,
             double *pslavsurf,double *pmastsurf,double *clearini,
             double *theta);
    
void FORTRAN(coriolissolve,(double *cc,int *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
             int *iter,double *d,double *temp));

void FORTRAN(createinterfacempcs,(int *imastnode,double *xmastnor,
	     int *nmastnode,int *ikmpc,int *ilmpc,int *nmpc,int *ipompc,
             int *nodempc,double *coefmpc,char *labmpc,int *mpcfree,
             int *ikboun,int *nboun));

void FORTRAN(createinum,(int *ipkon,int *inum,int *kon,char *lakon,int *nk,
             int *ne,char *cflag,int *nelemload,int *nload,int *nodeboun,
             int *nboun,int *ndirboun,int *ithermal,double *co,
             double *vold,int *mi));

void FORTRAN(createmddof,(int *imddof,int *nmddof,int *istartset,
       int *iendset,int *ialset,int *nactdof,int *ithermal,int *mi,
       int *imdnode,int *nmdnode,int *ikmpc,int *ilmpc,int *ipompc,
       int *nodempc,int *nmpc,int *imdmpc,
       int *nmdmpc,int *imdboun,int *nmdboun,int *ikboun,int *nboun,
       int *nset,int *ntie,char *tieset,char *set,char *lakon,int *kon,
       int *ipkon,char *labmpc,int *ilboun,char *filab,char *prlab,
       char *prset,int *nprint,int *ne,int *cyclicsymmetry));

void FORTRAN(createmdelem,(int *imdnode,int *nmdnode,double *xforc,
             int *ikmpc,int *ilmpc,int *ipompc,int *nodempc,int *nmpc,
             int *imddof,int *nmddof,int *nactdof,int *mi,int *imdmpc,
             int *nmdmpc,int *imdboun,int *nmdboun,int *ikboun,int *nboun,
             int *ilboun,int *ithermal,int *imdelem,int *nmdelem,
             int *iponoel,int *inoel,char *prlab,char *prset,int *nprint,
             char *lakon,char *set,int *nset,int *ialset,int *ipkon,
             int *kon,int *istartset,int *iendset,int *nforc,
             int *ikforc,int *ilforc));

void FORTRAN(createtiedsurfs,(int *nodface,int *ipoface,char *set,
             int *istartset,int *iendset,int *ialset,char *tieset,
             int *inomat,int *ne,int *ipkon,char *lakon,int *kon,
	     int *ntie,double *tietol,int *nalset,int *nk,int *nset,
             int *iactive));

void FORTRAN(dattime,(char *date,char *clock));

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

void FORTRAN(dnaupd,(int *ido,char *bmat,int *n,char *which,int *nev,
	     double *tol,double *resid,int *ncv,double *z,int *ldz,
	     int *iparam,int *ipntr,double *workd,double *workl,
	     int *lworkl,int *info));

void FORTRAN(dsaupd,(int *ido,char *bmat,int *n,char *which,int *nev,
	     double *tol,double *resid,int *ncv,double *z,int *ldz,
	     int *iparam,int *ipntr,double *workd,double *workl,
	     int *lworkl,int *info));

void FORTRAN(dneupd,(int *rvec,char *howmny,int *select,double *d,
	     double *di,double *z,int *ldz,double *sigma,double *sigmai,
             double *workev,char *bmat,int *neq,char *which, 
	     int *nev,double *tol,double *resid,int *ncv,double *v,
	     int *ldv,int *iparam,int *ipntr,double *workd,
	     double *workl,int *lworkl,int *info));

void FORTRAN(dseupd,(int *rvec,char *howmny,int *select,double *d,double *z,
	     int *ldz,double *sigma,char *bmat,int *neq,char *which, 
	     int *nev,double *tol,double *resid,int *ncv,double *v,
	     int *ldv,int *iparam,int *ipntr,double *workd,
	     double *workl,int *lworkl,int *info));

void FORTRAN(dsort,(double *dx,int *iy,int *n,int *kflag));

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
               double **xstatep,int *npmat_,char *matname,int *mi,
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
	       int *ikforc,int *ilforc,double *thicke,
               int *nslavs);

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
              int *itp,int *inext,
              int *imastop,int *nslavnode,int *islavnode,
              int *islavsurf,
              int *itiefac,double *areaslav,int *iponoels,int *inoels,
              double *springarea,int *izdof,int *nzdof,double *fn,
	      int *imastnode,int *nmastnode,double *xmastnor,
              double *xstateini,int *nslavs,
              int *cyclicsymmetry,double *xnoels,int *ielas);
 
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
	     int *nmdboun,int *imdmpc,int *nmdmpc,int *nmethod,double *time));

void electromagnetics(double **co,int *nk,int **konp,int **ipkonp,char **lakonp,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, int *nactdof, 
	     int **icolp,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int **ikmpcp,int **ilmpcp,int *ikboun, 
	     int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int **ielmatp,
	     int **ielorienp,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double **vold,int *iperturb,double *sti,int *nzs,  
	     int *kode, char *filab,int *idrct,
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
             char *output,
             double *shcon,int *nshcon,double *cocon,int *ncocon,
             double *physcon,int *nflow,double *ctrl, 
             char **setp,int *nset,int **istartsetp,
             int **iendsetp,int **ialsetp,int *nprint,char *prlab,
             char *prset,int *nener,int *ikforc,int *ilforc,double *trab, 
             int *inotr,int *ntrans,double **fmpcp,char *cbody,
             int *ibody,double *xbody,int *nbody,double *xbodyold,
             int *ielprop,double *prop,int *ntie,char **tiesetp,
	     int *itpamp,int *iviewfile,char *jobnamec,double **tietolp,
	     int *nslavs,double *thicke,int *ics,int *nalset,int *nmpc_);

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
                      int *nodempc,int *ipompc,char *labmpc,int *ikboun,
                      int *nasym));

void FORTRAN(equationcheck,(double *ac,int *nteq,int *nactdog,
                            int *itg,int *ntg,int *nacteq,int *network));

void FORTRAN(errorestimator,(double *yi,double *yn,int *ipkon,int *inum,
             int *kon,char *lakon,int *nk,int *ne,int *mi,int *ielmat,
	     double *thicke,int *nterms));

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
	     char* tieset,int* ntie,int *imddof,int *nmddof,
	     int *imdnode,int *nmdnode,int *imdboun,int *nmdboun,
             int *imdmpc,int *nmdmpc,int **izdofp,int *nzdof,int *nherm,
	     double *xmr,double *xmi);

void FORTRAN(extrapolate,(double *sti,double *stn,int *ipkon,int *inum,
             int *kon,char *lakon,int *nfield,int *nk,int *ne,int *mi,
             int *ndim,double *orab,int *ielorien,double *co,int *iorienglob,
	     char *cflag,int *nelemload,int *nload,int *nodeboun,int *nboun,
             int *ndirboun,double *vold,int *ithermal,int *force,
	     int *cfd,int *ielmat,double *thicke,char *filab));

void FORTRAN(extrapolate_gradv,(int *nface,int *ielfa,double *xrlfa,
			       double *gradv,double *gradvfa));

void FORTRAN(extrapolate_v,(int *nface,int *ielfa,double *xrlfa,double *vel,
             double *vfa,int *ifabou,double *xboun,double *sumfix,
             double *sumfree,double *xxn,double *area));

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
	      double *xloadact,int *nactdog,int *network,int *mi,
	      int *ne,int *ipkon,char *lakon,int *kon));

void FORTRAN(forcesolve,(double *zc,int *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
	     int *iter,double *d,int *neq,double *z,int *istartnmd,
	     int *iendnmd,int *nmd,int *cyclicsymmetry,int *neqact,
	     int *igeneralizedforce));

void frd(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne0,
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
	 double *thicke,char *jobnamec, char *output,double *qfx,
         double *cdn,int *mortar,double *cdnr,double *cdni);

void frdcyc(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne,double *v,
	    double *stn,int *inum,int *nmethod,int *kode,char *filab,
	    double *een,double *t1,double *fn,double *time,double *epn,
	    int *ielmat,char *matname,double *cs,int *mcs,int *nkon,
	    double *enern,double *xstaten,int *nstate_,int *istep,
            int *iinc,int *iperturb,double *ener,int *mi,char *output,
            int *ithermal,double *qfn,int *ialset,int *istartset,
            int *iendset,double *trab,int *inotr,int *ntrans,double *orab,
	    int *ielorien,int *norien,double *sti,double *veold,int *noddiam,
            char *set,int *nset,double *emn, double *thicke,char *jobnamec,
            int *ne0,double *cdn,int *mortar);

void FORTRAN(frdfluid,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *vold,int *kode,double *time,
             int *ielmat,char *matname,int *nnstep,double *vtu,
	     double *vcontu,double *voldaux,double *physcon,char *filab,
	     int *inomat,int *ntrans,int *inotr,double *trab,int *mi,
	     double *stn,double *qfn,int *istep));

void frdheader(int *icounter,double *oner,double *time,double *pi,
	       int *noddiam,double *cs,int *null,int *mode,
	       int *noutloc,char *description,int *kode,int *nmethod,
               FILE *f1,char *output,int *istep,int *iinc);

void FORTRAN(frditeration,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *time,int *ielmat,char *matname,
	     int *mi,int *istep,int *iinc,int *ithermal));

void frdselect(double *field1,double *field2,int *iset,int *nkcoords,int *inum,
     char *m1,int *istartset,int *iendset,int *ialset,int *ngraph,int *ncomp,
     int *ifield,int *icomp,int *nfield,int *iselect,char *m2,FILE *f1,
     char *output, char *m3);

void frdset(char *filabl,char *set,int *iset,int *istartset,int *iendset,
	    int *ialset,int *inum,int *noutloc,int *nout,int *nset,
	    int *noutmin,int *noutplus,int *iselect,int *ngraph);

void frdvector(double *v,int *iset,int *ntrans,char * filabl,int *nkcoords,
               int *inum,char *m1,int *inotr,double *trab,double *co,
               int *istartset,int *iendset,int *ialset,int *mi,int *ngraph,
               FILE *f1,char *output,char *m3);

void FORTRAN(fsub,(double *time,double *tend,double *aai,double *bbi,
		   double *ddj,double *h1,double *h2,double *h3,double *h4,
                   double *func,double *funcp));

void FORTRAN(fsuper,(double *time,double *tend,double *aai,double *bbi,
		       double *h1,double *h2,double *h3,double *h4,
		       double *h5,double *h6,double *func,double *funcp));

void FORTRAN(gasmechbc,(double *vold,int *nload,char *sideload,
			int *nelemload,double *xload,int *mi));

void FORTRAN(genadvecelem,(int *inodesd,int *ipkon,int *ne,char *lakon,
             int *kon,int *nload,char *sideload,int *nelemload,int *nkon));

void FORTRAN(gencontelem_f2f,(char *tieset,int *ntie,int *itietri,int *ne,
             int *ipkon,int *kon,char *lakon,double *cg,double *straight,
             int *ifree,int *koncont,double *co,double *vold,double *xo,
             double *yo,double *zo,double *x,double *y,double *z,int *nx,
             int *ny,int *nz,int *ielmat,double *elcon,int *istep,int *iinc,
             int *iit,int *ncmat_,int *ntmat_,int *mi,int *imastop,
             int *islavsurf,int *itiefac,double *springarea,double *tietol,
             double *reltime,char *filab,int *nasym,
	     double *pslavsurf,double *pmastsurf,double *clearini,
             double *theta));

void FORTRAN(gencontelem_n2f,(char *tieset,int *ntie,int *itietri,int *ne,
     int *ipkon,int *kon,char *lakon,
     double *cg,double *straight,int *ifree,int *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,int *nx,int *ny,int *nz,
     int *ielmat,double *elcon,int *istep,int *iinc,int *iit,
     int *ncmat_,int *ntmat_,
     int *nmethod,int *mi,int *imastop,int *nslavnode,
     int *islavnode,int *islavsurf,int *itiefac,double *areaslav,
     int *iponoels,int *inoels,double *springarea,
     char *set,int *nset,int *istartset,int *iendset,int *ialset,
     double *tietol,double *reltime,
     int *imastnode,int *nmastnode,char* filab,int *nasym,double *xnoels));

void FORTRAN(generateeminterfaces,(int *istartset,int *iendset,
	     int *ialset,int *iactive,int *ipkon,char *lakon,int *kon,
	     int *ikmpc,int *nmpc,int *nafaces));

void FORTRAN(generatetet,(int *kontet,int *ifatet,int *netet,
             int *inodfa,int *ifreefa,double *planfa,int *ipofa,
             int *nodes,double *cotet));

void  FORTRAN(gennactdofinv,(int *nactdof,int *nactdofinv,int *nk,
       int *mi,int *nodorig,int *ipkon,char *lakon,int *kon,int *ne));

void FORTRAN(gentiedmpc,(char *tieset,int *ntie,int *itietri,
          int *ipkon,int *kon,char *lakon,char *set,int *istartset,
          int *iendset,int *ialset,double *cg,double *straight,
	  int *koncont,double *co,double *xo,double *yo,double *zo,
          double *x,double *y,double *z,int *nx,int *ny,int *nz,int *nset,
          int *ifaceslave,int *istartfield,int *iendfield,int *ifield,
          int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *nmpc_,
          int *mpcfree,int *ikmpc,int *ilmpc,char *labmpc,int *ithermal,
	  double *tietol,int *icfd,int *ncont,int *imastop,int *ikboun,
	  int *nboun,char *kind));

void FORTRAN(geomview,(double *vold,double *co,double *pmid,double *e1,
             double *e2,double *e3,int *kontri,double *area,double *cs,
             int *mcs,int *inocs,int *ntrit,int *nk,int *mi,double *sidemean));

void getglobalresults (char *jobnamec,int **integerglobp,double **doubleglobp,
                       int *nboun,int *iamboun,double *xboun, int *nload,
                       char *sideload,int *iamload,int *iglob);

int getSystemCPUs();;

void FORTRAN(identamta,(double *amta,double *reftime,int *istart,int *iend,
               int *id));

void FORTRAN(identifytiedface,(char *tieset,int *ntie,char *set,int *nset,
			       int *faceslave,char *kind));

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
             int *iperturb,int *ikboun,int *nboun,double *co,int *istep,
             double **xnoelsp);

void FORTRAN(init,(int *nktet,int *inodfa,int *ipofa,int *netet_));

void FORTRAN(initialcfd,(int *ne,int *ipkon,int *kon,char *lakon,
             double *co,double *coel,double *cofa,int *nface,int *ielfa,
             double *area,int *ipnei,int *neiel,double *xxn,double *xxi,
             double *xle,double *xlen,double *xlet,double *xrlfa,double *cosa));

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

void insertrad(int *ipointer,int **mast1p,int **mast2p,int *i1,
	    int *i2,int *ifree,int *nzs_);

void FORTRAN(integral_boundary,(double *sumfix,double *sumfree,int *ifaext,
             int *nfaext,int *ielfa,int *ifabou,double *vfa));

void FORTRAN(interpolatestate,(int *ne,int *ipkon,int *kon,char *lakon,
             int *ne0,int *mi,double *xstate,
             double *pslavsurf,int *nstate_,double *xstateini,
             int *islavsurf,int *islavsurfold,
             double *pslavsurfold));

void FORTRAN(islavactive,(char *tieset,int *ntie,int *itietri,double *cg,
              double *straight,double *co,double *vold,double *xo,
              double *yo,double *zo,double *x,double *y,double *z,
              int *nx,int *ny,int *nz,int *mi,int *imastop,int *nslavnode,
              int *islavnode,int *islavact));

void FORTRAN(isortid,(int *ix,double *dy,int *n,int *kflag));

void FORTRAN(isortii,(int *ix,int *iy,int *n,int *kflag));

void FORTRAN(isortiid,(int *ix,int *iy,double *dy,int *n,int *kflag));

void FORTRAN(isortiddc,(int *ix,double *dy1,double *dy2,char *cy,int *n, 
                         int *kflag));

void FORTRAN(isortiiddc,(int *ix1,int *ix2,double *dy1,double *dy2, 
                         char *cy,int *n,int *kflag));

void FORTRAN(keystart,(int *ifreeinp,int *ipoinp,int *inp,char *name,
           int *iline,int *ikey));
  
void linstatic(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int *ipompc,int *nodempc,double *coefmpc,char *labmpc,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload, int *nactdof, 
	     int **icolp,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int *ikmpc,int *ilmpc,int *ikboun, 
	     int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int *ielmat,
	     int *ielorien,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double *vold,int *iperturb,double *sti,int *nzs, 
	     int *kode,char *filab,double *eme,
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
	     double *thicke, char *jobnamec,char *tieset,int *ntie,
             int *istep);

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

void FORTRAN(mafillem,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
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
	       double *xstate,double *thicke,
               int *integerglob,double *doubleglob,char *tieset,
	       int *istartset,int *iendset,int *ialset,int *ntie,
	       int *nasym,int *iactive,double *h0,double *pslavsurf,
               double *pmastsurf,int *mortar,double *clearini));

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
	       double *xstate,double *thicke,
               int *integerglob,double *doubleglob,char *tieset,
	       int *istartset,int *iendset,int *ialset,int *ntie,
	       int *nasym,double *pslavsurf,double *pmastsurf,int *mortar,
               double *clearini));

void FORTRAN(mafillsmas,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
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
	       double *xstate,double *thicke,
               int *integerglob,double *doubleglob,char *tieset,
	       int *istartset,int *iendset,int *ialset,int *ntie,
	       int *nasym,double *pslavsurf,double *pmastsurf,int *mortar,
               double *clearini));

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
	       double *springarea,double *thicke,
               int *integerglob,double *doubleglob,char *tieset,
	       int *istartset,int *iendset,int *ialset,int *ntie,
	       int *nasym,double *pslavsurf,double *pmastsurf,int *mortar,
               double *clearini));

void FORTRAN(mafillsmcsas,(double *co,int *nk,int *kon,int *ipkon, 
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
	       double *springarea,double *thicke,
               int *integerglob,double *doubleglob,char *tieset,
	       int *istartset,int *iendset,int *ialset,int *ntie,
	       int *nasym,int *nstate_,double *xstateini,double *xstate,
	       double *pslavsurf,double *pmastsurf,int *mortar,
               double *clearini));

void FORTRAN(mafillv,(int *ne,int *nactdoh,int *ipnei,int *neifa,int *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             int *jq,int *irow,int *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
             double *xxi,double *body,double *volume,int *compressible,
             int *ielfa,char *lakon,int *ifabou));

void mastruct(int *nk,int *kon,int *ipkon,char*lakon,int *ne,
	      int *nodeboun,int *ndirboun,int *nboun,int *ipompc,
	      int *nodempc,int *nmpc,int *nactdof,int *icol,
	      int *jq,int **mast1p,int **irowp,int *isolver,int *neq,
	      int *ikmpc,int *ilmpc,int *ipointer,int *nzs,int *nmethod,
              int *ithermal,int *ikboun,int *ilboun,int *iperturb,
              int *mi,int *mortar);

void mastructcs(int *nk,int *kon,int *ipkon,char *lakon,
	       int *ne,int *nodeboun,
	       int *ndirboun,int *nboun,int *ipompc,int *nodempc,
	       int *nmpc,int *nactdof,int *icol,int *jq,int **mast1p,
	       int **irowp,int *isolver,int *neq,
	       int *ikmpc,int *ilmpc,int *ipointer,
	       int *nzs,int *nmethod,int *ics,double *cs,
	       char *labmpc,int *mcs,int *mi,int *mortar);

void mastructem(int *nk, int *kon, int *ipkon, char *lakon, int *ne,
	      int *nodeboun, int *ndirboun, int *nboun, int *ipompc,
	      int *nodempc, int *nmpc, int *nactdof, int *icol,
	      int *jq, int **mast1p, int **irowp, int *isolver, int *neq,
	      int *ikmpc, int *ilmpc,int *ipointer, int *nzs, 
	      int *ithermal,int *mi,int *ielmat, double *elcon, int *ncmat_, 
	      int *ntmat_,int *inomat);

void mastructf(int *nk,int *kon,int *ipkon,char *lakon,int *ne,
	       int *nactdoh,int *icol,int *jq, int **mast1p, int **irowp,
	       int *isolver, int *neq,int *ipointer, int *nzs,
               int *ipnei,int *ineiel,int *mi);

void mastructrad(int *ntr,int *nloadtr,char *sideload,int *ipointerrad,
              int **mast1radp,int **irowradp,int *nzsrad,
	      int *jqrad,int *icolrad);

void FORTRAN(mpcrem,(int *i,int *mpcfree,int *nodempc,int *nmpc,int *ikmpc,
                     int *ilmpc,char *labmpc,double *coefmpc,int *ipompc));

void FORTRAN(mult,(double *matrix,double *trans,int *n));

void FORTRAN(networkinum,(int *ipkon,int *inum,int *kon,char *lakon,
       int *ne,int *itg,int *ntg));

void FORTRAN(nident,(int *x,int *px,int *n,int *id));

void FORTRAN(nidentll,(long long *x,long long *px,int *n,int *id));

void FORTRAN(nodestiedface,(char *tieset,int *ntie,int *ipkon,int *kon,
       char *lakon,char *set,int *istartset,int *iendset,int *ialset,
       int *nset,int *faceslave,int *istartfield,int *iendfield,
       int *ifield,int *nconf,int *ncone,char *kind));

void nonlingeo(double **co,int *nk,int **konp,int **ipkonp,char **lakonp,
	     int *ne, 
	     int *nodeboun,int *ndirboun,double *xboun,int *nboun, 
	     int **ipompcp,int **nodempcp,double **coefmpcp,char **labmpcp,
             int *nmpc, 
	     int *nodeforc,int *ndirforc,double *xforc,int *nforc, 
	     int *nelemload,char *sideload,double *xload,
	     int *nload,int *nactdof, 
	     int **icolp,int *jq,int **irowp,int *neq,int *nzl, 
	     int *nmethod,int **ikmpcp,int **ilmpcp,int *ikboun, 
	     int *ilboun,
	     double *elcon,int *nelcon,double *rhcon,int *nrhcon,
	     double *alcon,int *nalcon,double *alzero,int **ielmatp,
	     int **ielorienp,int *norien,double *orab,int *ntmat_,
	     double *t0,double *t1,double *t1old, 
	     int *ithermal,double *prestr,int *iprestr, 
	     double **vold,int *iperturb,double *sti,int *nzs,  
	     int *kode,char *filab,int *idrct,
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
             char *output,
             double *shcon,int *nshcon,double *cocon,int *ncocon,
             double *physcon,int *nflow,double *ctrl, 
             char *set,int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,int *nener,int *ikforc,int *ilforc,double *trab, 
             int *inotr,int *ntrans,double **fmpcp,char *cbody,
             int *ibody,double *xbody,int *nbody,double *xbodyold,
             int *ielprop,double *prop,int *ntie,char *tieset,
	     int *itpamp,int *iviewfile,char *jobnamec,double *tietol,
	     int *nslavs,double *thicke,int *ics,
	     int *nintpoint,int *mortar,int *ifacecount,char *typeboun,
	     int **islavsurfp,double **pslavsurfp,double **clearinip);

void FORTRAN(nonlinmpc,(double *co,double *vold,int *ipompc,int *nodempc,
		   double *coefmpc,char *labmpc,int *nmpc,int *ikboun,
		   int *ilboun,int *nboun,double *xbounact,double *aux,
		   int *iaux,int *maxlenmpc,int *ikmpc,int *ilmpc,
                   int *icascade,int *kon,int *ipkon,char *lakon,
		   int *ne,double *reltime,int *newstep,double *xboun,
		   double *fmpc,int *newinc,int *idiscon,int *ncont,
		   double *trab,int *ntrans,int *ithermal,int *mi));

void FORTRAN(normalsoninterface,(int *istartset,int *iendset,
	     int *ialset,int *imast,int *ipkon,int *kon,char *lakon,
             int *imastnode,int *nmastnode,double *xmastnor,double *co));

void FORTRAN(op,(int *n,double *x,double *y,double *ad,double *au,int *jq,int *irow));

void FORTRAN(opas,(int *n,double *x,double *y,double *ad,double *au,int *jq,
		   int *irow,int *nzs));

void FORTRAN(op_corio,(int *n,double *x,double *y,double *ad,double *au,
		       int *jq,int *irow));

void FORTRAN(openfile,(char *jobname,char *output));

void FORTRAN(openfilefluid,(char *jobname));

void FORTRAN(postview,(int *ntr,char *sideload,int *nelemload,int *kontri,
             int *ntri,int *nloadtr,double *tenv,double *adview,double *auview,
             double *area,double *fenv,int *jqrad,int *irowrad,int *nzsrad));

void FORTRAN(precfd,(int *ne,int *ipkon,int *kon,char *lakon,int *ipnei,
                     int *neifa,int *neiel,int *ipoface,int *nodface,
                     int *ielfa,int *nkonnei,int *nface,int *ifaext,
                     int *nfaext,int *isolidsurf,int *nsolidsurf,char *set,
                     int *nset,int *istartset,int *iendset,int *ialset,
                     double *vel,double *vfa,double *vold,int *mi));

void precontact(int *ncont, int *ntie, char *tieset, int *nset, char *set,
        int *istartset, int *iendset, int *ialset, int *itietri,
        char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
        double *cg, double *straight, double *co,double *vold,
        int *istep,int *iinc,int *iit,int *itiefac,
        int *islavsurf, int *islavnode, int *imastnode,
        int *nslavnode, int *nmastnode,int *imastop,int *mi,
	int *ipe, int *ime,double *tietol,int *iflagact,
	int *nintpoint,double **pslavsurfp,double *xmastnor,double *cs,
	int *mcs,int *ics,double *clearini,int *nslavs);

void prediction(double *uam,int *nmethod,double *bet,double *gam,double *dtime,
               int *ithermal,int *nk,double *veold,double *accold,double *v,
	       int *iinc,int *idiscon,double *vold,int *nactdof,int *mi);

void prediction_em(double *uam,int *nmethod,double *bet,double *gam,double *dtime,
               int *ithermal,int *nk,double *veold,double *v,
	       int *iinc,int *idiscon,double *vold,int *nactdof,int *mi);

void preiter(double *ad,double **aup,double *b,int **icolp,int **irowp, 
	     int *neq,int *nzs,int *isolver,int *iperturb);

void FORTRAN(printout,(char *set,int *nset,int *istartset,int *iendset,
             int *ialset,int *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,int *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,int *mi,int *nstate_,int *ithermal,
             double *co,int *kon,double *qfx,double *ttime,double *trab,
             int *inotr,int *ntrans,double *orab,int *ielorien,
	     int *norien,int *nk,int *ne,int *inum,char *filab,double *vold,
             int *ikin,int *ielmat,double *thicke,double *eme,int *islavsurf,
             int *mortar));

void FORTRAN(printoutfluid,(char *set,int *nset,int *istartset,int *iendset,
             int *ialset,int *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,int *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,int *mi,int *nstate_,int *ithermal,
             double *co,int *kon,double *qfx,double *ttime,double *trab,
             int *inotr,int *ntrans,double *orab,int *ielorien,
	     int *norien,int *nk,int *ne,int *inum,char *filab,double *vold,
             int *ikin,int *ielmat,double *thicke,double *eme,
             double *vcontu,double *physcon));

void FORTRAN(printoutface,(double *co,double *rhcon,int *nrhcon,int *ntmat_,
            double *vold,double *shcon,int *nshcon,double *cocon,
            int *ncocon,int *icompressible,int *istartset,int *iendset,
            int *ipkon,char *lakon,int *kon,int *ialset,char *prset,
	    double *timef,int *nset,char *set,int *nprint,char *prlab,
	    int *ielmat,int *mi));

int pthread_create (pthread_t *thread_id, const pthread_attr_t *attributes,
                    void *(*thread_function)(void *), void *arguments);

int pthread_join (pthread_t thread, void **status_ptr);

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
       int *iamload,int *jqrad,int *irowrad,int *nzsrad,int *icolrad,
       int *ne);

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

void FORTRAN(readforce,(double *zc,int *neq,int *nev,int *nactdof,
	     int *ikmpc,int *nmpc,int *ipompc,int *nodempc,int *mi, 
	     double *coefmpc,char *jobnamec,double *aa,
             int *igeneralizedforce));

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

void FORTRAN(rectcyltrfm,(int *node,double *co,double *cs,int *cntrl,
             double *fin,double *fout));

void FORTRAN(rectcylvi,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,int *nk, 
		      int *icntrl,double *t,char *filab,int *imag,int *mi,
                      double *emn));

void remastruct(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
              int *mpcfree,int *nodeboun,int *ndirboun,int *nboun,
              int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
              char *labmpc,int *nk,
              int *memmpc_,int *icascade,int *maxlenmpc,
              int *kon,int *ipkon,char *lakon,int *ne,
              int *nactdof,int *icol,int *jq,int **irowp,int *isolver,
              int *neq,int *nzs,int *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,int *ithermal,
	      int *iperturb,int *mass,int *mi,int *iexpl,int *mortar);

void remastructar(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
              int *mpcfree,int *nodeboun,int *ndirboun,int *nboun,
              int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
              char *labmpc,int *nk,
              int *memmpc_,int *icascade,int *maxlenmpc,
              int *kon,int *ipkon,char *lakon,int *ne,
              int *nactdof,int *icol,int *jq,int **irowp,int *isolver,
              int *neq,int *nzs,int *nmethod,int *ithermal,
	      int *iperturb,int *mass,int *mi,int *ics,double *cs,
	      int *mcs,int *mortar);

void remastructem(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
              int *mpcfree,int *nodeboun,int *ndirboun,int *nboun,
              int *ikmpc,int *ilmpc,int *ikboun,int *ilboun,
              char *labmpc,int *nk,
              int *memmpc_,int *icascade,int *maxlenmpc,
              int *kon,int *ipkon,char *lakon,int *ne,
              int *nactdof,int *icol,int *jq,int **irowp,int *isolver,
              int *neq,int *nzs,int *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,int *ithermal,
	      int *iperturb,int *mass,int *mi,int *ielmat,double *elcon,
	      int *ncmat_,int *ntmat_,int *inomat);

void FORTRAN(restartshort,(int *nset,int *nload,int *nbody,int *nforc,
    int *nboun,
    int *nk,int *ne,int *nmpc,int *nalset,int *nmat,int *ntmat,int *npmat,
    int *norien,int *nam,int *nprint,int *mint,int *ntrans,int *ncs,
    int *namtot,int *ncmat,int *memmpc,int *ne1d,int *ne2d,int *nflow,
    char *set,int *meminset,int *rmeminset,char *jobnamec,int *irestartstep,
    int *icntrl,int *ithermal,int *nener,int *nstate_,int *ntie,int *nslavs,
    int *nkon,int *mcs,int *nprop,int *mortar,int *ifacecount,int *nintpoint));

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
  double *xnor,int *knor,double *thicke,double *offset, 
  int *iponoel,int *inoel,int *rig, 
  double *shcon,int *nshcon,double *cocon,int *ncocon, 
  int *ics,double *sti,double *ener,double *xstate, 
  char *jobnamec,int *infree,double *prestr,int *iprestr,
  char *cbody,int *ibody,double *xbody,int *nbody,double *xbodyold,
  double *ttime,double *qaold,double *cs,
  int *mcs,char *output,double *physcon,double *ctrl,char *typeboun,
  double *fmpc,char *tieset,int *ntie,double *tietol,int *nslavs,
  double *t0g,double *t1g,int *nprop,int *ielprop,double *prop,int *mortar,
  int *nintpoint,int *ifacecount,int *islavsurf,double *pslavsurf,
  double *clearini));

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
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,int *ncocon,char *set, 
             int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             int *inotr,int *ntrans,double *fmpc,int *nelemload,
	     int *nload,int *ikmpc,int *ilmpc,int *istep,int *iinc,
	     double *springarea,double *reltime,int *ne0,double *xforc,
             int *nforc,double *thicke,
             double *shcon,int *nshcon,char *sideload,double *xload,
             double *xloadold,int *icfd,int *inomat,double *pslavsurf,
             double *pmastsurf,int *mortar,int *islavact,double *cdn,
             int *islavnode,int *nslavnode,int *ntie,double *clearini,
             int *islavsurf);

void FORTRAN(resultsem,(double *co,int *kon,int *ipkon,char *lakon,
             double *v,double *elcon,int *nelcon,int *ielmat,int *ntmat_,
             double *vold,double *dtime,char *matname,int *mi,int *ncmat_,
             int *nea,int *neb,double *sti,double *alcon,
             int *nalcon,double *h0));

void *resultsemmt(int *i);

void  FORTRAN(resultsforc,(int *nk,double *f,double *fn,int *nactdof,
       int *ipompc,int *nodempc,double *coefmpc,char *labmpc,int *nmpc,
       int *mi,double *fmpc,int *calcul_fn,
       int *calcul_f));

void  FORTRAN(resultsforcem,(int *nk,double *f,double *fn,int *nactdof,
       int *ipompc,int *nodempc,double *coefmpc,char *labmpc,int *nmpc,
       int *mi,double *fmpc,int *calcul_fn,int *calcul_f,int *inomat));

void FORTRAN(resultsini,(int *nk,double *v,int *ithermal,char *filab,
       int *iperturb,double *f,double *fn,int *nactdof,int *iout,
       double *qa,double *vold,double *b,int *nodeboun,int *ndirboun,
       double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
       char *labmpc,int *nmpc,int *nmethod,double *cam,int *neq,
       double *veold,double *accold,double *bet,double *gam,double *dtime,
       int *mi,double *vini,int *nprint,char *prlab,int *intpointvar,
       int *calcul_fn,int *calcul_f,int *calcul_qa,int *calcul_cauchy,
       int *iener,int *ikin,int *intpointvart,double *xforc,int *nforc));

void FORTRAN(resultsini_em,(int *nk,double *v,int *ithermal,char *filab,
       int *iperturb,double *f,double *fn,int *nactdof,int *iout,
       double *qa,double *vold,double *b,int *nodeboun,int *ndirboun,
       double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
       char *labmpc,int *nmpc,int *nmethod,double *cam,int *neq,
       double *veold,double *dtime,
       int *mi,double *vini,int *nprint,char *prlab,int *intpointvar,
       int *calcul_fn,int *calcul_f,int *calcul_qa,int *calcul_cauchy,
       int *iener,int *ikin,int *intpointvart,double *xforc,int *nforc));

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
	  double *emeini,double *pslavsurf,
	  double *pmastsurf,int *mortar,double *clearini,int *nea,int *neb));

void *resultsmechmt(int *i);

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
       double *thicke,double *eme,double *emn,double *rhcon,int *nrhcon,
       double *shcon,int *nshcon,double *cocon,int *ncocon,int *ntmat_,
       char *sideload,int *icfd,int *inomat,double *pslavsurf,
       int *islavact,double *cdn,int *mortar,int *islavnode,int *nslavnode,
       int *ntie,int *islavsurf,double *time));

void FORTRAN(resultstherm,(double *co,int *kon,int *ipkon,
       char *lakon,double *v,
       double *elcon,int *nelcon,double *rhcon,int *nrhcon,int *ielmat,
       int *ielorien,int *norien,double *orab,int *ntmat_,double *t0,
       int *iperturb,double *fn,double *shcon,int *nshcon,int *iout,
       double *qa,double *vold,int *ipompc,int *nodempc,
       double *coefmpc,int *nmpc,double *dtime,
       double *time,double *ttime,double *plkcon,int *nplkcon,double *xstateini,
       double *xstiff,double *xstate,int *npmat_,char *matname,
       int *mi,int *ncmat_,int *nstate_,double *cocon,int *ncocon,
       double *qfx,int *ikmpc,int *ilmpc,int *istep,
       int *iinc,double *springarea,int *calcul_fn,int *calcul_qa,int *nal,
       int *nea,int *neb,int *ithermal,int *nelemload,int *nload,
       int *nmethod,double *reltime,char *sideload,double *xload,
       double *xloadold,double *pslavsurf,double *pmastsurf,int *mortar,
       double *clearini,double *plicon,int *nplicon));

void *resultsthermemmt(int *i);

void *resultsthermmt(int *i);

void resultsinduction(double *co,int *nk,int *kon,int *ipkon,char *lakon,
	     int *ne,double *v,double *stn,int *inum, 
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
	     int *icmd,int *ncmat_,int *nstate_,double *sti,
	     double *vini,int *ikboun,int *ilboun,double *ener,
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,int *ncocon,char *set, 
             int *nset,int *istartset,
             int *iendset,int *ialset,int *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             int *inotr,int *ntrans,double *fmpc,int *nelemload,
	     int *nload,int *ikmpc,int *ilmpc,int *istep,int *iinc,
	     double *springarea,double *reltime,int *ne0,double *xforc,
             int *nforc,double *thicke,
             double *shcon,int *nshcon,char *sideload,double *xload,
	     double *xloadold,int *icfd,int *inomat,double *h0,
	     int *islavnode,int *nslavnode,int *ntie);

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
       
void FORTRAN(slavintpoints,(int *ntie,int *itietri,int *ipkon,
        int *kon,char *lakon,double *straight,
        int *nintpoint,int *koncont,double *co,double *vold,double *xo,
        double *yo,double *zo,double *x,double *y,double *z,int *nx,
        int *ny,int *nz,int *islavsurf,
        int *islavnode,int *nslavnode,int *imastop,
	int *mi,int *ncont,int *ipe,int *ime,double *pslavsurf,
        int *i,int *l,int *ntri));

void FORTRAN(sortev,(int *nev,int *nmd,double *eigxx,int *cyclicsymmetry,
		     double *xx,double *eigxr,int *pev,
		     int *istartnmd,int *iendnmd,double *aa,double *bb));

void FORTRAN(spcmatch,(double *xboun,int *nodeboun,int *ndirboun,int *nboun,
	       double *xbounold,int *nodebounold,int *ndirbounold,
	       int *nbounold,int *ikboun,int *ilboun,double *vold,
	       double *reorder,int *nreorder,int *mi));

void FORTRAN(splitline,(char *text,char *textpart,int *n));

void spooles(double *ad,double *au,double *adb,double *aub,
             double *sigma,double *b,
	     int *icol,int *irow,int *neq,int *nzs,int *symmtryflag,
             int *inputformat,int *nzs3);

void FORTRAN(springforc_n2f,(double *xl,int *konl,double *vl,int *imat,
             double *elcon,int *nelcon,double *elas,double *fnl,int *ncmat_,
             int *ntmat_,int *nope,char *lakonl,double *t1l,int *kode,
             double *elconloc,double *plicon,int *nplicon,int *npmat_,
             double *senergy,int *iener,double *cstr,int *mi,
             double *springarea,int *nmethod,int *ne0,int *nstate_,
	     double *xstateini,double *xstate,double *reltime,int *ielas));

void FORTRAN(springstiff_n2f,(double *xl,double *elas,int *konl,double *voldl,
             double *s,int *imat,double *elcon,int *nelcon,int *ncmat_,
             int *ntmat_,int *nope,char *lakonl,double *t1l,int *kode,
             double *elconloc,double *plicon,int *nplicon,int *npmat_,
             int *iperturb,double *springarea,int *nmethod,int *mi,int *ne0,
             int *nstate_,double *xstateini,double *xstate,double *reltime,
             int *nasym));

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
	  double **veoldp,char *amname,double *amta,
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
	  int *ics,double *cs,int *mpcend,double *ctrl,
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
             double *co,int *itg,int *ntg,double *vold,
	     int *ielmat, double *thicke));

int strcmp1(const char *s1, const char *s2);

int strcmp2(const char *s1, const char *s2,int length);

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
               int *inotr,double *veold,int *integerglob,
               double *doubleglob,char *tieset,int *istartset,
               int *iendset,int *ialset,int *ntie,int *nmpc,int *ipompc,
               int *ikmpc,int *ilmpc,int *nodempc,double *coefmpc));

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
       int *mortar,int *ipoface,int *nodface,int *nk,double *xnoels));   

void tiedcontact(int *ntie,char *tieset,int *nset,char *set,
               int *istartset,int *iendset,int *ialset,
               char *lakon,int *ipkon,int *kon,double *tietol,
               int *nmpc,int *mpcfree,int *memmpc_,
               int **ipompcp,char **labmpcp,int **ikmpcp,int **ilmpcp,
               double **fmpcp,int **nodempcp,double **coefmpcp,
	       int *ithermal,double *co,double *vold,int *cfd,
	       int *nmpc_,int *mi,int *nk,int *istep,int *ikboun,
	       int *nboun,char *kind1,char *kind2);
	
void FORTRAN(transformatrix,(double *xab,double *p,double *a));

void FORTRAN(trianeighbor,(int *ipe,int *ime,int *imastop,int *ncont,
               int *koncont,int *ifreeme));

void FORTRAN(triangucont,(int *ncont,int *ntie,char *tieset,int *nset,
          char *set,int *istartset,int *iendset,int *ialset,int *itietri,
	  char *lakon,int *ipkon,int *kon,int *koncont,char *kind1,
	  char *kind2,double *co,int *nk));

#ifdef BAM
void FORTRAN(uexternaldb,(int *lop,int *lrestart,double *time,double *dtime,
                          int *kstep,int *kinc));
#endif

void FORTRAN(ufaceload,(double *co,int *ipkon,int *kon,char *lakon,
			int *nboun, int *nodeboun,
                        int *nelemload,char *sideload,int *nload,
                        int *ne,int *nk));

void FORTRAN(uout,(double *v,int *mi,int *ithermal));

void FORTRAN(updatecont,(int *koncont,int *ncont,double *co,double *vold,
			 double *cg,double *straight,int *mi));

void FORTRAN(updatecontpen,(int *koncont,int *ncont,double *co,double *vold,
			 double *cg,double *straight,int *mi,int *imastnode,
                         int *nmastnode,double *xmastnor,int *ntie,
                         char *tieset,int *nset,char *set,int *istartset,
                         int *iendset,int *ialset,int *ipkon,char *lakon,
			 int *kon,double *cs,int *mcs,int *ics));

void *u_calloc(size_t num,size_t size);

void writeBasisParameter(FILE *f);

void FORTRAN(writeboun,(int *nodeboun,int *ndirboun,double *xboun,
      char *typeboun,int *nboun));

void FORTRAN(writebv,(double *,int *));

void FORTRAN(writeev,(double *,int *,double *,double *));

void FORTRAN(writeevcomplex,(double *eigxx,int *nev,double *fmin,double *fmax));

void FORTRAN(writeevcs,(double *,int *,int *,double *,double *));

void FORTRAN(writeevcscomplex,(double *eigxx,int *nev,int *nm,double *fmin,
            double *fmax));

void FORTRAN(writehe,(int *));

void FORTRAN(writeim,());

void FORTRAN(writeinput,(char *inpc,int *ipoinp,int *inp,int *nline,int *ninp,
                         int *ipoinpc));

void FORTRAN(writemac,(double *mac,int *nev));

void FORTRAN(writemaccs,(double *mac,int *nev,int* nm));

void FORTRAN(writempc,(int *,int *,double *,char *,int *));

void FORTRAN(writepf,(double *d,double *bjr,double *bji,double *freq , 
		      int *nev,int *mode,int *nherm));

void FORTRAN(writere,());

void FORTRAN(writesummary,(int *istep,int *j,int *icutb,int *l,double *ttime,
		   double *time,double *dtime));

void FORTRAN(writesummarydiv,(int *istep,int *j,int *icutb,int *l,double *ttime,
		   double *time,double *dtime));

void FORTRAN(writetetmesh,(int *kontet,int *netet_,double *cotet,
     int *nktet, double *field, int *nfield));
			
void FORTRAN(writeview,(int *ntr,double *adview,double *auview,double *fenv,
            int *nzsrad,char *jobnamef));

void FORTRAN(zienzhu,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
		      int *ne,double *stn,int *ipneigh,int *neigh,
		      double *sti,int *mi));

void FORTRAN(znaupd,(int *ido,char *bmat,int *n,char *which,int *nev,
	     double *tol,double *resid,int *ncv,double *z,int *ldz,
	     int *iparam,int *ipntr,double *workd,double *workl,
	     int *lworkl,double *rwork,int *info));

void FORTRAN(zneupd,(int *rvec,char *howmny,int *select,double *d,
	     double *z,int *ldz,double *sigma,
             double *workev,char *bmat,int *neq,char *which, 
	     int *nev,double *tol,double *resid,int *ncv,double *v,
	     int *ldv,int *iparam,int *ipntr,double *workd,
	     double *workl,int *lworkl,double *rwork,int *info));
