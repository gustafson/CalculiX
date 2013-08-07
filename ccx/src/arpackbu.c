/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2013 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef ARPACK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif

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
             double *xstate, int *npmat_, char *matname, int *mi,
             int *ncmat_, int *nstate_, double *ener, char *output, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime,double *fmpc,
	     char *cbody, int *ibody,double *xbody, int *nbody, 
	     double *thicke,char *jobnamec){
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A",
      description[13]="            ",*tieset=NULL;

  int *inum=NULL,k,ido,dz,iparam[11],ipntr[11],lworkl,im,nasym=0,
    info,rvec=1,*select=NULL,lfin,j,lint,iout,iconverged=0,ielas,icmd=0,
    iinc=1,istep=1,*ncocon=NULL,*nshcon=NULL,nev,ncv,mxiter,jrow,
    *ipobody=NULL,inewton=0,coriolis=0,ifreebody,symmetryflag=0,
    inputformat=0,ngraph=1,mt=mi[1]+1,mass[2]={0,0}, stiffness=1, buckling=0, 
    rhsi=1, intscheme=0, noddiam=-1,*ipneigh=NULL,*neigh=NULL,ne0,
    *integerglob=NULL,ntie,icfd=0,*inomat=NULL;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,
    *workl=NULL,*aux=NULL,*d=NULL,sigma,*temp_array=NULL,
    *een=NULL,cam[5],*f=NULL,*fn=NULL,qa[3],*fext=NULL,
    time=0.,*epn=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,
    *xstateini=NULL,*xstiff=NULL,*stiini=NULL,*vini=NULL,*stx=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*cocon=NULL,
    *shcon=NULL,*physcon=NULL,*qfx=NULL,*qfn=NULL,tol, *cgr=NULL,
    *xloadold=NULL,reltime,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*cs=NULL,*springarea=NULL,*eenmax=NULL,
    *xnormastface=NULL,*emeini=NULL,*doubleglob=NULL;

  /* buckling routine; only for mechanical applications */

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime;

#ifdef SGI
  int token;
#endif
 
  /* copying the frequency parameters */

  nev=mei[0];
  ncv=mei[1];
  mxiter=mei[2];
  tol=fei[0];

  /* calculating the stresses due to the buckling load; this is a second
     order calculation if iperturb != 0 */

  *nmethod=1;
  
  /* assigning the body forces to the elements */ 

  if(*nbody>0){
      ifreebody=*ne+1;
      ipobody=NNEW(int,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,int,2*(*ne+ifreebody));
      }
      RENEW(ipobody,int,2*(ifreebody-1));
  }

  /* determining the internal forces and the stiffness coefficients */

  f=NNEW(double,neq[0]);

  /* allocating a field for the stiffness matrix */

  xstiff=NNEW(double,(long long)27*mi[0]**ne);

//  iout=-1;
  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stx=NNEW(double,6*mi[0]**ne);

  iout=-1;
  inum=NNEW(int,*nk);
  if(*iperturb==0){
     results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t0,ithermal,
	       prestr,iprestr,filab,eme,emn,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	       fmpc,nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
	       sideload,xload,xloadold,&icfd,inomat);
  }else{
     results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	       prestr,iprestr,filab,eme,emn,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	       fmpc,nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xload,xloadold,&icfd,inomat);
  }

  free(v);free(fn);free(stx);free(inum);
  iout=1;

  /* determining the system matrix and the external forces */

  ad=NNEW(double,neq[0]);
  au=NNEW(double,nzs[0]);
  fext=NNEW(double,neq[0]);

  if(*iperturb==0){
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      &nzs[0],stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
	      ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
              shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,&coriolis,
	      ibody,xloadold,&reltime,veold,springarea,nstate_,
	      xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,&ntie,&nasym));
  }
  else{
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
	      &nzs[0],stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
              shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,&coriolis,
	      ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,&ntie,&nasym));
  }

  /* determining the right hand side */

  b=NNEW(double,neq[0]);
  for(k=0;k<neq[0];++k){
      b[k]=fext[k]-f[k];
  }
  free(fext);free(f);

  if(*nmethod==0){

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;
    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }

    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
	    nstate_,&istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
            thicke,jobnamec,output,qfx);
    
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    free(inum);FORTRAN(stop,());

  }

  *nmethod=3;
  buckling=1;rhsi=0;

  sigma=0.;
  if(*isolver==0){
#ifdef SPOOLES
    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],&symmetryflag,
            &inputformat,&nzs[2]);
#else
    printf("*ERROR in arpackbu: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
#else
    printf("*ERROR in arpackbu: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[0],&nzs[0]);
#else
    printf("*ERROR in arpackbu: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		 &symmetryflag,&inputformat,jq,&nzs[2]);
#else
    printf("*ERROR in arpackbu: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stn=NNEW(double,6**nk);
  inum=NNEW(int,*nk);
  stx=NNEW(double,6*mi[0]**ne);
  
  if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);

  eei=NNEW(double,6*mi[0]**ne);
  if(*nener==1){
      stiini=NNEW(double,6*mi[0]**ne);
      enerini=NNEW(double,mi[0]**ne);}

  if(*iperturb==0){
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	  stx,elcon,nelcon,
	  rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
	  ntmat_,t0,t0,ithermal,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
          f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xboun,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,&bet,
          &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
          xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
          ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	  nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,&reltime,
          &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
          sideload,xload,xloadold,&icfd,inomat);}
  else{
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
          ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
          f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xboun,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,&bet,
          &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
          xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
          ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	  nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,&reltime,
          &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
          sideload,xload,xloadold,&icfd,inomat);
  }

  for(k=0;k<mt**nk;++k){
    vold[k]=v[k];
  }

  ++*kode;
  if(strcmp1(&filab[1044],"ZZS")==0){
      neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
  }

  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
	    nstate_,&istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
            thicke,jobnamec,output,qfx);

  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
  free(v);free(fn);free(stn);free(inum);

  if(strcmp1(&filab[261],"E   ")==0) free(een);
  if(strcmp1(&filab[2697],"ME  ")==0) free(emn);
  if(strcmp1(&filab[522],"ENER")==0) free(enern);

  /* in buckling mode stx and sti are kept */


  /* calculation of the left hand matrix (ad and au) and the right
     hand matrix (adb and aub); stx are the stresses due to the buckling
     load, sti due to previous loads, if any */

  aub=NNEW(double,nzs[0]);
  adb=NNEW(double,neq[0]);

  fext=NNEW(double,neq[0]);

  if(*iperturb==0){
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      &nzs[0],stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
              shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,&coriolis,
	      ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,&ntie,&nasym));
  }
  else{
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
	      &nzs[0],stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
              shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,&coriolis,
	      ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,&ntie,&nasym));
  }

  free(stx);free(fext);if(*nbody>0) free(ipobody);

  if(*nmethod==1){return;}

  /* loop checking the plausibility of the buckling factor
     if (5*sigma<buckling factor<50000*sigma) the solution is accepted,
     else sigma is set to buckling factor/500, and a new iteration is
     started */

  sigma=1.;

  do{


  /* LU decomposition of the left hand matrix */

//  sigma=1.;
  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
                   &symmetryflag,&inputformat,&nzs[2]);
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=2;
    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],token);
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[0],&nzs[0]);
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#endif
  }

  /* calculating the bucking factors and buckling modes */

  printf(" Calculating the buckling factors and buckling modes:\n\n");

  ido=0;
  dz=neq[0];
  iparam[0]=1;
  iparam[2]=mxiter;
  iparam[3]=1;
  iparam[6]=4;

  lworkl=ncv*(8+ncv);
  info=0;

  resid=NNEW(double,neq[0]);
  z=NNEW(double,ncv*neq[0]);
  workd=NNEW(double,3*neq[0]);
  workl=NNEW(double,lworkl);

  FORTRAN(dsaupd,(&ido,bmat,&neq[0],which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
	  workl,&lworkl,&info));

  temp_array=NNEW(double,neq[0]);

  while((ido==-1)||(ido==1)||(ido==2)){
    if(ido==-1){
      FORTRAN(op,(&neq[0],aux,&workd[ipntr[0]-1],temp_array,ad,au,icol,irow,nzl));
    }
    if((ido==-1)||(ido==1)){

      /* solve the linear equation system  */

      if(ido==-1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(temp_array,&neq[0]);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
	  token=2;
          sgi_solve(temp_array,token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(temp_array,&neq[0]);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(temp_array,&neq[0],&symmetryflag);
#endif
        }
        for(jrow=0;jrow<neq[0];jrow++){
          workd[ipntr[1]-1+jrow]=temp_array[jrow];
        }
      }
      else if(ido==1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(&workd[ipntr[2]-1],&neq[0]);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
	  token=2;
          sgi_solve(&workd[ipntr[2]-1],token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(&workd[ipntr[2]-1],&neq[0]);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(&workd[ipntr[2]-1],&neq[0],
               &symmetryflag);
#endif
        }
        for(jrow=0;jrow<neq[0];jrow++){
          workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
        }
      }

    }

    if(ido==2){
      FORTRAN(op,(&neq[0],aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],ad,au,icol,irow,nzl));
    }

    FORTRAN(dsaupd,(&ido,bmat,&neq[0],which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
	    workl,&lworkl,&info));
  }

/*--------------------------------------------------------------------*/
/*
   -----------
   free memory
   -----------
*/
  if(*isolver==0){
#ifdef SPOOLES
    spooles_cleanup();
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=2;
    sgi_cleanup(token);
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_cleanup();
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_cleanup(&neq[0],&symmetryflag);
#endif
  }

  if(info!=0){
    printf("*ERROR in arpackbu: info=%d\n",info);
    printf("       # of converged eigenvalues=%d\n\n",iparam[4]);
  }         

  select=NNEW(int,ncv);
  d=NNEW(double,nev);

  FORTRAN(dseupd,(&rvec,howmny,select,d,z,&dz,&sigma,bmat,&neq[0],which,&nev,&tol,resid,
	  &ncv,z,&dz,iparam,ipntr,workd,workl,&lworkl,&info));

  printf("sigma=%f,d[0]=%f\n\n",sigma,d[0]);
  if((5.>d[0]/sigma)||(50000.<d[0]/sigma)){
    if(iconverged<-4) {
      printf("no convergence for the buckling factor; maybe no buckling occurs");
      FORTRAN(stop,());
    }
    sigma=d[0]/500.;
    printf("no convergence; new iteration\n\n");
    --iconverged;
    free(z);free(d);
  }
  else{iconverged=0;}
     
  free(resid);free(workd);free(workl);free(select);free(temp_array);

  } while(iconverged<0);

  free(aub);free(adb);free(au);free(ad);

  FORTRAN(writebv,(d,&nev));

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stn=NNEW(double,6**nk);
  inum=NNEW(int,*nk);
  stx=NNEW(double,6*mi[0]**ne);
  
  if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);

  lfin=0;
  for(j=0;j<nev;++j){

    for(k=0;k<6*mi[0]**ne;k++){eme[k]=0.;}

    lint=lfin;
    lfin=lfin+neq[0];

    if(*nprint>0) FORTRAN(writehe,(&j));

//    memset(&v[0],0.,sizeof(double)*mt**nk);
    DMEMSET(v,0,mt**nk,0.);
    if(*iperturb==0){
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	    stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t0,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,&z[lint],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,&bet,
            &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);}
    else{
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	    stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,&z[lint],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[0],veold,accold,&bet,
            &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);
    }

    ++*kode;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }

    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&d[j],epn,ielmat,matname,enern,xstaten,
	    nstate_,&istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
   	    thicke,jobnamec,output,qfx);

    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
  }

  free(v);free(fn);free(stn);free(inum);free(stx);free(z);free(d);free(eei);
  if(*nener==1){
      free(stiini);free(enerini);}

  if(strcmp1(&filab[261],"E   ")==0) free(een);
  if(strcmp1(&filab[2697],"ME  ")==0) free(emn);
  if(strcmp1(&filab[522],"ENER")==0) free(enern);

  return;
}

#endif
