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
#include <string.h>
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
#ifdef MATRIXSTORAGE
   #include "matrixstorage.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif

void arpack(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc, 
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int *icol, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *shcon, int *nshcon, double *cocon, int *ncocon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,
	     int *mei, double *fei,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double **xstatep, int *npmat_, char *matname, int *mi,
             int *ncmat_, int *nstate_, double **enerp, char *jobnamec,
             char *output, char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     char *cbody, int *ibody,double *xbody, int *nbody,
	     double *thicke, int *nslavs, double *tietol, int *nkon,
	     int *mpcinfo,int *ntie,int *istep,int *mcs,int *ics,
	     int *nnn,char *tieset,double *cs){

  /* calls the Arnoldi Package (ARPACK) */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A", fneig[132]="",
      description[13]="            ";

  int *inum=NULL,k,ido,ldz,iparam[11],ipntr[14],lworkl,ngraph=1,im,
    info,rvec=1,*select=NULL,lfin,j,lint,iout,ielas=0,icmd=0,mt=mi[1]+1,
    iinc=1,nev,ncv,mxiter,jrow,*ipobody=NULL,inewton=0,ifreebody,
    mass[2]={1,1}, stiffness=1, buckling=0, rhsi=0, intscheme=0,noddiam=-1,
    coriolis=0,symmetryflag=0,inputformat=0,*ipneigh=NULL,*neigh=NULL,ne0,
    *integerglob=NULL,nasym=0,zero=0,irenewxstate,ncont,*itietri=NULL,
    *koncont=NULL,ismallsliding=0,*itiefac=NULL,*islavsurf=NULL,
    *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,mortar=0,
    *imastop=NULL,*iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,ifacecount,
    mpcfree,memmpc_,icascade,maxlenmpc,nkon0,iit=-1,*irow=NULL,nherm=1,
    icfd=0,*inomat=NULL;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,
    *workl=NULL,*aux=NULL,*d=NULL,sigma=1,*temp_array=NULL,
    *een=NULL,sum,cam[5],*f=NULL,*fn=NULL,qa[3],*fext=NULL,
    *epn=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*emeini=NULL,
    *xstateini=NULL,*xstiff=NULL,*stiini=NULL,*vini=NULL,freq,*stx=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,
    *physcon=NULL,*qfx=NULL,*qfn=NULL,tol,fmin,fmax,pi,*cgr=NULL,
    *xloadold=NULL,reltime,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*springarea=NULL,*eenmax=NULL,
    *xnormastface=NULL,*doubleglob=NULL,*cg=NULL,*straight=NULL,
    *xmastnor=NULL,*areaslav=NULL,*xnoels=NULL,
    *di=NULL,sigmai=0,*workev=NULL,*ener=NULL,*xstate=NULL,*dc=NULL;

  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

#ifdef SGI
  int token;
#endif

  irow=*irowp;ener=*enerp;xstate=*xstatep;

  if((strcmp1(&filab[870],"PU  ")==0)||
     (strcmp1(&filab[1479],"PHS ")==0)||
     (strcmp1(&filab[1566],"MAXU")==0)||
     (strcmp1(&filab[1653],"MAXS")==0)){
      printf("*WARNING in arpack: PU, PHS, MAXU or MAX was selected in a frequency calculation without cyclic symmetry;\n this is not correct; output request is removed;\n");
      strcpy1(&filab[870],"    ",4);
      strcpy1(&filab[1479],"    ",4);
      strcpy1(&filab[1566],"    ",4);
      strcpy1(&filab[1653],"    ",4);
  }

  /* copying the frequency parameters */

  pi=4.*atan(1.);

  nev=mei[0];
  ncv=mei[1];
  mxiter=mei[2];
  tol=fei[0];
  fmin=2*pi*fei[1];
  fmax=2*pi*fei[2];
  
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
      if(inewton==1){
	  printf("*ERROR in arpackcs: generalized gravity loading is not allowed in frequency calculations");
      FORTRAN(stop,());
      }
  }

  ne0=*ne;nkon0=*nkon;

  /* contact conditions */
  
  if(*iperturb!=0){

      memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
      maxlenmpc=mpcinfo[3];

      if(*nslavs==0){irenewxstate=1;}else{irenewxstate=0;}
      inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	  lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
          &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
          &mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,&ifacecount,
          nmpc,&mpcfree,&memmpc_,
	  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	  iperturb,ikboun,nboun,co,istep,&xnoels);

      if(ncont!=0){

	  if(*nener==1){RENEW(ener,double,mi[0]*(*ne+*nslavs)*2);}
	  RENEW(ipkon,int,*ne+*nslavs);
	  RENEW(lakon,char,8*(*ne+*nslavs));
	  
	  if(*norien>0){
	      RENEW(ielorien,int,mi[2]*(*ne+*nslavs));
	      for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
	  }
	  RENEW(ielmat,int,mi[2]*(*ne+*nslavs));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;
	  cg=NNEW(double,3*ncont);
	  straight=NNEW(double,16*ncont);
      
    /* 11 instead of 10: last position is reserved for the
       local contact spring element number; needed as
       pointer into springarea */

	  RENEW(kon,int,*nkon+11**nslavs);
	  if((irenewxstate==1)&&(*nslavs!=0)){
	      RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	      for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
		  xstate[k]=0.;
	      }
	  }
	  xmastnor=NNEW(double,3*nmastnode[*ntie]);
	  xnormastface=NNEW(double,3*9**nslavs);
	  areaslav=NNEW(double,ifacecount);
	  springarea=NNEW(double,2**nslavs);
      
          /* generating contact spring elements */

	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
	     ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	     co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
	     &ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
	     &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	     iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
	     itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
	     imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
             xnoels);
	  
          /* determining the structure of the stiffness/mass matrix */
	  
	  remastructar(ipompc,&coefmpc,&nodempc,nmpc,
		 &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		 labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		 kon,ipkon,lakon,ne,nnn,nactdof,icol,jq,&irow,isolver,
		 neq,nzs,nmethod,ithermal,iperturb,mass,mi,ics,cs,
                 mcs);
      }
  }

  /* field for initial values of state variables (needed if
     previous static step was viscoplastic and for contact */

  if(*nstate_!=0){
      xstateini=NNEW(double,*nstate_*mi[0]*(ne0+*nslavs));
      for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
      xstateini[k]=xstate[k];
    }
  }

  /* determining the internal forces and the stiffness coefficients */

  f=NNEW(double,neq[1]);

  /* allocating a field for the stiffness matrix */

  xstiff=NNEW(double,(long long)27*mi[0]**ne);

  iout=-1;
  v=NNEW(double,mt**nk);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  fn=NNEW(double,mt**nk);
  stx=NNEW(double,6*mi[0]**ne);
  if(*ithermal>1){
      qfx=NNEW(double,3*mi[0]**ne);
  }
  inum=NNEW(int,*nk);
  if(*iperturb==0){
     results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t0,ithermal,
  	       prestr,iprestr,filab,eme,emn,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	       fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xload,xloadold,&icfd,inomat);
  }else{
     results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	       prestr,iprestr,filab,eme,emn,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	       fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xload,xloadold,&icfd,inomat);
  }
  free(f);free(v);free(fn);free(stx);if(*ithermal>1)free(qfx);free(inum);
  iout=1;

  /* for the frequency analysis linear strain and elastic properties
     are used */

  iperturb[1]=0;ielas=1;

  /* filling in the matrix */

  ad=NNEW(double,neq[1]);
  au=NNEW(double,nzs[2]);

  adb=NNEW(double,neq[1]);
  aub=NNEW(double,nzs[1]);

  fext=NNEW(double,neq[1]);

  if(*iperturb==0){
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
	      physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
	      &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,ntie,&nasym));
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
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
              physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
	      &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,ntie,&nasym));

      if(nasym==1){
	  RENEW(au,double,nzs[2]+nzs[1]);
	  RENEW(aub,double,nzs[2]+nzs[1]);
	  symmetryflag=2;
	  inputformat=1;
	  
	  FORTRAN(mafillsmas,(co,nk,kon,ipkon,lakon,ne,nodeboun,
                  ndirboun,xboun,nboun,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		  nforc,nelemload,sideload,xload,nload,xbody,ipobody,
		  nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		  nmethod,ikmpc,ilmpc,ikboun,ilboun,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
		  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		  xstiff,npmat_,&dtime,matname,mi,
                  ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
                  physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
                  &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
                  xstateini,xstate,thicke,
                  xnormastface,integerglob,doubleglob,tieset,istartset,iendset,
		  ialset,ntie,&nasym));
      }
  }

  free(fext);

  if(*nmethod==0){

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;time=0.;
    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }

    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);
    
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    free(inum);FORTRAN(stop,());

  }

  /* LU decomposition of the left hand matrix */

  if(nasym==1){sigma=0.;}else{sigma=1.;}

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
                   &symmetryflag,&inputformat,&nzs[2]);
#else
    printf("*ERROR in arpack: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
    printf("*ERROR in arpack: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
    printf("*ERROR in arpack: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==6){
#ifdef MATRIXSTORAGE
    matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
		  ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
                  lakon,kon,ne,mei,nboun,nmpc,cs,mcs);
#else
    printf("*ERROR in arpack: the MATRIXSTORAGE library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
    printf("*ERROR in arpack: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

/*  free(au);free(ad);*/

/* calculating the eigenvalues and eigenmodes */
  
  printf(" Calculating the eigenvalues and the eigenmodes\n\n");

  ido=0;
  ldz=neq[1];
  iparam[0]=1;
  iparam[2]=mxiter;
  iparam[3]=1;
  iparam[6]=3;

  info=0;

  resid=NNEW(double,neq[1]);
  z=NNEW(double,(long long)ncv*neq[1]);
  workd=NNEW(double,3*neq[1]);

  if(nasym==1){
      lworkl=3*ncv*(2+ncv);
      workl=NNEW(double,lworkl);
      FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
	  workl,&lworkl,&info));
  }else{
      lworkl=ncv*(8+ncv);
      workl=NNEW(double,lworkl);
      FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
	  workl,&lworkl,&info));
  }

  temp_array=NNEW(double,neq[1]);

  while((ido==-1)||(ido==1)||(ido==2)){
    if(ido==-1){
	if(nasym==1){
	    FORTRAN(opas,(&neq[1],aux,&workd[ipntr[0]-1],temp_array,adb,aub,icol,irow,nzl,nzs));
	}else{
	    FORTRAN(op,(&neq[1],aux,&workd[ipntr[0]-1],temp_array,adb,aub,icol,irow,nzl));
	}
    }
    if((ido==-1)||(ido==1)){

      /* solve the linear equation system  */

      if(ido==-1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(temp_array,&neq[1]);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(temp_array,token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(temp_array,&neq[1]);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(temp_array,&neq[1],&symmetryflag);
#endif
        }
        for(jrow=0;jrow<neq[1];jrow++){
          workd[ipntr[1]-1+jrow]=temp_array[jrow];
        }
      }
      else if(ido==1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(&workd[ipntr[2]-1],token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
        }
        if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag);
#endif
        }
        for(jrow=0;jrow<neq[1];jrow++){
          workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
        }
      }

    }

    if(ido==2){
	if(nasym==1){
	    FORTRAN(opas,(&neq[1],aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
                    adb,aub,icol,irow,nzl,nzs));
	}else{
	    FORTRAN(op,(&neq[1],aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
                    adb,aub,icol,irow,nzl));
	}
    }

    if(nasym==1){
	FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
                        iparam,ipntr,workd,workl,&lworkl,&info));
    }else{
	FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
                iparam,ipntr,workd,workl,&lworkl,&info));
    }
  }

/*--------------------------------------------------------------------*/
/*
   -----------
   free memory
   -----------
*/
  free(temp_array);
  if(*isolver==0){
#ifdef SPOOLES
    spooles_cleanup();
#endif
  }
  else if(*isolver==4){
#ifdef SGI
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
      pardiso_cleanup(&neq[1],&symmetryflag);
#endif
  }

  if(info!=0){
    printf("*ERROR in arpack: info=%d\n",info);
    printf("       # of converged eigenvalues=%d\n\n",iparam[4]);
  }         

  select=NNEW(int,ncv);

  if(nasym==1){
      d=NNEW(double,nev+1);
      di=NNEW(double,nev+1);
      workev=NNEW(double,3*ncv);
      FORTRAN(dneupd,(&rvec,howmny,select,d,di,z,&ldz,&sigma,&sigmai,
          workev,bmat,&neq[1],which,&nev,&tol,resid,
	  &ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
      free(workev);
      dc=NNEW(double,2*nev);

      /* storing as complex number and taking the square root */

      for(j=0;j<nev;j++){
	  dc[2*j]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])+d[j])/sqrt(2.);
	  dc[2*j+1]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])-d[j])/sqrt(2.);
	  if(di[j]<0.) dc[2*j+1]=-dc[2*j+1];
      }
      FORTRAN(writeevcomplex,(dc,&nev,&fmin,&fmax));
      free(di);free(dc);
  }else{
      d=NNEW(double,nev);
      FORTRAN(dseupd,(&rvec,howmny,select,d,z,&ldz,&sigma,bmat,&neq[1],which,
          &nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
      FORTRAN(writeev,(d,&nev,&fmin,&fmax));
  }
  free(select);free(workd);free(workl);free(resid);

  /* writing the eigenvalues and mass matrix to a binary file */

  if(mei[3]==1){

      strcpy(fneig,jobnamec);
      strcat(fneig,".eig");
      
      if((f1=fopen(fneig,"wb"))==NULL){
	  printf("*ERROR in arpack: cannot open eigenvalue file for writing...");

	  exit(0);
      }

      /* storing a zero as indication that this was not a
         cyclic symmetry calculation */

      if(fwrite(&zero,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving the cyclic symmetry flag to the eigenvalue file...");
	  exit(0);
      }

      /* Hermitian */

      if(fwrite(&nherm,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving the Hermitian flag to the eigenvalue file...");
	  exit(0);
      }

      /* storing the number of eigenvalues */

      if(fwrite(&nev,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving the number of eigenvalues to the eigenvalue file...");
	  exit(0);
      }

      /* the eigenfrequencies are stores as radians/time */

      if(fwrite(d,sizeof(double),nev,f1)!=nev){
	  printf("*ERROR saving the eigenfrequencies to the eigenvalue file...");
	  exit(0);
      }

      /* storing the stiffness matrix */

      if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR saving the diagonal of the stiffness matrix to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf("*ERROR saving the off-diagonal terms of the stiffness matrix to the eigenvalue file...");
	  exit(0);
      }

      /* storing the mass matrix */

      if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR saving the diagonal of the mass matrix to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf("*ERROR saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
	  exit(0);
      }
  }

  free(au);free(ad);

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stn=NNEW(double,6**nk);
  inum=NNEW(int,*nk);
  stx=NNEW(double,6*mi[0]**ne);
  if(*ithermal>1){
      qfn=NNEW(double,3**nk);
      qfx=NNEW(double,3*mi[0]**ne);
  }
  
  if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);

  temp_array=NNEW(double,neq[1]);

  lfin=0;
  for(j=0;j<nev;++j){
    lint=lfin;
    lfin=lfin+neq[1];

    for(k=0;k<6*mi[0]*ne0;k++){eme[k]=0.;}

    sum=0.;
    for(k=0;k<neq[1];++k)
      temp_array[k]=0.;
    if(nasym==1){
	FORTRAN(opas,(&neq[1],aux,&z[lint],temp_array,adb,aub,icol,irow,nzl,nzs));
    }else{
	FORTRAN(op,(&neq[1],aux,&z[lint],temp_array,adb,aub,icol,irow,nzl));
    }
    for(k=0;k<neq[1];++k)
      sum+=z[lint+k]*temp_array[k];
    for(k=0;k<neq[1];++k)
      z[lint+k]=z[lint+k]/sqrt(sum);

    if(mei[3]==1){
	if(fwrite(&z[lint],sizeof(double),neq[1],f1)!=neq[1]){
	    printf("*ERROR saving data to the eigenvalue file...");
	    exit(0);
	}
    }

    /* check whether the frequency belongs to the requested
       interval */

    if(fmin>-0.5){
	if(fmin*fmin>d[j]) continue;
    }
    if(fmax>-0.5){
	if(fmax*fmax<d[j]) continue;
    }

    if(*nprint>0) FORTRAN(writehe,(&j));

    eei=NNEW(double,6*mi[0]*ne0);
    if(*nener==1){
	stiini=NNEW(double,6*mi[0]*ne0);
	enerini=NNEW(double,mi[0]*ne0);}

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
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
            emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	    &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);}
    else{
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	    stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,&z[lint],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);
    }
    free(eei);
    if(*nener==1){
	free(stiini);free(enerini);}

    ++*kode;
    if(d[j]>=0.){
	freq=sqrt(d[j])/6.283185308;
    }else{
	freq=0.;
    }

    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }

    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&freq,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);

    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
  }

  if((fmax>-0.5)&&(fmax*fmax>d[nev-1])){
    printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
  }

  if(mei[3]==1){
      fclose(f1);
  }

  if(*iperturb!=0){
      if(ncont!=0){
	  *ne=ne0;*nkon=nkon0;
	  if(*nener==1){
	      RENEW(ener,double,mi[0]**ne*2);
	  }
	  RENEW(ipkon,int,*ne);
	  RENEW(lakon,char,8**ne);
	  RENEW(kon,int,*nkon);
	  if(*norien>0){
	      RENEW(ielorien,int,mi[2]**ne);
	  }
	  RENEW(ielmat,int,mi[2]**ne);
	  free(cg);free(straight);
	  free(imastop);free(itiefac);free(islavsurf);free(islavnode);
	  free(nslavnode);free(iponoels);free(inoels);free(imastnode);
	  free(nmastnode);free(itietri);free(koncont);
	  free(areaslav);free(springarea);free(xmastnor);free(xnormastface);
      }
  }
  
  free(adb);free(aub);free(temp_array);

  free(v);free(fn);free(stn);free(inum);free(stx);
  free(z);free(d);free(xstiff);free(ipobody);

  if(*ithermal>1){free(qfn);free(qfx);}

  if(*nstate_!=0){free(xstateini);}

  if(strcmp1(&filab[261],"E   ")==0) free(een);
  if(strcmp1(&filab[2697],"ME  ")==0) free(emn);
  if(strcmp1(&filab[522],"ENER")==0) free(enern);

  for(k=0;k<6*mi[0]*ne0;k++){eme[k]=0.;}

  if(*iperturb!=0){
      mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
      mpcinfo[3]=maxlenmpc;
  }

  *irowp=irow;*enerp=ener;*xstatep=xstate;

  return;
}

#endif
