/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

char *lakon1,*sideload1, *matname1, *sideface1;

int *nk1,*kon1,*ipkon1,*ne1,*nodeboun1,*ndirboun1,*nboun1,*ipompc1,
    *nodempc1,*nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,*nload1,
    *ipobody1,*nbody1,*nactdoh1,*icolv1,*jqv1,*irowv1,neqv1,nzlv1,*nmethod1,
    *ikmpc1,*ilmpc1,*ikboun1,*ilboun1,*nrhcon1,*ielmat1,*ntmat_1,*ithermal1,
    nzsv1,*mi1,*ncmat_1,*nshcon1,*istep1,*iinc1,*ibody1,*turbulent1,
    *nelemface1,*nface1,compressible1,num_cpus,*icolp1,*jqp1,*irowp1,
    neqp1,nzlp1,nzsp1,iexplicit1,*ncocon1,neqt1,nzst1;

double *co1,*xboun1,*coefmpc1,*xforc1,*xload1,*xbody1,*rhcon1,*t01,
    *vold1,*voldaux1,dtimef1,*physcon1,*shcon1,*ttime1,timef1,*xloadold1,
    *voldtu1,*yy1,*b=NULL,*xbounact1,theta11,*v1,theta21,*cocon1,
    reltimef1;

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
    int *mi, int *ncmat_, double *physcon, int *istep, int *iinc,
    int *ibody, double *xloadold, double *xboun,
    double *coefmpc, int *nmethod, double *xforcold, double *xforcact,
    int *iamforc,int *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, int *iamt1, double *amta,
    int *namta, int *nam, double *ampli, double *xbounold, double *xbounact,
    int *iamboun, int *itg, int *ntg, char *amname, double *t0, int *nelemface,
    int *nface, double *cocon, int *ncocon, double *xloadact, double *tper,
    int *jmax, int *jout, char *set, int *nset, int *istartset,
    int *iendset, int *ialset, char *prset, char *prlab, int *nprint,
    double *trab, int *inotr, int *ntrans, char *filab, char *labmpc){

    /* main computational fluid dynamics routine */

    /* References:

       Zienkiewicz, O.C., Taylor, R.L. and Nithiarasu, P., "The Finite
       Element Method for Fluid Dynamics", 6th Edition, Elsevier (2006)

       Menter, F.R., "Two-Equation Eddy-Viscosity Turbulence Models
       for Engineering Applications", AIAA Journal(1994), 32(8), 
       1598-1605                                                       */
  
  int *ipointer=NULL, *mast1=NULL, *irowt=NULL, *irowv=NULL, *irowp=NULL,
      *irowk=NULL, *icolt=NULL, *icolv=NULL, *icolp=NULL, *icolk=NULL,
      *jqt=NULL, *jqv=NULL, *jqp=NULL, *jqk=NULL, *nactdoh=NULL,i,j, 
      *nactdok=NULL, *nx=NULL, *ny=NULL, *nz=NULL,nzs,neqt,neqv,neqp,
      neqk,nzst,nzsv,nzsp,nzsk,iexplicit,nzlt,nzlv,nzlp,nzlk,kode,nnstep,
      convergence,iout,iit,symmetryflag=0,inputformat=0,compressible,
      nmethodd,nstate_,*ielorien=NULL,norien,*inum=NULL,ismooth=0,
      *inomat=NULL,ikin=0,mt=mi[1]+1;

  double *yy=NULL, *xsolidsurf=NULL, *dt=NULL, *voldaux=NULL, *x=NULL,
      *y=NULL, *z=NULL, *xo=NULL, *yo=NULL, *zo=NULL, *adbt=NULL,
      *aubt=NULL, *adbv=NULL, *aubv=NULL, *adbp=NULL, *aubp=NULL,
      *adbk=NULL, *aubk=NULL,*v=NULL, *vtu=NULL,timef,ttimef,
      dtimef,*addiv=NULL,*sol=NULL, *aux=NULL,shockscale,
      *bk=NULL,*bt=NULL,*solk=NULL,*solt=NULL,theta1,theta2,*adb=NULL,
      *aub=NULL,sigma=0.,*dh=NULL,reltimef,*fn=NULL,*stx=NULL,
      *eme=NULL,*qfx=NULL,*orab=NULL,*xstate=NULL,*ener=NULL,
      csmooth=0.,shockcoefref=2.,*sa=NULL,*sav=NULL,shockcoef=2.,
      *adlt=NULL,*adlv=NULL,*adlp=NULL,*adlk=NULL;

  /* standard: shockcoef=2 */

#ifdef SGI
  int token;
#endif

  /* variables for multithreading procedure */

  int sys_cpus;
  char *env;
  
  num_cpus = 0;
#ifdef _SC_NPROCESSORS_CONF
  sys_cpus = sysconf(_SC_NPROCESSORS_CONF);
  if (sys_cpus <= 0)
      sys_cpus = 1;
#else
  sys_cpus = 1;
#endif
  env = getenv("CCX_NPROC");
  if (env)
      num_cpus = atoi(env);
  if (num_cpus > 0) {
//      if (num_cpus > sys_cpus)
//	  num_cpus = sys_cpus;
  } else if (num_cpus == -1) {
      num_cpus = sys_cpus;
  } else {
      num_cpus = 1;
  }
  printf("Using up to %d cpu(s) for spooles.\n", num_cpus);
  
  pthread_t tid[num_cpus];

  kode=0;

  /*  *iexpl==0:  structure:implicit, fluid:semi-implicit
      *iexpl==1:  structure:implicit, fluid:explicit
      *iexpl==2:  structure:explicit, fluid:semi-implicit
      *iexpl==3:  structure:explicit, fluid:explicit */

  if((*iexpl==1)||(*iexpl==3)){
      iexplicit=1;theta1=0.5;theta2=0.;compressible=1;
  }else{
      iexplicit=0;
      theta1=1.0;theta2=1.0;compressible=0;
  }

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  if(*ithermal==1) *ithermal=2;

  /* determining the matrix structure */

  nzs=1000000;
  
  ipointer=NNEW(int,3**nk);
  mast1=NNEW(int,nzs);
  irowv=NNEW(int,nzs);
  irowp=NNEW(int,nzs);
  icolv=NNEW(int,3**nk);
  icolp=NNEW(int,*nk);
  jqv=NNEW(int,3**nk+1);
  jqp=NNEW(int,*nk+1);
  nactdoh=NNEW(int,mt**nk);
  inomat=NNEW(int,*nk);

  if(*ithermal>1){
      irowt=NNEW(int,nzs);
      icolt=NNEW(int,*nk);
      jqt=NNEW(int,*nk+1);
  }

  if(*turbulent!=0){
      irowk=NNEW(int,nzs);
      icolk=NNEW(int,*nk);
      jqk=NNEW(int,*nk+1);
      nactdok=NNEW(int,*nk);
  }

  mastructf(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	    nodempc,nmpc,nactdoh,icolt,icolv,icolp,icolk,jqt,jqv,jqp,
	    jqk,&mast1,&irowt,&irowv,&irowp,&irowk,isolver,&neqt,&neqv,
            &neqp,&neqk,ikmpc,ilmpc,ipointer,&nzst,&nzsv,&nzsp,&nzsk,
            ithermal,ikboun,ilboun,turbulent,nactdok,ifreestream,nfreestream,
	    isolidsurf,nsolidsurf,&nzs,&iexplicit,ielmat,inomat,labmpc);

  free(ipointer);free(mast1);

  /* initialization */

  yy=NNEW(double,*nk);
  xsolidsurf=NNEW(double,*nsolidsurf);
  dh=NNEW(double,*nk);
  voldaux=NNEW(double,mt**nk);
  x=NNEW(double,*nsolidsurf);
  y=NNEW(double,*nsolidsurf);
  z=NNEW(double,*nsolidsurf);
  xo=NNEW(double,*nsolidsurf);
  yo=NNEW(double,*nsolidsurf);
  zo=NNEW(double,*nsolidsurf);
  nx=NNEW(int,*nsolidsurf);
  ny=NNEW(int,*nsolidsurf);
  nz=NNEW(int,*nsolidsurf);
  
  FORTRAN(initialcfd,(yy,nk,co,ne,ipkon,kon,lakon,x,y,z,xo,yo,zo,
       nx,ny,nz,isolidsurf,neighsolidsurf,xsolidsurf,dh,nshcon,shcon,
       nrhcon,rhcon,vold,voldaux,ntmat_,iponoel,inoel,
       &iexplicit,ielmat,nsolidsurf,turbulent,physcon,&compressible,
       matname,inomat,voldtu,mi));
  
  free(x);free(y);free(z);free(xo);free(yo);free(zo);free(nx);free(ny);
  free(nz);

  /* composing those left hand sides which do not depend on the increment */

  /* lhs for the energy */

  if(*ithermal>1){
      adbt=NNEW(double,neqt);
      aubt=NNEW(double,nzst);

      FORTRAN(mafilltlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	      xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
	      nactdoh,icolt,jqt,irowt,&neqt,&nzlt,
	      ikmpc,ilmpc,ikboun,ilboun,&nzst,adbt,aubt));

      adlt=NNEW(double,neqt);
      FORTRAN(lump,(adbt,aubt,adlt,irowt,jqt,&neqt));
  }

  /* lhs for the velocity */

  adbv=NNEW(double,neqv);
  aubv=NNEW(double,nzsv);

  FORTRAN(mafillvlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
       nactdoh,icolv,jqv,irowv,&neqv,&nzlv,
       ikmpc,ilmpc,ikboun,ilboun,&nzsv,adbv,aubv));

  adlv=NNEW(double,neqv);
  FORTRAN(lump,(adbv,aubv,adlv,irowv,jqv,&neqv));

  /* lhs for the pressure  */

  adbp=NNEW(double,neqp);
  aubp=NNEW(double,nzsp);
      
  FORTRAN(mafillplhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nactdoh,icolp,jqp,
          irowp,&neqp,&nzlp,ikmpc,ilmpc,ikboun,ilboun,&nzsp,adbp,aubp,
          nmethod,&iexplicit));

  if(iexplicit==1){
      adlp=NNEW(double,neqp);
      FORTRAN(lump,(adbp,aubp,adlp,irowp,jqp,&neqp));
  }

  if((iexplicit!=1)&&(neqp>0)){

  /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
			 &symmetryflag,&inputformat);
#else
	  printf("*ERROR in compfluid: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,token);
#else
	  printf("*ERROR in compfluid: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(adbp,&aubp,adb,aub,&sigma,icolp,&irowp,&neqp,&nzsp);
#else
	  printf("*ERROR in compfluid: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp);
#else
	  printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      
  }

  /* lhs for the turbulent */

  if(*turbulent!=0){
      adbk=NNEW(double,neqk);
      aubk=NNEW(double,nzsk);
      FORTRAN(mafillklhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
         xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
         nactdok,icolk,jqk,irowk,&neqk,&nzlk,
         ikmpc,ilmpc,ikboun,ilboun,&nzsk,adbk,aubk));

      adlk=NNEW(double,neqk);
      FORTRAN(lump,(adbk,aubk,adlk,irowk,jqk,&neqk));
  }

  /* starting the main loop */

  v=NNEW(double,mt**nk);
  vtu=NNEW(double,2**nk);

  /* ttimef is the total time up to the start of the present increment
     timef is the step time up to the end of the present increment 
     dtimef is the present increment size */

  ttimef=*ttime;
  timef=*time-*dtime;
  dt=NNEW(double,*nk);

  iit=0;

  do{

      iit++;

      /* determining a new time increment */
  
      FORTRAN(compdt,(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,iponoel,
	      inoel,&dtimef,&iexplicit,ielmat,physcon,dh,cocon,ncocon,ithermal,
              mi));

      /* fixed time */

      /*if(iexplicit==1) dtimef=1.e-4;*/
      /*if(iexplicit==1) dtimef/=3.;*/

      timef+=dtimef;
      if((*dtime<timef)&&(*nmethod==4)){
	  dtimef-=timef-*dtime;
	  timef=*dtime;
      }
      reltimef=timef/(*tper);
      if(reltimef>1.) reltimef=1.;
      /*printf("timef=%e,dtimef=%e\n",timef,dtimef);*/

      /* determining the instantaneous load */

      if(*nmethod==1){
	  nmethodd=4;
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltimef,ttime,dtime,ithermal,&nmethodd,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi));
/*	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
             co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));*/
	     }else if(*nmethod==4){
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi));
      }

      /*    if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
	      &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
	      }*/

      /* STEP 1: velocity correction */

/*      printf("STEP1: velocity correction\n\n");*/

      b=NNEW(double,num_cpus*neqv);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc,coefmpc1=coefmpc;nmpc1=nmpc;
      nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
      nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
      xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
      icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;nzlv1=nzlv;
      nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;voldaux1=voldaux;
      nzsv1=nzsv;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
      physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
      timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
      turbulent1=turbulent;voldtu1=voldtu;yy1=yy;nelemface1=nelemface;
      sideface1=sideface;nface1=nface;compressible1=compressible;
  
  /* create threads and wait */
  
      for(i=0; i<num_cpus; i++)  {
//	  printf("cpu=%d\n",i);
	  pthread_create(&tid[i], NULL, mafillv1rhsmt, (void *)i);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqv;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqv];
	  }
      }
      RENEW(b,double,neqv);

      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
      free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv1,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
      free(sol);

/*          if((iit/jout[1])*jout[1]==iit){
	  nnstep=1;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
	     &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/
      
      /* inserting the velocity boundary conditions */

	   FORTRAN(applybounv,(nodeboun,ndirboun,nboun,xbounact,
	     ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
	     nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
             voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
	     &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,mi));

      /* STEP 2: pressure correction */

	   /*   b=NNEW(double,neqp);

      int nea=1;
      int neb=*ne;

      FORTRAN(mafillprhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xbounact,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
       nface,b,nactdoh,icolp,jqp,irowp,&neqp,&nzlp,nmethod,ikmpc,ilmpc,
       ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,&nzsp,
       &dtimef,matname,mi,ncmat_,shcon,nshcon,v,&theta1,
       &iexplicit,physcon,&nea,&neb));*/

      b=NNEW(double,num_cpus*neqp);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xbounact1=xbounact;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc,coefmpc1=coefmpc;nmpc1=nmpc;
      nelemface1=nelemface;sideface1=sideface;nface1=nface;
      nactdoh1=nactdoh;icolp1=icolp;jqp1=jqp;irowp1=irowp;neqp1=neqp;
      nzlp1=nzlp;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;vold1=vold;voldaux1=voldaux;nzsp1=nzsp;dtimef1=dtimef;
      matname1=matname;mi1=mi;ncmat_1=ncmat_;shcon1=shcon;nshcon1=nshcon;
      v1=v;theta11=theta1;iexplicit1=iexplicit;physcon1=physcon;
  
      for(i=0; i<num_cpus; i++)  {
	  pthread_create(&tid[i], NULL, mafillprhsmt, (void *)i);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqp;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqp];
	  }
      }
      RENEW(b,double,neqp);

      sol=NNEW(double,neqp);
      if((iexplicit==1)&&(neqp>0)){
	  aux=NNEW(double,neqp);
	  FORTRAN(solveeq,(adbp,aubp,adlp,addiv,b,sol,aux,icolp,irowp,jqp,
			   &neqp,&nzsp,&nzlp));
	  free(b);free(aux);
      }else if(neqp>0){

          /* solving the system of equations (only for liquids) */

	  if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(b,&neqp);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	      sgi_solve(b,token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	      tau_solve(b,&neqp);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(b,&neqp);
#endif
	  }

          /* copying the solution into field sol */

	  for(i=0;i<neqp;i++){
	      sol[i]=b[i]/(theta1*theta2*dtimef*dtimef);
	  }
	  free(b);

      }

      /* storing the pressure correction in v */

      FORTRAN(resultsp,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,
              mi));
      free(sol);

      if(iexplicit==0){
      
      /* inserting the pressure boundary conditions for liquids */

          FORTRAN(applybounp,(nodeboun,ndirboun,nboun,xbounact,
       ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
       nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
       voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
       ipompc,nodempc,coefmpc,nmpc,inomat,mi));
      }

      /*    if((iit/jout[1])*jout[1]==iit){
	  nnstep=2;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
	     &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/
      
      /* STEP 3: velocity correction */

/*      printf("STEP3: velocity correction\n\n");*/

      b=NNEW(double,num_cpus*neqv);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc,coefmpc1=coefmpc;nmpc1=nmpc;
      nactdoh1=nactdoh;icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;
      nzlv1=nzlv;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;vold1=vold;nzsv1=nzsv;dtimef1=dtimef;v1=v;
      theta21=theta2;iexplicit1=iexplicit;mi1=mi;
  
  /* create threads and wait */
  
      for(i=0; i<num_cpus; i++)  {
	  pthread_create(&tid[i], NULL, mafillv2rhsmt, (void *)i);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqv;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqv];
	  }
      }
      RENEW(b,double,neqv);

      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
      free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv2,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
      free(sol);

      /*     if((iit/jout[1])*jout[1]==iit){
	  nnstep=3;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		&kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/
      
      /* inserting the velocity boundary conditions */

       FORTRAN(applybounv,(nodeboun,ndirboun,nboun,xbounact,
       ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
       nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
       voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
       &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,mi));

      /* STEP 4: energy correction */

/*      printf("STEP4: energy correction\n\n");*/

      if(*ithermal>1){

/*	  b=NNEW(double,neqt);
	  int nea=1;
	  int neb=*ne;

	  FORTRAN(mafilltrhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	     xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
             nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,
             b,nactdoh,&neqt,nmethod,ikmpc,ilmpc,ikboun,
             ilboun,rhcon,nrhcon,ielmat,ntmat_,t0,ithermal,vold,voldaux,&nzst,
             &dtimef,matname,mi,ncmat_,physcon,shcon,nshcon,ttime,&timef,
             istep,iinc,ibody,xloadold,&reltimef,cocon,ncocon,nelemface,
	     sideface,nface,&compressible,v,voldtu,yy,turbulent,&nea,&neb));*/

	   b=NNEW(double,num_cpus*neqt);

	  co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
	  nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
	  ipompc1=ipompc;nodempc1=nodempc,coefmpc1=coefmpc;nmpc1=nmpc;
	  nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
	  nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
	  xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
	  neqt1=neqt;
	  nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
	  ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
	  ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;voldaux1=voldaux;
	  nzst1=nzst;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
	  physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
	  timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
	  reltimef1=reltimef;cocon1=cocon;ncocon1=ncocon;nelemface1=nelemface;
          sideface1=sideface;nface1=nface;compressible1=compressible;v1=v;
          voldtu1=voldtu;yy1=yy;turbulent1=turbulent;
	  
	  for(i=0; i<num_cpus; i++)  {
	      pthread_create(&tid[i], NULL, mafilltrhsmt, (void *)i);
	  }
	  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	  
	  for(i=0;i<neqt;i++){
	      for(j=1;j<num_cpus;j++){
		  b[i]+=b[i+j*neqt];
	      }
	  }
	  RENEW(b,double,neqt);
	  
	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  FORTRAN(solveeq,(adbt,aubt,adlt,addiv,b,sol,aux,icolt,irowt,jqt,
			   &neqt,&nzst,&nzlt));
	  free(b);free(aux);
	  
      /* storing the temperature correction in v */

	  FORTRAN(resultst,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
	  free(sol);
      }

      /*    if((iit/jout[1])*jout[1]==iit){
	  nnstep=4;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/

      /* STEP 5: turbulent correction */

      if(*turbulent!=0){
/*	        if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/

	  bk=NNEW(double,neqk);
	  bt=NNEW(double,neqk);
	  
	  FORTRAN(mafillkrhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	    xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
	    nface,nactdok,&neqk,nmethod,ikmpc,ilmpc,
	    ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,&nzsk,
	    &dtimef,matname,mi,ncmat_,shcon,nshcon,v,&theta1,
	    bk,bt,voldtu,isolidsurf,nsolidsurf,ifreestream,nfreestream,
	    xsolidsurf,yy,&compressible,turbulent));
	  
	  solk=NNEW(double,neqk);
	  aux=NNEW(double,neqk);
	  FORTRAN(solveeq,(adbv,aubv,adlk,addiv,bk,solk,aux,icolk,irowk,jqk,
			   &neqk,&nzsk,&nzlk));
	  for(i=0;i<neqk;i++){printf("i=%d,bk=%e,solk=%e\n",i,bk[i],solk[i]);}
	  free(bk);free(aux);
	  
	  solt=NNEW(double,neqk);
	  aux=NNEW(double,neqk);
	  FORTRAN(solveeq,(adbv,aubv,adlk,addiv,bt,solt,aux,icolk,irowk,jqk,
			   &neqk,&nzsk,&nzlk));
	  free(bt);free(aux);
	  
	  /* storing the turbulence correction in vtu */
	  
	  FORTRAN(resultsk,(nk,nactdok,vtu,solk,solt,ipompc,nodempc,
			       coefmpc,nmpc));
	  free(solk);free(solt);

/*	        if((iit/jout[1])*jout[1]==iit){
	  nnstep=5;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
		  physcon,filab,inomat,ntrans,inotr,trab,mi));
		  }*/
      }
      
      /* adding v to voldaux and vtu to voldtu; 
         for thermal incompressible fluids: determine the static
                           temperature  */

      FORTRAN(updatecfd,(vold,voldaux,v,nk,
           ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	   nmethod,&convergence,physcon,iponoel,inoel,ithermal,
	   nactdoh,&iit,&compressible,&ismooth,voldtu,vtu,turbulent,
	   inomat,nodeboun,ndirboun,nboun,mi,&shockscale));

      /* inserting the boundary conditions for the turbulence
         parameters */

      /*     if(*turbulent!=0){
	  FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
	    iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
	    ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	    inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
	    rhcon,voldtu,ntmat_,labmpc,inomat));
	    }*/
      
      /* smoothing the solution (only for compressible fluids) */

      if(compressible){

	  ismooth=1;

	  /* smoothing rho * total energy density */

	  /*  adl=NNEW(double,neqt);
	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  for(i=0;i<neqt;i++){sol[i]=voldaux[5*i];}
	  FORTRAN(smooth,(adbt,aubt,adl,sol,aux,icolt,irowt,jqt,
			  &neqt,&nzlt,&csmooth));
	  for(i=0;i<neqt;i++){voldaux[5*i]=sol[i];}
	  free(adl);free(sol);free(aux);*/

	  /* smoothing rho * velocity */

	  /* adl=NNEW(double,neqv);
	  sol=NNEW(double,neqv);
	  aux=NNEW(double,neqv);
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  sol[3*i+j]=voldaux[5*i+j+1];
	      }
	  }
	  FORTRAN(smooth,(adbv,aubv,adl,sol,aux,icolv,irowv,jqv,
			  &neqv,&nzlv,&csmooth));
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  voldaux[5*i+j+1]=sol[3*i+j];
	      }
	  }
	  free(adl);free(sol);free(aux);*/

	  /* smoothing rho */

	  /* adl=NNEW(double,neqp);
	  sol=NNEW(double,neqp);
	  aux=NNEW(double,neqp);
	  for(i=0;i<neqp;i++){sol[i]=voldaux[5*i+4];}
	  FORTRAN(smooth,(adbp,aubp,adl,sol,aux,icolp,irowp,jqp,
			  &neqp,&nzlp,&csmooth));
	  for(i=0;i<neqp;i++){voldaux[5*i+4]=sol[i];}
	  free(adl);free(sol);free(aux);*/
	 
	  /* inserting the velocity boundary conditions */
 
/*	  FORTRAN(applybounv,(nodeboun,ndirboun,nboun,xbounact,
	    ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
	    nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
	    voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
	    &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat));*/
	  
	  /* determine the static temperature and the static pressure */
	  
	  FORTRAN(updatecfd,(vold,voldaux,v,nk,
	    ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	    nmethod,&convergence,physcon,iponoel,inoel,ithermal,
	    nactdoh,&iit,&compressible,&ismooth,voldtu,vtu,turbulent,
	    inomat,nodeboun,ndirboun,nboun,mi,&shockscale));

	  /* determining a measure for the pressure gradients */

//	  shockcoef=shockscale*shockcoefref;
	  sa=NNEW(double,neqt);
	  sav=NNEW(double,neqv);
	  FORTRAN(presgradient,(iponoel,inoel,sa,sav,nk,dt,&shockcoef,
				&dtimef,ipkon,kon,lakon,vold,mi));

	  /* shocksmoothing rho * total energy density */

	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  for(i=0;i<neqt;i++){sol[i]=voldaux[mt*i];}
	  FORTRAN(smoothshock,(adbt,aubt,adlt,addiv,sol,aux,icolt,irowt,jqt,
			  &neqt,&nzlt,sa));
	  for(i=0;i<neqt;i++){voldaux[mt*i]=sol[i];}
	  free(sol);free(aux);

	  /* shocksmoothing rho * velocity */

	  sol=NNEW(double,neqv);
	  aux=NNEW(double,neqv);
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  sol[3*i+j]=voldaux[mt*i+j+1];
	      }
	  }
	  FORTRAN(smoothshock,(adbv,aubv,adlv,addiv,sol,aux,icolv,irowv,jqv,
			  &neqv,&nzlv,sav));
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  voldaux[mt*i+j+1]=sol[3*i+j];
	      }
	  }
	  free(sol);free(aux);

	  /* shocksmoothing rho */

	  sol=NNEW(double,neqp);
	  aux=NNEW(double,neqp);
	  for(i=0;i<neqp;i++){sol[i]=voldaux[mt*i+4];}
	  FORTRAN(smoothshock,(adbp,aubp,adlp,addiv,sol,aux,icolp,irowp,jqp,
			  &neqp,&nzlp,sa));
	  for(i=0;i<neqp;i++){voldaux[mt*i+4]=sol[i];}
	  free(sol);free(aux);

	  free(sa);free(sav);
	 
	  /* inserting the velocity boundary conditions */
 
	  FORTRAN(applybounv,(nodeboun,ndirboun,nboun,xbounact,
	    ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
	    nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
	    voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
	    &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,mi));
	  
	  /* determine the static temperature and the static pressure */

	  if((iit/jout[1])*jout[1]==iit){ismooth=2;}	  

	  FORTRAN(updatecfd,(vold,voldaux,v,nk,
	    ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	    nmethod,&convergence,physcon,iponoel,inoel,ithermal,
	    nactdoh,&iit,&compressible,&ismooth,voldtu,vtu,turbulent,
	    inomat,nodeboun,ndirboun,nboun,mi,&shockscale));
	  
	  /* inserting the static pressure boundary conditions */
	  
	  FORTRAN(applybounpgas,(nodeboun,ndirboun,nboun,xbounact,
		  iponoel,vold,ipompc,nodempc,coefmpc,nmpc,inomat,matname,
                  nshcon,shcon,nrhcon,rhcon,physcon,ntmat_,
		  voldaux,mi));
	  
	  ismooth=0;
	  
      }
      
      /* inserting the static temperature boundary conditions */
      
      FORTRAN(applybount,(nodeboun,ndirboun,nboun,xbounact,
	iponoel,vold,ipompc,nodempc,coefmpc,nmpc,inomat,matname,
        nshcon,shcon,nrhcon,rhcon,physcon,&compressible,ntmat_,
	voldaux,mi));
      
      /* inserting the boundary conditions for the turbulence
	 parameters */
      
      if(*turbulent!=0){
	  FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
		 iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
		 ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	         inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
		 rhcon,voldtu,ntmat_,labmpc,inomat,mi));
      }
      
      if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
			    physcon,filab,inomat,ntrans,inotr,trab,mi));
      }
      
      if((iit/jout[1])*jout[1]==iit){
	  FORTRAN(printout,(set,nset,istartset,iendset,ialset,nprint,
	    prlab,prset,vold,t1,fn,ipkon,lakon,stx,eme,xstate,ener,
	    mi,&nstate_,ithermal,co,kon,qfx,&timef,trab,inotr,ntrans,
	    orab,ielorien,&norien,nk,ne,inum,filab,vold,&ikin));

	  /* lift and drag force */

	  FORTRAN(printoutface,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
	    cocon,ncocon,&compressible,istartset,iendset,ipkon,lakon,kon,
	    ialset,prset,&timef,nset,set,nprint,prlab,ielmat,mi));
      }
      
      if(convergence==1){
	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux,
			    physcon,filab,inomat,ntrans,inotr,trab,mi));
	  break;
      }
      if(iit==jmax[1]) FORTRAN(stop,());
      
      ttimef+=dtimef;
  }while(1);
  
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
    pardiso_cleanup(&neqp);
#endif
  }

  free(yy);free(xsolidsurf);free(dt);free(dh);free(voldaux);

  free(irowv);free(irowp);
  free(icolv);free(icolp);
  free(jqv);free(jqp);
  free(nactdoh);free(inomat);

  free(adbv);free(adbp);
  free(aubv);free(aubp);
  free(adlv);if(iexplicit==1) free(adlp);

  if(*ithermal>1){
      free(irowt);free(icolt);free(jqt);free(adbt);free(aubt);free(adlt);
  }

  if(*turbulent!=0){
      free(irowk);free(icolk);free(jqk);free(nactdok);
      free(adbk);free(aubk);free(adlk);
  }

  free(v);free(vtu);
  
  return;
  
} 

/* subroutine for multithreading of mafillv1rhs */

void *mafillv1rhsmt(void *i){

    int index,nea,neb,nedelta;

    index=((int)i)*neqv1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=((int)i)*nedelta+1;
    neb=(((int)i)+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillv1rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
	 xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,ndirforc1,xforc1,
	 nforc1,nelemload1,sideload1,xload1,nload1,xbody1,ipobody1,nbody1,
	 &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,nmethod1,ikmpc1,ilmpc1,ikboun1,
         ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,vold1,voldaux1,&nzsv1,
         &dtimef1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,&timef1,
         istep1,iinc1,ibody1,xloadold1,turbulent1,voldtu1,yy1,
	 nelemface1,sideface1,nface1,&compressible1,&nea,&neb));

    return NULL;
}

/* subroutine for multithreading of mafillprhs */

void *mafillprhsmt(void *i){

    int index,nea,neb,nedelta;

    index=((int)i)*neqp1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=((int)i)*nedelta+1;
    neb=(((int)i)+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillprhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
       xbounact1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nelemface1,sideface1,
       nface1,&b[index],nactdoh1,icolp1,jqp1,irowp1,&neqp1,&nzlp1,nmethod1,ikmpc1,ilmpc1,
       ikboun1,ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,vold1,voldaux1,&nzsp1,
       &dtimef1,matname1,mi1,ncmat_1,shcon1,nshcon1,v1,&theta11,
       &iexplicit1,physcon1,&nea,&neb));

    return NULL;
}

/* subroutine for multithreading of mafillv2rhs */

void *mafillv2rhsmt(void *i){

    int index,nea,neb,nedelta;

    index=((int)i)*neqv1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=((int)i)*nedelta+1;
    neb=(((int)i)+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillv2rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
       xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,
       &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,nmethod1,ikmpc1,ilmpc1,ikboun1,
       ilboun1,vold1,&nzsv1,&dtimef1,v1,&theta21,&iexplicit1,&nea,&neb,mi1));

    return NULL;
}

/* subroutine for multithreading of mafilltrhs */

void *mafilltrhsmt(void *i){

    int index,nea,neb,nedelta;

    index=((int)i)*neqt1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=((int)i)*nedelta+1;
    neb=(((int)i)+1)*nedelta;
    if(neb>*ne1) neb=*ne1;
	  
    FORTRAN(mafilltrhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
	     xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,ndirforc1,xforc1,
             nforc1,nelemload1,sideload1,xload1,nload1,xbody1,ipobody1,nbody1,
             &b[index],nactdoh1,&neqt1,nmethod1,ikmpc1,ilmpc1,ikboun1,
             ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,vold1,voldaux1,&nzst1,
             &dtimef1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,&timef1,
             istep1,iinc1,ibody1,xloadold1,&reltimef1,cocon1,ncocon1,nelemface1,
	     sideface1,nface1,&compressible1,v1,voldtu1,yy1,turbulent1,&nea,&neb));

    return NULL;
}







