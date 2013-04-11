
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
    int *nface, double *cocon, int *ncocon, double *xloadact, double *tper){

    /* main computational fluid dynamics routine */

    /* References:

       Zienkiewicz, O.C., Taylor, R.L. and Nithiarasu, P., "The Finite
       Element Method for Fluid Dynamics", 6th Edition, Elsevier (2006)

       Menter, F.R., "Two-Equation Eddy-Viscosity Turbulence Models
       for Engineering Applications", AIAA Journal(1994), 32(8), 
       1598-1605                                                       */
  
  int *ipointer=NULL, *mast1=NULL, *irowt=NULL, *irowv=NULL, *irowp=NULL,
      *irowk=NULL, *icolt=NULL, *icolv=NULL, *icolp=NULL, *icolk=NULL,
      *jqt=NULL, *jqv=NULL, *jqp=NULL, *jqk=NULL, *nactdoh=NULL,i, 
      *nactdok=NULL, *nx=NULL, *ny=NULL, *nz=NULL,nzs,neqt,neqv,neqp,
      neqk,nzst,nzsv,nzsp,nzsk,iexplicit,nzlt,nzlv,nzlp,nzlk,kode,nnstep,
      convergence,iout,iit,ifreq=10,symmetryflag=0,inputformat=0;
  double *yy=NULL, *xsolidsurf=NULL, *dt=NULL, *voldaux=NULL, *x=NULL,
      *y=NULL, *z=NULL, *xo=NULL, *yo=NULL, *zo=NULL, *adbt=NULL,
      *aubt=NULL, *adbv=NULL, *aubv=NULL, *adbp=NULL, *aubp=NULL,
      *adbk=NULL, *aubk=NULL,*v=NULL, *vtu=NULL,timef,ttimef,
      dtimef,*b=NULL,*adl=NULL,*addiv=NULL,*sol=NULL, *aux=NULL,
      *bk=NULL,*bt=NULL,*solk=NULL,*solt=NULL,theta1,theta2,*adb=NULL,
      *aub=NULL,sigma=0.,*dh=NULL,dtimefend,reltimef;

#ifdef SGI
  int token;
#endif

  kode=0;

  /*  *iexpl==0:  structure:implicit, fluid:semi-implicit
      *iexpl==1:  structure:implicit, fluid:explicit
      *iexpl==2:  structure:explicit, fluid:semi-implicit
      *iexpl==3:  structure:explicit, fluid:explicit */

  if((*iexpl==1)||(*iexpl==3)){
      iexplicit=1;theta2=0.;
  }else{
      iexplicit=0;
      theta2=0.5;
  }
  theta1=0.5;

  /* initialization */

  yy=NNEW(double,*nk);
  xsolidsurf=NNEW(double,*nsolidsurf);
  dt=NNEW(double,*nk);
  dh=NNEW(double,*nk);
  voldaux=NNEW(double,5**nk);
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
       &iexplicit,ielmat,nsolidsurf,turbulent,physcon));
  
  free(x);free(y);free(z);free(xo);free(yo);free(zo);free(nx);free(ny);
  free(nz);

  /* determining the matrix structure */

  nzs=1000000;
  
  ipointer=NNEW(int,3**nk);
  mast1=NNEW(int,nzs);
  irowv=NNEW(int,nzs);
  irowp=NNEW(int,nzs);
  icolv=NNEW(int,3**nk);
  icolp=NNEW(int,*nk);
  jqv=NNEW(int,3**nk);
  jqp=NNEW(int,*nk);
  nactdoh=NNEW(int,5**nk);

  if(*ithermal>1){
      irowt=NNEW(int,nzs);
      icolt=NNEW(int,*nk);
      jqt=NNEW(int,*nk);
  }

  if(*turbulent==1){
      irowk=NNEW(int,nzs);
      icolk=NNEW(int,*nk);
      jqk=NNEW(int,*nk);
      nactdok=NNEW(int,2**nk);
  }

  mastructf(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	    nodempc,nmpc,nactdoh,icolt,icolv,icolp,icolk,jqt,jqv,jqp,
	    jqk,&mast1,&irowt,&irowv,&irowp,&irowk,isolver,&neqt,&neqv,
            &neqp,&neqk,ikmpc,ilmpc,ipointer,&nzst,&nzsv,&nzsp,&nzsk,
            ithermal,ikboun,ilboun,turbulent,nactdok,ifreestream,nfreestream,
	    isolidsurf,nsolidsurf,&nzs);

  free(ipointer);free(mast1);

  /* composing those left hand sides which do not depend on the increment */

  /* lhs for the energy */

  if(*ithermal>0){
      adbt=NNEW(double,neqt);
      aubt=NNEW(double,nzst);
      FORTRAN(mafilltlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	      xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
	      nactdoh,icolt,jqt,irowt,&neqt,&nzlt,
	      ikmpc,ilmpc,ikboun,ilboun,&nzst,adbt,aubt));
  }

  /* lhs for the velocity */

  adbv=NNEW(double,neqv);
  aubv=NNEW(double,nzsv);

  FORTRAN(mafillvlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
       nactdoh,icolv,jqv,irowv,&neqv,&nzlv,
       ikmpc,ilmpc,ikboun,ilboun,&nzsv,adbv,aubv));

  /* lhs for the pressure  */

  adbp=NNEW(double,neqp);
  aubp=NNEW(double,nzsp);
      
  FORTRAN(mafillplhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nactdoh,icolp,jqp,
          irowp,&neqp,&nzlp,ikmpc,ilmpc,ikboun,ilboun,&nzsp,adbp,aubp,
          nmethod,&iexplicit));

  if((iexplicit!=1)&&(neqp>0)){

  /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
			 &symmetryflag,&inputformat);
#else
	  printf("*ERROR in arpack: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,token);
#else
	  printf("*ERROR in arpack: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(adbp,&aubp,adb,aub,&sigma,icolp,&irowp,&neqp,&nzsp);
#else
	  printf("*ERROR in arpack: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      
  }

  /* lhs for the turbulent */

  if(*turbulent==1){
      adbk=NNEW(double,neqk);
      aubk=NNEW(double,nzsk);
      FORTRAN(mafillklhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
         xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
         nactdok,icolk,jqk,irowk,&neqk,&nzlk,
         ikmpc,ilmpc,ikboun,ilboun,&nzsk,adbk,aubk));
  }

  /* starting the main loop */

  v=NNEW(double,5**nk);
  vtu=NNEW(double,2**nk);

  /* ttimef is the total time up to the start of the present increment
     timef is the step time up to the end of the present increment 
     dtimef is the present increment size */

  ttimef=*ttime;
  timef=*time-*dtime;

  iit=0;

  do{

      iit++;

      /* determining a new time increment */
  
      FORTRAN(compdt,(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,iponoel,
       inoel,&dtimef,&iexplicit,ielmat,physcon,dh));

      timef+=dtimef;
      if((*dtime<timef)&&(*nmethod==4)){
	  dtimef-=timef-*dtime;
	  timef=*dtime;
      }
      reltimef=timef/(*tper);
      if(reltimef>1.) reltimef=1.;
      printf("timef=%e,dtimef=%e\n",timef,dtimef);

      /* determining the instantaneous load */

      if(*nmethod==1){
/*	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltimef,ttime,dtime,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
             co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));*/
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
             co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));
	     }else if(*nmethod==4){
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
             co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));
      }

      /* insert boundary conditions */
      
      FORTRAN(applyboun,(nodeboun,ndirboun,nboun,xbounact,
       ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
       nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
       voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon));

      /* check whether time increment is also OK for the
         velocity boundary conditions at the end of the increment
         This is especially important in case that there is a
         jump in the velocity, e.g. because the boundary conditions
         are different from the initial conditions */
  
      FORTRAN(compdt,(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,iponoel,
       inoel,&dtimefend,&iexplicit,ielmat,physcon,dh));

      if(dtimefend<0.9*dtimef){

	  /* velocity boundary conditions at the end of the
             increment require a smaller time increment for
             stability */

	  timef-=dtimef-dtimefend;
	  dtimef=dtimefend;
	  reltimef=timef/(*tper);
	  if(reltimef>1.) reltimef=1.;

	  printf("correction: timef=%e,dtimef=%e\n",timef,dtimef);

	  /* determining the new instantaneous load */

	  if(*nmethod==1){
/*	      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
		xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
		xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
		namta,nam,ampli,time,&reltimef,ttime,dtime,ithermal,nmethod,
		xbounold,xboun,xbounact,iamboun,nboun,
		nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
		co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));*/
	      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
		xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
		xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
		namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
		xbounold,xboun,xbounact,iamboun,nboun,
		nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
		co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));
	  }else if(*nmethod==4){
	      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
             co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload));
	  }

	  /* insert boundary conditions */
      
	  FORTRAN(applyboun,(nodeboun,ndirboun,nboun,xbounact,
            ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
            nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
            voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon));
      }

      /* STEP 1: velocity correction */

      printf("STEP1: velocity correction\n\n");

      b=NNEW(double,neqv);

      FORTRAN(mafillv1rhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	 xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	 nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,
         b,nactdoh,icolv,jqv,irowv,&neqv,&nzlv,nmethod,ikmpc,ilmpc,ikboun,
         ilboun,rhcon,nrhcon,ielmat,ntmat_,t0,ithermal,vold,voldaux,&nzsv,
         &dtimef,matname,mint_,ncmat_,physcon,shcon,nshcon,ttime,&timef,
         istep,iinc,ibody,xloadold,turbulent,voldtu,yy,
         nelemface,sideface,nface));

      adl=NNEW(double,neqv);
      addiv=NNEW(double,neqv);
      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adl,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
      free(adl);free(addiv);free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv1,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc));
      free(sol);

       if((iit/ifreq)*ifreq==iit){
	  nnstep=1;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
			    &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
       }

      /* STEP 2: pressure correction */

      printf("STEP2: pressure correction\n\n");

      b=NNEW(double,neqp);

      FORTRAN(mafillprhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
       nface,b,nactdoh,icolp,jqp,irowp,&neqp,&nzlp,nmethod,ikmpc,ilmpc,
       ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,&nzsp,
       &dtimef,matname,mint_,ncmat_,shcon,nshcon,v,&theta1,
       &iexplicit,physcon));

      sol=NNEW(double,neqp);
      if((iexplicit==1)&&(neqp>0)){
	  adl=NNEW(double,neqp);
	  addiv=NNEW(double,neqp);
	  aux=NNEW(double,neqp);
	  FORTRAN(solveeq,(adbp,aubp,adl,addiv,b,sol,aux,icolp,irowp,jqp,
			   &neqp,&nzsp,&nzlp));
	  free(adl);free(addiv);free(b);free(aux);
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

          /* copying the solution into field sol */

	  for(i=0;i<neqp;i++){
	      sol[i]=b[i]/(theta1*theta2*dtimef);
	  }
	  free(b);

      }

      /* storing the pressure correction in v */

      FORTRAN(resultsp,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc));
      free(sol);

      if((iit/ifreq)*ifreq==iit){
	  nnstep=2;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
			    &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
      }
      
      /* STEP 3: velocity correction */

      printf("STEP3: velocity correction\n\n");

      b=NNEW(double,neqv);

      FORTRAN(mafillv2rhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
       b,nactdoh,icolv,jqv,irowv,&neqv,&nzlv,nmethod,ikmpc,ilmpc,ikboun,
       ilboun,vold,&nzsv,&dtimef,v,&theta2,&iexplicit));

      adl=NNEW(double,neqv);
      addiv=NNEW(double,neqv);
      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adl,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
      free(adl);free(addiv);free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv2,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc));
      free(sol);

       if((iit/ifreq)*ifreq==iit){
	  nnstep=3;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
			    &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
       }

      /* STEP 4: energy correction */

      printf("STEP4: energy correction\n\n");

      if(*ithermal>1){

	  b=NNEW(double,neqt);
	  
	  FORTRAN(mafilltrhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	     xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
             nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,
             b,nactdoh,&neqt,nmethod,ikmpc,ilmpc,ikboun,
             ilboun,rhcon,nrhcon,ielmat,ntmat_,t0,ithermal,vold,voldaux,&nzst,
             &dtimef,matname,mint_,ncmat_,physcon,shcon,nshcon,ttime,&timef,
             istep,iinc,ibody,xloadold,&reltimef,cocon,ncocon,nelemface,
             sideface,nface));
	  
	  adl=NNEW(double,neqt);
	  addiv=NNEW(double,neqt);
	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  FORTRAN(solveeq,(adbt,aubt,adl,addiv,b,sol,aux,icolt,irowt,jqt,
			   &neqt,&nzst,&nzlt));
	  free(adl);free(addiv);free(b);free(aux);
	  
      /* storing the temperature correction in v */

	  FORTRAN(resultst,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc));
	  free(sol);
      }

       if((iit/ifreq)*ifreq==iit){
	  nnstep=4;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
			    &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
       }

      /* STEP 5: turbulent correction */

      if(*turbulent==1){

	  bk=NNEW(double,neqk);
	  bt=NNEW(double,neqk);
	  
	  FORTRAN(mafillkrhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	    xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
	    nface,nactdok,&neqk,nmethod,ikmpc,ilmpc,
	    ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,&nzsk,
	    &dtimef,matname,mint_,ncmat_,shcon,nshcon,v,&theta1,
	    bk,bt,voldtu,isolidsurf,nsolidsurf,ifreestream,nfreestream,
	    xsolidsurf));
	  
	  adl=NNEW(double,neqk);
	  addiv=NNEW(double,neqk);
	  solk=NNEW(double,neqk);
	  aux=NNEW(double,neqk);
	  FORTRAN(solveeq,(adbv,aubv,adl,addiv,bk,solk,aux,icolk,irowk,jqk,
			   &neqk,&nzsk,&nzlk));
	  free(adl);free(addiv);free(bk);free(aux);
	  
	  adl=NNEW(double,neqk);
	  addiv=NNEW(double,neqk);
	  solt=NNEW(double,neqk);
	  aux=NNEW(double,neqk);
	  FORTRAN(solveeq,(adbv,aubv,adl,addiv,bt,solt,aux,icolk,irowk,jqk,
			   &neqk,&nzsk,&nzlk));
	  free(adl);free(addiv);free(bt);free(aux);
	  
	  /* storing the turbulence correction in vtu */
	  
	  FORTRAN(resultsk,(nk,nactdoh,vtu,solk,solt,ipompc,nodempc,
			       coefmpc,nmpc));
	  free(solk);free(solt);
      }
      
      /* adding v to vold and vtu to voldtu; updating voldaux */

      FORTRAN(updatecfd,(vold,voldaux,v,nk,
           ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	   nmethod,&convergence,physcon,iponoel,inoel,ithermal,
           nactdoh));

       if((iit/ifreq)*ifreq==iit){
	  nnstep=6;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
			    &kode,&timef,ielmat,matname,&nnstep,vtu,voldtu,voldaux));
       }

      if(convergence==1) break;
      if(iit==500) FORTRAN(stop,());
      
      ttimef+=dtimef;
  }while(1);

  free(yy);free(xsolidsurf);free(dt);free(dh);free(voldaux);

  free(irowv);free(irowp);
  free(icolv);free(icolp);
  free(jqv);free(jqp);
  free(nactdoh);

  free(adbv);free(adbp);
  free(aubv);free(aubp);

  if(*ithermal>1){
      free(irowt);free(icolt);free(jqt);free(adbt);free(aubt);
  }

  if(*turbulent==1){
      free(irowk);free(icolk);free(jqk);free(nactdok);
      free(adbk);free(aubk);
  }

  free(v);free(vtu);
  
  return;
  
} 

