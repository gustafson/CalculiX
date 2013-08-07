/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
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

static char *lakon1,*sideload1, *matname1, *sideface1;

static int *nk1,*kon1,*ipkon1,*ne1,*nodeboun1,*ndirboun1,*nboun1,*ipompc1,
    *nodempc1,*nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,*nload1,
    *ipobody1,*nbody1,*nactdoh1,*icolv1,*jqv1,*irowv1,neqv1,nzlv1,*nmethod1,
    *ikmpc1,*ilmpc1,*ikboun1,*ilboun1,*nrhcon1,*ielmat1,*ntmat_1,*ithermal1,
    nzsv1,*mi1,*ncmat_1,*nshcon1,*istep1,*iinc1,*ibody1,*iturbulent1,
    *nelemface1,*nface1,compressible1,num_cpus,*icolp1,*jqp1,*irowp1,
    neqp1,nzlp1,nzsp1,iexplicit1,*ncocon1,neqt1,nzst1,*ipvar1,*ipvarf1,
    *nactdok1,neqk1,nzsk1,*isolidsurf1,*nsolidsurf1,*ifreestream1,
    *nfreestream1;

static double *co1,*xboun1,*coefmpc1,*xforc1,*xload1,*xbody1,*rhcon1,*t01,
    *vold1,*vcon1,dtimef1,*physcon1,*shcon1,*ttime1,timef1,*xloadold1,
    *vcontu1,*yy1,*b=NULL,*xbounact1,theta11,*v1,theta21,*cocon1,
    reltimef1,*dt1,*var1,*varf1,*sti1,*bk=NULL,*bt=NULL,*xsolidsurf1;

void compfluid(double **cop, int *nk, int **ipkonp, int **konp, char **lakonp,
    int *ne, char **sidefacep, int *ifreestream, 
    int *nfreestream, int *isolidsurf, int *neighsolidsurf,
    int *nsolidsurf, int **iponoelp, int **inoelp, int *nshcon, double *shcon,
    int *nrhcon, double *rhcon, double **voldp, int *ntmat_,int *nodeboun, 
    int *ndirboun, int *nboun, int **ipompcp,int **nodempcp, int *nmpc,
    int **ikmpcp, int **ilmpcp, int *ithermal, int *ikboun, int *ilboun,
    int *iturbulent, int *isolver, int *iexpl, double *vcontu, double *ttime,
    double *time, double *dtime, int *nodeforc,int *ndirforc,double *xforc,
    int *nforc, int *nelemload, char *sideload, double *xload,int *nload,
    double *xbody,int *ipobody,int *nbody, int **ielmatp, char *matname,
    int *mi, int *ncmat_, double *physcon, int *istep, int *iinc,
    int *ibody, double *xloadold, double *xboun,
    double **coefmpcp, int *nmethod, double *xforcold, double *xforcact,
    int *iamforc,int *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, int *iamt1, double *amta,
    int *namta, int *nam, double *ampli, double *xbounold, double *xbounact,
    int *iamboun, int *itg, int *ntg, char *amname, double *t0, 
    int **nelemfacep,
    int *nface, double *cocon, int *ncocon, double *xloadact, double *tper,
    int *jmax, int *jout, char *set, int *nset, int *istartset,
    int *iendset, int *ialset, char *prset, char *prlab, int *nprint,
    double *trab, int *inotr, int *ntrans, char *filab, char **labmpcp, 
    double *sti, int *norien, double *orab, char *jobnamef,char *tieset,
    int *ntie, int *mcs, int *ics, double *cs, int *nkon, int *mpcfree,
    int *memmpc_,double **fmpcp,int *nef,int **inomatp,double *qfx){

    /* main computational fluid dynamics routine */

    /* References:

       Zienkiewicz, O.C., Taylor, R.L. and Nithiarasu, P., "The Finite
       Element Method for Fluid Dynamics", 6th Edition, Elsevier (2006)

       Menter, F.R., "Two-Equation Eddy-Viscosity Turbulence Models
       for Engineering Applications", AIAA Journal(1994), 32(8), 
       1598-1605                                                       */
  
  char cflag[1],*labmpc=NULL,*lakon=NULL,*sideface=NULL;

  int *ipointer=NULL, *mast1=NULL, *irowt=NULL, *irowv=NULL, *irowp=NULL,
      *irowk=NULL, *icolt=NULL, *icolv=NULL, *icolp=NULL, *icolk=NULL,
    *jqt=NULL, *jqv=NULL, *jqp=NULL, *jqk=NULL, *nactdoh=NULL,i,j,k,
      *nactdok=NULL, *nx=NULL, *ny=NULL, *nz=NULL,nzs,neqt,neqv,neqp,
      neqk,nzst,nzsv,nzsp,nzsk,iexplicit,nzlt,nzlv,nzlp,nzlk,kode,nnstep,
      iconvergence,iout,iit,symmetryflag=0,inputformat=0,compressible,
      nmethodd,nstate_=0,*ielorien=NULL,*inum=NULL,ismooth=0,iqfx=0,isti=0,
      ikin=0,mt=mi[1]+1,*ipvar=NULL,*ipvarf=NULL,nvar_,nvarf_,
      nfield,ndim,iorienglob,icfdout=1,force=0,euler=1,*ithread=NULL,
      *integerglob=NULL,nslav,*islav=NULL,ncs,nmast,*imast=NULL,
      nkref,neref,*nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,
      *ipkon=NULL,*kon=NULL,*ielmat=NULL,*nelemface=NULL,*inoel=NULL,
      *iponoel=NULL,*ipoface=NULL,*nodface=NULL,inoelfree,*inoslav=NULL,
      *inomast=NULL,*ielslav=NULL,*ielmast=NULL,*nr=NULL,*inomat=NULL,
      nstart=501,memmpcref_,mpcfreeref,nmpcref,nefref,*nodempcref=NULL,
      ithermalref;

  double *yy=NULL, *xsolidsurf=NULL, *dt=NULL, *vcon=NULL, *x=NULL,
      *y=NULL, *z=NULL, *xo=NULL, *yo=NULL, *zo=NULL, *adbt=NULL,
      *aubt=NULL, *adbv=NULL, *aubv=NULL, *adbp=NULL, *aubp=NULL,
      *adbk=NULL, *aubk=NULL,*v=NULL, *vtu=NULL,timef,ttimef,
      dtimef,*addiv=NULL,*sol=NULL, *aux=NULL,shockscale,*stn=NULL,
      *solk=NULL,*solt=NULL,theta1,theta2,*adb=NULL,*qfn=NULL,
      *aub=NULL,sigma=0.,*dh=NULL,reltimef,*fn=NULL,*thicke=NULL,
      *eme=NULL,*xstate=NULL,*ener=NULL,
      csmooth=0.,shockcoef,*sa=NULL,*sav=NULL,*varf=NULL,
      *adlt=NULL,*adlv=NULL,*adlp=NULL,*adlk=NULL,factor=1.,*var=NULL,
      *vconini=NULL,*doubleglob=NULL,*coefmpc=NULL,*fmpc=NULL,
      *co=NULL,*vold=NULL,*rcs=NULL,*zcs=NULL,*rcs0=NULL,*zcs0=NULL,
      sum=0.,sumx=0.,sumxx=0.,sumy[7]={0.,0.,0.,0.,0.,0.,0.},
      sumxy[7]={0.,0.,0.,0.,0.,0.,0.},*del=NULL,reltime,*coefmpcref=NULL;

  nodempc=*nodempcp;ipompc=*ipompcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;fmpc=*fmpcp;co=*cop;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielmat=*ielmatp;
  nelemface=*nelemfacep;sideface=*sidefacep;inoel=*inoelp;
  iponoel=*iponoelp;vold=*voldp;inomat=*inomatp;

  /* standard: shockcoef=0 */

#ifdef SGI
  int token;
#endif

  /* open frd-file for fluids */

  FORTRAN(openfilefluid,(jobnamef));

  /* variables for multithreading procedure */

  int sys_cpus;
  char *env,*envloc,*envsys;
      
  num_cpus = 0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
      sys_cpus=atoi(envsys);
      if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
      sys_cpus = getSystemCPUs();
      if(sys_cpus<1) sys_cpus=1;
  }
  
  /* local declaration prevails, if strictly positive */
  
  envloc = getenv("CCX_NPROC_CFD");
  if(envloc){
      num_cpus=atoi(envloc);
      if(num_cpus<0){
	  num_cpus=0;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
      if (env)
	  num_cpus = atoi(env);
      if (num_cpus < 1) {
	  num_cpus=1;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
// next line is to be inserted in a similar way for all other paralell parts
  
  if(*ne<num_cpus) num_cpus=*ne;
  
  printf(" Using up to %d cpu(s) for CFD.\n", num_cpus);
  
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

  ithermalref=*ithermal;
  if(*ithermal==1){
    *ithermal=2;
  }

  /* generating additional elements and MPC's for cyclic symmetric
     structures */

  if(*mcs==1){

      /* find the cyclic symmetry slave nodes */

      /* ncs is the number of fluid and solid master nodes,
         nslav is the number of fluid master nodes only */

      ncs=cs[3];

      islav=NNEW(int,ncs);
      imast=NNEW(int,ncs);

      nr=NNEW(int,ncs);
      nz=NNEW(int,ncs);
      rcs=NNEW(double,ncs);
      zcs=NNEW(double,ncs);
      rcs0=NNEW(double,ncs);
      zcs0=NNEW(double,ncs);

      inoslav=NNEW(int,*nk);
      inomast=NNEW(int,*nk);

      FORTRAN(findslavcfd,(nmpc,labmpc,ipompc,nodempc,islav,&nslav,
			   inoslav,inomast,ics,cs,imast,&nmast,co,inomat,
                           nr,nz,rcs,zcs,rcs0,zcs0,&ncs));
      RENEW(islav,int,nslav);
      RENEW(imast,int,nmast);

      free(nr);free(nz);free(rcs);free(zcs);free(rcs0);free(zcs0);

      /* generate new nodes and elements in a layer on the slave
         and master side */

      nkref=*nk;neref=*ne;nefref=*nef;

      ielslav=NNEW(int,*ne);
      ielmast=NNEW(int,*ne);

      RENEW(co,double,3*(3**nk));
      RENEW(vold,double,mt*(3**nk));
      RENEW(ipkon,int,3**ne);
      RENEW(lakon,char,8*(3**ne));
      RENEW(kon,int,*nkon+8*2**ne);
      RENEW(ielmat,int,mi[2]*3**ne);
      RENEW(inomat,int,3**nk);

      FORTRAN(gencycsymelemcfd,(cs,islav,
	 &nslav,imast,&nmast,inomat,
	 nk,co,ne,ipkon,lakon,kon,nkon,ielmat,mi,vold,
         ielslav,ielmast,inoslav,inomast,iponoel,inoel));

      free(ielslav);free(ielmast);

      RENEW(co,double,3**nk);
      RENEW(vold,double,mt**nk);
      RENEW(ipkon,int,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,int,*nkon);
      RENEW(ielmat,int,mi[2]**ne);
      RENEW(inomat,int,*nk);

      /* generate new MPC's */

      if(*ithermal>1){
	  RENEW(ipompc,int,*nmpc+5*(2**nk));
	  RENEW(ikmpc,int,*nmpc+5*(2**nk));
	  RENEW(ilmpc,int,*nmpc+5*(2**nk));
	  RENEW(labmpc,char,20*(*nmpc+5*(2**nk))+1);
      }else{
	  RENEW(ipompc,int,*nmpc+4*(2**nk));
	  RENEW(ikmpc,int,*nmpc+4*(2**nk));
	  RENEW(ilmpc,int,*nmpc+4*(2**nk));
	  RENEW(labmpc,char,20*(*nmpc+4*(2**nk))+1);
      }

      memmpcref_=*memmpc_;mpcfreeref=*mpcfree;nmpcref=*nmpc;
      nodempcref=NNEW(int,3**memmpc_);
      for(k=0;k<3**memmpc_;k++){nodempcref[k]=nodempc[k];}
      coefmpcref=NNEW(double,*memmpc_);
      for(k=0;k<*memmpc_;k++){coefmpcref[k]=coefmpc[k];}

      interpolcycsymcfd(&nkref,co,&neref,ipkon,kon,&nodempc,ipompc,nmpc,
                        ikmpc,ilmpc,&coefmpc,labmpc,mpcfree,memmpc_,lakon,
                        &nmast,&nslav,ithermal,cs,inoslav,inomast,imast,islav);

      RENEW(ipompc,int,*nmpc);
      RENEW(ikmpc,int,*nmpc);
      RENEW(ilmpc,int,*nmpc);
      RENEW(labmpc,char,20**nmpc);

      free(inoslav);free(inomast);

      /* decascading the new MPC's */

      int mpcend,mpcmult,maxlenmpc;
      int icascade=0;
      int callfrommain=0;
      int iperturb=0;
      cascade(ipompc,&coefmpc,&nodempc,nmpc,mpcfree,nodeboun,ndirboun,
              nboun,ikmpc,ilmpc,ikboun,ilboun,&mpcend,&mpcmult,labmpc,
              nk,memmpc_,&icascade,&maxlenmpc,&callfrommain,&iperturb,
              ithermal);

      /* update the field with external faces
	 and node-to-element dependence */

      free(iponoel);free(inoel);free(sideface);free(nelemface);
      *nef+=(*ne-neref);
      sideface=NNEW(char,6**nef);
      nelemface=NNEW(int,6**nef);
      ipoface=NNEW(int,*nk);
      nodface=NNEW(int,5*6**nef);
      iponoel=NNEW(int,*nk);
      inoel=NNEW(int,3*20**nef);

      FORTRAN(precfdcyc,(nelemface,sideface,nface,ipoface,nodface,
       ne,ipkon,kon,lakon,ikboun,ilboun,xboun,nboun,nk,isolidsurf,
       nsolidsurf,ifreestream,nfreestream,neighsolidsurf,iponoel,inoel,
       &inoelfree,nef,co,ipompc,nodempc,ikmpc,ilmpc,nmpc,set,istartset,
       iendset,ialset,nset,iturbulent));

      RENEW(sideface,char,*nface);
      RENEW(nelemface,int,*nface);
      free(ipoface);free(nodface);
      RENEW(inoel,int,3*inoelfree);

  }else if(*mcs>1){
      printf(" *ERROR in compfluid: for CFD only one cyclic symmetry\n");
      printf("        conditions is allowed\n");
      FORTRAN(stop,());
  }

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

  if(*ithermal>1){
      irowt=NNEW(int,nzs);
      icolt=NNEW(int,*nk);
      jqt=NNEW(int,*nk+1);
  }

  if(*iturbulent!=0){
      irowk=NNEW(int,nzs);
      icolk=NNEW(int,*nk);
      jqk=NNEW(int,*nk+1);
      nactdok=NNEW(int,*nk);
  }

  mastructf(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	    nodempc,nmpc,nactdoh,icolt,icolv,icolp,icolk,jqt,jqv,jqp,
	    jqk,&mast1,&irowt,&irowv,&irowp,&irowk,isolver,&neqt,&neqv,
            &neqp,&neqk,ikmpc,ilmpc,ipointer,&nzst,&nzsv,&nzsp,&nzsk,
            ithermal,ikboun,ilboun,iturbulent,nactdok,ifreestream,nfreestream,
	    isolidsurf,nsolidsurf,&nzs,&iexplicit,ielmat,inomat,labmpc);

  free(ipointer);free(mast1);

  /* initialization */

  yy=NNEW(double,*nk);
  xsolidsurf=NNEW(double,*nsolidsurf);
  dh=NNEW(double,*nk);
  vcon=NNEW(double,mt**nk);
  vconini=NNEW(double,mt**nk);
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
       nrhcon,rhcon,vold,vcon,ntmat_,iponoel,inoel,
       &iexplicit,ielmat,nsolidsurf,iturbulent,physcon,&compressible,
		      matname,inomat,vcontu,mi,&euler,ithermal));
  
  free(x);free(y);free(z);free(xo);free(yo);free(zo);free(nx);free(ny);
  free(nz);

  /* calculating the shape functions, their derivatives and the
     Jacobian determinant in the integration points of the elements */

  nvar_=35**ne;
  ipvar=NNEW(int,*ne);
  var=NNEW(double,nvar_);

  nvarf_=8**ne;
  ipvarf=NNEW(int,*ne);
  varf=NNEW(double,nvarf_);

  calcshapef(&nvar_,ipvar,&var,ne,lakon,co,ipkon,kon,
             nelemface,sideface,nface,&nvarf_,ipvarf,
             &varf);

  /* composing those left hand sides which do not depend on the increment */

  /* lhs for the energy */

  if(*ithermal>1){
      adbt=NNEW(double,neqt);
      aubt=NNEW(double,nzst);

      FORTRAN(mafilltlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	      xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
	      nactdoh,icolt,jqt,irowt,&neqt,&nzlt,
	      ikmpc,ilmpc,ikboun,ilboun,&nzst,adbt,aubt,ipvar,var));

      adlt=NNEW(double,neqt);
      FORTRAN(lump,(adbt,aubt,adlt,irowt,jqt,&neqt));
  }

  /* lhs for the velocity */

  adbv=NNEW(double,neqv);
  aubv=NNEW(double,nzsv);

  FORTRAN(mafillvlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
       nactdoh,icolv,jqv,irowv,&neqv,&nzlv,
       ikmpc,ilmpc,ikboun,ilboun,&nzsv,adbv,aubv,ipvar,var));

  adlv=NNEW(double,neqv);
  FORTRAN(lump,(adbv,aubv,adlv,irowv,jqv,&neqv));

  /* lhs for the pressure  */

  adbp=NNEW(double,neqp);
  aubp=NNEW(double,nzsp);
      
  FORTRAN(mafillplhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
	  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nactdoh,icolp,jqp,
          irowp,&neqp,&nzlp,ikmpc,ilmpc,ikboun,ilboun,&nzsp,adbp,aubp,
	  nmethod,&iexplicit,ipvar,var));

  /* lumping is only applied to compressible fluids */

  if(iexplicit==1){
      adlp=NNEW(double,neqp);
      FORTRAN(lump,(adbp,aubp,adlp,irowp,jqp,&neqp));
  }

  if((iexplicit!=1)&&(neqp>0)){

  /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
			 &symmetryflag,&inputformat,&nzsp);
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
	  pardiso_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
                         &symmetryflag,&inputformat,jqp,&nzsp);
#else
	  printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      
  }

  /* lhs for the turbulent parameters */

  if(*iturbulent!=0){
      adbk=NNEW(double,neqk);
      aubk=NNEW(double,nzsk);
      FORTRAN(mafillklhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
         xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
         nactdok,icolk,jqk,irowk,&neqk,&nzlk,
	 ikmpc,ilmpc,ikboun,ilboun,&nzsk,adbk,aubk,ipvar,var));

      adlk=NNEW(double,neqk);
      FORTRAN(lump,(adbk,aubk,adlk,irowk,jqk,&neqk));
  }

  /* starting the main loop */

  v=NNEW(double,mt**nk);
  vtu=NNEW(double,2**nk);
      
  /* inserting the velocity and temperature conditions 
     for incompressible materials*/
    
  if(compressible==0){
      FORTRAN(applyboun,(nodeboun,ndirboun,nboun,xbounact,
	     ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
	     nsolidsurf,xsolidsurf,nfreestream,ifreestream,iturbulent,
	     vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
	     &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,mi,
	     ikboun,ilboun,ilmpc,labmpc));
  }
  
  /* ttime is the total time up to the start of the present increment
     time is the step time up to the end of the present increment 
     dtime is the present increment size */

  reltime=*time/(*tper);
//  printf("reltime= %e %e %e\n",*time,*tper,reltime);

  ttimef=*ttime;
  timef=*time-*dtime;
  dt=NNEW(double,*nk);

  if(compressible){
      sa=NNEW(double,neqt);
      sav=NNEW(double,neqv);
      del=NNEW(double,7*jmax[1]);
      shockcoef=physcon[9];
  }

  iit=0;

  do{

      iit++;

      /* determining a new time increment */

      FORTRAN(compdt,(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,iponoel,
	      inoel,&dtimef,&iexplicit,ielmat,physcon,dh,cocon,ncocon,ithermal,
	      mi,ipkon,kon,lakon,ne,v,co,iturbulent,vcontu,vcon));

      /* correction for too large changes */

//      if(iit>2){
//ccc      if(iit>=2){
//ccc	  if((iexplicit==1)&&(*nmethod==1)){dtimef*=factor;}
//	  if(*nmethod==1){dtimef*=factor;}
//      }

      timef+=dtimef;
/*      if((*dtime<timef)&&(*nmethod==4)){
	  dtimef-=timef-*dtime;
	  timef=*dtime;
	  }*/
      if((*time<timef)&&(*nmethod==4)){
	  dtimef-=timef-*time;
	  timef=*time;
	  iconvergence=1;
      }
      reltimef=timef/(*tper);
      if(reltimef>1.) reltimef=1.;

      /* determining the instantaneous load */

      if(*nmethod==1){
/*	  nmethodd=4;
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltimef,ttime,dtime,ithermal,&nmethodd,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	     ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));*/

          /* boundary conditions at end of mechanical increment */

	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltime,ttime,dtime,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	     ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
      }else if(*nmethod==4){

          /* boundary conditions at end of fluid increment */

	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
             ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
      }

      /*    if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
	      &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon));
	      }*/

      /* STEP 1: velocity correction */

      b=NNEW(double,num_cpus*neqv);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
      nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
      nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
      xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
      icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;nzlv1=nzlv;
      nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;vcon1=vcon;
      nzsv1=nzsv;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
      physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
      timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
      iturbulent1=iturbulent;vcontu1=vcontu;yy1=yy;nelemface1=nelemface;
      sideface1=sideface;nface1=nface;compressible1=compressible;
      dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;sti1=sti;
  
  /* create threads and wait */
  
      ithread=NNEW(int,num_cpus);
      for(i=0; i<num_cpus; i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)mafillv1rhsmt, (void *)&ithread[i]);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqv;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqv];
	  }
      }
      RENEW(b,double,neqv);free(ithread);

      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
	      
      free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv1,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
      free(sol);

      /*  if((iit/jout[1])*jout[1]==iit){
	  nnstep=1;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
	     &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
			    physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn));
			    }*/

      /* STEP 2: pressure correction */

      b=NNEW(double,num_cpus*neqp);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xbounact1=xbounact;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
      nelemface1=nelemface;sideface1=sideface;nface1=nface;
      nactdoh1=nactdoh;icolp1=icolp;jqp1=jqp;irowp1=irowp;neqp1=neqp;
      nzlp1=nzlp;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;vold1=vold;vcon1=vcon;nzsp1=nzsp;dtimef1=dtimef;
      matname1=matname;mi1=mi;ncmat_1=ncmat_;shcon1=shcon;nshcon1=nshcon;
      v1=v;theta11=theta1;iexplicit1=iexplicit;physcon1=physcon;
      dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
  
      ithread=NNEW(int,num_cpus);
      for(i=0; i<num_cpus; i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)mafillprhsmt, (void *)&ithread[i]);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqp;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqp];
	  }
      }
      RENEW(b,double,neqp);free(ithread);

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
	      pardiso_solve(b,&neqp,&symmetryflag);
#endif
	  }

          /* copying the solution into field sol */

	  for(i=0;i<neqp;i++){
	      sol[i]=b[i]/(theta1*theta2*dtimef*dtimef);
	  }
	  free(b);

      }

      /* storing the pressure (incompressible) or density (compressible)
         correction in v */

      FORTRAN(resultsp,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,
              mi));
      free(sol);

      if(iexplicit==0){
      
      /* inserting the pressure boundary conditions for liquids */

          FORTRAN(applybounp,(nodeboun,ndirboun,nboun,xbounact,
       ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
       nsolidsurf,xsolidsurf,nfreestream,ifreestream,iturbulent,
       vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
       ipompc,nodempc,coefmpc,nmpc,inomat,mi));
      }

      /*     if((iit/jout[1])*jout[1]==iit){
	  nnstep=2;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
	     &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn));
		  }*/
      
      /* STEP 3: velocity correction */

/*      printf("STEP3: velocity correction\n\n");*/

      b=NNEW(double,num_cpus*neqv);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
      nactdoh1=nactdoh;icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;
      nzlv1=nzlv;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;vold1=vold;nzsv1=nzsv;dtimef1=dtimef;v1=v;
      theta21=theta2;iexplicit1=iexplicit;mi1=mi;
      dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
  
  /* create threads and wait */
  
      ithread=NNEW(int,num_cpus);
      for(i=0; i<num_cpus; i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)mafillv2rhsmt, (void *)&ithread[i]);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

      for(i=0;i<neqv;i++){
	  for(j=1;j<num_cpus;j++){
	      b[i]+=b[i+j*neqv];
	  }
      }
      RENEW(b,double,neqv);free(ithread);

      sol=NNEW(double,neqv);
      aux=NNEW(double,neqv);
      FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
	 &neqv,&nzsv,&nzlv));
      free(b);free(aux);

      /* storing the velocity correction in v */

      FORTRAN(resultsv2,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
      free(sol);

      /*   if((iit/jout[1])*jout[1]==iit){
	  nnstep=3;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		&kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn));
		  }*/

      /* STEP 4: energy correction */

      if(*ithermal>1){


	   b=NNEW(double,num_cpus*neqt);

	  co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
	  nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
	  ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
	  nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
	  nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
	  xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
	  neqt1=neqt;
	  nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
	  ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
	  ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;vcon1=vcon;
	  nzst1=nzst;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
	  physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
	  timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
	  reltimef1=reltimef;cocon1=cocon;ncocon1=ncocon;nelemface1=nelemface;
          sideface1=sideface;nface1=nface;compressible1=compressible;
          vcontu1=vcontu;yy1=yy;iturbulent1=iturbulent;
	  dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
	  
	  ithread=NNEW(int,num_cpus);
	  for(i=0; i<num_cpus; i++)  {
	      ithread[i]=i;
	      pthread_create(&tid[i], NULL, (void *)mafilltrhsmt, (void *)&ithread[i]);
	  }
	  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	  
	  for(i=0;i<neqt;i++){
	      for(j=1;j<num_cpus;j++){
		  b[i]+=b[i+j*neqt];
	      }
	  }
	  RENEW(b,double,neqt);free(ithread);
	  
	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  FORTRAN(solveeq,(adbt,aubt,adlt,addiv,b,sol,aux,icolt,irowt,jqt,
			   &neqt,&nzst,&nzlt));
	  free(b);free(aux);
	  
      /* storing the temperature correction in v */

	  FORTRAN(resultst,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
	  free(sol);
      }

      /* if((iit/jout[1])*jout[1]==iit){
	  nnstep=4;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn));
		  }*/

      /* STEP 5: turbulent correction */

      if(*iturbulent!=0){
/*	        if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn));
		  }*/

	  bk=NNEW(double,num_cpus*neqk);
	  bt=NNEW(double,num_cpus*neqk);

	  co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
	  nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
	  ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
	  nelemface1=nelemface;sideface1=sideface;nface1=nface;
          nactdok1=nactdok;neqk1=neqk;
	  nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
	  ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
	  ntmat_1=ntmat_;vold1=vold;vcon1=vcon;
	  nzsk1=nzsk;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
	  shcon1=shcon;nshcon1=nshcon;theta11=theta1;
          vcontu1=vcontu;isolidsurf1=isolidsurf;nsolidsurf1=nsolidsurf;
          ifreestream1=ifreestream;nfreestream1=nfreestream;
          xsolidsurf1=xsolidsurf;yy1=yy;compressible1=compressible;
          iturbulent1=iturbulent;ithermal1=ithermal;ipvar1=ipvar;var1=var;
          ipvarf1=ipvarf;varf1=varf;dt1=dt;
	  
	  ithread=NNEW(int,num_cpus);
	  for(i=0; i<num_cpus; i++)  {
	      ithread[i]=i;
	      pthread_create(&tid[i], NULL, (void *)mafillkrhsmt, (void *)&ithread[i]);
	  }
	  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	  
	  for(i=0;i<neqk;i++){
	      for(j=1;j<num_cpus;j++){
		  bk[i]+=bk[i+j*neqk];
		  bt[i]+=bt[i+j*neqk];
	      }
	  }
	  RENEW(bk,double,neqk);
	  RENEW(bt,double,neqk);free(ithread);
	  
	  solk=NNEW(double,neqk);
	  aux=NNEW(double,neqk);
	  FORTRAN(solveeq,(adbv,aubv,adlk,addiv,bk,solk,aux,icolk,irowk,jqk,
			   &neqk,&nzsk,&nzlk));
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

	  /*       if((iit/jout[1])*jout[1]==iit){
	  nnstep=5;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn));
		  }*/
      }
      
      /* update the conservative variables
         (for incompressible fluids: pressure instead of density  */

      FORTRAN(updatecon,(vold,vcon,v,nk,
	      ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	      nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
	      nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
	      inomat,nodeboun,ndirboun,nboun,mi,co,&factor));

      /* inserting the boundary conditions for the turbulence
         parameters */

      /*   if(*iturbulent!=0){
	  FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
	    iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
	    ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	    inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
	    rhcon,vcontu,ntmat_,labmpc,inomat));
	    }*/
 
          /* extrapolating the stresses (for debugging purposes) */

/*	      nfield=6;
	      ndim=6;
	      cfd=1;
	      if((*norien>0)&&(strcmp1(&filab[179],"L")==0)){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[178],1);
	      stn=NNEW(double,6**nk);
	      inum=NNEW(int,*nk);
	      FORTRAN(extrapolate,(sti,stn,ipkon,inum,kon,lakon,
		      &nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
		      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
		      vold,ithermal,&force,&cfd));*/
     
      /* smoothing the solution (only for compressible fluids) */

      if((compressible)&&(shockcoef>0.)){

	       /*    ismooth=1;*/

	  /* shocksmoothing rho * total energy density */

	  sol=NNEW(double,neqt);
	  aux=NNEW(double,neqt);
	  for(i=0;i<neqt;i++){sol[i]=vcon[mt*i];}
	  FORTRAN(smoothshock,(adbt,aubt,adlt,addiv,sol,aux,icolt,irowt,jqt,
			  &neqt,&nzlt,sa));
	  for(i=0;i<neqt;i++){vcon[mt*i]=sol[i];}
	  free(sol);free(aux);

	  /* shocksmoothing rho * velocity */

	  sol=NNEW(double,neqv);
	  aux=NNEW(double,neqv);
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  sol[3*i+j]=vcon[mt*i+j+1];
	      }
	  }
	  FORTRAN(smoothshock,(adbv,aubv,adlv,addiv,sol,aux,icolv,irowv,jqv,
			  &neqv,&nzlv,sav));
	  for(i=0;i<neqv/3;i++){
	      for(j=0;j<3;j++){
		  vcon[mt*i+j+1]=sol[3*i+j];
	      }
	  }
	  free(sol);free(aux);

	  /* shocksmoothing rho */

	  sol=NNEW(double,neqp);
	  aux=NNEW(double,neqp);
	  for(i=0;i<neqp;i++){sol[i]=vcon[mt*i+4];}
	  FORTRAN(smoothshock,(adbp,aubp,adlp,addiv,sol,aux,icolp,irowp,jqp,
			  &neqp,&nzlp,sa));
	  for(i=0;i<neqp;i++){vcon[mt*i+4]=sol[i];}
	  free(sol);free(aux);
	  
      }

      /* deriving the physical variables from the conservative
         variables */

      FORTRAN(con2phys,(vold,vcon,v,nk,
	      ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	      nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
	      nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
	      inomat,nodeboun,ndirboun,nboun,mi,co,&factor));
      
      /* inserting boundary conditions to the physical variables */
      
      FORTRAN(applyboun,(nodeboun,ndirboun,nboun,xbounact,
	     ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
	     nsolidsurf,xsolidsurf,nfreestream,ifreestream,iturbulent,
	     vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
	     &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,mi,
	     ikboun,ilboun,ilmpc,labmpc));

      /* calculating the pressure gradient for the shock smoothing in
         the next iteration */

      if((compressible)&&(shockcoef>0.)){
	  FORTRAN(presgradient,(iponoel,inoel,sa,sav,&neqt,dt,&shockcoef,
				&dtimef,ipkon,kon,lakon,vold,mi,
                                &compressible,nmethod,dt,isolidsurf,
                                nsolidsurf,co,&euler));
      }
      
      /* inserting the boundary conditions for the turbulence
	 parameters */
      
      /*    if(*iturbulent!=0){
	  FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
		 iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
		 ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	         inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
		 rhcon,vcontu,ntmat_,labmpc,inomat,mi,ithermal));
		 }*/

      /* check iconvergence */

      FORTRAN(cfdconv,(vold,vcon,v,nk,
	      ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
	      nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
	      nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
	      inomat,nodeboun,ndirboun,nboun,mi,co,&factor,
	      vconini,&dtimef,del,&sum,&sumx,&sumxx,sumy,sumxy,&nstart,
              &shockcoef));
      
      if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)||
	 (iit==jmax[1])){

	  /* calculating the stress and the heat flow at the
             integration points, if requested */

	  if(strcmp1(&filab[3306],"SF  ")==0)isti=1;
          if(strcmp1(&filab[3393],"HFLF")==0)iqfx=1;
	  for(i=0;i<*nprint;i++){
	      if(strcmp1(&prlab[6*i],"SF")==0) isti=1;
	      if(strcmp1(&prlab[6*i],"HFLF")==0)iqfx=1;
	  }
	  if((isti==1)||(iqfx==1)){
	      FORTRAN(calcstressheatflux,(kon,lakon,ipkon,ielmat,ntmat_,
		      vold,matname,mi,shcon,nshcon,iturbulent,&compressible,
                      ipvar,var,sti,qfx,cocon,ncocon,ne,&isti,&iqfx));
	  }
 
          /* extrapolating the stresses */

	  if(strcmp1(&filab[3306],"SF  ")==0){
	      nfield=6;
	      ndim=6;
	      if((*norien>0)&&(strcmp1(&filab[2962],"L")==0)){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[2962],1);
	      stn=NNEW(double,6**nk);
	      inum=NNEW(int,*nk);
	      FORTRAN(extrapolate,(sti,stn,ipkon,inum,kon,lakon,
		      &nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
		      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
		      vold,ithermal,&force,&icfdout,ielmat,thicke,filab));
	      free(inum);
	  }

	  /* extrapolating the heat flow */

	  
	  if(strcmp1(&filab[3393],"HFLF")==0){
	      nfield=3;
	      ndim=3;
	      if((*norien>0)&&(strcmp1(&filab[3049],"L")==0)){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[3049],1);
	      qfn=NNEW(double,3**nk);
	      inum=NNEW(int,*nk);
	      FORTRAN(extrapolate,(qfx,qfn,ipkon,inum,kon,lakon,
		      &nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
		      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
		      vold,ithermal,&force,&icfdout,ielmat,thicke,filab));
	      free(inum);
	  }
	 
	  /* check whether the Mach number is requested */

	  if((strcmp1(&filab[1914],"MACH")==0)|| 
             (strcmp1(&filab[3219],"TTF")==0)||
             (strcmp1(&filab[3132],"PTF")==0)){
	      FORTRAN(calcmach,(vold,vcon,v,nk,
		      ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
		      nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
		      nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
		      inomat,nodeboun,ndirboun,nboun,mi,co,&factor));
	  }

          /* print output */

	  if(iconvergence==1) timef=*time;
	  FORTRAN(printoutfluid,(set,nset,istartset,iendset,ialset,nprint,
	    prlab,prset,vold,t1,fn,ipkon,lakon,sti,eme,xstate,ener,
	    mi,&nstate_,ithermal,co,kon,qfx,&timef,trab,inotr,ntrans,
	    orab,ielorien,norien,nk,ne,inum,filab,vold,&ikin,ielmat,thicke,
            eme,vcontu,physcon));

	  /* lift and drag force */

	  /*  FORTRAN(printoutface,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
	    cocon,ncocon,&compressible,istartset,iendset,ipkon,lakon,kon,
	    ialset,prset,&timef,nset,set,nprint,prlab,ielmat,mi));*/

          /* frd output */

	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));

	  if(strcmp1(&filab[3306],"SF  ")==0){free(stn);}
	  if(strcmp1(&filab[3393],"HFLF")==0){free(qfn);}

      }
      
      /*  if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)){
	  FORTRAN(printoutfluid,(set,nset,istartset,iendset,ialset,nprint,
	    prlab,prset,vold,t1,fn,ipkon,lakon,sti,eme,xstate,ener,
	    mi,&nstate_,ithermal,co,kon,qfx,&timef,trab,inotr,ntrans,
	    orab,ielorien,norien,nk,ne,inum,filab,vold,&ikin,ielmat,thicke,
            eme,vcontu,physcon));*/

	  /* lift and drag force */
	  
/*	  FORTRAN(printoutface,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
	    cocon,ncocon,&compressible,istartset,iendset,ipkon,lakon,kon,
	    ialset,prset,&timef,nset,set,nprint,prlab,ielmat,mi));
	    }*/
      
      /*   if(iconvergence==1){
	  nnstep=6;
	  FORTRAN(frdfluid,(co,nk,kon,ipkon,lakon,ne,v,vold,
		  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
		  physcon,filab,inomat,ntrans,inotr,trab,mi,stn));
	  break;
      }*/

//      if((iit==jmax[1])||(iconvergence==1)) FORTRAN(stop,());
      if((iit==jmax[1])||(iconvergence==1)) break;
      
      ttimef+=dtimef;
  }while(1);
  
  if((iexplicit!=1)&&(neqp>0)){
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
	  pardiso_cleanup(&neqp,&symmetryflag);
#endif
      }
  }

  if(compressible){free(sa);free(sav);free(del);}

  free(yy);free(xsolidsurf);free(dt);free(dh);free(vcon);free(vconini);

  free(irowv);free(irowp);
  free(icolv);free(icolp);
  free(jqv);free(jqp);
  free(nactdoh);

  free(adbv);free(adbp);
  free(aubv);free(aubp);
  free(adlv);if(iexplicit==1) free(adlp);

  if(*ithermal>1){
      free(irowt);free(icolt);free(jqt);free(adbt);free(aubt);free(adlt);
  }

  if(*iturbulent!=0){
      free(irowk);free(icolk);free(jqk);free(nactdok);
      free(adbk);free(aubk);free(adlk);
  }

  free(v);free(vtu);free(var);free(ipvar);free(varf);free(ipvarf);

  if(*mcs==1){
    *memmpc_=memmpcref_;*mpcfree=mpcfreeref;*nmpc=nmpcref;
    *nk=nkref;*ne=neref;*nef=nefref;

    RENEW(nodempc,int,3*memmpcref_);
    for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
    RENEW(coefmpc,double,memmpcref_);
    for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
    
    RENEW(ipompc,int,*nmpc);
    RENEW(ikmpc,int,*nmpc);
    RENEW(ilmpc,int,*nmpc);
    RENEW(labmpc,char,20**nmpc);
    
    RENEW(co,double,3**nk);
    RENEW(vold,double,mt**nk);
    RENEW(ipkon,int,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(kon,int,*nkon);
    RENEW(ielmat,int,mi[2]**ne);
    RENEW(inomat,int,*nk);

    sideface=NNEW(char,6**nef);
    nelemface=NNEW(int,6**nef);
    iponoel=NNEW(int,*nk);
    inoel=NNEW(int,3*20**nef);
  }

  *ithermal=ithermalref;

  *nodempcp=nodempc;*ipompcp=ipompc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*fmpcp=fmpc;*cop=co;
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielmatp=ielmat;
  *nelemfacep=nelemface;*sidefacep=sideface;*inoelp=inoel;
  *iponoelp=iponoel;*voldp=vold;*inomatp=inomat;
  
  //  FORTRAN(stop,());
  FORTRAN(closefilefluid,());

  return;
  
} 

/* subroutine for multithreading of mafillv1rhs */

void *mafillv1rhsmt(int *i){

    int index,nea,neb,nedelta;

    index=*i*neqv1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillv1rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
	 xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,ndirforc1,xforc1,
	 nforc1,nelemload1,sideload1,xload1,nload1,xbody1,ipobody1,nbody1,
	 &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,nmethod1,ikmpc1,ilmpc1,ikboun1,
         ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,vold1,vcon1,&nzsv1,
         dt1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,&timef1,
         istep1,iinc1,ibody1,xloadold1,iturbulent1,vcontu1,yy1,
	 nelemface1,sideface1,nface1,&compressible1,&nea,&neb,&dtimef1,
	 ipvar1,var1,ipvarf1,varf1,sti1));

    return NULL;
}

/* subroutine for multithreading of mafillprhs */

void *mafillprhsmt(int *i){

    int index,nea,neb,nedelta;

    index=*i*neqp1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillprhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
       xbounact1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nelemface1,sideface1,
       nface1,&b[index],nactdoh1,icolp1,jqp1,irowp1,&neqp1,&nzlp1,nmethod1,ikmpc1,ilmpc1,
       ikboun1,ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,vold1,vcon1,&nzsp1,
       dt1,matname1,mi1,ncmat_1,shcon1,nshcon1,v1,&theta11,
       &iexplicit1,physcon1,&nea,&neb,&dtimef1,ipvar1,var1,ipvarf1,varf1));

    return NULL;
}

/* subroutine for multithreading of mafillv2rhs */

void *mafillv2rhsmt(int *i){

    int index,nea,neb,nedelta;

    index=*i*neqv1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if(neb>*ne1) neb=*ne1;

    FORTRAN(mafillv2rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
       xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,
       &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,nmethod1,ikmpc1,ilmpc1,ikboun1,
       ilboun1,vold1,&nzsv1,dt1,v1,&theta21,&iexplicit1,&nea,&neb,mi1,&dtimef1,
       ipvar1,var1,ipvarf1,varf1));

    return NULL;
}

/* subroutine for multithreading of mafilltrhs */

void *mafilltrhsmt(int *i){

    int index,nea,neb,nedelta;

    index=*i*neqt1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if(neb>*ne1) neb=*ne1;
	  
    FORTRAN(mafilltrhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
	     xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,ndirforc1,
             xforc1,
             nforc1,nelemload1,sideload1,xload1,nload1,xbody1,ipobody1,nbody1,
             &b[index],nactdoh1,&neqt1,nmethod1,ikmpc1,ilmpc1,ikboun1,
             ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,vold1,
             vcon1,&nzst1,
             dt1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,&timef1,
             istep1,iinc1,ibody1,xloadold1,&reltimef1,cocon1,ncocon1,nelemface1,
	     sideface1,nface1,&compressible1,vcontu1,yy1,iturbulent1,&nea,
	     &neb,&dtimef1,ipvar1,var1,ipvarf1,varf1));

    return NULL;
}

/* subroutine for multithreading of mafillkrhs */

void *mafillkrhsmt(int *i){

    int index,nea,neb,nedelta;

    index=*i*neqk1;
    
    nedelta=(int)ceil(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if(neb>*ne1) neb=*ne1;
	  
    FORTRAN(mafillkrhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
	    xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nelemface1,sideface1,
	    nface1,nactdok1,&neqk1,nmethod1,ikmpc1,ilmpc1,
	    ikboun1,ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,vold1,vcon1,
            &nzsk1,
	    &dtimef1,matname1,mi1,ncmat_1,shcon1,nshcon1,&theta11,
	    &bk[index],&bt[index],vcontu1,isolidsurf1,nsolidsurf1,
            ifreestream1,nfreestream1,
	    xsolidsurf1,yy1,&compressible1,iturbulent1,ithermal1,ipvar1,var1,
	    ipvarf1,varf1,&nea,&neb,dt1));

    return NULL;
}







