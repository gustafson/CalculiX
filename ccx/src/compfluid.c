/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                     */

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

static int num_cpus;

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
    int *memmpc_,double **fmpcp,int *nef,int **inomatp,double *qfx,
    int *neifa,int *neiel,int *ielfa,int *ifaext,double *vfa,double *vel,
    int *ipnei,int *nflnei,int *nfaext,char *typeboun){

    /* main computational fluid dynamics routine */
  
  char *labmpc=NULL,*lakon=NULL,*sideface=NULL;

  int *ipointer=NULL,*mast1=NULL,*irow=NULL,*icol=NULL,*jq=NULL,
      *nactdoh=NULL,nzs,neq,kode,compressible,mt=mi[1]+1,*ifabou=NULL,
      *nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,nfabou,
      *ipkon=NULL,*kon=NULL,*ielmat=NULL,*nelemface=NULL,*inoel=NULL,
      *iponoel=NULL,*inomat=NULL,ithermalref,*integerglob=NULL,iit,
      iconvergence;

  double *coefmpc=NULL,*fmpc=NULL,*umfa=NULL,reltime,*doubleglob=NULL,
      *co=NULL,*vold=NULL,*coel=NULL,*cosa=NULL,*gradv=NULL,*gradvfa=NULL,
      *xxn=NULL,*xxi=NULL,*xle=NULL,*xlen=NULL,*xlet=NULL,timef,dtimef,
      *cofa=NULL,*area=NULL,*xrlfa=NULL,reltimef,ttimef,sumfix,sumfree,
      *au=NULL,*ad=NULL,*b=NULL,*volume=NULL,body[3];

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
  
  /*  *iexpl==0:  structure:implicit, fluid:incompressible
      *iexpl==1:  structure:implicit, fluid:compressible
      *iexpl==2:  structure:explicit, fluid:incompressible
      *iexpl==3:  structure:explicit, fluid:compressible */

  if((*iexpl==1)||(*iexpl==3)){
      compressible=1;
  }else{
      compressible=0;
  }

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  ithermalref=*ithermal;
  if(*ithermal==1){
    *ithermal=2;
  }

  /* determining the matrix structure */
  
  ipointer=NNEW(int,3**nk);
  mast1=NNEW(int,nzs);
  irow=NNEW(int,nzs);
  icol=NNEW(int,3**nk);
  jq=NNEW(int,3**nk+1);
  nactdoh=NNEW(int,mt**nk);

  mastructf(nk,kon,ipkon,lakon,ne,nactdoh,icol,jq,&mast1,&irow,
	    isolver,&neq,ipointer,&nzs,ipnei,neiel,mi);

  free(ipointer);free(mast1);
 
  /* calculation geometric data */

  coel=NNEW(double,3**ne);
  volume=NNEW(double,*ne);
  cosa=NNEW(double,3**nflnei);
  xxn=NNEW(double,3**nflnei);
  xxi=NNEW(double,3**nflnei);
  xle=NNEW(double,*nflnei);
  xlen=NNEW(double,*nflnei);
  xlet=NNEW(double,*nflnei);
  cofa=NNEW(double,3**nface);
  area=NNEW(double,*nface);
  xrlfa=NNEW(double,3**nface);

  FORTRAN(initialcfd,(ne,ipkon,kon,lakon,co,coel,cofa,nface,
	  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa));

  /* storing pointers to the boundary conditions in ielfa */

  ifabou=NNEW(int,7**nfaext);
  FORTRAN(applyboun,(ifaext,nfaext,ielfa,ikboun,ilboun,
       nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
       ifabou,&nfabou));
  RENEW(ifabou,int,nfabou);

  /* material properties for athermal calculations 
     = calculation for which no initial thermal conditions
     were defined */

  umfa=NNEW(double,*nface);

  if(*ithermal==0){
      
      /* calculating the density a the element centers */

      FORTRAN(calcrhoel,(ne,lakon,vel,rhcon,nrhcon,ielmat,ntmat_,
			 ithermal,mi));
      
      /* calculating the density at the face centers */

      FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmat,ntmat_,
			 ithermal,mi,ielfa));
      
      /* calculating the dynamic viscosity at the face centers */

      FORTRAN(calcdvifa,(nface,vfa,shcon,nshcon,ielmat,ntmat_,
			 ithermal,mi,ielfa,umfa));

  }

  /* next section is for stationary calculations */
  
  if(*nmethod==1){
      
      /* boundary conditions at the end of the mechanical
	 increment */
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	     xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltime,ttime,dtime,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	     ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
  }

  /* extrapolating the velocity from the elements centers to the face
     centers, thereby taking the boundary conditions into account */

  FORTRAN(extrapolate_v,(nface,ielfa,xrlfa,vel,vfa,
			 ifabou,xboun,&sumfix,&sumfree,xxn,area));

  /*   modifying the velocity values along the boundary which are not fixed
       by boundary conditions such that the mass conservation holds */


  FORTRAN(integral_boundary,(&sumfix,&sumfree,ifaext,nfaext,ielfa,
			     ifabou,vfa));

  /* start of the major loop */

  gradv=NNEW(double,9**ne);
  gradvfa=NNEW(double,9**nface);

  au=NNEW(double,nzs);
  ad=NNEW(double,*nef);
  b=NNEW(double,3**nef);

  iit=0;

  do{

      iit++;

      /* determining dtimef (to do ) */

      timef+=dtimef;
      if((*time<timef)&&(*nmethod==4)){
	  dtimef-=timef-*time;
	  timef=*time;
	  iconvergence=1;
      }
      reltimef=timef/(*tper);
      if(reltimef>1.) reltimef=1.;

      /* conditions for transient calculations */
     
      if(*nmethod==4){

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

      /* calculating the gradient of the velocity at the element
         centers */

      FORTRAN(calcgradv,(ne,lakon,ipnei,vfa,area,xxn,gradv,neifa));

      /* extrapolating the gradient of the velocity from the element
         centers to the face centers */

      FORTRAN(extrapolate_gradv,(nface,ielfa,xrlfa,gradv,gradvfa));

      /* filling the lhs and rhs's for the balance of momentum
         equations */

      FORTRAN(mafillv,(ne,nactdoh,ipnei,neifa,neiel,vfa,xxn,area,
          au,ad,jq,irow,&nzs,b,vel,cosa,umfa,xlet,xle,gradvfa,xxi,
	  body,volume,&compressible,ielfa,lakon,ifabou));






  }while(1);
  
  FORTRAN(closefilefluid,());

  free(irow);free(icol);free(jq);free(nactdoh);
  
  free(coel);free(cosa);free(xxn);free(xxi);free(xle);free(xlen);
  free(xlet);free(cofa);free(area);free(xrlfa);free(volume);

  free(ifabou);free(umfa);

  free(gradv);free(gradvfa);free(au);free(ad);free(b);

  return;
  
} 
