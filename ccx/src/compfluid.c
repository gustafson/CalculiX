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

static ITG num_cpus;

void compfluid(double **cop, ITG *nk, ITG **ipkonp, ITG **konp, char **lakonp,
    ITG *ne, char **sidefacep, ITG *ifreestream, 
    ITG *nfreestream, ITG *isolidsurf, ITG *neighsolidsurf,
    ITG *nsolidsurf, ITG **iponoelp, ITG **inoelp, ITG *nshcon, double *shcon,
    ITG *nrhcon, double *rhcon, double **voldp, ITG *ntmat_,ITG *nodeboun, 
    ITG *ndirboun, ITG *nboun, ITG **ipompcp,ITG **nodempcp, ITG *nmpc,
    ITG **ikmpcp, ITG **ilmpcp, ITG *ithermal, ITG *ikboun, ITG *ilboun,
    ITG *iturbulent, ITG *isolver, ITG *iexpl, double *vcontu, double *ttime,
    double *time, double *dtime, ITG *nodeforc,ITG *ndirforc,double *xforc,
    ITG *nforc, ITG *nelemload, char *sideload, double *xload,ITG *nload,
    double *xbody,ITG *ipobody,ITG *nbody, ITG **ielmatp, char *matname,
    ITG *mi, ITG *ncmat_, double *physcon, ITG *istep, ITG *iinc,
    ITG *ibody, double *xloadold, double *xboun,
    double **coefmpcp, ITG *nmethod, double *xforcold, double *xforcact,
    ITG *iamforc,ITG *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, ITG *iamt1, double *amta,
    ITG *namta, ITG *nam, double *ampli, double *xbounold, double *xbounact,
    ITG *iamboun, ITG *itg, ITG *ntg, char *amname, double *t0, 
    ITG **nelemfacep,
    ITG *nface, double *cocon, ITG *ncocon, double *xloadact, double *tper,
    ITG *jmax, ITG *jout, char *set, ITG *nset, ITG *istartset,
    ITG *iendset, ITG *ialset, char *prset, char *prlab, ITG *nprint,
    double *trab, ITG *inotr, ITG *ntrans, char *filab, char **labmpcp, 
    double *sti, ITG *norien, double *orab, char *jobnamef,char *tieset,
    ITG *ntie, ITG *mcs, ITG *ics, double *cs, ITG *nkon, ITG *mpcfree,
    ITG *memmpc_,double **fmpcp,ITG *nef,ITG **inomatp,double *qfx,
    ITG *neifa,ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
    ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun){

    /* main computational fluid dynamics routine */
  
  char *labmpc=NULL,*lakon=NULL,*sideface=NULL;

  ITG *ipointer=NULL,*mast1=NULL,*irow=NULL,*icol=NULL,*jq=NULL,
      *nactdoh=NULL,nzs=20000000,neq,kode,compressible,mt=mi[1]+1,*ifabou=NULL,
      *nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,nfabou,
      *ipkon=NULL,*kon=NULL,*ielmat=NULL,*nelemface=NULL,*inoel=NULL,
      *iponoel=NULL,*inomat=NULL,ithermalref,*integerglob=NULL,iit,
      iconvergence,symmetryflag=2,inputformat=1,i,*neij=NULL;

  double *coefmpc=NULL,*fmpc=NULL,*umfa=NULL,reltime,*doubleglob=NULL,
      *co=NULL,*vold=NULL,*coel=NULL,*cosa=NULL,*gradv=NULL,*gradvfa=NULL,
      *xxn=NULL,*xxi=NULL,*xle=NULL,*xlen=NULL,*xlet=NULL,timef,dtimef,
      *cofa=NULL,*area=NULL,*xrlfa=NULL,reltimef,ttimef,sumfix,sumfree,
      *au=NULL,*ad=NULL,*b=NULL,*volume=NULL,*body=NULL,sigma=0.,
      *adb=NULL,*aub=NULL,*adfa=NULL,*ap=NULL,*ppel=NULL,*ppfa=NULL,
      *gradpp=NULL;

  nodempc=*nodempcp;ipompc=*ipompcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;fmpc=*fmpcp;co=*cop;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielmat=*ielmatp;
  nelemface=*nelemfacep;sideface=*sidefacep;inoel=*inoelp;
  iponoel=*iponoelp;vold=*voldp;inomat=*inomatp;

  /* standard: shockcoef=0 */

#ifdef SGI
  ITG token;
#endif

  /* open frd-file for fluids */

  FORTRAN(openfilefluid,(jobnamef));

  /* variables for multithreading procedure */

  ITG sys_cpus;
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
  
  printf(" Using up to %" ITGFORMAT " cpu(s) for CFD.\n", num_cpus);
  
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
  
  ipointer=NNEW(ITG,3**nk);
  mast1=NNEW(ITG,nzs);
  irow=NNEW(ITG,nzs);
  icol=NNEW(ITG,3**nk);
  jq=NNEW(ITG,3**nk+1);
  nactdoh=NNEW(ITG,mt**nk);

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

  ifabou=NNEW(ITG,7**nfaext);
  FORTRAN(applyboun,(ifaext,nfaext,ielfa,ikboun,ilboun,
       nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
       ifabou,&nfabou));
  RENEW(ifabou,ITG,nfabou);

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

  if(*nbody>0) body=NNEW(double,3**ne);

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

      /* body forces (gravity, centrifugal and Coriolis forces */

      if(*nbody>0){
	  FORTRAN(calcbody,(ne,body,ipobody,ibody,xbody,coel,vel,lakon));
      }
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
  adfa=NNEW(double,*nface);
  ap=NNEW(double,*nface);

  au=NNEW(double,nzs);
  ad=NNEW(double,neq);
  b=NNEW(double,3*neq);

  ppel=NNEW(double,*ne);
  ppfa=NNEW(double,*nface);
  gradpp=NNEW(double,3**ne);

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

	  /* body forces (gravity, centrifugal and Coriolis forces */
	  
	  if(*nbody>0){
	      FORTRAN(calcbody,(ne,body,ipobody,ibody,xbody,coel,vel,lakon));
	  }
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
	  body,volume,&compressible,ielfa,lakon,ifabou,nbody));

      /* LU decomposition */

      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,&nzs);
#else
	  printf("*ERROR in arpack: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,jq,&nzs);
#else
	  printf("*ERROR in arpack: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }

      /* solving the system of equations (for x, y and z
         separately */

      for(i=0;i<3;i++){

        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(&b[i*neq],&neq);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(&b[i*neq],&neq,&symmetryflag);
#endif
        }
      }

      /* free memory */
      
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_cleanup(&neq,&symmetryflag);
#endif
      }

      /* storing the solution into vel */

      FORTRAN(calcvel,(nk,nactdoh,vel,b,&neq));

      /* generating 1/ad at the face centers */

      FORTRAN(extrapolate_ad,(nface,ielfa,xrlfa,ad,adfa,nactdoh));

      /* extrapolating the velocity from the elements centers to the face
	 centers, thereby taking the boundary conditions into account */
      
      FORTRAN(extrapolate_v,(nface,ielfa,xrlfa,vel,vfa,
			     ifabou,xboun,&sumfix,&sumfree,xxn,area));
      
      /*   modifying the velocity values along the boundary which are not fixed
	   by boundary conditions such that the mass conservation holds */
      
      
      FORTRAN(integral_boundary,(&sumfix,&sumfree,ifaext,nfaext,ielfa,
				 ifabou,vfa));

      /* calculating the lhs and rhs of the equation system to determine
         p' (balance of mass) */

      FORTRAN(mafillp,(ne,lakon,nactdoh,ipnei,neifa,neiel,vfa,area,
              adfa,xlet,cosa,volume,au,ad,jq,irow,&nzs,ap,ielfa,ifabou,xle,
	      b,xxn,&compressible));

      /* LU decomposition of the p' system*/

      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,&nzs);
#else
	  printf("*ERROR in arpack: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,jq,&nzs);
#else
	  printf("*ERROR in arpack: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      
      /* solving the system of equations for p'  */
      
      if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(b,&neq);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(b,&neq,&symmetryflag);
#endif
      }

      /* storing the solution p' into ppel */

      FORTRAN(calcppel,(ne,nactdoh,ppel,b));

      /* extrapolating the p' from the elements centers to the face
	 centers  */
      
      FORTRAN(extrapolate_pp,(nface,ielfa,xrlfa,ppel,ppfa,ifabou));

      /* calculation of the gradient of p' at the center
	 of the elements from the p' values at the neighboring
	 faces */

      FORTRAN(calcgradpp,(ne,lakon,ipnei,ppfa,area,xxn,gradpp,neifa));

      /* calculating the right hand side of the equations to determine p'' */

      FORTRAN(mafillppprhs,(ne,lakon,nactdoh,ipnei,neifa,neiel,neij,
	      gradpp,xxi,cosa,xxn,ap,xle,xlen,b,ielfa,ifabou));




  }while(1);
  
  FORTRAN(closefilefluid,());

  free(irow);free(icol);free(jq);free(nactdoh);
  
  free(coel);free(cosa);free(xxn);free(xxi);free(xle);free(xlen);
  free(xlet);free(cofa);free(area);free(xrlfa);free(volume);

  free(ifabou);free(umfa);

  free(gradv);free(gradvfa);free(au);free(ad);free(b);free(adfa);
  free(ap);free(ppel);free(ppfa);free(gradpp);

  if(*nbody>0) free(body);

  return;
  
} 
