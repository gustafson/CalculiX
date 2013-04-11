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

void radflowload(int *itg,int *ieg,int *ntg,int *ntr,int *ntm,
                 double *ac,double *bc,int *nload,char *sideload,
                 int *nelemload,double *xloadact,char *lakon,int *ipiv,
                 int *ntmat_,double *vold,double *shcon,
                 int *nshcon,int *ipkon,int *kon,double *co,double *pmid,
                 double *e1,double *e2,double *e3,int *iptri,int *kontri,
                 int *ntri,int *nloadtr,double *tarea,double *tenv,
                 double *physcon,double *erad,double *f,double *dist,
                 int *idist,double *area,int *nflow,int *ikboun,
                 double *xbounact,int *nboun,int *ithermal,
                 int *iinc,int *iit,double *cs, int *mcs, int *inocs, 
                 int *ntrit,int *nk, double *fenv,int *istep,double *dtime,
                 double *ttime,double *time,int *ilboun,int *ikforc,
                 int *ilforc,double *xforcact,int *nforc,double *cam,
                 int *ielmat,int *nteq,double *prop,int *ielprop,int *nactdog,
                 int *nacteq,int *nodeboun,int *ndirboun,
                 int *network, double *rhcon, int *nrhcon, int *ipobody,
                 int *ibody, double *xbodyact, int *nbody,int *iviewfile,
                 char *jobnamef, double *ctrl, double *xloadold,
                 double *reltime, int *nmethod, char *set){
  
  /* network=0: purely thermal
     network=1: general case (temperatures, fluxes and pressures unknown)
     network=2: purely aerodynamic, i.e. only fluxes and pressures unknown */
  
  int nhrs=1,info=0,i,iin=0,symmetryflag=2,inputformat=2,*icol=NULL,
      *irow=NULL,icntrl,icutb=0,iin_abs=0;
  double uamt=0,uamf=0,uamp=0,camt,camf,camp,*au=NULL,*adb=NULL,
    *aub=NULL,sigma=0.,ramt=0.,ramf=0.,ramp=0.,ram1t=0.,ram1f=0.,ram1p=0.,
    ram2t=0.,ram2f=0.,ram2p=0.,dtheta=1.,*v=NULL;
  
  /* check whether there are any gas temperature nodes; this check should
     NOT be done on nteq, since also for zero equations the temperature
     of the gas nodes with boundary conditions must be stored in v
     (in initialgas) */ 

  v=NNEW(double,5**nk);

  if(*ntg!=0) {
      icntrl=0;
      while(icntrl==0) {
	  
	  if(iin==0){
	      
	      for(i=0;i<5**nk;i++) v[i]=vold[i];

	      FORTRAN(initialgas,(itg,ieg,ntm,ntg,ac,bc,lakon,v,
                           ipkon,kon,nflow,
			   ikboun,nboun,prop,ielprop,nactdog,ndirboun,
			   nodeboun,xbounact,ielmat,ntmat_,shcon,nshcon,
			   physcon,ipiv,nteq,rhcon,nrhcon,ipobody,ibody,
			   xbodyact,co,nbody,network,&iin_abs,vold,set,
			   istep,iit));
      
	      FORTRAN(resultgas,(itg,ieg,ntg,ntm,bc,nload,sideload,
			  nelemload,xloadact,
			  lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,nflow,
			  iinc,istep,dtime,ttime,time,
			  ikforc,ilforc,xforcact,
                          nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,&iin,
			  physcon,&camt,&camf,&camp,rhcon,nrhcon,ipobody,
			  ibody,xbodyact,nbody,&dtheta,vold,xloadold,
			  reltime,nmethod,set));
	  }
	  
	  iin++;
	  iin_abs++;
	  printf("      gas iteration %d \n \n",iin);
	  
	  FORTRAN(mafillgas,(itg,ieg,ntg,ntm,ac,nload,sideload,
			     nelemload,xloadact,lakon,ntmat_,v,
			     shcon,nshcon,ipkon,kon,co,nflow,iinc,
			     istep,dtime,ttime,time,
			     ielmat,nteq,prop,ielprop,nactdog,nacteq,
			     physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
			     nbody,vold,xloadold,reltime,nmethod,set));
	  
	  if(*ntm>0){
	      FORTRAN(dgesv,(nteq,&nhrs,ac,ntm,ipiv,bc,ntm,&info)); 
	  }

	    /*spooles(ac,au,adb,aub,&sigma,bc,icol,irow,nteq,ntm,
	      &symmetryflag,&inputformat);*/
	  
	  if (info!=0) {
	      printf(" *WARNING in radflowload: singular matrix\n");
	    
	      FORTRAN(mafillgas,(itg,ieg,ntg,ntm,ac,nload,sideload,
				 nelemload,xloadact,lakon,ntmat_,v,
				 shcon,nshcon,ipkon,kon,co,nflow,iinc,
				 istep,dtime,ttime,time,
				 ielmat,nteq,prop,ielprop,nactdog,nacteq,
				 physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
				 nbody,vold,xloadold,reltime,nmethod,set));
	    
	      FORTRAN(equationcheck,(ac,ntm,nteq,nactdog,itg,ntg,nacteq,network));
	    
	      iin=0;

	  }
	  else {
	      FORTRAN(resultgas,(itg,ieg,ntg,ntm,bc,nload,sideload,nelemload,
	       xloadact,lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,
	       nflow,iinc,istep,dtime,ttime,time,
	       ikforc,ilforc,xforcact,
	       nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,
	       &iin,physcon,&camt,&camf,&camp,rhcon,nrhcon,ipobody,
	       ibody,xbodyact,nbody,&dtheta,vold,xloadold,
	       reltime,nmethod,set));
	    
	      if(*network!=2){ 
		  ram2t=ram1t;
		  ram1t=ramt;
		  ramt=camt;
		  if (camt>uamt) {uamt=camt;}
		  printf
                    ("      largest increment of gas temperature=%e\n",uamt);
		  printf
                    ("      largest correction to gas temperature=%e\n",camt);
	      }
	      
	      if(*network!=0){
		  ram2f=ram1f;
		  ram1f=ramf;
		  ramf=camf;
		  if (camf>uamf) {uamf=camf;}
		  printf("      largest increment of gas massflow=%e\n",uamf);
		  printf("      largest correction to gas massflow=%e\n",camf);
		  
		  ram2p=ram1p;
		  ram1p=ramp;
		  ramp=camp;
		  if (camp>uamp) {uamp=camp;}
		  printf("      largest increment of gas pressure=%e\n",uamp);
		  printf("      largest correction to gas pressure=%e\n",camp);
	      }
	      
	  }
	  
	  printf("\n");
	  
	  /* for purely thermal calculations no iterations are
	     deemed necessary */
	  
	  if(*network==0) {icntrl=1;}
	  else {
	      checkconvgas (&icutb,&iin,&uamt,&uamf,&uamp,
		 &ram1t,&ram1f,&ram1p,&ram2t,&ram2f,&ram2p,&ramt,&ramf,
		 &ramp,&icntrl,&dtheta,ctrl);
	  }
      }
	
      FORTRAN(flowresult,(ntg,itg,cam,vold,v,nload,sideload,
			    nelemload,xloadact,nactdog,network));
#ifdef NETWORKOUT

      FORTRAN(flowoutput,(itg,ieg,ntg,ntm,bc,lakon,ntmat_,
			  v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
			  ielmat,prop,ielprop,nactdog,nacteq,&iin,physcon,
			  &camt,&camf,&camp,&uamt,&uamf,&uamp,rhcon,nrhcon,
			  vold,jobnamef,set));
#endif
  }
      
  if(*ntr>0){
	
      FORTRAN(radmatrix, (ntr,ntm,ac,bc,sideload,nelemload,xloadact,lakon,
			    vold,ipkon,kon,co,pmid,e1,e2,e3,iptri,kontri,ntri,
			    nloadtr,tarea,tenv,physcon,erad,f,dist,idist,area,
			    ithermal,iinc,iit,cs,mcs,inocs,ntrit,nk,fenv,istep,
			    dtime,ttime,time,iviewfile,jobnamef,xloadold,
                            reltime,nmethod));
		
#ifdef SPOOLES
      spooles(ac,au,adb,aub,&sigma,bc,icol,irow,ntr,ntm,
	      &symmetryflag,&inputformat);
#else
      FORTRAN(dgesv,(ntr,&nhrs,ac,ntm,ipiv,bc,ntm,&info));
#endif
	
      if (info!=0){
	  printf("*ERROR IN RADFLOWLOAD: SINGULAR MATRIX*\n");}   
      
      else{ FORTRAN(radresult, (ntr,xloadact,ntm,bc,nloadtr,tarea,
				tenv,physcon,erad,f,fenv));}
  }

  free(v);

  return;

} 

