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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
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

void dynacont(double *co, int *nk, int *kon, int *ipkon, char *lakon, int *ne, 
	      int *nodeboun, int *ndirboun, double *xboun, int *nboun,
	      int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
	      int *nmpc, int *nodeforc,int *ndirforc,double *xforc, 
	      int *nforc,int *nelemload, char *sideload,double *xload,
	      int *nload, 
	      int *nactdof,int *neq, int *nzl,int *icol, int *irow, 
	      int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	      int *ilboun,double *elcon, int *nelcon, double *rhcon, 
	      int *nrhcon,double *cocon, int *ncocon,
	      double *alcon, int *nalcon, double *alzero, 
	      int *ielmat,int *ielorien, int *norien, double *orab, 
	      int *ntmat_,double *t0, 
	      double *t1,int *ithermal,double *prestr, int *iprestr, 
	      double *vold,int *iperturb, double *sti, int *nzs, 
	      double *tinc, double *tper, double *xmodalsteady,
	      double *veold, char *amname, double *amta,
	      int *namta, int *nam, int *iamforc, int *iamload,
	      int *iamt1,int *jout,char *filab,double *eme, double *xforcold, 
	      double *xloadold,
	      double *t1old, int *iamboun, double *xbounold, int *iexpl,
	      double *plicon, int *nplicon, double *plkcon,int *nplkcon,
	      double *xstate, int *npmat_, char *matname, int *mi,
	      int *ncmat_, int *nstate_, double *ener, char *jobnamec,
	      double *ttime, char *set, int *nset, int *istartset,
	      int *iendset, int *ialset, int *nprint, char *prlab,
	      char *prset, int *nener, double *trab, 
	      int *inotr, int *ntrans, double *fmpc, char *cbody, int *ibody,
	      double *xbody, int *nbody, double *xbodyold, int *istep,
	      int *isolver,int *jq, char *output, int *mcs, int *nkon,
	      int *mpcend, int *ics, double *cs, int *ntie, char *tieset,
	      int *idrct, int *jmax, double *tmin, double *tmax,
	      double *ctrl, int *itpamp, double *tietol,int *iit,
	      int *ncont,int *ne0, double *reltime, double *dtime,
	      double *bcontini, double *bj, double *aux, int *iaux,
	      double *bcont, int *nev, double *v,
              int *nkon0, double *deltmx, double *dtheta, double *theta,
              int *iprescribedboundary, int *mpcfree, int *memmpc_,
              int *itietri, int *koncont, double *cg, double *straight,
              int *iinc, double *vini,
              double *aa, double *bb, double *aanew, double *d, 
	      double *z, double *zeta,double *b, double *time0,double *time, 
	      int *ipobody,
              double *xforcact, double *xloadact, double *t1act, 
              double *xbounact, double *xbodyact, double *cd, double *cv,
              double *ampli, double *dthetaref, double *bjp, double *bp,
              double *cstr,int *imddof, int *nmddof, 
              int **ikactcontp, int *nactcont,int *nactcont_,
              double *aamech, double *bprev, int *iprev, int *inonlinmpc,
              int **ikactmechp, int *nactmech,int *imdnode,int *nmdnode,
              int *imdboun,int *nmdboun,int *imdmpc,int *nmdmpc,
              int *itp, int *inext,int *ifricdamp,double *aafric,
              double *bfric, int *imastop,int *nslavnode,int *islavnode,
              int *islavsurf,
              int *itiefac,double *areaslav,int *iponoels,int *inoels,
              double *springarea,int *izdof,int *nzdof,double *fn,
	      int *imastnode,int *nmastnode,double *xmastnor,
              double *xnormastface, double *xstateini, int *nslavs,
              int *cyclicsymmetry,double *xnoels,int *ielas){
    
  char lakonl[9]="        \0";

  int i,j,k,l,init,*itg=NULL,ntg=0,maxlenmpc,icascade=0,loop,
      konl[20],imat,nope,kodem,indexe,j1,jdof,kmin,kmax,
      id,newstep=0,idiscon,*ipiv=NULL,info,nrhs=1,kode,iener=0,
      *ikactcont=NULL,*ilactcont=NULL,*ikactcont1=NULL,nactcont1=0,
      i1,icutb=0,iconvergence=0,idivergence=0,mt=mi[1]+1,
      nactcont1_=100,*ikactmech=NULL,iabsload=0,nactfric_,nactfric,
      *ikactfric=NULL,im,nasym=0;

  long long i2;

  double *adb=NULL,*aub=NULL,*cgr=NULL, *au=NULL,fexp,fcos,fsin,fexm,
      physcon[1],zetaj,dj,ddj,h1,h2,h3,h4,h5,h6,sum,aai,bbi,tstart,tend,
      *ad=NULL,sigma=0.,alpham,betam,*bact=NULL,*bmin=NULL,*bv=NULL,
      xl[27],voldl[mt*9],elas[21],fnl[27],t0l,t1l,elconloc[21],veoldl[mt*9],
      bbmax,s[3600],*aaa=NULL,*bbb=NULL,func,funcp,*bjbasp=NULL,
      *bjbas=NULL, *bjinc=NULL, *dbj=NULL, *lhs=NULL,dbjmax,bjmax,
      *bjincp=NULL,sump,h14,*dbjp=NULL,senergy=0.0,*xforcdiff=NULL,
      df,i0,ic,ia,dbjmaxOLD1,dbjmaxOLD2,*xloaddiff=NULL,*dbcont=NULL,
      zl=0.0,*xbodydiff=NULL,*t1diff=NULL,*xboundiff=NULL,*bdiff=NULL;

  ikactcont=*ikactcontp;ikactmech=*ikactmechp;

  /* check for cyclic symmetry */

//  if((*mcs==0)||(cs[1]<0)){cyclicsymmetry=0;}else{cyclicsymmetry=1;}

  if(*inonlinmpc==1) iabsload=2;

  if(ithermal[0]<=1){
      kmin=1;kmax=3;
  }else if(ithermal[0]==2){
      kmin=0;kmax=mi[1];if(kmax>2)kmax=2;
  }else{
      kmin=0;kmax=3;
  }

  xforcdiff=NNEW(double,*nforc);
  xloaddiff=NNEW(double,2**nload);
  xbodydiff=NNEW(double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodydiff[k]=xbody[k];}
  xboundiff=NNEW(double,*nboun);
  if(*ithermal==1) t1diff=NNEW(double,*nk);

  /* load the convergence constants from ctrl*/

  i0=ctrl[0];ic=ctrl[3];ia=ctrl[7];df=ctrl[10];

  /* set the convergence parameters*/

  dbjmaxOLD1=0.0;
  dbjmaxOLD2=0.0;

//  printf("\nstart dynacont\n");

  /* calculating the contact forces */
	
//  memset(&bcont[0],0,sizeof(double)*neq[1]);
  for(j=0;j<*nactcont;j++){bcont[ikactcont[j]]=0.;}
  
  *ne=*ne0;*nkon=*nkon0;
  
  contact(ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
	  straight,nkon,co,vold,ielmat,cs,elcon,istep,
	  iinc,iit,ncmat_,ntmat_,ne0,
	  vini,nmethod,nmpc,mpcfree,memmpc_,
	  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
          iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
          itiefac,areaslav,iponoels,inoels,springarea,tietol,reltime,
	  imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
          xnoels);

  ikactcont1=NNEW(int,nactcont1_);

  for(i=*ne0;i<*ne;i++){
      indexe=ipkon[i];
      imat=ielmat[mi[2]*i];
      kodem=nelcon[2*imat-2];
      for(j=0;j<8;j++){lakonl[j]=lakon[8*i+j];}
      nope=atoi(&lakonl[7])+1;
      for(j=0;j<nope;j++){
	  konl[j]=kon[indexe+j];
	  for(j1=0;j1<3;j1++){
	      xl[j*3+j1]=co[3*(konl[j]-1)+j1];
	      voldl[mt*j+j1+1]=vold[mt*(konl[j]-1)+j1+1];
	      veoldl[mt*j+j1+1]=veold[mt*(konl[j]-1)+j1+1];
	  }
      }
      konl[nope]=kon[indexe+nope];

      FORTRAN(springforc,(xl,konl,voldl,&imat,elcon,nelcon,elas,
			  fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,
			  &t1l,&kodem,elconloc,plicon,nplicon,npmat_,
			  veoldl,&senergy,&iener,cstr,mi,
                          &springarea[2*(konl[nope]-1)],nmethod,ne0,iperturb,
                          nstate_,xstateini,xstate,reltime,
                          &xnormastface[27*(konl[nope]-1)],ielas));

      storecontactdof(&nope,nactdof,&mt,konl,&ikactcont1,&nactcont1,
		      &nactcont1_,bcont,fnl,ikmpc,nmpc,ilmpc,ipompc,nodempc, 
		      coefmpc);

  }
  RENEW(ikactcont1,int,nactcont1);

  /* merging ikactcont with ikactcont1; the result ist
     stored in ikactcont */

  for(i=0;i<nactcont1;i++){
      jdof=ikactcont1[i];
      FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
      do{
	  if(id>0){
	      if(ikactcont[id-1]==jdof){
		  break;
	      }
	  }
	  (*nactcont)++;
	  if(*nactcont>*nactcont_){
	      *nactcont_=(int)(1.1**nactcont_);
	      RENEW(ikactcont,int,*nactcont_);
	  }
	  k=*nactcont-1;
	  l=k-1;
	  while(k>id){
	      ikactcont[k--]=ikactcont[l--];
	  }
	  ikactcont[id]=jdof;
	  break;
      }while(1);
  }
//  free(ikactcont1);
    
    /* calculate the change in contact force */
    
  bbmax=0.;
  if(icutb==0){
      for(i=0;i<*nactcont;i++){
	  jdof=ikactcont[i];
	  if(fabs(bcont[jdof]-bcontini[jdof])>bbmax){
	      bbmax=fabs(bcont[jdof]-bcontini[jdof]);
	  }
      }
  }

  /* removing entries in bcont */

  for(j=0;j<nactcont1;j++){bcont[ikactcont1[j]]=0.;}
  free(ikactcont1);
  *nactcont=0;
  
  /* major loop to calculate the correction of bj due to contact */

  ilactcont=NNEW(int,*nactcont_);
  dbcont=NNEW(double,*nactcont_**nev);

  icutb=0;

  do{

      /* restoring initial values */

      if(*nmdnode>0){
	  for(i=0;i<*nmdnode;i++){
	      i1=mt*(imdnode[i]-1);
	      for(j=kmin;j<=kmax;j++){
		  vold[i1+j]=vini[i1+j];
	      }
	  }
      }else{
	  memcpy(&vold[0],&vini[0],sizeof(double)*mt**nk);
      }

      if(*nstate_!=0){
	for(k=0;k<*nstate_*mi[0]*(*ne0+*nslavs);++k){
	  xstate[k]=xstateini[k];
	}	
      }
    
    /* restoring aa[(iinc-1)*nev+i] (before change of *dtime) */

    for(i=0;i<*nev;i++){
      aa[i]+=bb[i]*(*time-*dtime);
    }
    
    /* increment size is reduced if:
            - the contact force change is too large (only in first iteration)
            - or the increment did not converge  */
    
    if((bbmax>*deltmx || icutb>0)&&(((*itp==1)&&(*dtheta>*tmin))||(*itp==0))){
      
      /* force increase too big: increment size is decreased */

      if(icutb>0){
	*dtheta=*dtheta*df;
//	printf("*INFORMATION: increment size is decreased to %e\nthe increment is reattempted\n\n",*dtheta**tper);
      }
      else{
	*dtheta=*dtheta**deltmx/bbmax;
      }
//      printf("correction of dtime due to contact: %e\n",*dtheta**tper);
      *dthetaref=*dtheta;
      if(*itp==1){
	  (*inext)--;
	  *itp=0;
      }
      
      /* check whether the new increment size is not too small */
      
      if(*dtheta<*tmin){
//	  printf("\n *WARNING: increment size %e smaller than minimum %e\n",*dtheta**tper,*tmin**tper);
//	  printf("             minimum is taken\n");
	  *dtheta=*tmin;
	  *dthetaref=*dtheta;
      }
      
      *reltime=*theta+(*dtheta);
      *time=*reltime**tper;
      *dtime=*dtheta**tper;
      
      /* calculating the instantaneous loads (forces, surface loading, 
	 centrifugal and gravity loading or temperature) */
      
      FORTRAN(temploaddiff,(xforcold,xforc,xforcact,iamforc,nforc,
		xloadold,xload,xloadact,iamload,nload,ibody,xbody,
		nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
		namta,nam,ampli,time,reltime,ttime,dtime,ithermal,
	        nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
	        ndirboun,nodeforc,
	        ndirforc,istep,iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
		nelemload,sideload,mi,
		xforcdiff,xloaddiff,xbodydiff,t1diff,xboundiff,&iabsload,
		iprescribedboundary,ntrans,trab,inotr,veold,nactdof,bcont,
                fn));
	      
      /* calculating the instantaneous loading vector */
      
      if(iabsload!=2){
	  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		   ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcdiff,
		   nforc,nelemload,sideload,xloaddiff,nload,xbodydiff,
		   ipobody,nbody,cgr,b,nactdof,&neq[1],nmethod,
		   ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		   alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		   t0,t1diff,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		   nplicon,plkcon,nplkcon,
		   npmat_,ttime,time,istep,iinc,dtime,physcon,ibody,
		   xbodyold,reltime,veold,matname,mi,ikactmech,nactmech));
      }else{
	  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		   ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		   nforc,nelemload,sideload,xloadact,nload,xbodyact,
		   ipobody,nbody,cgr,b,nactdof,&neq[1],nmethod,
		   ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		   alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		   t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		   nplicon,plkcon,nplkcon,
		   npmat_,ttime,time,istep,iinc,dtime,physcon,ibody,
		   xbodyold,reltime,veold,matname,mi,ikactmech,nactmech));
      }
      
      /* correction for nonzero SPC's */
      
      if(*iprescribedboundary){
	  dynboun(amta,namta,nam,ampli,time,ttime,dtime,
		  xbounold,xboun,
		  xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		  aub,icol,irow,neq,nzs,&sigma,b,isolver,
		  &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv,
  	          bprev,bdiff,nactmech,&iabsload,iprev);
      }

      /* correcting aamech */

      if(!(*cyclicsymmetry)){
	  for(i=0;i<*nev;i++){
	      i2=(long long)i*neq[1];
	      
	      if(iabsload==2){aamech[i]=0.;}
	      if(*nactmech<neq[1]/2){
		  for(j=0;j<*nactmech;j++){
		      aamech[i]+=z[i2+ikactmech[j]]*b[ikactmech[j]];
		  }
	      }else{
		  for(j=0;j<neq[1];j++){
		      aamech[i]+=z[i2+j]*b[j];
		  }
	      }
	  }
      }else{
	  for(i=0;i<*nev;i++){
	      if(iabsload==2){aamech[i]=0.;}
	  }
	  for(j=0;j<*nactmech;j++){
	      FORTRAN(nident,(izdof,&ikactmech[j],nzdof,&id));
	      if(id!=0){
		  if(izdof[id-1]==ikactmech[j]){
		      for(i=0;i<*nev;i++){
			  aamech[i]+=z[i**nzdof+id-1]*b[ikactmech[j]];
		      }
		  }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	      }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	  }
      }
      
    }

    bbmax=0.;
    
    /* calculating the linearized force function connecting the
       mechanical+friction+contact load at the start of the increment with
       the mechanical+friction load at the end of the increment
       = base load */

    for(i=0;i<*nev;i++){

      aanew[i]=aamech[i];
      if(*ifricdamp==1){aanew[i]+=aafric[i];}

      bb[i]=(aanew[i]-aa[i])/(*dtime);
      aa[i]=aanew[i]-bb[i]**time;
    }
    
    /* calculating the base response */
    bjbas=NNEW(double,*nev); /* basis response modal decomposition */
    bjbasp=NNEW(double,*nev);
    for(l=0;l<*nev;l++){
      zetaj=zeta[l];
      dj=d[l];
      
      /* zero eigenfrequency: rigid body mode */
      
      if(fabs(d[l])<=1.e-10){
	  aai=aa[l];
	  bbi=bb[l];
	  tstart=*time0;
	  tend=*time;
	  sum=tend*(aai**time+
		    tend*((bbi**time-aai)/2.-bbi*tend/3.))-
	      tstart*(aai**time+
		      tstart*((bbi**time-aai)/2.-bbi*tstart/3.));
	  sump=tend*(aai+bbi*tend/2.)-tstart*(aai+bbi*tstart/2.);
	  bjbas[l]=sum+cd[l]+*dtime*cv[l];
	  bjbasp[l]=sump+cv[l];
      }
      
      /*   subcritical damping */
      
      else if(zetaj<1.-1.e-6){
	ddj=dj*sqrt(1.-zetaj*zetaj);
	h1=zetaj*dj;
	h2=h1*h1+ddj*ddj;
	h3=h1*h1-ddj*ddj;
	h4=2.*h1*ddj/h2;
	h14=zetaj*dj/ddj;
	tstart=0;
	FORTRAN(fsub,(time,dtime,&aa[l],&bb[l],&ddj,
		      &h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsub,(time,&tstart,&aa[l],&bb[l],&ddj,
		      &h1,&h2,&h3,&h4,&func,&funcp));
	sum-=func;sump-=funcp;
	fexp=exp(-h1**dtime);
	fsin=sin(ddj**dtime);
	fcos=cos(ddj**dtime);
	
	bjbas[l]=sum/ddj+fexp*(fcos+zetaj/sqrt(1.-zetaj*zetaj)*fsin)*cd[l]+
	  fexp*fsin*cv[l]/ddj;
	bjbasp[l]=sump/ddj+fexp*((-h1+ddj*h14)*fcos+(-ddj-h1*h14)*fsin)*cd[l]
	  +fexp*(-h1*fsin+ddj*fcos)*cv[l]/ddj;
	
      }
      
      /*      supercritical damping */
      
      else if(zetaj>1.+1.e-6){
	ddj=dj*sqrt(zetaj*zetaj-1.);
	h1=ddj-zetaj*dj;
	h2=ddj+zetaj*dj;
	h3=1./h1;
	h4=1./h2;
	h5=h3*h3;
	h6=h4*h4;
	tstart=0;
	FORTRAN(fsuper,(time,dtime,&aa[l],&bb[l],
			&h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsuper,(time,&tstart,&aa[l],&bb[l],
			&h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	sum-=func;sump-=funcp;
	fexm=exp(h1**dtime);
	fexp=exp(-h2**dtime);
	h14=zetaj*dj/ddj;
	
	bjbas[l]=sum/(2.*ddj)+(fexm+fexp)*cd[l]/2.+zetaj*(fexm-fexp)/(2.*sqrt(zetaj*zetaj-1.))*cd[l]+(fexm-fexp)*cv[l]/(2.*ddj);
	bjbasp[l]=sump/(2.*ddj)+(h1*fexm-h2*fexp)*cd[l]/2.+(h14*cd[l]+cv[l]/ddj)*(h1*fexm+h2*fexp)/2.;
      }
      
      /* critical damping */
      
      else{
	h1=zetaj*dj;
	h2=1./h1;
	h3=h2*h2;
	h4=h2*h3;
	tstart=0;
	FORTRAN(fcrit,(time,dtime,&aa[l],&bb[l],&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fcrit,(time,&tstart,&aa[l],&bb[l],&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum-=func;sump+=funcp;
	fexp=exp(-h1**dtime);
	bjbas[l]=sum+fexp*((1.+h1**dtime)*cd[l]+*dtime*cv[l]);
	bjbasp[l]=sump+fexp*(-h1*h1**dtime*cd[l]+(1.-h1**dtime)*cv[l]);
      }
    }
    
    /* calculating the incremental response due to contact */
    
    aai=-(*time-*dtime)/(*dtime);
    bbi=1./(*dtime);
    
    bjinc=NNEW(double,*nev); /* incremental response modal decomposition */
    bjincp=NNEW(double,*nev);
    for(l=0;l<*nev;l++){
      zetaj=zeta[l];
      dj=d[l];
      
      /* zero eigenfrequency: rigid body mode */
      
      if(fabs(d[l])<=1.e-10){
	tstart=*time0;
	tend=*time;
	sum=tend*(aai**time+
		  tend*((bbi**time-aai)/2.-bbi*tend/3.))-
	  tstart*(aai**time+
		  tstart*((bbi**time-aai)/2.-bbi*tstart/3.));
	sump=tend*(aai+bbi*tend/2.)-tstart*(aai+bbi*tstart/2.);
	
	bjinc[l]=sum;
	bjincp[l]=sump;
      }
      
      /*   subcritical damping */
      
      else if(zetaj<1.-1.e-6){
	ddj=dj*sqrt(1.-zetaj*zetaj);
	h1=zetaj*dj;
	h2=h1*h1+ddj*ddj;
	h3=h1*h1-ddj*ddj;
	h4=2.*h1*ddj/h2;
	tstart=0.;
	FORTRAN(fsub,(time,dtime,&aai,&bbi,&ddj,
		      &h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsub,(time,&tstart,&aai,&bbi,&ddj,
		      &h1,&h2,&h3,&h4,&func,&funcp));
	sum-=func;sump-=funcp;
	
	bjinc[l]=sum/ddj;
	bjincp[l]=sump/ddj;
	
      }
      
      /*      supercritical damping */
      
      else if(zetaj>1.+1.e-6){
	ddj=dj*sqrt(zetaj*zetaj-1.);
	h1=ddj-zetaj*dj;
	h2=ddj+zetaj*dj;
	h3=1./h1;
	h4=1./h2;
	h5=h3*h3;
	h6=h4*h4;
	tstart=0.;
	FORTRAN(fsuper,(time,dtime,&aai,&bbi,
			&h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsuper,(time,&tstart,&aai,&bbi,
			&h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	sum-=func;sump-=funcp;
	
	bjinc[l]=sum/(2.*ddj);
	bjincp[l]=sump/(2.*ddj);
	
      }
      
      /* critical damping */
      
      else{
	h1=zetaj*dj;
	h2=1./h1;
	h3=h2*h2;
	h4=h2*h3;
	tstart=0.;
	FORTRAN(fcrit,(time,dtime,&aai,&bbi,&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fcrit,(time,&tstart,&aai,&bbi,&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum-=func;sump-=funcp;
	
	bjinc[l]=sum;
	bjincp[l]=sump;
	
      }
    }

    aaa=NNEW(double,*nev);
    bbb=NNEW(double,*nev**nev);
    lhs=NNEW(double,*nev**nev);
    ipiv=NNEW(int,*nev);
    dbj=NNEW(double,*nev); /* change in bj */
    dbjp=NNEW(double,*nev); /* change in djp */
    
    /* starting solution for the iteration loop = base solution */

    memcpy(&bj[0],&bjbas[0],sizeof(double)**nev);
    memcpy(&bjp[0],&bjbasp[0],sizeof(double)**nev);
    
    /* major iteration loop for the contact response */
    
    loop=0;
//    printf("Contact-Iteration\n");
    do{
      loop++;
//      printf("loop=%d\n",loop);
      
      /* composing the response */
      
      if(*iprescribedboundary){
	  if(*nmdnode==0){
	      memcpy(&b[0],&bmin[0],sizeof(double)*neq[1]);
	      memcpy(&bp[0],&bv[0],sizeof(double)*neq[1]);
	  }else{
	      for(i=0;i<*nmddof;i++){
		  b[imddof[i]]=bmin[imddof[i]];
		  bp[imddof[i]]=bv[imddof[i]];
	      }
	  }
      }
      else{
	  if(*nmdnode==0){
	      DMEMSET(b,0,neq[1],0.);
	      DMEMSET(bp,0,neq[1],0.);
	  }else{
	      for(i=0;i<*nmddof;i++){
		  b[imddof[i]]=0.;
		  bp[imddof[i]]=0.;
	      }
	  }
      }
      
      if(!(*cyclicsymmetry)){
	  if(*nmdnode==0){
	      for(i=0;i<neq[1];i++){
		  for(j=0;j<*nev;j++){
		      b[i]+=bj[j]*z[(long long)j*neq[1]+i];
		      bp[i]+=bjp[j]*z[(long long)j*neq[1]+i];
		  }
	      }
	  }else{
	      for(i=0;i<*nmddof;i++){
		  for(j=0;j<*nev;j++){
		      b[imddof[i]]+=bj[j]*z[(long long)j*neq[1]+imddof[i]];
		      bp[imddof[i]]+=bjp[j]*z[(long long)j*neq[1]+imddof[i]];
		  }
	      }
	  }
      }else{
	  for(i=0;i<*nmddof;i++){
	      FORTRAN(nident,(izdof,&imddof[i],nzdof,&id));
	      if(id!=0){
		  if(izdof[id-1]==imddof[i]){
		      for(j=0;j<*nev;j++){
			  b[imddof[i]]+=bj[j]*z[j**nzdof+id-1];
			  bp[imddof[i]]+=bjp[j]*z[j**nzdof+id-1];
		      }
		  }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	      }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	  }
      }
      
      /* update nonlinear MPC-coefficients (e.g. for rigid
	 body MPC's */

      if(*inonlinmpc==1){
	  FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,reltime,&newstep,xboun,fmpc,
			 iit,&idiscon,ncont,trab,ntrans,ithermal,mi));
      }
      
      /* calculating displacements/temperatures */
      
      FORTRAN(dynresults,(nk,v,ithermal,nactdof,vold,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
	      b,bp,veold,dtime,mi,imdnode,nmdnode,imdboun,nmdboun,
	      imdmpc,nmdmpc,nmethod,time));
      
      if(iconvergence==1){
	break;
      }
      
      /* creating contact elements and calculating the contact forces
	 based on the displacements at the end of the present increment */
      
	  for(j=0;j<*nactcont;j++){bcont[ikactcont[j]]=0.;}

      RENEW(dbcont,double,*nactcont_**nev);
      RENEW(ikactcont,int,*nactcont_);
      RENEW(ilactcont,int,*nactcont_);
      *nactcont=0;
      
      DMEMSET(dbcont,0,*nactcont_**nev,0.);
      DMEMSET(ikactcont,0,*nactcont_,0.);
      
      *ne=*ne0;*nkon=*nkon0;
      contact(ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
	      straight,nkon,co,vold,ielmat,cs,elcon,istep,
	      iinc,iit,ncmat_,ntmat_,ne0,
	      vini,nmethod,nmpc,mpcfree,memmpc_,
	      &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	      iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
              itiefac,areaslav,iponoels,inoels,springarea,tietol,reltime,
	      imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
              xnoels);

//      printf("number of contact springs = %d\n",*ne-*ne0);

      for(i=*ne0;i<*ne;i++){
	indexe=ipkon[i];
	imat=ielmat[mi[2]*i];
	kodem=nelcon[2*imat-2];
	for(j=0;j<8;j++){lakonl[j]=lakon[8*i+j];}
	nope=atoi(&lakonl[7])+1;
	for(j=0;j<nope;j++){
	  konl[j]=kon[indexe+j];
	  for(j1=0;j1<3;j1++){
	    xl[j*3+j1]=co[3*(konl[j]-1)+j1];
	    voldl[mt*j+j1+1]=vold[mt*(konl[j]-1)+j1+1];
	    veoldl[mt*j+j1+1]=veold[mt*(konl[j]-1)+j1+1];
	  }
	}
	konl[nope]=kon[indexe+nope];

	FORTRAN(springforc,(xl,konl,voldl,&imat,elcon,nelcon,elas,
			    fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,
			    &t1l,&kodem,elconloc,plicon,nplicon,npmat_,
			    veoldl,&senergy,&iener,cstr,mi,
                            &springarea[2*(konl[nope]-1)],nmethod,ne0,iperturb,
                            nstate_,xstateini,xstate,reltime,
                            &xnormastface[27*(konl[nope]-1)],ielas));
	
	FORTRAN(springstiff,(xl,elas,konl,voldl,s,&imat,elcon,nelcon,
	        ncmat_,ntmat_,&nope,lakonl,&t0l,&t1l,&kode,elconloc,
		plicon,nplicon,npmat_,iperturb,&springarea[2*(konl[nope]-1)],
		nmethod,mi,ne0,nstate_,xstateini,xstate,reltime,
			     &xnormastface[27*(konl[nope]-1)],&nasym));

	dfdbj(bcont,&dbcont,&neq[1],&nope,konl,nactdof,
	      s,z,ikmpc,ilmpc,ipompc,nodempc,nmpc,coefmpc,
	      fnl,nev,&ikactcont,&ilactcont,nactcont,nactcont_,mi,
              cyclicsymmetry,izdof,nzdof);
      }

      if(*nactcont>100){*nactcont_=*nactcont;}else{*nactcont_=100;}
      RENEW(ikactcont,int,*nactcont_);
      RENEW(ilactcont,int,*nactcont_);
      RENEW(dbcont,double,*nactcont_**nev);

      /* aaa(i) is the internal product of the contact force at the end of the
	 increment with eigenmode i
	 bbb(i,j) is the internal product of the change of the contact force with
	 respect to modal coordinate j with the eigenmode i */

      DMEMSET(bbb,0,*nev**nev,0.);
      DMEMSET(aaa,0,*nev,0.);

      if(!(*cyclicsymmetry)){
	  for(k=0; k<*nactcont; k++){
	      i1=ikactcont[k];
	      i2=(ilactcont[k]-1)**nev;
	      for(j=0; j<*nev; j++){
		  zl=z[(long long)j*neq[1]+i1];
		  aaa[j]+=zl*bcont[i1];
		  for(l=0; l<*nev; l++){
		      bbb[l**nev+j]+=zl*dbcont[i2+l];
		  }
	      }
	  }
      }else{
	  for(k=0; k<*nactcont; k++){
	      i1=ikactcont[k];
	      i2=(ilactcont[k]-1)**nev;
	      FORTRAN(nident,(izdof,&i1,nzdof,&id));
	      if(id!=0){
		  if(izdof[id-1]==i1){
		      for(j=0; j<*nev; j++){
			  zl=z[j**nzdof+id-1];
			  aaa[j]+=zl*bcont[i1];
			  for(l=0; l<*nev; l++){
			      bbb[l**nev+j]+=zl*dbcont[i2+l];
			  }
		      }
		  }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	      }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	  }
      }
      
      for(l=0;l<*nev;l++){
	i1=l**nev;
	for(j=0;j<*nev;j++){
	  if(j==l){lhs[i1+j]=1.;}else{lhs[i1+j]=0.;}
	  lhs[i1+j]-=bjinc[j]*bbb[i1+j];
	}
	dbj[l]=bjbas[l]+bjinc[l]*aaa[l]-bj[l];
      }

      /* solve the system of equations; determine dbj */
      
      FORTRAN(dgesv,(nev,&nrhs,lhs,nev,ipiv,dbj,nev,&info));
      
      /* check the size of dbj */
      
      bjmax=0.;
      dbjmaxOLD2=dbjmaxOLD1;
      dbjmaxOLD1=dbjmax;
      dbjmax=0.;
      for(i=0;i<*nev;i++){
	if(fabs(bj[i])>bjmax) bjmax=fabs(bj[i]);
	if(fabs(dbj[i])>dbjmax) dbjmax=fabs(dbj[i]);
      }

      iconvergence=0;
      idivergence=0;

      if(dbjmax<=0.005*bjmax){
	
	//calculate bjp: the derivative of bj w.r.t. time
	
	for(j=0; j<*nev; j++){
	  bjp[j]=bjbasp[j]+bjincp[j]*aaa[j];
	}
	FORTRAN(dgetrs,("No transpose",nev,&nrhs,lhs,nev,ipiv,bjp,nev,&info));
	iconvergence=1;	   
      }
      else{
	  if(loop>=i0 && loop<=ic){
	  /* check for divergence */
	  if((dbjmax>dbjmaxOLD1) && (dbjmax>dbjmaxOLD2)){
	    /* divergence --> cutback */ 
//	    printf("*INFORMATION: divergence --> cutback\n");
	    idivergence=1;
	    icutb++;
	    break;
	  }
	}
	else{
	  if(loop>ic){
	    /* cutback after ic iterations*/
//	    printf("*INFORMATION: too many iterations --> cutback\n");
	    idivergence=1;
	    icutb++;
	    break;
	  }
	}
      }

      /* add dbj to db */
      
      for(j=0;j<*nev;j++){
	bj[j]+=dbj[j];
      }
      
    }while(1);
  }while(idivergence==1 && icutb<10);

//  printf("Contact-Iteration Done\n");

  if(icutb>=10){
    //no convergence, stop all
    printf("*ERROR: Contact did not converge.\n");
    FORTRAN(stop,());
  }
  
  /* convergence has been reached */

  /* calculating the damping/friction contribution */

  if(*ifricdamp==1){
      nactfric_=*nactcont_;
      nactfric=0;
      ikactfric=NNEW(int,nactfric_);

//      memset(&ikactfric[0],0,sizeof(int)*nactfric_);
      DMEMSET(ikactfric,0,nactfric_,0.);
      
      *ne=*ne0;*nkon=*nkon0;
      contact(ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
	      straight,nkon,co,vold,ielmat,cs,elcon,istep,
	      iinc,iit,ncmat_,ntmat_,ne0,
	      vini,nmethod,nmpc,mpcfree,memmpc_,
	      &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	      iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
              itiefac,areaslav,iponoels,inoels,springarea,tietol,reltime,
	      imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
              xnoels);

//      printf("number of contact springs = %d\n",*ne-*ne0);

      for(i=*ne0;i<*ne;i++){
	indexe=ipkon[i];
	imat=ielmat[mi[2]*i];
	kodem=nelcon[2*imat-2];
	for(j=0;j<8;j++){lakonl[j]=lakon[8*i+j];}
	nope=atoi(&lakonl[7])+1;
	for(j=0;j<nope;j++){
	  konl[j]=kon[indexe+j];
	  for(j1=0;j1<3;j1++){
	    xl[j*3+j1]=co[3*(konl[j]-1)+j1];
	    voldl[mt*j+j1+1]=vold[mt*(konl[j]-1)+j1+1];
	    veoldl[mt*j+j1+1]=veold[mt*(konl[j]-1)+j1+1];
	  }
	}
	konl[nope]=kon[indexe+nope];

	FORTRAN(fridaforc,(xl,konl,voldl,&imat,elcon,nelcon,elas,
			  fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,
			  &t1l,&kodem,elconloc,plicon,nplicon,npmat_,
			  veoldl,&senergy,&iener,cstr,mi,
                          &springarea[konl[nope]-1]));
	
	storecontactdof(&nope,nactdof,&mt,konl,&ikactfric,&nactfric,
			  &nactfric_,bfric,fnl,ikmpc,nmpc,ilmpc,ipompc,nodempc, 
			  coefmpc);
      }

      /* calculating aafric: contains the force contributions due
         to damping and friction in the contact areas */

      if(!(*cyclicsymmetry)){
	  for(i=0;i<*nev;i++){
	      i2=(long long)i*neq[1];
	      aafric[i]=0.;
	      for(j=0;j<nactfric;j++){
		  aafric[i]+=z[i2+ikactfric[j]]*bfric[ikactfric[j]];
	      }
	  }
      }else{
	  for(i=0;i<*nev;i++){aafric[i]=0.;}
	  for(j=0;j<nactfric;j++){
	      FORTRAN(nident,(izdof,&ikactfric[j],nzdof,&id));
	      if(id!=0){
		  if(izdof[id-1]==ikactfric[j]){
		      for(i=0;i<*nev;i++){
			  aafric[i]+=z[i**nzdof+id-1]*bfric[ikactfric[j]];
		      }
		  }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	      }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	  }
      }
      
      /* setting bfric to zero */

      for(j=0;j<nactfric;j++){bfric[ikactfric[j]]=0.;}
      free(ikactfric);

  }

  /* restoring aa[(*iinc-1)*nev+i] */

  for(i=0;i<*nev;i++){
    aa[i]+=bb[i]*(*time-*dtime);
  }
  
  /* calculating the linearized force function connecting the
     mechanical+contact load at the start of the increment with
     the mechanical+contact load at the end of the increment */
  
  if(!(*cyclicsymmetry)){
      for(i=0;i<*nev;i++){
	  i2=(long long)i*neq[1];
	  
	  aanew[i]=aamech[i];
	  if(*ifricdamp==1){aanew[i]+=aafric[i];}
	  for(j=0;j<*nactcont;j++){
	      aanew[i]+=z[i2+ikactcont[j]]*bcont[ikactcont[j]];
	  }
	  
	  bb[i]=(aanew[i]-aa[i])/(*dtime);
	  aa[i]=aanew[i]-bb[i]**time;
      }
  }else{
      memcpy(&aanew[0],&aamech[0],sizeof(double)**nev);
      if(*ifricdamp==1){
	  for(i=0;i<*nev;i++){aanew[i]+=aafric[i];}
      }
      for(j=0;j<*nactcont;j++){
	  FORTRAN(nident,(izdof,&ikactcont[j],nzdof,&id));
	  if(id!=0){
	      if(izdof[id-1]==ikactcont[j]){
		  for(i=0;i<*nev;i++){
		      aanew[i]+=z[i**nzdof+id-1]*bcont[ikactcont[j]];
		  }
	      }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
	  }else{printf("*ERROR in dynacont\n");FORTRAN(stop,());}
      }
      for(i=0;i<*nev;i++){
      	  bb[i]=(aanew[i]-aa[i])/(*dtime);
	  aa[i]=aanew[i]-bb[i]**time;
      }
  }
  
  free(aaa);free(bbb);free(bjbas);free(bjinc);free(dbj);free(lhs);
  free(ipiv);free(bjbasp);free(bjincp);free(dbjp);free(ilactcont);
  free(dbcont);
  free(xforcdiff);free(xloaddiff);free(xboundiff),free(xbodydiff);

  if(*ithermal==1) free(t1diff);
  
  *ikactcontp=ikactcont;*ikactmechp=ikactmech;
  
  return;
}


