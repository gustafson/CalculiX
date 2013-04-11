/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

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
#include <omp.h>
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
	      double *tinc, double *tper, double *xmodal,
	      double *veold, char *amname, double *amta,
	      int *namta, int *nam, int *iamforc, int *iamload,
	      int *iamt1,int *jout,char *filab,double *eme, double *xforcold, 
	      double *xloadold,
	      double *t1old, int *iamboun, double *xbounold, int *iexpl,
	      double *plicon, int *nplicon, double *plkcon,int *nplkcon,
	      double *xstate, int *npmat_, char *matname, int *mint_,
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
	      double *bcont, double **dbcontp, int *nev, double *v,
              int *nkon0, double *deltmx, double *dtheta, double *theta,
              int *iprescribedboundary, int *mpcfree, int *memmpc_,
              int *itietri, int *koncont, double *cg, double *straight,
              int *iinc, int *ifcont1, int *ifcont2, double *vini,
              double *aa, double *bb, double *aanew, double *d, 
	      double *z, double *zeta,
              double *b, double *bmech, double *time0,double *time1, 
	      int *ipobody,
              double *xforcact, double *xloadact, double *t1act, 
              double *xbounact, double *xbodyact, double *cd, double *cv,
              double *ampli, double *dthetaref, double *bjp, double *bp,
              double *cstr,double **dbcontinip){
    
  char lakonl[9]="        \0";

  int i,j,k,l,init,*itg=NULL,ntg=0,maxlenmpc,icascade=0,loop,
      konl[20],imat,nope,kodem,indexe,j1,ist,jdof,index,
      id,node,ndir,itp=0,newstep=0,idiscon,*ipiv=NULL,info,nrhs=1,kode,iener=0,
      nactcont=0,*ikactcont=NULL,*ilactcont,nactcont_=0,
      i1,i2,i3,i4,i5,icutb=0,iconvergence=0,idivergence=0;
  time_t start,end;
  struct timeval starttime,endtime;
  double *adb=NULL,*aub=NULL,*cgr=NULL, *au=NULL,fexp,fcos,fsin,fexm,
      physcon[1],zetaj,dj,ddj,h1,h2,h3,h4,h5,h6,sum,aai,bbi,tstart,tend,
      *ad=NULL,sigma=0.,alpham,betam,*bact=NULL,*bmin=NULL,*bv=NULL,
      xl[27],voldl[36],elas[21],fnl[27],t0l,t1l,elconloc[21],veoldl[27],
      bbmax,s[3600],*aaa=NULL,*bbb=NULL,func,funcp,*bjbasp=NULL,
      *bjbas=NULL, *bjinc=NULL, *dbj=NULL, *lhs=NULL,dbjmax,bjmax,
      *bjincp=NULL,sump,h14,*dbjp=NULL,senergy=0.0,
      df,i0,ic,ia,dbjmaxOLD1,dbjmaxOLD2,dbdbcontmax=0,
      *dbcont=*dbcontp,*dbcontini=*dbcontinip,zl=0.0;

  /* load the convergence constants from ctrl*/

  i0=ctrl[0];ic=ctrl[3];ia=ctrl[7];df=ctrl[10];

  nactcont_=100;
  ikactcont=NNEW(int,nactcont_);
  ilactcont=NNEW(int,nactcont_);
  dbcont=NNEW(double,nactcont_**nev);

  /* set the convergence parameters*/

  dbjmaxOLD1=0.0;
  dbjmaxOLD2=0.0;

  printf("\nstart dynacont\n");

  /* calculating the contact forces */
	
  memset(&bcont[0],0,sizeof(double)*neq[1]);
  //  for(i=0;i<neq[1];i++){bcont[i]=0.;}
  
  *ne=*ne0;*nkon=*nkon0;
  
  contact(ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
	  straight,nkon,co,vold,ielmat,cs,elcon,istep,
	  iinc,iit,ncmat_,ntmat_,ifcont1,ifcont2,ne0,
	  vini,nmethod,nmpc,mpcfree,memmpc_,
	  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
          iperturb,ikboun,nboun);
  for(i=*ne0;i<*ne;i++){
      indexe=ipkon[i];
      imat=ielmat[i];
      kodem=nelcon[2*imat-2];
      for(j=0;j<8;j++){lakonl[j]=lakon[8*i+j];}
      nope=atoi(&lakonl[7]);
      for(j=0;j<nope;j++){
	  konl[j]=kon[indexe+j];
	  for(j1=0;j1<3;j1++){
	      xl[j*3+j1]=co[3*(konl[j]-1)+j1];
	      voldl[j*4+j1+1]=vold[5*(konl[j]-1)+j1+1];
	      veoldl[j*3+j1]=veold[4*(konl[j]-1)+j1+1];
	  }
      }
      konl[nope]=kon[indexe+nope];
      FORTRAN(springforc,(xl,konl,voldl,&imat,elcon,nelcon,elas,
			  fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,
			  &t1l,&kodem,elconloc,plicon,nplicon,npmat_,
			  veoldl,&senergy,&iener,cstr));

      for(j=0;j<nope;j++){
	  for(j1=0;j1<3;j1++){
	      jdof=nactdof[4*(konl[j]-1)+j1+1];
	      if(jdof!=0){
		  bcont[jdof-1]-=fnl[3*j+j1];
	      }else{
		  jdof=8*(konl[j]-1)+j1+1;
		  FORTRAN(nident,(ikmpc,&jdof,nmpc,&id));
		  if(id>0){
		      if(ikmpc[id-1]==jdof){
			  id=ilmpc[id-1];
			  ist=ipompc[id-1];
			  index=nodempc[3*ist-1];
			  if(index==0) continue;
			  do{
			      node=nodempc[3*index-3];
			      ndir=nodempc[3*index-2];
			      jdof=nactdof[4*(node-1)+ndir];
			      if(jdof!=0){
				  bcont[jdof-1]+=coefmpc[index-1]*
				      fnl[3*j+j1]/coefmpc[ist-1];
			      }
			      index=nodempc[3*index-1];
			      if(index==0) break;
			  }while(1);
		      }
		  }
	      }
	  }
      }
  }

  /* major loop to calculate the correction of bj due to contact */
  icutb=0;
  do{

    memcpy(&vini[0],&vold[0],sizeof(double)**nk);
    //for(i=0;i<*nk;i++){vold[i]=vini[i];}
    
    /* restoring aa[(iinc-1)*nev+i] (before change of *dtime) */
    i1=(*iinc-1)**nev;
    for(i=0;i<*nev;i++){
      aa[i]+=bb[i]*(*time1-*dtime);
    }
    
    /* reduce the size of the increment if the change in contact force
       is too large */
    
    bbmax=0.;
    if(icutb==0){
      for(i=0;i<neq[1];i++){
	if(fabs(bcont[i]-bcontini[i])>bbmax) {
	  bbmax=fabs(bcont[i]-bcontini[i]);
	}
      }
    }
    
    /* check for size of contact force */
    
    if(bbmax>*deltmx || icutb>0){
      
      /* force increase too big: increment size is decreased */
      if(icutb>0){
	*dtheta=*dtheta*df;
	printf("*INFORMATION: increment size is decreased to %e\nthe increment is reattempted\n\n",*dtheta**tper);
      }
      else{
	*dtheta=*dtheta**deltmx/bbmax;
      }
      printf("correction of dtheta due to contact: %e\n",*dtheta);
      *dthetaref=*dtheta;
      itp=0;
      
      /* check whether the new increment size is not too small */
      
      if(*dtheta<*tmin){
	  printf("\n *WARNING: increment size %e smaller than minimum %e\n",*dtheta**tper,*tmin**tper);
	  printf("             minimum is taken\n");
	  *dtheta=*tmin;
	  *dthetaref=*dtheta;
      }
      
      *reltime=*theta+(*dtheta);
      *time1=*reltime**tper;
      *dtime=*dtheta**tper;
      
      /* calculating the instantaneous loads (forces, surface loading, 
	 centrifugal and gravity loading or temperature) */
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
		xloadold,xload,xloadact,iamload,nload,ibody,xbody,
		nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
		namta,nam,ampli,time1,reltime,ttime,dtime,ithermal,
	        nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
	        ndirboun,nodeforc,
	        ndirforc,istep,iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
	        nelemload,sideload));
	      
      /* calculating the instantaneous loading vector */
      
      FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		   ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		   nforc,nelemload,sideload,xloadact,nload,xbodyact,
		   ipobody,nbody,cgr,b,nactdof,&neq[1],nmethod,
		   ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		   alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		   t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		   nplicon,plkcon,nplkcon,
		   npmat_,ttime,time1,istep,iinc,dtime,physcon,ibody,
		   xbodyold,reltime,veold,matname));
      
      /* correction for nonzero SPC's */
      
      if(*iprescribedboundary){
	  dynboun(amta,namta,nam,ampli,time1,ttime,dtime,
		  xbounold,xboun,
		  xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		  aub,icol,irow,neq,nzs,&sigma,b,isolver,
		  &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv);
      }
      memcpy(&bmech[0],&b[0],sizeof(double)*neq[1]);
      //      for(i=0;i<neq[1];i++){bmech[i]=b[i];}
      
    }
    
    /* calculating the linearized force function connecting the
       mechanical+contact load at the start of the increment with
       the mechanical load at the end of the increment
       = base load */
    i1=*iinc**nev;
    i3=(*iinc-1)**nev;
    for(i=0;i<*nev;i++){
      aa[i]=aanew[i];
      i2=i*neq[1];
      aanew[i]=0.;
      for(j=0;j<neq[1];j++){
	aanew[i]+=z[i2+j]*bmech[j];
      }
      bb[i]=(aanew[i]-aa[i])/(*dtime);
      aa[i]=aanew[i]-bb[i]**time1;
    }
    
    /* calculating the base response */
    bjbas=NNEW(double,*nev); /* basis response modal decomposition */
    bjbasp=NNEW(double,*nev);
    for(l=0;l<*nev;l++){
      zetaj=zeta[l];
      dj=d[l];
      
      /* zero eigenfrequency: rigid body mode */
      
      if(fabs(d[l])<=1.e-10){
	sum=0.;sump=0.;
	for(i=*iinc-1;i<*iinc;i++){
	  aai=aa[l];
	  bbi=bb[l];
	  tstart=*time0;
	  tend=*time1;
	  sum+=tend*(aai**time1+
		     tend*((bbi**time1-aai)/2.-bbi*tend/3.))-
	    tstart*(aai**time1+
		    tstart*((bbi**time1-aai)/2.-bbi*tstart/3.));
	  sump+=tend*(aai+bbi*tend/2.)-tstart*(aai+bbi*tstart/2.);
	}
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
	sum=0.;sump=0.;
	for(i=*iinc-1;i<*iinc;i++){
	  aai=aa[l];
	  bbi=bb[l];
	  //	  tstart=timea[*iinc]-timea[i+1];
	  tstart=0;
	  tend=*time1-*time0;
	  FORTRAN(fsub,(time1,&tend,&aai,&bbi,&ddj,
			&h1,&h2,&h3,&h4,&func,&funcp));
	  sum+=func;sump+=funcp;
	  FORTRAN(fsub,(time1,&tstart,&aai,&bbi,&ddj,
			&h1,&h2,&h3,&h4,&func,&funcp));
	  sum-=func;sump-=funcp;
	}
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
	sum=0.;sump=0.;
	for(i=*iinc-1;i<*iinc;i++){
	  aai=aa[l];
	  bbi=bb[l];
	  //	  tstart=timea[*iinc]-timea[i+1];
	  tstart=0;
	  tend=*time1-*time0;
	  FORTRAN(fsuper,(time1,&tend,&aai,&bbi,
			  &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	  sum+=func;sump+=funcp;
	  FORTRAN(fsuper,(time1,&tstart,&aai,&bbi,
			  &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	  sum-=func;sump-=funcp;
	}
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
	sum=0.;sump=0.;
	for(i=*iinc-1;i<*iinc;i++){
	  aai=aa[l];
	  bbi=bb[l];
	  //	  tstart=timea[*iinc]-timea[i+1];
	  tstart=0;
	  tend=*time1-*time0;
	  FORTRAN(fcrit,(time1,&tend,&aai,&bbi,&zetaj,&dj,
			 &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	  sum+=func;sump+=funcp;
	  FORTRAN(fcrit,(time1,&tstart,&aai,&bbi,&zetaj,&dj,
			 &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	  sum-=func;sump+=funcp;
	}
	fexp=exp(-h1**dtime);
	bjbas[l]=sum+fexp*((1.+h1**dtime)*cd[l]+*dtime*cv[l]);
	bjbasp[l]=sump+fexp*(-h1*h1**dtime*cd[l]+(1.-h1**dtime)*cv[l]);
      }
    }
    
    /* calculating the incremental response due to contact */
    
    aai=-(*time1-*dtime)/(*dtime);
    bbi=1./(*dtime);
    
    
    bjinc=NNEW(double,*nev); /* incremental response modal decomposition */
    bjincp=NNEW(double,*nev);
    for(l=0;l<*nev;l++){
      zetaj=zeta[l];
      dj=d[l];
      
      /* zero eigenfrequency: rigid body mode */
      
      if(fabs(d[l])<=1.e-10){
	tstart=*time0;
	tend=*time1;
	sum=tend*(aai**time1+
		  tend*((bbi**time1-aai)/2.-bbi*tend/3.))-
	  tstart*(aai**time1+
		  tstart*((bbi**time1-aai)/2.-bbi*tstart/3.));
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
	tend=*dtime;
	FORTRAN(fsub,(time1,&tend,&aai,&bbi,&ddj,
		      &h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsub,(time1,&tstart,&aai,&bbi,&ddj,
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
	tend=*dtime;
	FORTRAN(fsuper,(time1,&tend,&aai,&bbi,
			&h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fsuper,(time1,&tstart,&aai,&bbi,
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
	tend=*dtime;
	FORTRAN(fcrit,(time1,&tend,&aai,&bbi,&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum=func;sump=funcp;
	FORTRAN(fcrit,(time1,&tstart,&aai,&bbi,&zetaj,&dj,
		       &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	sum-=func;sump-=funcp;
	
	bjinc[l]=sum;
	bjincp[l]=sump;
	
      }
    }
    
    aaa=NNEW(double,neq[1]);
    bbb=NNEW(double,*nev*neq[1]);
    lhs=NNEW(double,*nev**nev);
    ipiv=NNEW(int,*nev);
    dbj=NNEW(double,*nev); /* change in bj */
    dbjp=NNEW(double,*nev); /* change in djp */
    
    memcpy(&bj[0],&bjbas[0],sizeof(double)**nev);
    memcpy(&bjp[0],&bjbasp[0],sizeof(double)**nev);
    /*    for(i=0; i<*nev; i++){
      bj[i]=bjbas[i];
      bjp[i]=bjbasp[i];
      }*/
    
    /* major iteration loop for the contact response */
    
    loop=0;
    printf("Contact-Iteration\n");
    do{
      start=time(NULL);
      loop++;
      printf("loop=%d\n",loop);
      
      /* composing the response */
      
      if(*iprescribedboundary){
	for(i=0;i<neq[1];i++){
	  b[i]=bmin[i];
	}
      }
      else{
	for(i=0;i<neq[1];i++){
	  b[i]=0.;
	  bp[i]=0.;
	}
      }
      
      for(i=0;i<neq[1];i++){
	for(j=0;j<*nev;j++){
	  b[i]+=bj[j]*z[j*neq[1]+i];
	  bp[i]+=bjp[j]*z[j*neq[1]+i];
	}
      }
      
      /* update nonlinear MPC-coefficients (e.g. for rigid
	 body MPC's */
      printf("Nonlinmpc_dynresults\n");
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,reltime,&newstep,xboun,fmpc,
			 iit,&idiscon,ncont,trab,ntrans,ithermal));
      
      /* calculating displacements/temperatures */
      
      FORTRAN(dynresults,(nk,v,ithermal,nactdof,vold,nodeboun,
			  ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			  b,bp,veold,dtime));
      
      memcpy(&vold[0],&v[0],sizeof(double)*5**nk);
      //      for(i=0;i<5**nk;i++){vold[i]=v[i];}
      
      if(iconvergence==1){
	break;
      }
      
      /* creating contact elements and calculating the contact forces
	 based on the displacements at the end of the present increment */
      
      memset(&bcont[0],0,sizeof(double)*neq[1]);
      //for(i=0;i<neq[1];i++){bcont[i]=0.;}
      //for(i=0;i<*nev*neq[1];i++){dbcont[i]=0.;}

      RENEW(dbcont,double,nactcont_**nev);
      RENEW(ikactcont,int,nactcont_);
      RENEW(ilactcont,int,nactcont_);
      
      //      memset(&dbcont[0],0,sizeof(double)*nactcont_**nev);
      //      memset(&iactcont[0],0,sizeof(double)*nactcont_);
      for(i=0; i<nactcont_**nev; i++){dbcont[i]=0.0;}
      for(i=0; i<nactcont_; i++){ikactcont[i]=0.0;ilactcont[i]=0.0;}
      
      *ne=*ne0;*nkon=*nkon0;
      end=time(NULL);
      printf("Zeit zum erstellen der L\"osung: %e\n",difftime(end,start));
      start=time(NULL);
      contact(ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
	      straight,nkon,co,vold,ielmat,cs,elcon,istep,
	      iinc,iit,ncmat_,ntmat_,ifcont1,ifcont2,ne0,
	      vini,nmethod,nmpc,mpcfree,memmpc_,
	      &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	      iperturb,ikboun,nboun);
      nactcont=0;
      printf("number of contact springs = %d\n",*ne-*ne0);
      for(i=*ne0;i<*ne;i++){
	indexe=ipkon[i];
	imat=ielmat[i];
	kodem=nelcon[2*imat-2];
	for(j=0;j<8;j++){lakonl[j]=lakon[8*i+j];}
	nope=atoi(&lakonl[7]);
	for(j=0;j<nope;j++){
	  konl[j]=kon[indexe+j];
	  for(j1=0;j1<3;j1++){
	    xl[j*3+j1]=co[3*(konl[j]-1)+j1];
	    voldl[j*4+j1+1]=vold[5*(konl[j]-1)+j1+1];
	    veoldl[j*3+j1]=veold[4*(konl[j]-1)+j1+1];
	  }
	}
	konl[nope]=kon[indexe+nope];
	FORTRAN(springforc,(xl,konl,voldl,&imat,elcon,nelcon,elas,
			    fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,
			    &t1l,&kodem,elconloc,plicon,nplicon,npmat_,
			    veoldl,&senergy,&iener,cstr));
	
	if(i==*ne0) printf("spring forc %e\n",fnl[3*nope-3]);
	if(i==*ne0) printf("spring disp %e\n",voldl[4*nope-3]);
	
	FORTRAN(springstiff,(xl,elas,konl,voldl,s,&imat,elcon,nelcon,
			     ncmat_,ntmat_,&nope,lakonl,&t0l,&t1l,&kode,elconloc,plicon,
			     nplicon,npmat_,iperturb));
	dfdbj(bcont,&dbcont,&neq[1],&nope,konl,nactdof,
	      s,z,ikmpc,ilmpc,ipompc,nodempc,nmpc,coefmpc,
	      fnl,nev,&ikactcont,&ilactcont,&nactcont,&nactcont_);
      }
      end=time(NULL);
      printf("Zeit zum berechnen der Kontaktelemente: %e\n",difftime(end,start));
      start=time(NULL);
      nactcont_=nactcont;
      RENEW(ikactcont,int,nactcont);
      RENEW(ilactcont,int,nactcont);
      RENEW(dbcont,double,nactcont**nev);
      /* aaa(i) is the internal product of the contact force at the end of the
	 increment with eigenmode i
	 bbb(i,j) is the internal product of the change of the contact force with
	 respect to modal coordinate j with the eigenmode i */

      //print iactcont,bbb for savety purposes

      /*      printf("nactcont: %d\n",nactcont);
            for(l=0; l<nactcont; l++){
	printf("%d -- %d\n",l,ikactcont[l]);
      }


      printf("dbcont-Ausgabe:\n");
      for(l=0; l<*nev; l++){
	for(j=0; j<nactcont; j++){
	  printf("%d -- %e\n ",j**nev+l,dbcont[(ilactcont[j]-1)**nev+l]);
	}
      }
      
            printf("bbb-Ausgabe:\n");
      
      */
      memset(&bbb[0],0,*nev**nev);
      memset(&aaa[0],0,neq[1]);
      /*#pragma omp parallel shared(bbb,z,dbcont,ikactcont,nev,neq) \
                     private(k,j,l,i1,i2,zl)
		       {
		       #pragma omp for schedule(static,100)*/
      for(k=0; k<nactcont; k++){
	i1=ikactcont[k]-1;
	i2=(ilactcont[k]-1)**nev;
	for(j=0; j<*nev; j++){
	  zl=z[j*neq[1]+i1];
	  aaa[j]+=zl*bcont[i1];
	  for(l=0; l<*nev; l++){
	    bbb[l**nev+j]+=zl*dbcont[i2+l];
	  }
	}
      }
      //	       }
      /*    for(l=0; l<*nev; l++){
	for(j=0; j<*nev; j++){
	  printf("%d -- %e\n",l**nev+j,bbb[l**nev+j]);
	}
	}*/
      

      end=time(NULL);
      printf("Zeit zum berechnen von bbb und aaa: %e\n",difftime(end,start));
      start=time(NULL);
      /* calculating the rhs and lhs of the system of equations; the unknowns are
	 the changes of the modal coordinates */
      
      /*      for (j=0;j<*nev;j++){
	dbj[j]=bjbas[j]+bjinc[j]*aaa[j]-bj[j];
	}*/

      end=time(NULL);
      printf("Zeit zum berechnen von RHS: %e\n",difftime(end,start));
      start=time(NULL);
      
      for(l=0;l<*nev;l++){
	i1=l**nev;
	for(j=0;j<*nev;j++){
	  if(j==l){lhs[i1+j]=1.;}else{lhs[i1+j]=0.;}
	  lhs[i1+j]-=bjinc[j]*bbb[i1+j];
	}
	dbj[l]=bjbas[l]+bjinc[l]*aaa[l]-bj[l];
      }

      end=time(NULL);
      printf("Zeit zum berechnen von LHS: %e\n",difftime(end,start));
      start=time(NULL);

      /* solve the system of equations; determine dbj */
      
      printf("Solving the system of equations\n");
      
      FORTRAN(dgesv,(nev,&nrhs,lhs,nev,ipiv,dbj,nev,&info));

      end=time(NULL);
      printf("Zeit zum Lösen des Gleichungssystems: %e\n",difftime(end,start));
      start=time(NULL);
      
      /* check the size of dbj */
      
      bjmax=0.;
      dbjmax=0.;
      dbjmaxOLD1=dbjmax;
      dbjmaxOLD2=dbjmaxOLD1;
      for(i=0;i<*nev;i++){
	if(fabs(bj[i])>bjmax) bjmax=fabs(bj[i]);
	if(fabs(dbj[i])>dbjmax) dbjmax=fabs(dbj[i]);
      }
      printf("bjmax= %e, dbjmax = %e\n",bjmax,dbjmax);
      
      //determine the change of the change of the contact force
      //compare dbcont with dbcontini, if fabs(dbcont-dbcontini)>0.5 then cutback
      
      /*dbdbcontmax=0.;
	for(j=0; j<neq[1]; j++){
	if(fabs(((bcont[j]-bcontini[j])/(*dtime))-dbcontini[j])>dbdbcontmax){
	dbdbcontmax=fabs((bcont[j]-bcontini[j])/(*dtime)-dbcontini[j]);
	}
	}
	printf("dbddbcontmax: %e\n",dbdbcontmax);*/
      

      iconvergence=0;
      idivergence=0;

      dbdbcontmax=0.0;
      if(dbjmax<=0.005*bjmax && dbdbcontmax<(0.2/(*dtime))){
	
	//calculate bjp: the derivative of bj w.r.t. time
	
	for(j=0; j<*nev; j++){
	  bjp[j]=bjbasp[j]+bjincp[j]*aaa[j];
	}
	FORTRAN(dgetrs,("No transpose",nev,&nrhs,lhs,nev,ipiv,bjp,nev,&info));
	//save dbcont for convergency purposes
	dbcontini=NNEW(double,nactcont);
	for(j=0; j<nactcont; j++){
	  dbcontini[j]=fabs((bcont[ikactcont[j]-1]-bcontini[ikactcont[j]-1])/(*dtime));
	}
	iconvergence=1;	   
      }
      else{
	if(dbdbcontmax>=(0.2/(*dtime))){
	  //cutback
	  printf("*INFORMATION: the change of the change of the contactforce is too big --> cutback\n");
	  idivergence=1;
	  icutb++;
	  break;
	}
	else if(loop>=i0 && loop<=ic){
	  /* check for divergence */
	  if((dbjmax>dbjmaxOLD1) && (dbjmax>dbjmaxOLD2)){
	    /* divergence --> cutback */ 
	    printf("*INFORMATION: divergence --> cutback\n");
	    idivergence=1;
	    icutb++;
	    break;
	  }
	}
	else{
	  if(loop>ic){
	    /* cutback after ic iterations*/
	    printf("*INFORMATION: to many iterations --> cutback\n");
	    idivergence=1;
	    icutb++;
	    break;
	  }
	}
      }

      end=time(NULL);
      printf("Zeit zur Konvergenzüberprüfung: %e\n",difftime(end,start));
      
      /*      if(loop>10){
	printf("*WARNING: nonconvergence at contact iteration -- increment: %d \n",*iinc);
      }
      //calculate bjp 
      for(j=0; j<*nev; j++){
	bjp[j]=bjbasp[j]+bjincp[j]*aaa[j];
      }
      FORTRAN(dgetrs,("No transpose",nev,&nrhs,lhs,nev,ipiv,bjp,nev,&info));
      //istiter=1;
      */
      /* add dbj to db */
      
      for(j=0;j<*nev;j++){
	bj[j]+=dbj[j];
      }
      
    }while(1);
  }while(idivergence==1 && icutb<10);

  printf("Contact-Iteration Done\n");

  if(icutb>=10){
    //non convergence, stop all
    printf("*ERROR: Contact did not converge.\n");
    FORTRAN(stop,());
  }
  
  /* restoring aa[(*iinc-1)*nev+i] */

  i1=(*iinc-1)**nev;
  for(i=0;i<*nev;i++){
    aa[i]+=bb[i]*(*time1-*dtime);
  }
  
  /* calculating the linearized force function connecting the
     mechanical+contact load at the start of the increment with
     the mechanical+contact load at the end of the increment
     = base load */
  
  i1=*iinc**nev;
  i3=(*iinc-1)**nev;
  for(i=0;i<*nev;i++){
    aa[i]=aanew[i];
    i2=i*neq[1];
    aanew[i]=0.;
    for(j=0;j<neq[1];j++){
      aanew[i]+=z[i2+j]*(bmech[j]+bcont[j]);
    }
    bb[i]=(aanew[i]-aa[i])/(*dtime);
    aa[i]=aanew[i]-bb[i]**time1;
  }

  *dbcontp=dbcont;
  *dbcontinip=dbcontini;
  
  free(aaa);free(bbb);free(bjbas);free(bjinc);free(dbj);free(lhs);
  free(ipiv);free(bjbasp);free(bjincp);free(dbjp);free(ikactcont);
  free(ilactcont);
  
  
  return;
}


