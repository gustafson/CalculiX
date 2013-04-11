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

void steadystate(double **cop, int *nk, int **konp, int **ipkonp, char **lakonp, int *ne, 
	       int **nodebounp, int **ndirbounp, double **xbounp, int *nboun,
	       int **ipompcp, int **nodempcp, double **coefmpcp, char **labmpcp,
               int *nmpc, int *nodeforc,int *ndirforc,double *xforc, 
               int *nforc,int *nelemload, char *sideload,double *xload,
	       int *nload, 
	       int **nactdofp,int *neq, int *nzl,int *icol, int *irow, 
	       int *nmethod, int **ikmpcp, int **ilmpcp, int **ikbounp, 
	       int **ilbounp,double *elcon, int *nelcon, double *rhcon, 
	       int *nrhcon,double *cocon, int *ncocon,
               double *alcon, int *nalcon, double *alzero, 
               int **ielmatp,int **ielorienp, int *norien, double *orab, 
               int *ntmat_,double **t0p, 
	       double **t1p,int *ithermal,double *prestr, int *iprestr, 
	       double **voldp,int *iperturb, double *sti, int *nzs, 
	       double *tinc, double *tper, double *xmodal,
	       double *veold, char *amname, double *amta,
	       int *namta, int *nam, int *iamforc, int *iamload,
	       int **iamt1p,int *jout,
	       int *kode, char *filab,double **emep, double *xforcold, 
	       double *xloadold,
               double **t1oldp, int **iambounp, double **xbounoldp, int *iexpl,
               double *plicon, int *nplicon, double *plkcon,int *nplkcon,
               double *xstate, int *npmat_, char *matname, int *mint_,
               int *ncmat_, int *nstate_, double **enerp, char *jobnamec,
               double *ttime, char *set, int *nset, int *istartset,
               int *iendset, int *ialset, int *nprint, char *prlab,
               char *prset, int *nener, double *trab, 
               int **inotrp, int *ntrans, double **fmpcp, char *cbody, int *ibody,
               double *xbody, int *nbody, double *xbodyold, int *istep,
               int *isolver, int *jq, char *output, int *mcs,int *nkon, 
               int *ics, double *cs, int *mpcend){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *labmpcold=NULL;

  int nev,i,j,k, *inum=NULL,*ipobody=NULL,inewton=0,nsectors,
    iinc,l,iout,ielas,icmd,iprescribedboundary,ndata,nmd,nevd,
    ndatatot,*iphaseforc=NULL,*iphaseload=NULL,*iphaseboun=NULL,
    *isave=NULL,nfour,ii,ir,ic,mode=-1,noddiam=-1,*nm=NULL,
    *kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,*inotr=NULL,
    *nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,*ipompcold=NULL,*nodempcold=NULL,*ikmpcold=NULL,
    *ilmpcold=NULL,nmpcold,mpcendold,kflag=2,*iamt1=NULL,ifreebody,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,dashpot,nrhs=1,
    *ipiv=NULL,info,nev2;

  double *d=NULL, *z=NULL,*stiini=NULL,*vini=NULL,*freqnh=NULL,
    *xforcact=NULL, *xloadact=NULL,y,*fr=NULL,*fi=NULL,*cc=NULL,
    *t1act=NULL, *ampli=NULL, *aa=NULL, *bb=NULL, *vr=NULL,*vi=NULL,
    *stn=NULL, *stx=NULL, *een=NULL, *adb=NULL,*xstiff=NULL,
    *aub=NULL, *aux=NULL, *bjr=NULL, *bji=NULL,*xbodyr=NULL,
    *f=NULL, *fn=NULL, *xbounact=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL, *xbodyact=NULL, *cgr=NULL, *au=NULL,*xbodyi=NULL,
    time,dtime,reltime,*co=NULL,*xboun=NULL,*xbounold=NULL,
    physcon[1],qa[2],cam[2],accold[1],bet,gam,*ad=NULL,sigma=0.,alpham,betam,
    fmin,fmax,bias,*freq=NULL,*xforcr=NULL,dd,pi,vreal,constant,
    *xforci=NULL,*xloadr=NULL,*xloadi=NULL,*xbounr=NULL,*xbouni=NULL,
    *br=NULL,*bi=NULL,*ubr=NULL,*ubi=NULL,*mubr=NULL,*mubi=NULL,
    *wsave=NULL,*r=NULL,*xbounacttime=NULL,*btot=NULL,breal,tmin,tmax,
    *vold=NULL,*eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,
    *coefmpcold=NULL,*t0=NULL,*t1=NULL,*t1old=NULL,*adc=NULL,*auc=NULL,
    *am=NULL,*bm=NULL,*zc=NULL;

  FILE *f1;

  int *ipneigh=NULL,*neigh=NULL;

#ifdef SGI
  int token;
#endif

  co=*cop;kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;
  ielorien=*ielorienp;inotr=*inotrp;nodeboun=*nodebounp;
  ndirboun=*ndirbounp;iamboun=*iambounp;xboun=*xbounp;
  xbounold=*xbounoldp;ikboun=*ikbounp;ilboun=*ilbounp;nactdof=*nactdofp;
  vold=*voldp;eme=*emep;ener=*enerp;ipompc=*ipompcp;nodempc=*nodempcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;

  xstiff=NNEW(double,27**mint_**ne);

  pi=4.*atan(1.);
  iout=1;

  alpham=xmodal[0];
  betam=xmodal[1];
  fmin=2.*pi*xmodal[2];
  fmax=2.*pi*xmodal[3];
  ndata=floor(xmodal[4]);
  bias=xmodal[5];
  nfour=floor(xmodal[6]);
  if(nfour>0){
      tmin=xmodal[7];
      tmax=xmodal[8];
  }

  /* determining nzl */

  *nzl=0;
  for(i=neq[1];i>0;i--){
      if(icol[i-1]>0){
	  *nzl=i;
	  break;
      }
  }

  strcpy(fneig,jobnamec);
  strcat(fneig,".eig");

  if((f1=fopen(fneig,"rb"))==NULL){
    printf("*ERROR: cannot open eigenvalue file for reading...");
    exit(0);
  }

  if(*mcs==0){
      if(fread(&nev,sizeof(int),1,f1)!=1){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      d=NNEW(double,nev);
      
      if(fread(d,sizeof(double),nev,f1)!=nev){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      ad=NNEW(double,neq[1]);
      adb=NNEW(double,neq[1]);
      au=NNEW(double,nzs[2]);
      aub=NNEW(double,nzs[1]);
      
      if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
      
      z=NNEW(double,neq[1]*nev);
      
      if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
	  printf("*ERROR reading the eigenvalue file...");
	  exit(0);
      }
  }
  else{
      nev=0;
      do{
	  if(fread(&nmd,sizeof(int),1,f1)!=1){
	      break;
	  }
	  if(fread(&nevd,sizeof(int),1,f1)!=1){
	      printf("*ERROR reading the eigenvalue file...");
	      exit(0);
	  }
	  if(nev==0){
	      d=NNEW(double,nevd);
	      nm=NNEW(int,nevd);
	  }else{
	      RENEW(d,double,nev+nevd);
	      RENEW(nm,int,nev+nevd);
	  }
	  
	  if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
	      printf("*ERROR reading the eigenvalue file...");
	      exit(0);
	  }
	  for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	  
	  if(nev==0){
	      adb=NNEW(double,neq[1]);
	      aub=NNEW(double,nzs[1]);

	      if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
		  printf("*ERROR reading the eigenvalue file...");
		  exit(0);
	      }
	      
	      if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
		  printf("*ERROR reading the eigenvalue file...");
		  exit(0);
	      }
	  }
	  
	  if(nev==0){
	      z=NNEW(double,neq[1]*nevd);
	  }else{
	      RENEW(z,double,neq[1]*(nev+nevd));
	  }
	  
	  if(fread(&z[neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	      printf("*ERROR reading the eigenvalue file...");
	      exit(0);
	  }
	  nev+=nevd;
      }while(1);

      /* determining the maximum amount of segments */

      nsectors=0;
      for(i=0;i<*mcs;i++){
	  if(cs[17*i]>nsectors) nsectors=cs[17*i];
      }

      /* allocating field for the expanded structure */

      RENEW(co,double,3**nk*nsectors);
      if(*ithermal!=0){
	  RENEW(t0,double,*nk*nsectors);
	  RENEW(t1old,double,*nk*nsectors);
	  RENEW(t1,double,*nk*nsectors);
	  if(*nam>0) RENEW(iamt1,int,*nk*nsectors);
      }
      RENEW(nactdof,int,4**nk*nsectors);
      if(*ntrans>0) RENEW(inotr,int,2**nk*nsectors);
      RENEW(kon,int,*nkon*nsectors);
      RENEW(ipkon,int,*ne*nsectors);
      RENEW(lakon,char,8**ne*nsectors);
      RENEW(ielmat,int,*ne*nsectors);
      if(*norien>0) RENEW(ielorien,int,*ne*nsectors);
      RENEW(z,double,neq[1]*nev*nsectors/2);

      RENEW(nodeboun,int,*nboun*nsectors);
      RENEW(ndirboun,int,*nboun*nsectors);
      if(*nam>0) RENEW(iamboun,int,*nboun*nsectors);
      RENEW(xboun,double,*nboun*nsectors);
      RENEW(xbounold,double,*nboun*nsectors);
      RENEW(ikboun,int,*nboun*nsectors);
      RENEW(ilboun,int,*nboun*nsectors);

      ipompcold=NNEW(int,*nmpc);
      nodempcold=NNEW(int,3**mpcend);
      coefmpcold=NNEW(double,*mpcend);
      labmpcold=NNEW(char,20**nmpc);
      ikmpcold=NNEW(int,*nmpc);
      ilmpcold=NNEW(int,*nmpc);

      for(i=0;i<*nmpc;i++){ipompcold[i]=ipompc[i];}
      for(i=0;i<3**mpcend;i++){nodempcold[i]=nodempc[i];}
      for(i=0;i<*mpcend;i++){coefmpcold[i]=coefmpc[i];}
      for(i=0;i<20**nmpc;i++){labmpcold[i]=labmpc[i];}
      for(i=0;i<*nmpc;i++){ikmpcold[i]=ikmpc[i];}
      for(i=0;i<*nmpc;i++){ilmpcold[i]=ilmpc[i];}
      nmpcold=*nmpc;
      mpcendold=*mpcend;

      RENEW(ipompc,int,*nmpc*nsectors);
      RENEW(nodempc,int,3**mpcend*nsectors);
      RENEW(coefmpc,double,*mpcend*nsectors);
      RENEW(labmpc,char,20**nmpc*nsectors);
      RENEW(ikmpc,int,*nmpc*nsectors);
      RENEW(ilmpc,int,*nmpc*nsectors);
      RENEW(fmpc,double,*nmpc*nsectors);

      expand(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	ipompc,nodempc,coefmpc,labmpc,nmpc,nodeforc,ndirforc,xforc,
        nforc,nelemload,sideload,xload,nload,nactdof,neq,
	nmethod,ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	t0,ithermal,prestr,iprestr,vold,iperturb,sti,nzs,
	adb,aub,filab,eme,plicon,nplicon,plkcon,nplkcon,
        xstate,npmat_,matname,mint_,ics,cs,mpcend,ncmat_,
        nstate_,mcs,nkon,ener,jobnamec,output,set,nset,istartset,
        iendset,ialset,nprint,prlab,prset,nener,trab,
        inotr,ntrans,ttime,fmpc,&nev,z,iamboun,xbounold,
	&nsectors,nm,icol,irow,nzl,nam,ipompcold,nodempcold,coefmpcold,
        labmpcold,&nmpcold,xloadold,iamload,t1old,t1,iamt1,xstiff);

      free(vold);vold=NNEW(double,4**nk);
      RENEW(eme,double,6**mint_**ne);
      RENEW(xstiff,double,27**mint_**ne);
      if(*nener==1) RENEW(ener,double,*mint_**ne);
  }

  fclose(f1);

  /* check whether there are dashpot elements */

  dashpot=0;
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[i*8],"ED")==0){
	  dashpot=1;break;}
  }
  if(dashpot){
      adc=NNEW(double,neq[1]);
      auc=NNEW(double,nzs[1]);
      FORTRAN(mafilldm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mint_,ncmat_,
	      ttime,&time,istep,&iinc,ibody));

      /*  zc = damping matrix * eigenmodes */

      zc=NNEW(double,neq[1]*nev);
      for(i=0;i<nev;i++){
	  FORTRAN(op,(&neq[1],aux,&z[i*neq[1]],&zc[i*neq[1]],adc,auc,
	  icol,irow,nzl));
      }

      /* cc is the reduced damping matrix (damping matrix mapped onto
         space spanned by eigenmodes) */

      cc=NNEW(double,nev*nev);
      for(i=0;i<nev;i++){
	  for(j=0;j<=i;j++){
	      for(k=0;k<neq[1];k++){
		  cc[i*nev+j]+=z[j*neq[1]+k]*zc[i*neq[1]+k];
	      }
	  }
      }

      /* symmetric part of cc matrix */

      for(i=0;i<nev;i++){
	  for(j=i;j<nev;j++){
	      cc[i*nev+j]=cc[j*nev+i];
	  }
      }
      free(zc);
      nev2=2*nev;
      am=NNEW(double,nev2*nev2);
      bm=NNEW(double,nev2);
      ipiv=NNEW(int,nev2);
  }

  if(nfour<=0){

      /* harmonic excitation */
      
      /* determining the frequency data points */
      
      freq=NNEW(double,ndata*(nev+1));
      
      ndatatot=0.;
      freq[0]=fmin;
      if(fabs(fmax-fmin)<1.e-10){
	  ndatatot=1;
      }else{
	  for(i=0;i<nev;i++){
	      if(d[i]>=fmin){
		  if(d[i]<=fmax){
		      for(j=1;j<ndata;j++){
			  y=-1.+2.*j/((double)(ndata-1));
			  if(fabs(y)<1.e-10){freq[ndatatot+j]=
                                 (freq[ndatatot]+d[i])/2.;}
			  else{
			      freq[ndatatot+j]=(freq[ndatatot]+d[i])/2.+
				  (d[i]-freq[ndatatot])*pow(fabs(y),1./bias)*
                                  y/(2.*fabs(y));
			  }
		      }
		      ndatatot+=ndata-1;
		  }
		  else{break;}
	      }
	  }
	  for(j=1;j<ndata;j++){
	      y=-1.+2.*j/((double)(ndata-1));
	      if(fabs(y)<1.e-10){freq[ndatatot+j]=(freq[ndatatot]+fmax)/2.;}
	      else{
		  freq[ndatatot+j]=(freq[ndatatot]+fmax)/2.+
		      (fmax-freq[ndatatot])*pow(fabs(y),1./bias)*
                      y/(2.*fabs(y));
	      }
	  }
	  ndatatot+=ndata;
      }
      RENEW(freq,double,ndatatot);
      
      /*  check for nonzero SPC's */
      
      iprescribedboundary=0;
      for(i=0;i<*nboun;i++){
	  if(fabs(xboun[i])>1.e-10){
	      iprescribedboundary=1;
	      break;
	  }
      }
      
      /* check whether the loading is real or imaginary */
      
      iphaseforc=NNEW(int,*nforc);
      for (i=0;i<*nforc;i++){
	  if(nodeforc[2*i]>*nk){
	      iphaseforc[i]=1;
	      nodeforc[2*i]=nodeforc[2*i]-*nk;
	  }
      }
      
      iphaseload=NNEW(int,*nload);
      for (i=0;i<*nload;i++){
	  if(nelemload[2*i]>*ne){
	      iphaseload[i]=1;
	      nelemload[2*i]=nelemload[2*i]-*ne;
	  }
      }
      
      if(iprescribedboundary){
	  iphaseboun=NNEW(int,*nboun);
	  for (i=0;i<*nboun;i++){
	      if(nodeboun[i]>*nk){
		  iphaseboun[i]=1;
		  nodeboun[i]=nodeboun[i]-*nk;
	      }
	  }
      }
      
      /* allocating actual loading fields */
      
      xforcact=NNEW(double,*nforc);
      xforcr=NNEW(double,*nforc);
      xforci=NNEW(double,*nforc);
      
      xloadact=NNEW(double,2**nload);
      xloadr=NNEW(double,2**nload);
      xloadi=NNEW(double,2**nload);
      
      xbodyact=NNEW(double,7**nbody);
      xbodyr=NNEW(double,7**nbody);
      xbodyi=NNEW(double,7**nbody);
      /* copying the rotation axis and/or acceleration vector */
      for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
      
      xbounact=NNEW(double,*nboun);
      
      if(*ithermal==1) t1act=NNEW(double,*nk);
      
      /* assigning the body forces to the elements */ 

      if(*nbody>0){
	  ifreebody=*ne+1;
	  ipobody=NNEW(int,2*ifreebody);
	  for(k=1;k<=*nbody;k++){
	      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
				 iendset,ialset,&inewton,nset,&ifreebody,&k));
	      RENEW(ipobody,int,2*(*ne+ifreebody));
	  }
	  RENEW(ipobody,int,2*(ifreebody-1));
      }
      
      br=NNEW(double,neq[1]); /* load rhs vector */
      bi=NNEW(double,neq[1]); /* load rhs vector */
      
      if(iprescribedboundary){
	  xbounr=NNEW(double,*nboun);
	  xbouni=NNEW(double,*nboun);
	  
	  fr=NNEW(double,neq[1]); /* force corresponding to real particular solution */
	  fi=NNEW(double,neq[1]); /* force corresponding to imaginary particular solution */
	  
	  ubr=NNEW(double,neq[1]); /* real particular solution */
	  ubi=NNEW(double,neq[1]); /* imaginary particular solution */
	  
	  mubr=NNEW(double,neq[1]); /* mass times real particular solution */
	  mubi=NNEW(double,neq[1]); /* mass times imaginary particular solution */
      }
      
      bjr=NNEW(double,nev); /* real response modal decomposition */
      bji=NNEW(double,nev); /* imaginary response modal decomposition */
      
      ampli=NNEW(double,*nam); /* instantaneous amplitude */
      
      aa=NNEW(double,nev); /* modal coefficients of the real loading */
      bb=NNEW(double,nev); /* modal coefficients of the imaginary loading */
      
      /* result fields */
      
      vr=NNEW(double,4**nk);
      vi=NNEW(double,4**nk);
      fn=NNEW(double,4**nk);
      stn=NNEW(double,6**nk);
      inum=NNEW(int,*nk);
      stx=NNEW(double,6**mint_**ne);
      if(*ithermal>1) {qfn=NNEW(double,3**nk);qfx=NNEW(double,3**mint_**ne);}
      
      if(strcmp1(&filab[18],"E   ")==0) een=NNEW(double,6**nk);
      if(strcmp1(&filab[36],"ENER")==0) enern=NNEW(double,*nk);
      
      eei=NNEW(double,6**mint_**ne);
      if(*nener==1){
	  stiini=NNEW(double,6**mint_**ne);
	  enerini=NNEW(double,*mint_**ne);}
      
      if(iprescribedboundary){
	  
	  /* LU decomposition of the stiffness matrix */
	  
	  if(*isolver==0){
#ifdef SPOOLES
	      spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
                      &symmetryflag,&inputformat);
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
      }
      
      for(l=0;l<ndatatot;l=l+*jout){
	  for(i=0;i<6**mint_**ne;i++){eme[i]=0.;}
	  time=freq[l]/(2.*pi);
	  *ttime=time;
	  
	  /* calculating the instantaneous loads (forces, surface loading, 
	     centrifugal and gravity loading or temperature) */
	  
	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	       xloadold,xload,xloadact,iamload,nload,ibody,xbody,
	       nbody,xbodyold,xbodyact,t1old,t1,t1act,
	       iamt1,nk,amta,namta,nam,ampli,&time,
	       &reltime,ttime,&dtime,ithermal,nmethod,
	       xbounold,xboun,xbounact,iamboun,nboun,
	       nodeboun,ndirboun,nodeforc,
	       ndirforc,istep,&iinc,co,vold,itg,&ntg));
	  
	  /* real part of forces */
	  
	  for (i=0;i<*nforc;i++){
	      xforcr[i]=xforcact[i]*(1-iphaseforc[i]);
	  }
	  
	  for (i=0;i<*nload;i++){
	      for(j=0;j<2;j++){
		  xloadr[2*i+j]=xloadact[2*i+j]*(1-iphaseload[i]);
	      }
	  }

	  for(i=0;i<*nbody;i++){
	      for(j=0;j<7;j++){
		  xbodyr[7*i+j]=xbodyact[7*i+j];
	      }
	      if(ibody[3*i+2]==2){
		  xbodyr[7*i]=0.;
	      }
	  }
	  
	  /* calculating the instantaneous loading vector */
	  
	  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcr,
		       nforc,nelemload,sideload,xloadr,nload,xbodyr,
		       ipobody,nbody,cgr,br,nactdof,&neq[1],nmethod,
		       ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		       alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		       t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		       nplicon,plkcon,nplkcon,
		       npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody));
	  
	  /* correction for nonzero SPC's */
	  
	  if(iprescribedboundary){
	      
	      /* real part of boundary conditions */
	      
	      for (i=0;i<*nboun;i++){
		  xbounr[i]=xbounact[i]*(1-iphaseboun[i]);
	      }
	      
	      for(j=0;j<neq[1];j++){fr[j]=0.;ubr[j]=0.;}
	      for(i=0;i<*nboun;i++){
		  ic=neq[1]+i;
		  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
		      ir=irow[j]-1;
		      fr[ir]=fr[ir]-au[j]*xbounr[i];
		      ubr[ir]=fr[ir];
		  }
	      }
	      if(*isolver==0){
#ifdef SPOOLES
		  spooles_solve(ubr,&neq[1]);
#endif
	      }
	      else if(*isolver==4){
#ifdef SGI
		  sgi_solve(ubr,token);
#endif
	      }
	      else if(*isolver==5){
#ifdef TAUCS
		  tau_solve(ubr,&neq[1]);
#endif
	      }
	      FORTRAN(op,(&neq[1],aux,ubr,mubr,adb,aub,
			  icol,irow,nzl));
	  }
	  
	  /* imaginary part of forces */
	  
	  for (i=0;i<*nforc;i++){
	      xforci[i]=xforcact[i]*iphaseforc[i];
	  }
	  
	  for (i=0;i<*nload;i++){
	      for(j=0;j<2;j++){
		  xloadi[2*i+j]=xloadact[2*i+j]*iphaseload[i];
	      }
	  }
	  
	  for(i=0;i<*nbody;i++){
	      for(j=0;j<7;j++){
		  xbodyi[7*i+j]=xbodyact[7*i+j];
	      }
	      if(ibody[3*i+2]==1){
		  xbodyi[7*i]=0.;
	      }
	  }
	  
	  /* calculating the instantaneous loading vector */
	  
	  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforci,
		       nforc,nelemload,sideload,xloadi,nload,xbodyi,
		       ipobody,nbody,cgr,bi,nactdof,&neq[1],nmethod,
		       ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		       alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		       t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		       nplicon,plkcon,nplkcon,
		       npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody));
	  
	  /* correction for nonzero SPC's */
	  
	  if(iprescribedboundary){
	      
	      /* imaginary part of boundary conditions */
	      
	      for (i=0;i<*nboun;i++){
		  xbouni[i]=xbounact[i]*iphaseboun[i];
	      }
	      
	      for(j=0;j<neq[1];j++){fi[j]=0.;ubi[j]=0.;}
	      for(i=0;i<*nboun;i++){
		  ic=neq[1]+i;
		  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
		      ir=irow[j]-1;
		      fi[ir]=fi[ir]-au[j]*xbouni[i];
		      ubi[ir]=fi[ir];
		  }
	      }
	      if(*isolver==0){
#ifdef SPOOLES
		  spooles_solve(ubi,&neq[1]);
#endif
	      }
	      else if(*isolver==4){
#ifdef SGI
		  sgi_solve(ubi,token);
#endif
	      }
	      else if(*isolver==5){
#ifdef TAUCS
		  tau_solve(ubi,&neq[1]);
#endif
	      }
	      FORTRAN(op,(&neq[1],aux,ubi,mubi,adb,aub,
			  icol,irow,nzl));
	  }
	  
	  /* correction for prescribed boundary conditions */
	  
	  if(iprescribedboundary){
	      for(i=0;i<neq[1];i++){
		  br[i]+=freq[l]*(freq[l]*mubr[i]+alpham*mubi[i]+betam*fi[i]);
		  bi[i]+=freq[l]*(freq[l]*mubi[i]-alpham*mubr[i]-betam*fr[i]);
	      }
	  }
	  
	  /* real and imaginary modal coefficients */
	  
	  for(i=0;i<nev;i++){
	      aa[i]=0.;
	      for(j=0;j<neq[1];j++){
		  aa[i]+=z[i*neq[1]+j]*br[j];
	      }
	  }
	  
	  for(i=0;i<nev;i++){
	      bb[i]=0.;
	      for(j=0;j<neq[1];j++){
		  bb[i]+=z[i*neq[1]+j]*bi[j];
	      }
	  }
	  
	  /* calculating the modal coefficients */
	  
	  if(dashpot==0){
	      for(i=0;i<nev;i++){
		  dd=pow(pow(d[i],2)-pow(freq[l],2),2)+
		      pow(alpham+betam*pow(d[i],2),2)*pow(freq[l],2);
		  bjr[i]=(aa[i]*(d[i]*d[i]-freq[l]*freq[l])+
			  bb[i]*(alpham+betam*d[i]*d[i])*freq[l])/dd;
		  bji[i]=(bb[i]*(d[i]*d[i]-freq[l]*freq[l])-
			  aa[i]*(alpham+betam*d[i]*d[i])*freq[l])/dd;
	      }
	      printf("old l=%d,bjr=%f,bji=%f\n",l,bjr[0],bji[0]);
	  }else{
	      for(i=0;i<nev2;i++){
		  for(j=0;j<nev2;j++){
		      am[i*nev2+j]=0.;
		  }
		  bm[i]=0.;
	      }
	      for(i=0;i<nev;i++){
		  am[i*nev2+i]=d[i]*d[i]-freq[l]*freq[l];
		  am[(i+nev)*nev2+i]=-(alpham+betam*d[i]*d[i])*freq[l];
		  bm[i]=aa[i];
		  am[i*nev2+nev+i]=-am[(i+nev)*nev2+i];
		  am[(i+nev)*nev2+nev+i]=am[i*nev2+i];
		  bm[nev+i]=bb[i];
		  for(j=0;j<nev;j++){
		      am[(j+nev)*nev2+i]=am[(j+nev)*nev2+i]
                           -cc[i*nev+j]*freq[l];
		      am[j*nev2+nev+i]=am[j*nev2+nev+i]
			   +cc[i*nev+j]*freq[l];
		  }
	      }

              /* solving the system of equations */

	      FORTRAN(dgesv,(&nev2,&nrhs,am,&nev2,ipiv,bm,&nev2,&info));
	      if(info!=0){
		  printf("*ERROR in steadystate: fatal termination of dgesv\n");
		  printf("       info=%d\n",info);
/*		  FORTRAN(stop,());*/
	      }
	      
	      /* storing the solution in bjr and bji */

	      for(i=0;i<nev;i++){
		  bjr[i]=bm[i];
		  bji[i]=bm[nev+i];
	      }

	  }
	  /*     printf("new l=%d,bjr=%f,bji=%f\n",l,bjr[0],bji[0]);*/
	  
	  /* calculating the real response */
	  
	  if(iprescribedboundary){
	      for(i=0;i<neq[1];i++){
		  br[i]=ubr[i];
	      }
	  }
	  else{
	      for(i=0;i<neq[1];i++){
		  br[i]=0.;
	      }
	  }
	  for(i=0;i<neq[1];i++){
	      for(j=0;j<nev;j++){
		  br[i]+=bjr[j]*z[j*neq[1]+i];
	      }
	  }
	  
	  if(iprescribedboundary){
	      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,
		  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		  ithermal,prestr,iprestr,filab,eme,een,
		  iperturb,f,fn,nactdof,&iout,qa,
		  vold,br,nodeboun,ndirboun,xbounr,nboun,
		  ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		  veold,accold,&bet,&gam,&dtime,&time,ttime,
		  plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
		  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		  enern,sti,xstaten,eei,enerini,cocon,ncocon,
		  set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		  qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));}
	  else{
	      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,
		  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		  ithermal,prestr,iprestr,filab,eme,een,
		  iperturb,f,fn,nactdof,&iout,qa,
		  vold,br,nodeboun,ndirboun,xbounact,nboun,
		  ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		  veold,accold,&bet,&gam,&dtime,&time,ttime,
		  plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
		  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		  enern,sti,xstaten,eei,enerini,cocon,ncocon,
		  set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		  qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));}
	  
	  (*kode)++;
	  

    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
	  FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,nmethod,kode,filab,
            een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,nstate_,istep,&iinc,
	    iperturb,ener,mint_,output,ithermal,qfn,&mode,&noddiam,
            trab,inotr,ntrans,orab,ielorien,norien,description,
                     ipneigh,neigh,sti));free(ipneigh);free(neigh);
	  
	  /* calculating the imaginary response */
	  
	  if(iprescribedboundary){
	      for(i=0;i<neq[1];i++){
		  bi[i]=ubi[i];
	      }
	  }
	  else{
	      for(i=0;i<neq[1];i++){
		  bi[i]=0.;
	      }
	  }
	  for(i=0;i<neq[1];i++){
	      for(j=0;j<nev;j++){
		  bi[i]+=bji[j]*z[j*neq[1]+i];
	      }
	  }
	  
	  if(iprescribedboundary){
	      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vi,stn,inum,
		  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		  ithermal,prestr,iprestr,filab,eme,een,
                  iperturb,f,fn,nactdof,&iout,qa,
                  vold,bi,nodeboun,ndirboun,xbouni,nboun,
                  ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
                  veold,accold,&bet,&gam,&dtime,&time,ttime,
                  plicon,nplicon,plkcon,nplkcon,
                  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
                  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
                  enern,sti,xstaten,eei,enerini,cocon,ncocon,
                  set,nset,istartset,iendset,ialset,nprint,prlab,prset,
	          qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));}
	  else{ 
	      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vi,stn,inum,
                  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
                  ielmat,ielorien,norien,orab,ntmat_,t0,t1,
                  ithermal,prestr,iprestr,filab,eme,een,
                  iperturb,f,fn,nactdof,&iout,qa,
                  vold,bi,nodeboun,ndirboun,xbounact,nboun,
                  ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
                  veold,accold,&bet,&gam,&dtime,&time,ttime,
                  plicon,nplicon,plkcon,nplkcon,
                  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
                  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
                  enern,sti,xstaten,eei,enerini,cocon,ncocon,
                  set,nset,istartset,iendset,ialset,nprint,prlab,prset,
	          qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));}
 
	  (*kode)++;
	  
    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
	  FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vi,stn,inum,nmethod,kode,filab,
            een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,nstate_,istep,&iinc,
	    iperturb,ener,mint_,output,ithermal,qfn,&mode,&noddiam,
            trab,inotr,ntrans,orab,ielorien,norien,description,
                     ipneigh,neigh,sti));free(ipneigh);free(neigh);
	  
	  /* calculating the magnitude and phase */
	  
	  if(strcmp1(&filab[60],"PU")==0){
	      
	      constant=180./pi;
	      
	      if(*ithermal<=1){
		  for(i=0;i<*nk;i++){
		      for(j=1;j<4;j++){
			  vreal=vr[4*i+j];
			  vr[4*i+j]=sqrt(vr[4*i+j]*vr[4*i+j]+vi[4*i+j]*vi[4*i+j]);
			  if(fabs(vreal)<1.e-10){
			      if(vi[4*i+j]>0.){vi[4*i+j]=90.;}
			      else{vi[4*i+j]=-90.;}
			  }
			  else{
			      vi[4*i+j]=atan(vi[4*i+j]/vreal)*constant;
			      if(vreal<0.) vi[4*i+j]+=180.;
			  }
		      }
		  }
	      }
	      else{
		  for(i=0;i<*nk;i++){
		      vreal=vr[4*i];
		      vr[4*i]=sqrt(vr[4*i]*vr[4*i]+vi[4*i]*vi[4*i]);
		      if(fabs(vreal)<1.e-10){
			  if(vi[4*i]>0){vi[4*i]=90.;}
			  else{vi[4*i]=-90.;}
		      }
		      else{
			  vi[4*i]=atan(vi[4*i]/vreal)*constant;
			  if(vreal<0.) vi[4*i]+=180.;
		      }
		  }
	      }
	      
	      (*kode)++;
	      
	      FORTRAN(frdphase,(kode,&time,nk,inum,vr,vi,filab));
	  }
	  
      }
      
      /* restoring the imaginary loading */
      
      for (i=0;i<*nforc;i++){
	  if(iphaseforc[i]==1){
	      nodeforc[2*i]=nodeforc[2*i]+*nk;
	  }
      }
      free(iphaseforc);free(xforcr);free(xforci);

      for (i=0;i<*nload;i++){
	  if(iphaseload[i]==1){
	      nelemload[2*i]=nelemload[2*i]+*ne;
	  }
      }
      free(iphaseload);free(xloadr);free(xloadi);
      
      free(xbodyr);free(xbodyi);
      
      if(iprescribedboundary){
	  for (i=0;i<*nboun;i++){
	      if(iphaseboun[i]==1){
		  nodeboun[i]=nodeboun[i]+*nk;
	      }
	  }
	  free(iphaseboun);
      }
      
  /* freeing the result fields */
      
      free(eei);
      if(*nener==1){free(stiini);free(enerini);}
      
      if(strcmp1(&filab[18],"E   ")==0) free(een);
      if(strcmp1(&filab[36],"ENER")==0) free(enern);
      if(*ithermal>1) {free(qfn);free(qfx);}
      
      /* updating the loading at the end of the step; 
	 important in case the amplitude at the end of the step
	 is not equal to one */
      
      for(k=0;k<*nboun;++k){xboun[k]=xbounact[k];}
      for(k=0;k<*nforc;++k){xforc[k]=xforcact[k];}
      for(k=0;k<2**nload;++k){xload[k]=xloadact[k];}
      for(k=0;k<7**nbody;k=k+7){xbody[k]=xbodyact[k];}
      if(*ithermal==1){
	  for(k=0;k<*nk;++k){t1[k]=t1act[k];}
      }
      
      free(fn);free(stn);free(inum);free(stx);
      free(br);free(bi);free(bjr);free(bji),free(freq);
      free(xforcact);free(xloadact);free(xbounact);free(aa);free(bb);
      free(ampli);free(xbodyact);free(vr);free(vi);if(*nbody>0) free(ipobody);
      
      if(*ithermal==1) free(t1act);
      
      if(iprescribedboundary){
	  free(xbounr);free(xbouni);free(fr);free(fi);free(ubr);free(ubi);
	  free(mubr);free(mubi);
      }
  }

  else{

      /* steady state response to a nonharmonic periodic loading */
      
      xforcact=NNEW(double,nfour**nforc);
      xloadact=NNEW(double,nfour*2**nload);
      xbodyact=NNEW(double,nfour*7**nbody);
      xbounact=NNEW(double,nfour**nboun);
      xbounacttime=NNEW(double,nfour**nboun);
      if(*ithermal==1) t1act=NNEW(double,*nk);

      r=NNEW(double,nfour);
      wsave=NNEW(double,2*nfour);
      isave=NNEW(int,15);
      
      /*  check for nonzero SPC's */
      
      iprescribedboundary=0;
      for(i=0;i<*nboun;i++){
	  if(fabs(xboun[i])>1.e-10){
	      iprescribedboundary=1;
	      break;
	  }
      }

      /* determining the load time history */
      
      ampli=NNEW(double,*nam); /* instantaneous amplitude */

      for(l=0;l<nfour;l++){

	  time=tmin+(tmax-tmin)*(double)l/(double)nfour;
	      
	  FORTRAN(tempload,(xforcold,xforc,&xforcact[l**nforc],iamforc,nforc,
	    xloadold,xload,&xloadact[l*2**nload],iamload,nload,ibody,xbody,
	    nbody,xbodyold,&xbodyact[l*7**nbody],t1old,t1,t1act,
	    iamt1,nk,amta,namta,nam,ampli,&time,
	    &reltime,ttime,&dtime,ithermal,nmethod,
	    xbounold,xboun,&xbounact[l**nboun],iamboun,nboun,
	    nodeboun,ndirboun,nodeforc,
	    ndirforc,istep,&iinc,co,vold,itg,&ntg));
	  
      }

      free(ampli);

      for(i=0;i<l**nboun;i++){xbounacttime[i]=xbounact[i];}

      /* determining the load frequency history:
         frequency transform of the load time history */

      for(i=0;i<*nforc;i++){
	  for(l=0;l<nfour;l++){
	      r[l]=xforcact[l**nforc+i];
	  }
	  FORTRAN(drffti,(&nfour,wsave,isave));
	  FORTRAN(drfftf,(&nfour,r,wsave,isave));
	  for(l=0;l<nfour;l++){
	      xforcact[l**nforc+i]=r[l]/nfour*2.;
	  }
	  xforcact[i]=xforcact[i]/2.;
      }

      for(i=0;i<*nload;i++){
	  for(l=0;l<nfour;l++){
	      r[l]=xloadact[l*2**nload+2*i];
	  }
	  FORTRAN(drffti,(&nfour,wsave,isave));
	  FORTRAN(drfftf,(&nfour,r,wsave,isave));
	  for(l=0;l<nfour;l++){
	      xloadact[l*2**nload+2*i]=r[l]/nfour*2.;
	  }
	  xloadact[2*i]=xloadact[2*i]/2.;
      }

      for(i=0;i<*nbody;i++){
	  for(l=0;l<nfour;l++){
	      r[l]=xbodyact[l**nbody+7*i];
	  }
	  FORTRAN(drffti,(&nfour,wsave,isave));
	  FORTRAN(drfftf,(&nfour,r,wsave,isave));
	  for(l=0;l<nfour;l++){
	      xbodyact[l**nbody+7*i]=r[l]/nfour*2.;
	  }
	  xbodyact[7*i]=xbodyact[7*i]/2.;
      }

      if(iprescribedboundary){
	  for(i=0;i<*nboun;i++){
	      for(l=0;l<nfour;l++){
		  r[l]=xbounact[l**nboun+i];
	      }
	      FORTRAN(drffti,(&nfour,wsave,isave));
	      FORTRAN(drfftf,(&nfour,r,wsave,isave));
	      for(l=0;l<nfour;l++){
		  xbounact[l**nboun+i]=r[l]/nfour*2.;
	      }
	      xbounact[i]=xbounact[i]/2.;
	  }
	  
	  /* LU decomposition of the stiffness matrix */
	  
	  if(*isolver==0){
#ifdef SPOOLES
	      spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
                        &symmetryflag,&inputformat);
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

      }

      free(r);free(wsave);free(isave);

      /* determining the frequency data points */
      
      freqnh=NNEW(double,ndata*(nev+1));
      
      ndatatot=0.;
      freqnh[0]=fmin;
      if(fabs(fmax-fmin)<1.e-10){
	  ndatatot=1;
      }else{
	  for(i=0;i<nev;i++){
	      if(d[i]>=fmin){
		  if(d[i]<=fmax){
		      for(j=1;j<ndata;j++){
			  y=-1.+2.*j/((double)(ndata-1));
			  if(fabs(y)<1.e-10){freqnh[ndatatot+j]=
                                 (freqnh[ndatatot]+d[i])/2.;}
			  else{
			      freqnh[ndatatot+j]=(freqnh[ndatatot]+d[i])/2.+
				  (d[i]-freqnh[ndatatot])*pow(fabs(y),1./bias)
                                  *y/(2.*fabs(y));
			  }
		      }
		      ndatatot+=ndata-1;
		  }
		  else{break;}
	      }
	  }
	  for(j=1;j<ndata;j++){
	      y=-1.+2.*j/((double)(ndata-1));
	      if(fabs(y)<1.e-10){freqnh[ndatatot+j]=
                         (freqnh[ndatatot]+fmax)/2.;}
	      else{
		  freqnh[ndatatot+j]=(freqnh[ndatatot]+fmax)/2.+
		      (fmax-freqnh[ndatatot])*pow(fabs(y),1./bias)*
                      y/(2.*fabs(y));
	      }
	  }
	  ndatatot+=ndata;
      }
      RENEW(freqnh,double,ndatatot);

      for(ii=0;ii<ndatatot;ii++){
	  for(i=0;i<6**mint_**ne;i++){eme[i]=0.;}

	  sprintf(description,"%12f",freqnh[ii]/(2.*pi));
	  
	  xforcr=NNEW(double,*nforc);
	  xloadr=NNEW(double,2**nload);
	  xbodyr=NNEW(double,7**nbody);
	  for(k=0;k<7**nbody;k++){xbodyr[k]=xbody[k];}
	  if(iprescribedboundary) xbounr=NNEW(double,*nboun);
	  
	  /* assigning the body forces to the elements */ 

	  if(*nbody>0){
	      ifreebody=*ne+1;
	      ipobody=NNEW(int,2*ifreebody);
	      for(k=1;k<=*nbody;k++){
		  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			    iendset,ialset,&inewton,nset,&ifreebody,&k));
		  RENEW(ipobody,int,2*(*ne+ifreebody));
	      }
	      RENEW(ipobody,int,2*(ifreebody-1));
	  }
	  
	  br=NNEW(double,neq[1]); /* load rhs vector (real part) */
	  bi=NNEW(double,neq[1]); /* load rhs vector (imaginary part) */
	  fr=NNEW(double,neq[1]); /* force corresponding to real particular solution */
	  ubr=NNEW(double,neq[1]); /* real particular solution */
	  mubr=NNEW(double,neq[1]); /* mass times real particular solution */
	  btot=NNEW(double,nfour*neq[1]);
	  
	  bjr=NNEW(double,nev); /* real response modal decomposition */
	  bji=NNEW(double,nev); /* imaginary response modal decomposition */
	  
	  aa=NNEW(double,nev); /* modal coefficients of the real loading */
	  bb=NNEW(double,nev); /* modal coefficients of the imaginary loading */
	  
	  /* loop over all Fourier frequencies */
	  
	  freq=NNEW(double,nfour);
	  
	  for(l=0;l<nfour;l++){
/*      for(l=1;l<2;l++){ */
	      
	      /* frequency */
	      
	      freq[l]=freqnh[ii]*floor((l+1.)/2.+0.1);
	      
	      /* loading for this frequency */
	      
	      for(i=0;i<*nforc;i++){
		  xforcr[i]=xforcact[l**nforc+i];
	      }
	      
	      for(i=0;i<*nload;i++){
		  xloadr[2*i]=xloadact[l*2**nload+2*i];
	      }
	      
	      for(i=0;i<*nbody;i++){
		  xbodyr[7*i]=xbodyact[l**nbody+7*i];
	      }
	      
	      /* calculating the instantaneous loading vector */
	      
	      FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcr,
		nforc,nelemload,sideload,xloadr,nload,xbodyr,
		ipobody,nbody,cgr,br,nactdof,&neq[1],nmethod,
		ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		nplicon,plkcon,nplkcon,
		npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody));
	      
	      for(i=0;i<neq[1];i++){bi[i]=0.;}
	      
	      if(iprescribedboundary){
		  
		  for(i=0;i<*nboun;i++){
		      xbounr[i]=xbounact[l**nboun+i];
		  }
		  
		  for(j=0;j<neq[1];j++){fr[j]=0.;ubr[j]=0.;}
		  for(i=0;i<*nboun;i++){
		      ic=neq[1]+i;
		      for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
			  ir=irow[j]-1;
			  fr[ir]=fr[ir]-au[j]*xbounr[i];
			  ubr[ir]=fr[ir];
		      }
		  }
		  if(*isolver==0){
#ifdef SPOOLES
		      spooles_solve(ubr,&neq[1]);
#endif
		  }
		  else if(*isolver==4){
#ifdef SGI
		      sgi_solve(ubr,token);
#endif
		  }
		  else if(*isolver==5){
#ifdef TAUCS
		      tau_solve(ubr,&neq[1]);
#endif
		  }
		  FORTRAN(op,(&neq[1],aux,ubr,mubr,adb,aub,
			      icol,irow,nzl));
		  
		  for(i=0;i<neq[1];i++){
		      br[i]+=freq[l]*(freq[l]*mubr[i]);
		      bi[i]+=freq[l]*(-alpham*mubr[i]-betam*fr[i]);
		  }
		  
	      }
	      
	      /* real and imaginary modal coefficients */
	      
	      for(i=0;i<nev;i++){
		  aa[i]=0.;
		  for(j=0;j<neq[1];j++){
		      aa[i]+=z[i*neq[1]+j]*br[j];
		  }
	      }
	      
	      for(i=0;i<nev;i++){
		  bb[i]=0.;
		  for(j=0;j<neq[1];j++){
		      bb[i]+=z[i*neq[1]+j]*bi[j];
		  }
	      }
	      
	      /* calculating the modal coefficients */
	      
	      if(dashpot==0){
		  for(i=0;i<nev;i++){
		      dd=pow(pow(d[i],2)-pow(freq[l],2),2)+
			  pow(alpham+betam*pow(d[i],2),2)*pow(freq[l],2);
		      bjr[i]=(aa[i]*(d[i]*d[i]-freq[l]*freq[l])+
			      bb[i]*(alpham+betam*d[i]*d[i])*freq[l])/dd;
		      bji[i]=(bb[i]*(d[i]*d[i]-freq[l]*freq[l])-
			      aa[i]*(alpham+betam*d[i]*d[i])*freq[l])/dd;
		  }
	      }else{
		  for(i=0;i<nev2;i++){
		      for(j=0;j<nev2;j++){
			  am[i*nev2+j]=0.;
		      }
		      bm[i]=0.;
		  }
		  for(i=0;i<nev;i++){
		      am[i*nev2+i]=d[i]*d[i]-freq[l]*freq[l];
		      am[(i+nev)*nev2+i]=-(alpham+betam*d[i]*d[i])*freq[l];
		      bm[i]=aa[i];
		      am[i*nev2+nev+i]=-am[(i+nev)*nev2+i];
		      am[(i+nev)*nev2+nev+i]=am[i*nev2+i];
		      bm[nev+i]=bb[i];
		      for(j=0;j<nev;j++){
			  am[(j+nev)*nev2+i]=am[(j+nev)*nev2+i]
			      -cc[i*nev+j]*freq[l];
			  am[j*nev2+nev+i]=am[j*nev2+nev+i]
			      +cc[i*nev+j]*freq[l];
		      }
		  }
		  
		  /* solving the system of equations */
		  
		  FORTRAN(dgesv,(&nev2,&nrhs,am,&nev2,ipiv,bm,&nev2,&info));
		  if(info!=0){
		      printf("*ERROR in steadystate: fatal termination of dgesv\n");
		      printf("       info=%d\n",info);
/*		  FORTRAN(stop,());*/
		  }
		  
		  /* storing the solution in bjr and bji */
		  
		  for(i=0;i<nev;i++){
		      bjr[i]=bm[i];
		      bji[i]=bm[nev+i];
		  }
		  
	      }
	      
	      /* calculating the real response */
	      
	      if(iprescribedboundary){
		  for(i=0;i<neq[1];i++){
		      br[i]=ubr[i];
		  }
	      }
	      else{
		  for(i=0;i<neq[1];i++){
		      br[i]=0.;
		  }
	      }
	      for(i=0;i<neq[1];i++){
		  for(j=0;j<nev;j++){
		      br[i]+=bjr[j]*z[j*neq[1]+i];
		  }
	      }
	      
	      /* calculating the imaginary response */
	      
	      for(i=0;i<neq[1];i++){
		  bi[i]=0.;
	      }
	      for(i=0;i<neq[1];i++){
		  for(j=0;j<nev;j++){
		      bi[i]+=bji[j]*z[j*neq[1]+i];
		  }
	      }
	      
	      /* magnitude and phase of the response */
	      
	      for(i=0;i<neq[1];i++){
		  breal=br[i];
		  br[i]=sqrt(br[i]*br[i]+bi[i]*bi[i]);
		  if(fabs(breal)<1.e-10){
		      if(bi[i]>0.){bi[i]=pi/2.;}
		      else{bi[i]=-pi/2.;}
		  }
		  else{
		      bi[i]=atan(bi[i]/breal);
		      if(breal<0.){bi[i]+=pi;}
		  }
	      }
	      
	      /* correction for the sinus terms */
	      
	      if((l!=0)&&(2*(int)floor(l/2.+0.1)==l)){
		  for(i=0;i<neq[1];i++){
		      bi[i]-=pi/2.;}
	      }
	      
	      /* contribution to the time response */
	      
	      for(j=0;j<nfour;j++){
		  time=tmin+2.*pi/freqnh[ii]*(double)j/(double)nfour;
		  for(i=0;i<neq[1];i++){
		      btot[j*neq[1]+i]+=br[i]*cos(freq[l]*time+bi[i]);
		  }
	      }
	      
	  }
	  
	  free(mubr);
	  free(xforcr);free(xloadr);free(xbodyr);free(br);free(bi);free(freq);
	  free(fr);free(ubr);free(bjr);free(bji);free(aa);free(bb);
          if(*nbody>0) free(ipobody);
	  if(iprescribedboundary) free(xbounr);
	  
	  /* result fields */
	  
	  vr=NNEW(double,4**nk);
	  fn=NNEW(double,4**nk);
	  stn=NNEW(double,6**nk);
	  inum=NNEW(int,*nk);
	  stx=NNEW(double,6**mint_**ne);
	  if(*ithermal>1) {qfn=NNEW(double,3**nk);qfx=NNEW(double,3**mint_**ne);}
	  
	  if(strcmp1(&filab[18],"E   ")==0) een=NNEW(double,6**nk);
	  if(strcmp1(&filab[36],"ENER")==0) enern=NNEW(double,*nk);
	  
	  eei=NNEW(double,6**mint_**ne);
	  if(*nener==1){
	      stiini=NNEW(double,6**mint_**ne);
	      enerini=NNEW(double,*mint_**ne);}
	  
	  /* storing the results */
	  
	  for(l=0;l<nfour;l++){
	      time=tmin+2.*pi/freqnh[ii]*(double)l/(double)nfour;
	      *ttime=time;
	      
	      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,
		  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		  ithermal,prestr,iprestr,filab,eme,een,
		  iperturb,f,fn,nactdof,&iout,qa,
		  vold,&btot[l*neq[1]],nodeboun,ndirboun,&xbounacttime[l**nboun],nboun,
		  ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		  veold,accold,&bet,&gam,&dtime,&time,ttime,
		  plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
		  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		  enern,sti,xstaten,eei,enerini,cocon,ncocon,
		  set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		  qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));
	  
	      (*kode)++;
	      
    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
	      FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,nmethod,kode,
		   filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
                   nstate_,istep,&iinc,iperturb,ener,mint_,output,ithermal,
                   qfn,&mode,&noddiam,trab,inotr,ntrans,orab,ielorien,norien,
                   description,
                     ipneigh,neigh,sti));free(ipneigh);free(neigh);

	  }
	  
	  free(vr);free(fn);free(stn);free(inum);free(stx);free(eei);free(btot);
	  if(*ithermal>1) {free(qfn);free(qfx);}
	  
	  if(strcmp1(&filab[18],"E   ")==0) free(een);
	  if(strcmp1(&filab[36],"ENER")==0) free(enern);
	  
	  if(*nener==1){free(stiini);free(enerini);}
	  
      }
      free(xforcact);free(xloadact);free(xbodyact);free(xbounact);
      free(xbounacttime);free(freqnh);
      if(*ithermal==1) free(t1act);

  }

  free(adb);free(aub);free(z);free(d);

  if(*mcs==0){
      free(ad);free(au);
  }else{

      *nk/=nsectors;
      *ne/=nsectors;
      *nboun/=nsectors;
      neq[1]=neq[1]*2/nsectors;

      RENEW(co,double,3**nk);
      if(*ithermal!=0){
	  RENEW(t0,double,*nk);
	  RENEW(t1old,double,*nk);
	  RENEW(t1,double,*nk);
	  if(*nam>0) RENEW(iamt1,int,*nk);
      }
      RENEW(nactdof,int,4**nk);
      if(*ntrans>0) RENEW(inotr,int,2**nk);
      RENEW(kon,int,*nkon);
      RENEW(ipkon,int,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(ielmat,int,*ne);
      if(*norien>0) RENEW(ielorien,int,*ne);
      RENEW(nodeboun,int,*nboun);
      RENEW(ndirboun,int,*nboun);
      if(*nam>0) RENEW(iamboun,int,*nboun);
      RENEW(xboun,double,*nboun);
      RENEW(xbounold,double,*nboun);
      RENEW(ikboun,int,*nboun);
      RENEW(ilboun,int,*nboun);

      /* recovering the original multiple point constraints */

      RENEW(ipompc,int,*nmpc);
      RENEW(nodempc,int,3**mpcend);
      RENEW(coefmpc,double,*mpcend);
      RENEW(labmpc,char,20**nmpc);
      RENEW(ikmpc,int,*nmpc);
      RENEW(ilmpc,int,*nmpc);
      RENEW(fmpc,double,*nmpc);

      *nmpc=nmpcold;
      *mpcend=mpcendold;
      for(i=0;i<*nmpc;i++){ipompc[i]=ipompcold[i];}
      for(i=0;i<3**mpcend;i++){nodempc[i]=nodempcold[i];}
      for(i=0;i<*mpcend;i++){coefmpc[i]=coefmpcold[i];}
      for(i=0;i<20**nmpc;i++){labmpc[i]=labmpcold[i];}
      for(i=0;i<*nmpc;i++){ikmpc[i]=ikmpcold[i];}
      for(i=0;i<*nmpc;i++){ilmpc[i]=ilmpcold[i];}
      free(ipompcold);free(nodempcold);free(coefmpcold);
      free(labmpcold);free(ikmpcold);free(ilmpcold);

      RENEW(vold,double,4**nk);
      RENEW(eme,double,6**mint_**ne);
      if(*nener==1)RENEW(ener,double,*mint_**ne);

/* distributed loads */

      for(i=0;i<*nload;i++){
	  if(nelemload[2*i]<=*ne*nsectors){
	      nelemload[2*i]-=*ne*nelemload[2*i+1];
	  }else{
	      nelemload[2*i]-=*ne*(nsectors+nelemload[2*i+1]-1);
	  }
      }

  /*  sorting the elements with distributed loads */

      if(*nload>0){
	  if(*nam>0){
	      FORTRAN(isortiddc2,(nelemload,iamload,xload,xloadold,sideload,nload,&kflag));
	  }else{
	      FORTRAN(isortiddc1,(nelemload,xload,xloadold,sideload,nload,&kflag));
	  }
      }
      
/* point loads */
      
      for(i=0;i<*nforc;i++){
	  if(nodeforc[2*i]<=*nk*nsectors){
	      nodeforc[2*i]-=*nk*nodeforc[2*i+1];
	  }else{
	      nodeforc[2*i]-=*nk*(nsectors+nodeforc[2*i+1]-1);
	  }
      }
  }

  free(xstiff);

  if(dashpot){free(adc);free(auc);free(cc);free(am);free(bm);free(ipiv);}


  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;

  return;
}
