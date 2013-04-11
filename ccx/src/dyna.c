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

void dyna(double **cop, int *nk, int **konp, int **ipkonp, char **lakonp, int *ne, 
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
	       double **veoldp, char *amname, double *amta,
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
               int *isolver,int *jq, char *output, int *mcs, int *nkon,
               int *mpcend, int *ics, double *cs){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *labmpcold=NULL;

  int nev,i,j,k,m,idof,jinc, *inum=NULL,*ipobody=NULL,inewton=0,
    iinc,l,iout,ielas,icmd,iprescribedboundary,init,ifreebody,
    mode=-1,noddiam=-1,*kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,
    *inotr=NULL,*nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,nsectors,nmpcold,mpcendold,*ipompcold=NULL,*nodempcold=NULL,
    *ikmpcold=NULL,*ilmpcold=NULL,kflag=2,nmd,nevd,*nm=NULL,*iamt1=NULL,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,dashpot,lrw,liw,iddebdf=0,
    *iwork=NULL;

  double *d=NULL, *z=NULL, *b=NULL, *zeta=NULL,*stiini=NULL,*vini=NULL,
    *cd=NULL, *cv=NULL, *xforcact=NULL, *xloadact=NULL,*cc=NULL,
    *t1act=NULL, *ampli=NULL, *aa=NULL, *bb=NULL, *bj=NULL, *v=NULL,
    *stn=NULL, *stx=NULL, *een=NULL, *adb=NULL,*xstiff=NULL,
    *aub=NULL, *temp_array1=NULL, *temp_array2=NULL, *aux=NULL,
    *f=NULL, *fn=NULL, *xbounact=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL, *xbodyact=NULL, *cgr=NULL, *au=NULL, *vbounact=NULL,
    *abounact=NULL,time,dtime,reltime,*t0=NULL,*t1=NULL,*t1old=NULL,
    physcon[1],zetaj,dj,ddj,h1,h2,h3,h4,h5,h6,sum,aai,bbi,tstart,tend,
    qa[2],cam[2],accold[1],bet,gam,*ad=NULL,sigma=0.,alpham,betam,
    *bact=NULL,*bmin=NULL,*co=NULL,*xboun=NULL,*xbounold=NULL,*vold=NULL,
    *eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,*coefmpcold,*veold=NULL,
    *xini=NULL,*rwork=NULL,*adc=NULL,*auc=NULL,*zc=NULL, *rpar=NULL;

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
  fmpc=*fmpcp;veold=*veoldp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;

  xstiff=NNEW(double,27**mint_**ne);

  alpham=xmodal[0];
  betam=xmodal[1];

  /* determining nzl */

  *nzl=0;
  for(i=neq[1];i>0;i--){
      if(icol[i-1]>0){
	  *nzl=i;
	  break;
      }
  }

  /* reading the eigenvalue and eigenmode information */

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
      free(veold);veold=NNEW(double,4**nk);
      RENEW(eme,double,6**mint_**ne);
      RENEW(xstiff,double,27**mint_**ne);
      if(*nener==1) RENEW(ener,double,*mint_**ne);
  }

  fclose(f1);
	  
/*  jinc=(*tper/(*tinc)/(*jout))**jout;
    if(jinc**tinc<*tper) jinc++;*/
  jinc=(*tper/(*tinc)+0.5);
  if((jinc/(*jout))*(*jout)!=jinc) jinc+=(*jout);

  /* check whether there are dashpot elements */

  dashpot=0;
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[i*8],"ED")==0){
	  dashpot=1;break;}
  }
  if(dashpot){
      liw=51;
      iwork=NNEW(int,liw);
      lrw=130+42*nev;
      rwork=NNEW(double,lrw);
      xini=NNEW(double,2*nev);
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
      rpar=NNEW(double,3+nev*(3+2*jinc+nev));
  }

  zeta=NNEW(double,nev);

  /* calculating the damping coefficients*/

  for(i=0;i<nev;i++){
    if(fabs(d[i])>(1.e-10)){
      zeta[i]=(alpham+betam*d[i]*d[i])/(2.*d[i]);
    }
    else {
      printf("*ERROR in dyna: one of the frequencies is zero");
      exit(0);
    }
  }

  /* modal decomposition of the initial conditions */

  cd=NNEW(double,nev);
  cv=NNEW(double,nev);

  for(i=0;i<nev;i++){cd[i]=0;cv[i]=0;}

  temp_array1=NNEW(double,neq[1]);
  temp_array2=NNEW(double,neq[1]);
  for(i=0;i<neq[1];i++){temp_array1[i]=0;temp_array2[i]=0;}

  /* displacement initial conditions */

  for(i=0;i<*nk;i++){
    for(j=0;j<4;j++){
      if(nactdof[4*i+j]!=0){
	idof=nactdof[4*i+j]-1;
	temp_array1[idof]=vold[4*i+j];
      }
    }
  }

  FORTRAN(op,(&neq[1],aux,temp_array1,temp_array2,adb,aub,icol,irow,nzl));

  for(i=0;i<neq[1];i++){
    for(k=0;k<nev;k++){
      cd[k]+=z[k*neq[1]+i]*temp_array2[i];
    }
  }

  /* velocity initial conditions */

  for(i=0;i<neq[1];i++){temp_array1[i]=0;temp_array2[i]=0;}
  for(i=0;i<*nk;i++){
    for(j=0;j<4;j++){
      if(nactdof[4*i+j]!=0){
	idof=nactdof[4*i+j]-1;
	temp_array1[idof]=veold[4*i+j];
      }
    }
  }
  
  FORTRAN(op,(&neq[1],aux,temp_array1,temp_array2,adb,aub,icol,irow,nzl));
  
  for(i=0;i<neq[1];i++){
    for(k=0;k<nev;k++){
      cv[k]+=z[k*neq[1]+i]*temp_array2[i];
    }
  }

  free(temp_array1);free(temp_array2);

  xforcact=NNEW(double,*nforc);
  xloadact=NNEW(double,2**nload);
  xbodyact=NNEW(double,7**nbody);
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
	       
  b=NNEW(double,neq[1]); /* load rhs vector */
  bj=NNEW(double,nev); /* response modal decomposition */
  ampli=NNEW(double,*nam); /* instantaneous amplitude */

  /* constant coefficient of the linear amplitude function */
  aa=NNEW(double,(jinc+1)*nev); 
  /* linear coefficient of the linear amplitude function */
  bb=NNEW(double,(jinc+1)*nev);
  
  v=NNEW(double,4**nk);
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

  time=0.;
  dtime=*tinc;

/*     calculating the instantaneous loads (forces, surface loading, 
       centrifugal and gravity loading or temperature) at time 0 */

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
     xload,xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
     xbodyact,t1old,t1,t1act,iamt1,nk,
     amta,namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
     xbounold,xboun,xbounact,iamboun,nboun,
     nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
     co,vold,itg,&ntg));

  /*  calculating the instantaneous loading vector at time 0 */

  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
       nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,nbody,
       cgr,b,nactdof,&neq[1],nmethod,
       ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,alcon,
       nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,t0,t1act,
       ithermal,iprestr,vold,iperturb,iexpl,plicon,
       nplicon,plkcon,nplkcon,npmat_,ttime,&time,istep,&iinc,&dtime,
       physcon,ibody));

  /*  check for nonzero SPC's */

  iprescribedboundary=0;
  for(i=0;i<*nboun;i++){
      if(fabs(xboun[i])>1.e-10){
	  iprescribedboundary=1;
	  break;
      }
  }

  /*  correction for nonzero SPC's */

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

      bact=NNEW(double,neq[1]);
      bmin=NNEW(double,neq[1]);

      init=1;
      dynboun(amta,namta,nam,ampli,&time,ttime,&dtime,xbounold,xboun,
	      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
	      aub,icol,irow,neq,nzs,&sigma,b,isolver,
	      &alpham,&betam,nzl,&init,bact,bmin,jq);
      init=0;
  }

  for(i=0;i<nev;i++){
      aa[i]=0.;
      for(j=0;j<neq[1];j++){
	  aa[i]+=z[i*neq[1]+j]*b[j];
      }
  }

  for(k=0;k<=jinc-*jout;k+=*jout){
      for(l=1;l<=*jout;l++){
	  m=k+l;
	  time=m**tinc;

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

	  /* calculating the instantaneous loading vector */

          FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
                ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
                nforc,nelemload,sideload,xloadact,nload,xbodyact,
                ipobody,nbody,cgr,b,nactdof,&neq[1],nmethod,
                ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
                alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
                t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
                nplicon,plkcon,nplkcon,
		npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody));

	  /* correction for nonzero SPC's */

	  if(iprescribedboundary){
	      dynboun(amta,namta,nam,ampli,&time,ttime,&dtime,xbounold,xboun,
		      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		      aub,icol,irow,neq,nzs,&sigma,b,isolver,
		      &alpham,&betam,nzl,&init,bact,bmin,jq);
	  }

	  for(i=0;i<nev;i++){
	      aa[m*nev+i]=0.;
	      for(j=0;j<neq[1];j++){
		  aa[m*nev+i]+=z[i*neq[1]+j]*b[j];
	      }
	      bb[(m-1)*nev+i]=(aa[m*nev+i]-aa[(m-1)*nev+i])/(*tinc);
	      aa[(m-1)*nev+i]=aa[m*nev+i]-bb[(m-1)*nev+i]*time;
	  }
      }

      if(dashpot){
	  FORTRAN(subspace,(d,aa,bb,cc,tinc,&alpham,&betam,&nev,xini,
			    cd,cv,&time,rwork,&lrw,&k,jout,rpar,bj,
                            iwork,&liw,&iddebdf));
	  if(iddebdf==2){
	      liw=56+2*nev;
	      RENEW(iwork,int,liw);
	      for(i=0;i<liw;i++){iwork[i]=0;}
	      lrw=250+20*nev+4*nev*nev;
	      RENEW(rwork,double,lrw);
	      for(i=0;i<lrw;i++){rwork[i]=0.;}
	      iddebdf=1;
	      FORTRAN(subspace,(d,aa,bb,cc,tinc,&alpham,&betam,&nev,xini,
			    cd,cv,&time,rwork,&lrw,&k,jout,rpar,bj,
                            iwork,&liw,&iddebdf));
	  }
	      
/*	  printf("bjsubspace=%le\n",bj[0]);
	  printf("bjsubspace=%le\n",bj[1]);
	  printf("bjsubspace=%le\n",bj[2]);
	  printf("bjsubspace=%le\n",bj[3]);
	  printf("bjsubspace=%le\n",bj[4]);*/
      }
      else{
	  for(l=0;l<nev;l++){
	      zetaj=zeta[l];
	      dj=d[l];
	      
	      /*   subcritical damping */
	      
	      if(zetaj<1.-1.e-6){
		  ddj=dj*sqrt(1.-zetaj*zetaj);
		  h1=zetaj*dj;
		  h2=h1*h1+ddj*ddj;
		  h3=h1*h1-ddj*ddj;
		  h4=2.*h1*ddj/h2;
		  sum=0.;
		  for(i=0;i<m;i++){
		      aai=aa[i*nev+l];
		      bbi=bb[i*nev+l];
		      tstart=time-(i+1)**tinc;
		      tend=time-i**tinc;
		      sum+=FORTRAN(fsub,(&time,&tend,&aai,&bbi,&ddj,
					 &h1,&h2,&h3,&h4));
		      sum-=FORTRAN(fsub,(&time,&tstart,&aai,&bbi,&ddj,
					 &h1,&h2,&h3,&h4));
		  }
		  bj[l]=sum/ddj+exp(-h1*time)*(cos(ddj*time)+zetaj/
		      sqrt(1.-zetaj*zetaj)*sin(ddj*time))*cd[l]+
		      exp(-h1*time)*sin(ddj*time)*cv[l]/ddj;
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
		  sum=0.;
		  for(i=0;i<m;i++){
		      aai=aa[i*nev+l];
		      bbi=bb[i*nev+l];
		      tstart=time-(i+1)**tinc;
		      tend=time-i**tinc;
		      sum+=FORTRAN(fsuper,(&time,&tend,&aai,&bbi,
					   &h1,&h2,&h3,&h4,&h5,&h6));
		      sum-=FORTRAN(fsuper,(&time,&tstart,&aai,&bbi,
					   &h1,&h2,&h3,&h4,&h5,&h6));
		  }
		  bj[l]=sum/(2.*ddj)+(exp(h1*time)+exp(-h2*time))*cd[l]
		      /2.+zetaj*(exp(h1*time)-exp(-h2*time))/(2.*
		      sqrt(zetaj*zetaj-1.))*cd[l]+(exp(h1*time)-
		      exp(-h2*time))*cv[l]/(2.*ddj);
	      }
	      
	      /* critical damping */
	      
	      else{
		  h1=zetaj*dj;
		  h2=1./h1;
		  h3=h2*h2;
		  h4=h2*h3;
		  sum=0.;
		  for(i=0;i<m;i++){
		      aai=aa[i*nev+l];
		      bbi=bb[i*nev+l];
		      tstart=time-(i+1)**tinc;
		      tend=time-i**tinc;
		      sum+=FORTRAN(fcrit,(&time,&tend,&aai,&bbi,&zetaj,&dj,
					  &ddj,&h1,&h2,&h3,&h4));
		      sum-=FORTRAN(fcrit,(&time,&tstart,&aai,&bbi,&zetaj,&dj,
					  &ddj,&h1,&h2,&h3,&h4));
		  }
		  bj[l]=sum+exp(-h1*time)*
		      ((1.+h1*time)*cd[l]+time*cv[l]);
	      }
	  }
/*	  printf("bj=%le\n",bj[0]);
	  printf("bj=%le\n",bj[1]);
	  printf("bj=%le\n",bj[2]);
	  printf("bj=%le\n",bj[3]);
	  printf("bj=%le\n",bj[4]);*/
      }

      /* composing the response */
       
      if(iprescribedboundary){
	  for(i=0;i<neq[1];i++){
	      b[i]=bmin[i];
	  }
      }
      else{
	  for(i=0;i<neq[1];i++){
	      b[i]=0.;
	  }
      }
 
/*	  for(i=0;i<neq[1];i++){
	      printf("particular solution: i=%d,b[i]=%f\n",i,b[i]);
	      }*/
      
      for(i=0;i<neq[1];i++){
	  for(j=0;j<nev;j++){
	      b[i]+=bj[j]*z[j*neq[1]+i];
	  }
      }
 
/*	  for(i=0;i<neq[1];i++){
	      printf("total solution: i=%d,b[i]=%f\n",i,b[i]);
	      }*/

      /*  storing the displacement field of the one but last
        calculation: needed for the velocity field at the start
        of the next step */

      if(k+*jout>jinc-*jout){
	  for(i=0;i<4**nk;i++){
	      vold[i]=v[i];
	  }
      }

      iout=1;

      *ttime+=*jout**tinc;

      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
             stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
             ielmat,ielorien,norien,orab,ntmat_,t0,t1,
             ithermal,prestr,iprestr,filab,eme,een,
             iperturb,f,fn,nactdof,&iout,qa,
             vold,b,nodeboun,ndirboun,xbounact,nboun,
             ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
             veold,accold,&bet,&gam,&dtime,&time,ttime,
             plicon,nplicon,plkcon,nplkcon,
             xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
             &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
             enern,sti,xstaten,eei,enerini,cocon,ncocon,
             set,nset,istartset,iendset,ialset,nprint,prlab,prset,
	     qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload));

      if((dashpot)||(k+*jout>jinc-*jout)){
	  for(i=0;i<4**nk;i++){
	      veold[i]=(v[i]-vold[i])/(*jout*dtime);
	      vold[i]=v[i];
	  }
      }

      (*kode)++;

    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
      FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,
          een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,nstate_,istep,&iinc,
	  iperturb,ener,mint_,output,ithermal,qfn,&mode,&noddiam,
          trab,inotr,ntrans,orab,ielorien,norien,description,
                     ipneigh,neigh,sti));free(ipneigh);free(neigh);

  }

  /* determining the initial displacement and velocity field for
     the next step */

  for(i=0;i<4**nk;i++){
      veold[i]=(v[i]-vold[i])/dtime;
      vold[i]=v[i];
  }

  free(eei);free(vbounact);free(abounact);
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
  
  free(v);free(fn);free(stn);free(inum);free(stx);free(adb);
  free(aub);free(z);free(b);free(zeta);free(bj);free(cd);free(cv);
  free(xforcact);free(xloadact);free(xbounact);free(aa);free(bb);
  free(ampli);free(xbodyact);

  if(*ithermal==1) free(t1act);

  if(iprescribedboundary){free(bact);free(bmin);}

  if(*mcs==0){
      free(ad);free(au);
  }else{

      *nk/=nsectors;
      *ne/=nsectors;
      *nboun/=nsectors;
      neq[1]=neq[1]*2/nsectors;

      RENEW(co,double,3**nk);
      if((*ithermal!=0)&&(*nam>0)) RENEW(iamt1,int,*nk);
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
      RENEW(veold,double,4**nk);
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

  free(xstiff);if(*nbody>0) free(ipobody);

  if(dashpot){
      free(xini);free(rwork);free(adc);free(auc);free(cc);
      free(rpar);free(iwork);}

  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*veoldp=veold;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;

  return;
}
