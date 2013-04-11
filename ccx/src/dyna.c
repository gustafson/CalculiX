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
#ifdef PARDISO
   #include "pardiso.h"
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
	       double **voldp,int *iperturb, double **stip, int *nzs, 
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
               int *iendset, int **ialsetp, int *nprint, char *prlab,
               char *prset, int *nener, double *trab, 
               int **inotrp, int *ntrans, double **fmpcp, char *cbody, int *ibody,
               double *xbody, int *nbody, double *xbodyold, int *istep,
               int *isolver,int *jq, char *output, int *mcs, int *nkon,
               int *mpcend, int *ics, double *cs, int *ntie, char *tieset,
               int *idrct, int *jmax, double *tmin, double *tmax,
	       double *ctrl, int *itpamp, double *tietol,int *nalset,
               int **nnnp){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *labmpcold=NULL,lakonl[9]="        \0",*tchar1=NULL,*tchar2=NULL,*tchar3;

  int nev,i,j,k,idof,*inum=NULL,*ipobody=NULL,inewton=0,
    iinc=0,jprint=0,l,iout,ielas,icmd,iprescribedboundary,init,ifreebody,
    mode=-1,noddiam=-1,*kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,
    *inotr=NULL,*nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,nsectors,nmpcold,mpcendold,*ipompcold=NULL,*nodempcold=NULL,
    *ikmpcold=NULL,*ilmpcold=NULL,kflag=2,nmd,nevd,*nm=NULL,*iamt1=NULL,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,dashpot,lrw,liw,iddebdf=0,
    *iwork=NULL,ngraph,nkg,neg,ncont,ncone,ne0,nkon0, *itietri=NULL,
    *koncont=NULL,konl[20],imat,nope,kodem,indexe,j1,ist,jdof,index,
    id,node,ndir,*ipneigh=NULL,*neigh=NULL,niter,inext,itp=0,icutb=0,
    ismallsliding=0,isteadystate,*ifcont1=NULL,*ifcont2=NULL,mpcfree,
    memmpc_,imax,iener=0,*icole=NULL,*irowe=NULL,*jqe=NULL,nzse[3],
    nalset_=*nalset,*ialset=*ialsetp,*istartset_=NULL,*iendset_=NULL,
    *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,
    *imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,mortar=0,*imastop=NULL,
    *iponoels=NULL,*inoels=NULL,*nnn=*nnnp;

  double *d=NULL, *z=NULL, *b=NULL, *zeta=NULL,*stiini=NULL,
    *cd=NULL, *cv=NULL, *xforcact=NULL, *xloadact=NULL,*cc=NULL,
    *t1act=NULL, *ampli=NULL, *aa=NULL, *bb=NULL, *aanew=NULL, *bj=NULL, 
    *v=NULL,
    *stn=NULL, *stx=NULL, *een=NULL, *adb=NULL,*xstiff=NULL,*bjp=NULL,
    *aub=NULL, *temp_array1=NULL, *temp_array2=NULL, *aux=NULL,
    *f=NULL, *fn=NULL, *xbounact=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL, *xbodyact=NULL, *cgr=NULL, *au=NULL, *vbounact=NULL,
    *abounact=NULL,*time=NULL,dtime,reltime,*t0=NULL,*t1=NULL,*t1old=NULL,
    physcon[1],zetaj,dj,ddj,h1,h2,h3,h4,h5,h6,sum,aai,bbi,tstart,tend,
    qa[3],cam[3],accold[1],bet,gam,*ad=NULL,sigma=0.,alpham,betam,
    *bact=NULL,*bmin=NULL,*co=NULL,*xboun=NULL,*xbounold=NULL,*vold=NULL,
    *eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,*coefmpcold,*veold=NULL,
    *xini=NULL,*rwork=NULL,*adc=NULL,*auc=NULL,*zc=NULL, *rpar=NULL,
    *cg=NULL,*straight=NULL,xl[27],voldl[36],elas[21],fnl[27],t0l,t1l,
    elconloc[21],veoldl[27],setnull,deltmx,bbmax,dd,dtheta,dthetaref,
    theta,*vini=NULL,dthetaold,*bcont=NULL,*aatrial=NULL,*bmech=NULL,
    *bbtrial=NULL,*vr=NULL,*vi=NULL,*bcontini=NULL,*dbcont=NULL,*dbcontini=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,precision,resultmaxprev,
    resultmax,func,funcp,fexp,fexm,fcos,fsin,sump,*bp=NULL,h14,senergy=0.0,
    *bv=NULL,*cstr=NULL,*aube=NULL,*adbe=NULL,*sti=*stip,time0=0.0,
    time1=0.0,dtemp=0.0;

  FILE *f1;

  /* dummy variables for nonlinmpc */

  int *iaux=NULL,maxlenmpc,icascade=0,newstep=0,iit=1,idiscon;

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

  /*    time=NNEW(double,1);
	time[0]=0.;*/
//  time0=0.;
//  time1=0.;
  dtime=*tinc;

  alpham=xmodal[0];
  betam=xmodal[1];

  dd=ctrl[16];deltmx=ctrl[26];
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
  nsectors=1;

  if(*mcs==0){

      nkg=*nk;
      neg=*ne;

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

      for(i=0;i<*mcs;i++){
	  if(cs[17*i]>nsectors) nsectors=cs[17*i];
      }

        /* determining the maximum number of sectors to be plotted */

      ngraph=1;
      for(j=0;j<*mcs;j++){
	  if(cs[17*j+4]>ngraph) ngraph=cs[17*j+4];
      }
      nkg=*nk*ngraph;
      neg=*ne*ngraph;
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

      /* determining the space needed to expand the
         contact surfaces */

      tchar1=NNEW(char,81);
      tchar2=NNEW(char,81);
      tchar3=NNEW(char,81);
      for(i=0; i<*ntie; i++){
	if(tieset[i*(81*3)+80]=='C'){
	  //a contact constrain was found, so increase nalset
	  memcpy(tchar2,&tieset[i*(81*3)+81],81);
	  tchar2[80]='\0';
	  memcpy(tchar3,&tieset[i*(81*3)+81+81],81);
	  tchar3[80]='\0';
	  for(j=0; j<*nset; j++){
	    memcpy(tchar1,&set[j*81],81);
	    tchar1[80]='\0';
	    if(strcmp(tchar1,tchar2)==0){
	      //dependent nodal surface was found
	      (*nalset)+=(iendset[j]-istartset[j]+1)*(nsectors);
	    }
	    else if(strcmp(tchar1,tchar3)==0){
	      //independent element face surface was found
	      (*nalset)+=(iendset[j]-istartset[j]+1)*(nsectors);
	    }
	  }
	}
      }
      free(tchar1);
      free(tchar2);
      free(tchar3);

      RENEW(ialset,int,*nalset);

      /* save the information in istarset and isendset */
      istartset_=NNEW(int,*nset);
      iendset_=NNEW(int,*nset);
      for(j=0; j<*nset; j++){
	istartset_[j]=istartset[j];
	iendset_[j]=iendset[j];
      }

      RENEW(xstiff,double,27**mint_**ne*nsectors);
      
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
        labmpcold,&nmpcold,xloadold,iamload,t1old,t1,iamt1,xstiff,
        &icole,&jqe,&irowe,isolver,nzse,&adbe,&aube,iexpl,
	ibody,xbody,nbody,cocon,ncocon,tieset,ntie,&nnn);

      free(vold);vold=NNEW(double,5**nk);
      free(veold);veold=NNEW(double,4**nk);
      RENEW(eme,double,6**mint_**ne);
      RENEW(sti,double,6**mint_**ne);

      if(*nener==1) RENEW(ener,double,*mint_**ne*2);
  }

  fclose(f1);
	  
  /* checking for steadystate calculations */

  if(*tper<0){
      precision=-*tper;
      *tper=1.e10;
      isteadystate=1;
  }else{
      isteadystate=0;
  }

  /* normalizing the time */

  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,&inext,&itp));
  dtheta=(*tinc)/(*tper);
  dthetaref=dtheta;
  dthetaold=dtheta;

  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;

  /* check for rigid body modes 
     if there is a jump of 1.e4 in two subsequent eigenvalues
     all eigenvalues preceding the jump are considered to
     be rigid body modes and their frequency is set to zero */

  setnull=1.;
  for(i=nev-2;i>-1;i--){
      if(fabs(d[i])<0.0001*fabs(d[i+1])) setnull=0.;
      d[i]*=setnull;
  }

  /* check whether there are dashpot elements */

  dashpot=0;
  for(i=0;i<*ne;i++){
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[i*8],"ED")==0){
	  dashpot=1;break;}
  }

  if(dashpot){

      if(*mcs!=0){
	  printf("*ERROR in dyna: dashpots are not allowed in combination with cyclic symmetry\n");
	  FORTRAN(stop,());
      }

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
	      ttime,&time0,istep,&iinc,ibody));

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
      rpar=NNEW(double,1);
  }

  if(*idrct==0){bmech=NNEW(double,neq[1]);}

  /* contact conditions */

  inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	  lakon,ipkon,kon,&koncont,&ncone,tietol,&ismallsliding,&itiefac,
          &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
          &mortar,&imastop,nkon,&iponoels,&inoels);

  if(ncont!=0){

      if(*idrct==1){
	  printf("*ERROR in dyna: contact is not allowed in combination with fixed increments\n");
	  FORTRAN(stop,());
      }
      if(dashpot){
	  printf("*ERROR in dyna: contact is not allowed in combination with dashpots\n");
	  FORTRAN(stop,());
      }
      RENEW(ipkon,int,*ne+ncone);
      RENEW(lakon,char,8*(*ne+ncone));
      if(*nener==1){
	RENEW(ener,double,*mint_*(*ne+ncone)*2);
      }

      /* 10 instead of 9: last position is reserved for how
         many dependent nodes are paired to this face */
    
      RENEW(kon,int,*nkon+10*ncone);
      if(*norien>0){
	  RENEW(ielorien,int,*ne+ncone);
	  for(k=*ne;k<*ne+ncone;k++) ielorien[k]=0;
      }
      RENEW(ielmat,int,*ne+ncone);
      for(k=*ne;k<*ne+ncone;k++) ielmat[k]=1;
      cg=NNEW(double,3*ncont);
      straight=NNEW(double,16*ncont);
      ifcont1=NNEW(int,ncone);
      ifcont2=NNEW(int,ncone);
      vini=NNEW(double,5**nk);
      //      for(i=0;i<*nk;i++){vini[i]=vold[i];}
      bcontini=NNEW(double,neq[1]);
      bcont=NNEW(double,neq[1]);
      aatrial=NNEW(double,nev);
      bbtrial=NNEW(double,nev);
  }

  /* storing the element and topology information before introducing 
     contact elements */

  ne0=*ne;nkon0=*nkon;

  zeta=NNEW(double,nev);
  cstr=NNEW(double,6);
  //  stx=NNEW(double,6**mint_**ne);

  /* calculating the damping coefficients*/
  if(xmodal[9]<0){
      for(i=0;i<nev;i++){
	if(fabs(d[i])>(1.e-10)){
	  zeta[i]=(alpham+betam*d[i]*d[i])/(2.*d[i]);
	}
	else {
	  printf("*WARNING in dyna: one of the frequencies is zero\n");
	  printf("         no Rayleigh mass damping allowed\n");
	  zeta[i]=0.;
	}
      }
  }
  else{
    /*copy the damping coefficients for every eigenfrequencie from xmodal[10....] */
    if(nev<(int)xmodal[9]){
      imax=nev;
      printf("*WARNING in dyna: too many modal damping coefficients applied\n");
      printf("         damping coefficients corresponding to nonexisting eigenvalues are ignored\n");
    }
    else{
      imax=(int)xmodal[9];
    }
    for(i=0; i<imax; i++){
      zeta[i]=xmodal[10+i];     
    }
    
  }

  /* modal decomposition of the initial conditions */
  
  cd=NNEW(double,nev);
  cv=NNEW(double,nev);

  if(*mcs==0){
    temp_array1=NNEW(double,neq[1]);
    temp_array2=NNEW(double,neq[1]);
    for(i=0;i<neq[1];i++){temp_array1[i]=0;temp_array2[i]=0;}
    
    /* displacement initial conditions */
    
    for(i=0;i<*nk;i++){
      for(j=0;j<4;j++){
	if(nactdof[4*i+j]!=0){
	  idof=nactdof[4*i+j]-1;
	  temp_array1[idof]=vold[5*i+j];
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
		      
  }
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
/*      ipobody=NNEW(int,2*ifreebody**nbody);*/
      ipobody=NNEW(int,2**ne);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,int,2*(*ne+ifreebody));
      }
      RENEW(ipobody,int,2*(ifreebody-1));
  }
	       
  b=NNEW(double,neq[1]); /* load rhs vector and displacement solution vector */
  bp=NNEW(double,neq[1]); /* velocity solution vector */
  bj=NNEW(double,nev); /* response modal decomposition */
  bjp=NNEW(double,nev); /* derivative of the response modal decomposition */
  ampli=NNEW(double,*nam); /* instantaneous amplitude */

  /* constant coefficient of the linear amplitude function */
  aa=NNEW(double,nev); 
  aanew=NNEW(double,nev);
  /* linear coefficient of the linear amplitude function */
  bb=NNEW(double,nev);
  
  v=NNEW(double,5**nk);
  fn=NNEW(double,4**nk);
  stn=NNEW(double,6**nk);
  inum=NNEW(int,*nk);

  if(*ithermal>1) {qfn=NNEW(double,3**nk);qfx=NNEW(double,3**mint_**ne);}

  if(strcmp1(&filab[18],"E   ")==0) een=NNEW(double,6**nk);
  if(strcmp1(&filab[36],"ENER")==0) enern=NNEW(double,*nk);

  eei=NNEW(double,6**mint_**ne);
  if(*nener==1){
    stiini=NNEW(double,6**mint_**ne);
    enerini=NNEW(double,*mint_**ne);}

/*     calculating the instantaneous loads (forces, surface loading, 
       centrifugal and gravity loading or temperature) at time 0 */

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
     xload,xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
     xbodyact,t1old,t1,t1act,iamt1,nk,
     amta,namta,nam,ampli,&time0,&reltime,ttime,&dtime,ithermal,nmethod,
     xbounold,xboun,xbounact,iamboun,nboun,
     nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
     co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload));

  /*  calculating the instantaneous loading vector at time 0 */

  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
       nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,nbody,
       cgr,b,nactdof,&neq[1],nmethod,
       ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,alcon,
       nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,t0,t1act,
       ithermal,iprestr,vold,iperturb,iexpl,plicon,
       nplicon,plkcon,nplkcon,npmat_,ttime,&time0,istep,&iinc,&dtime,
       physcon,ibody,xbodyold,&reltime,veold,matname));

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

      if(*mcs!=0){
	  printf("*ERROR in dyna: prescribed boundaries are not allowed in combination with cyclic symmetry\n");
	  FORTRAN(stop,());
      }

      if(*idrct!=1){
	  printf("*ERROR in dyna: variable increment length is not allwed in combination with prescribed boundaries\n");
	  FORTRAN(stop,());
      }
      
      /* LU decomposition of the stiffness matrix */
      
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
                         &symmetryflag,&inputformat);
#else
	  printf("*ERROR in dyna: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
	  printf("*ERROR in dyna: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
	  printf("*ERROR in dyna: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1]);
#else
	  printf("*ERROR in dyna: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }

      bact=NNEW(double,neq[1]);
      bmin=NNEW(double,neq[1]);
      bv=NNEW(double,neq[1]);

      init=1;
      dynboun(amta,namta,nam,ampli,&time0,ttime,&dtime,xbounold,xboun,
	      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
	      aub,icol,irow,neq,nzs,&sigma,b,isolver,
	      &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv);
      init=0;
  }

  if(*idrct==0){
      for(i=0;i<neq[1];i++){bmech[i]=b[i];}
  }

/* creating contact elements and calculating the contact forces
   (normal and shear) */

  if(ncont!=0){
      for(i=0;i<neq[1];i++){bcont[i]=0.;}
      contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	      co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
              ifcont1,ifcont2,&ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
              &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,iperturb,
              ikboun,nboun);
      for(i=ne0;i<*ne;i++){
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
	      fnl,ncmat_,ntmat_,&nope,lakonl,&t0l,&t1l,&kodem,elconloc,
	      plicon,nplicon,npmat_,veoldl,&senergy,&iener,cstr));
	  if(i==ne0) printf("spring start step %e\n",fnl[3*nope-3]);
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
      for(i=0;i<neq[1];i++){b[i]+=bcont[i];}
  }

  for(i=0;i<nev;i++){
      aanew[i]=0.;
      for(j=0;j<neq[1];j++){
	  aanew[i]+=z[i*neq[1]+j]*b[j];
      }
  }
      
  if(*idrct==1){
      niter=1;
  }else{
      niter=2;
  }

  /* major loop */

  resultmaxprev=0.;
  resultmax=0.;

  while(1.-theta>1.e-6){
    
//    dtemp=time1;
//    time1=time0;
//    time0=dtemp;
    time0=time1;
    
    printf("\nnew increment\n");
    
    if(*nener==1){
      for(k=0; k<*mint_*ne0; ++k){
	enerini[k]=ener[k];
      }
      if(*ithermal!=2){
	for(k=0; k<6**mint_*ne0; ++k){
	  stiini[k]=sti[k];
	}
      }
    }

    if(ncont!=0){
      for(i=0;i<*nk;i++){vini[i]=vold[i];}
    }
    iinc++;
    jprint++;
    //      RENEW(time,double,iinc+1);
    //      RENEW(rpar,double,3+iinc+nev*(1+2*iinc+nev));
    RENEW(rpar,double,4+nev*(3+nev));
    //      RENEW(aa,double,(iinc+1)*nev); 
    //RENEW(bb,double,(iinc+1)*nev);
      
    /* check for max. # of increments */
    
    if(iinc>*jmax){
      printf(" *ERROR: max. # of increments reached\n\n");
      FORTRAN(stop,());
    }
    
      if(iinc>1){
	for(i=0; i<nev; i++){
	  cd[i]=bj[i];
	  cv[i]=bjp[i];
	}
      }

      
      if((*idrct!=1)&&(iinc!=1)){
	
	/* increasing the increment size */
	
        dthetaold=dtheta;
        dtheta=dthetaref*dd;
	
	/* check increment length whether
	   - it does not exceed tmax
	   - the step length is not exceeded
	   - a time point is not exceeded  */
	
	checkinclength(&time0,ttime,&theta,&dtheta,idrct,tper,tmax,
		       tmin,ctrl, amta,namta,itpamp,&inext,&dthetaref,&itp,
		       &jprint,jout);
	
	dthetaref=dtheta;
    }
      
      reltime=theta+dtheta;
      //      time0=time1;
      time1=reltime**tper;
      dtime=dtheta**tper;

      //      printf("dtime=%e\n time=%e\n reltime=%e\n theta=%e\n dtheta=%e\n time0=%e\n tper=%e\n",dtime,time1,reltime,theta,dtheta,time0,*tper);

      printf("dtime=%e\n time=%e\n",dtime,time1);

      
      /* calculating the instantaneous loads (forces, surface loading, 
	 centrifugal and gravity loading or temperature) */
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	      xloadold,xload,xloadact,iamload,nload,ibody,xbody,
	      nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time1,&reltime,ttime,&dtime,ithermal,
	      nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
	      ndirboun,nodeforc,
	      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
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
		    npmat_,ttime,&time1,istep,&iinc,&dtime,physcon,ibody,
		    xbodyold,&reltime,veold,matname));
	      
      /* correction for nonzero SPC's */
      
      if(iprescribedboundary){
	  dynboun(amta,namta,nam,ampli,&time1,ttime,&dtime,
                      xbounold,xboun,
		      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		      aub,icol,irow,neq,nzs,&sigma,b,isolver,
		      &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv);
      }

      if(*idrct==0){
	  bbmax=0.;
	  for(i=0;i<neq[1];i++){
	      if(fabs(b[i]-bmech[i])>bbmax) bbmax=fabs(b[i]-bmech[i]);
	  }

	  /* check for size of mechanical force */

	  if(bbmax>deltmx){

	      /* force increase too big: increment size is decreased */
	      
	      dtheta=dtheta*deltmx/bbmax;
	      printf("correction of dtheta due to force increase: %e\n",dtheta);
	      dthetaref=dtheta;
	      itp=0;
	      
	      /* check whether the new increment size is not too small */
	      
	      if(dtheta<*tmin){
		  printf("\n *WARNING: increment size %e smaller than minimum %e\n",dtheta**tper,*tmin**tper);
		  printf("             minimum is taken\n");
		  dtheta=*tmin;
	      }

	      reltime=theta+dtheta;
	      time1=reltime**tper;
	      dtime=dtheta**tper;
	      
	      /* calculating the instantaneous loads (forces, surface loading, 
		 centrifugal and gravity loading or temperature) */
	      
	      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	        xloadold,xload,xloadact,iamload,nload,ibody,xbody,
	        nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
	        namta,nam,ampli,&time1,&reltime,ttime,&dtime,ithermal,
	        nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
	        ndirboun,nodeforc,
	        ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
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
		    npmat_,ttime,&time1,istep,&iinc,&dtime,physcon,ibody,
		    xbodyold,&reltime,veold,matname));
	      
	      /* correction for nonzero SPC's */
	      
	      if(iprescribedboundary){
		  dynboun(amta,namta,nam,ampli,&time1,ttime,&dtime,
                      xbounold,xboun,
		      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		      aub,icol,irow,neq,nzs,&sigma,b,isolver,
		      &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv);
	      }
	  }

	  for(i=0;i<neq[1];i++){bmech[i]=b[i];}
	  
	  if(ncont!=0){
	      for(i=0;i<neq[1];i++){b[i]+=bcont[i];}
	      for(i=0;i<neq[1];i++){bcontini[i]=bcont[i];}
	  }

      }

      /* step length is OK for mechanical load
	 calculating equation for linearized loading */

      for(i=0;i<nev;i++){
	  aa[i]=aanew[i];
	  aanew[i]=0.;
	  for(j=0;j<neq[1];j++){
	      aanew[i]+=z[i*neq[1]+j]*b[j];
	  }
	  bb[i]=(aanew[i]-aa[i])/dtime;
	  aa[i]=aanew[i]-bb[i]*time1;
      }
      /* calculating the response due to unchanged contact force during
         the increment */

      if(dashpot){
	  FORTRAN(subspace,(d,aa,bb,cc,&alpham,&betam,&nev,xini,
			    cd,cv,&time1,rwork,&lrw,&iinc,jout,rpar,bj,
			    iwork,&liw,&iddebdf,bjp));
	  if(iddebdf==2){
	      liw=56+2*nev;
	      RENEW(iwork,int,liw);
	      for(i=0;i<liw;i++){iwork[i]=0;}
	      lrw=250+20*nev+4*nev*nev;
	      RENEW(rwork,double,lrw);
	      for(i=0;i<lrw;i++){rwork[i]=0.;}
	      iddebdf=1;
	      FORTRAN(subspace,(d,aa,bb,cc,&alpham,&betam,&nev,xini,
				cd,cv,&time1,rwork,&lrw,&iinc,jout,rpar,bj,
				iwork,&liw,&iddebdf,bjp));
	  }
      }
      else{
	  for(l=0;l<nev;l++){
	      zetaj=zeta[l];
	      dj=d[l];
	      
	      /* zero eigenfrequency: rigid body mode */
	      
	      if(fabs(d[l])<=1.e-10){
		  sum=0.;sump=0.;
		  for(i=iinc-1;i<iinc;i++){
		      aai=aa[l];
		      bbi=bb[l];
		      tstart=time0;
		      tend=time1;
		      sum+=tend*(aai*time1+
				 tend*((bbi*time1-aai)/2.-bbi*tend/3.))-
			  tstart*(aai*time1+
				  tstart*((bbi*time1-aai)/2.-bbi*tstart/3.));
		      sump+=tend*(aai+bbi*tend/2.)-tstart*(aai+bbi*tstart/2.);
		  }
		  /*		  bj[l]=sum+cd[l]+time[iinc]*cv[l];
				  bjp[l]=sump+cv[l];*/
		  bj[l]=sum+cd[l]+dtime*cv[l];
		  bjp[l]=sump+cv[l];
	      }
	      
	      /*   subcritical damping */
	      
	      else if(zetaj<1.-1.e-6){
		  ddj=dj*sqrt(1.-zetaj*zetaj);
		  h1=zetaj*dj;
		  h2=h1*h1+ddj*ddj;
		  h3=h1*h1-ddj*ddj;
		  h4=2.*h1*ddj/h2;
		  h14=h1/ddj;
		  sum=0.;sump=0.;
		  for(i=iinc-1;i<iinc;i++){
		      aai=aa[l];
		      bbi=bb[l];
		      tstart=0;		      
		      tend=time1-time0;
		      FORTRAN(fsub,(&time1,&tend,&aai,&bbi,&ddj,
					 &h1,&h2,&h3,&h4,&func,&funcp));
		      sum+=func;sump+=funcp;
		      FORTRAN(fsub,(&time1,&tstart,&aai,&bbi,&ddj,
					 &h1,&h2,&h3,&h4,&func,&funcp));
		      sum-=func;sump-=funcp;
		  }
		  /* fexp=exp(-h1*time[iinc]);
		  fsin=sin(ddj*time[iinc]);
		  fcos=cos(ddj*time[iinc]);*/
		  fexp=exp(-h1*dtime);
		  fsin=sin(ddj*dtime);
		  fcos=cos(ddj*dtime);

		  bj[l]=sum/ddj+fexp*(fcos+zetaj/sqrt(1.-zetaj*zetaj)*fsin)*cd[l]+
		      fexp*fsin*cv[l]/ddj;
		  bjp[l]=sump/ddj+fexp*((-h1+ddj*h14)*fcos+(-ddj-h1*h14)*fsin)*cd[l]
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
		  for(i=iinc-1;i<iinc;i++){
		      aai=aa[l];
		      bbi=bb[l];
		      tstart=0;		      
		      tend=time1-time0;
		      FORTRAN(fsuper,(&time1,&tend,&aai,&bbi,
					   &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
		      sum+=func;sump+=funcp;
		      FORTRAN(fsuper,(&time1,&tstart,&aai,&bbi,
					   &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
		      sum-=func;sump-=funcp;
	      
		  }
		  /*		  fexm=exp(h1*time[iinc]);
				  fexp=exp(-h2*time[iinc]);*/
		  fexm=exp(h1*dtime);
		  fexp=exp(-h2*dtime);
		  h14=zetaj*dj/ddj;
		  bj[l]=sum/(2.*ddj)+(fexm+fexp)*cd[l]/2.+zetaj*(fexm-fexp)/(2.*
		      sqrt(zetaj*zetaj-1.))*cd[l]+(fexm-fexp)*cv[l]/(2.*ddj);
		  bjp[l]=sump/(2.*ddj)+(h1*fexm-h2*fexp)*cd[l]/2.
		    +(h14*cd[l]+cv[l]/ddj)*(h1*fexm+h2*fexp)/2.;
	      }
	      
	      /* critical damping */
	      
	      else{
		  h1=zetaj*dj;
		  h2=1./h1;
		  h3=h2*h2;
		  h4=h2*h3;
		  sum=0.;sump=0.;
		  for(i=iinc-1;i<iinc;i++){
		      aai=aa[l];
		      bbi=bb[l];
		      tstart=0;
		      tend=time1-time0;
		      FORTRAN(fcrit,(&time1,&tend,&aai,&bbi,&zetaj,&dj,
					  &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
		      sum+=func;sump+=funcp;
		      FORTRAN(fcrit,(&time1,&tstart,&aai,&bbi,&zetaj,&dj,
					  &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
		      sum-=func;sump-=funcp;
		  }
		  /*		  fexp=exp(-h1*time[iinc]);
		  bj[l]=sum+fexp*((1.+h1*time[iinc])*cd[l]+time[iinc]*cv[l]);
		  bjp[l]=sump+fexp*(-h1*h1*time[iinc]*cd[l]+
		  (1.-h1*time[iinc])*cv[l]);*/
		  fexp=exp(-h1*dtime);
		  bj[l]=sum+fexp*((1.+h1*dtime)*cd[l]+dtime*cv[l]);
		  bjp[l]=sump+fexp*(-h1*h1*dtime*cd[l]+
                                    (1.-h1*dtime)*cv[l]);
	      }
	  }
      }
      
      /* composing the response */
      
      if(iprescribedboundary){
	  for(i=0;i<neq[1];i++){
	      b[i]=bmin[i];
	      bp[i]=bv[i];
	  }
      }
      else{
	for(i=0;i<neq[1];i++){
	  b[i]=0.;
	     bp[i]=0.;
	}
      }
      
      for(i=0;i<neq[1];i++){
	  for(j=0;j<nev;j++){
	    b[i]+=bj[j]*z[j*neq[1]+i];
	    bp[i]+=bjp[j]*z[j*neq[1]+i];
	  }
      }
      
      /* update nonlinear MPC-coefficients (e.g. for rigid
	 body MPC's */
      
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal));
      
      /* calculating displacements/temperatures */
      
      FORTRAN(dynresults,(nk,v,ithermal,nactdof,vold,nodeboun,
			  ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			  b,bp,veold,&dtime));
      
      /* calculating the velocity field if there are dashpots
	 or if the step is finished */
      
      /* if((dashpot)||(1.-theta<=1.e-6)||(ncont!=0)){
         for(i=0;i<*nk;i++){
	      for(j=1;j<4;j++){
		  veold[4*i+j]=(v[5*i+j]-vold[5*i+j])/(*jout*dtime);
	      }
	  }
	   }*/
	  
      //     if(ncont!=0){
      //	  for(i=0;i<5**nk;i++){vini[i]=vold[i];}
      //      }
      for(i=0;i<5**nk;i++){vold[i]=v[i];}
      
      /* creating contact elements and calculating the contact forces
	 based on the displacements at the end of the present increment */
      
      if(ncont!=0){
	  dynacont(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,labmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,nactdof,neq,nzl,icol,
              irow,
	      nmethod,ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
	      nrhcon,cocon,ncocon,alcon,nalcon,alzero,ielmat,ielorien,
              norien,orab,ntmat_,t0,t1,ithermal,prestr,iprestr,
	      vold,iperturb,sti,nzs,tinc,tper,xmodal,veold,amname,amta,
	      namta,nam,iamforc,iamload,iamt1,jout,filab,eme,xforcold,
	      xloadold,t1old,iamboun,xbounold,iexpl,plicon,nplicon,plkcon,
              nplkcon,xstate,npmat_,matname,mint_,ncmat_,nstate_,ener,jobnamec,
	      ttime,set,nset,istartset,iendset,ialset,nprint,prlab,
	      prset,nener,trab,inotr,ntrans,fmpc,cbody,ibody,xbody,nbody,
              xbodyold,istep,isolver,jq,output,mcs,nkon,mpcend,ics,cs,ntie,
              tieset,idrct,jmax,tmin,tmax,ctrl,itpamp,tietol,&iit,
	      &ncont,&ne0,&reltime,&dtime,bcontini,bj,aux,iaux,bcont,&dbcont,
	      &nev,v,&nkon0,&deltmx,&dtheta,&theta,&iprescribedboundary,
	      &mpcfree,&memmpc_,itietri,koncont,cg,straight,&iinc,
	      ifcont1,ifcont2,vini,aa,bb,aanew,d,z,zeta,b,bmech,&time0,&time1,ipobody,
	      xforcact,xloadact,t1act,xbounact,xbodyact,cd,cv,ampli,
              &dthetaref,bjp,bp,cstr,&dbcontini);
      }   
	  
      theta+=dtheta;
      (*ttime)+=dtime;
      
      /* check whether a time point was reached */
      
      if((*itpamp>0)&&(*idrct==0)){
	  if(itp==1){
	      jprint=*jout;
	  }else{
	      jprint=*jout+1;
	  }
      }
      
      /* check whether output is needed */
      
      if((*jout==jprint)||(1.-theta<=1.e-6)){
	  iout=2;
	  jprint=0;
      }else if(*nener==1){
	  iout=-2;
      }else{
	  iout=0;
      }
      
      if((iout==2)||(iout==-2)){
	stx=NNEW(double,6**mint_**ne);
	FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	         stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		 ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		 ithermal,prestr,iprestr,filab,eme,een,
		 iperturb,f,fn,nactdof,&iout,qa,
		 vold,b,nodeboun,ndirboun,xbounact,nboun,
		 ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		 veold,accold,&bet,&gam,&dtime,&time1,ttime,
		 plicon,nplicon,plkcon,nplkcon,
		 xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
		 &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		 enern,sti,xstaten,eei,enerini,cocon,ncocon,
		 set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		 qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,
		 ilmpc,istep,&iinc));
	
	
	if(*ithermal!=2){
	  for(k=0;k<6**mint_*ne0;++k){
	    sti[k]=stx[k];
	  }
	}
      }
      if(iout==2){
    	(*kode)++;
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	FORTRAN(out,(co,&nkg,kon,ipkon,lakon,&neg,v,stn,inum,nmethod,kode,filab,
		     een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,nstate_,istep,&iinc,
		     iperturb,ener,mint_,output,ithermal,qfn,&mode,&noddiam,
		     trab,inotr,ntrans,orab,ielorien,norien,description,
		     ipneigh,neigh,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ne,cs));
	
	free(ipneigh);free(neigh);
      }
      
      if((iout==2)||(iout==-2)){
	free(stx);
      }
      
      
      FORTRAN(writesummary,(istep,&iinc,&icutb,&iit,ttime,&time1,&dtime));
      
      if(isteadystate==1){
	  
	  /* calculate maximum displacement/temperature */
	  
	  resultmax=0.;
	  if(*ithermal<2){
	      for(i=1;i<5**nk;i=i+5){
		  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	      for(i=2;i<5**nk;i=i+5){
		  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	      for(i=3;i<5**nk;i=i+5){
		  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	  }else if(*ithermal==2){
	      for(i=0;i<5**nk;i=i+5){
		  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	  }else{
	      printf("*ERROR in dyna: coupled temperature-displacement calculations are not allowed\n");
	  }
	  if(fabs((resultmax-resultmaxprev)/resultmax)<precision){
	      break;
	  }else{resultmaxprev=resultmax;}
      }
     
  }
  
  free(eei);
  free(vbounact);
  free(abounact);

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
  
  free(v);free(fn);free(stn);free(inum);free(adb);
  free(aub);free(z);free(b);free(zeta);free(bj);free(cd);free(cv);
  free(xforcact);free(xloadact);free(xbounact);free(aa);free(bb);free(aanew);
  free(ampli);free(xbodyact);free(bjp);free(bp);

  if(*ithermal==1) free(t1act);

  if(iprescribedboundary){
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
	  pardiso_cleanup(&neq[1]);
#endif
      }
      free(bact);free(bmin);free(bv);
  }

  /* deleting the contact information */
  *ne=ne0; *nkon=nkon0;
  if(ncont!=0){
      RENEW(ipkon,int,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,int,*nkon);
      if(*nener==1){
	RENEW(ener,double,*mint_**ne*2);
      }
      if(*norien>0){
	  RENEW(ielorien,int,*ne);
      }
      RENEW(ielmat,int,*ne);
      free(cg);free(straight);free(vini);free(bcont);free(bmech);
      free(aatrial);free(bbtrial);free(ifcont1);free(ifcont2);
  }

  if(*mcs==0){
      free(ad);free(au);
  }else{
      free(adbe); free(aube);free(icole); free(irowe); free(jqe);

      *nk/=nsectors;
      *ne/=nsectors;
      *nkon/=nsectors;
      *nboun/=nsectors;
      neq[1]=neq[1]*2/nsectors;

      RENEW(nnn,int,*nk);

      RENEW(ialset,int,nalset_);
      /* restore the infomration in istartset and iendset */
      for(j=0; j<*nset; j++){
	istartset[j]=istartset_[j];
	iendset[j]=iendset_[j];
      }
      free(istartset_);
      free(iendset_);

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

      RENEW(vold,double,5**nk);
      RENEW(veold,double,4**nk);
      RENEW(eme,double,6**mint_**ne);
      RENEW(sti,double,6**mint_**ne);
      if(*nener==1)RENEW(ener,double,*mint_**ne*2);

/* distributed loads */

      for(i=0;i<*nload;i++){
	  if(nelemload[2*i]<nsectors){
	      nelemload[2*i]-=*ne*nelemload[2*i+1];
	  }else{
	      nelemload[2*i]-=*ne*(nelemload[2*i+1]-nsectors);
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
	  if(nodeforc[2*i+1]<nsectors){
	      nodeforc[2*i]-=*nk*nodeforc[2*i+1];
	  }else{
	      nodeforc[2*i]-=*nk*(nodeforc[2*i+1]-nsectors);
	  }
      }
  }

  free(xstiff);if(*nbody>0) free(ipobody);

  if(dashpot){
      free(xini);free(rwork);free(adc);free(auc);free(cc);
      free(rpar);free(iwork);}

  free(cstr);
  //  free(time);

  *ialsetp=ialset;
  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*veoldp=veold;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;
  *nnnp=nnn;*stip=sti;

  return;
}

