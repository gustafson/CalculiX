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

void complexfreq(double **cop, int *nk, int **konp, int **ipkonp, char **lakonp, int *ne, 
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
               double *xstate, int *npmat_, char *matname, int *mi,
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
	       int **nnnp, int *ikforc, int *ilforc, double *thicke,
	       char *jobnamef,int *mei){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *lakont=NULL;

  int nev,i,j,k,idof,*inum=NULL,*ipobody=NULL,inewton=0,id,
    iinc=0,l,iout=1,ielas,icmd,ifreebody,mode,m,nherm,
    *kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,
    *inotr=NULL,*nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,nsectors,nmd,nevd,*nm=NULL,*iamt1=NULL,
    ngraph=1,nkg,neg,ne0,ij,lprev,nope,indexe,ilength,
    *ipneigh=NULL,*neigh=NULL,index,im,cyclicsymmetry,inode,
    *ialset=*ialsetp,mt=mi[1]+1,*nnn=*nnnp,kmin,kmax,i1,
    *iter=NULL,lint,lfin,kk,kkv,kk6,kkx,icomplex,igeneralizedforce,
    idir,*inumt=NULL,icntrl,imag,jj,is,l1,*inocs=NULL,ml1,l2,nkt,net,
    *ipkont=NULL,*ielmatt=NULL,*inotrt=NULL,*kont=NULL,node,iel,*ielcs=NULL,
    ielset,*istartnmd=NULL,*iendnmd=NULL,inmd,neqact,*nshcon=NULL,
    *ipev=NULL,icfd=0,*inomat=NULL;

  long long i2;

  double *d=NULL, *z=NULL,*stiini=NULL,*cc=NULL,*v=NULL,*zz=NULL,*emn=NULL,
    *stn=NULL, *stx=NULL, *een=NULL, *adb=NULL,*xstiff=NULL,
    *aub=NULL,*aux=NULL,*f=NULL, *fn=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL, *cgr=NULL, *au=NULL,dtime,reltime,*t0=NULL,*t1=NULL,*t1old=NULL,
    sum,qa[3],cam[5],accold[1],bet,gam,*ad=NULL,alpham,betam,
    *co=NULL,*xboun=NULL,*xbounold=NULL,*vold=NULL,*emeini=NULL,
    *eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,*veold=NULL,
    *adc=NULL,*auc=NULL,*zc=NULL,*fnr=NULL,*fni=NULL,setnull,deltmx,dd,
    theta,*vini=NULL,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,*vmax=NULL,
    *stnmax=NULL,*cstr=NULL,*sti=*stip,time0=0.0,time=0.0,zero=0.0,
    *springarea=NULL,*eenmax=NULL,*aa=NULL,*bb=NULL,*xx=NULL,
    *eiga=NULL,*eigb=NULL,*eigxx=NULL,*temp=NULL,*coefmpcnew=NULL,xreal,
    ximag,t[3],*vt=NULL,*t1t=NULL,*stnt=NULL,*eent=NULL,*fnt=NULL,*enernt=NULL,
    *stxt=NULL,pi,ctl,stl,*cot=NULL,*qfnt=NULL,vreal,vimag,constant,stnreal,
    stnimag,freq,*emnt=NULL,*xnormastface=NULL,*shcon=NULL,*eig=NULL,
    *eigxr=NULL,*eigxi=NULL,*xmac=NULL,*bett=NULL,*betm=NULL,*xmaccpx=NULL,
    fmin=0.,fmax=1.e30,*xmr=NULL,*xmi=NULL,*zi=NULL,*eigx=NULL;

  FILE *f1;

#ifdef SGI
  int token;
#endif

  pi=4.*atan(1.);
  constant=180./pi;

  co=*cop;kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;
  ielorien=*ielorienp;inotr=*inotrp;nodeboun=*nodebounp;
  ndirboun=*ndirbounp;iamboun=*iambounp;xboun=*xbounp;
  xbounold=*xbounoldp;ikboun=*ikbounp;ilboun=*ilbounp;nactdof=*nactdofp;
  vold=*voldp;eme=*emep;ener=*enerp;ipompc=*ipompcp;nodempc=*nodempcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;veold=*veoldp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;

  if(ithermal[0]<=1){
      kmin=1;kmax=3;
  }else if(ithermal[0]==2){
      kmin=0;kmax=mi[1];if(kmax>2)kmax=2;
  }else{
      kmin=0;kmax=3;
  }

  xstiff=NNEW(double,(long long)27*mi[0]**ne);

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

  /* opening the eigenvalue file and checking for cyclic symmetry */

  strcpy(fneig,jobnamec);
  strcat(fneig,".eig");

  if((f1=fopen(fneig,"rb"))==NULL){
    printf("*ERROR in complexfreq: cannot open eigenvalue file for reading");
    exit(0);
  }

  printf(" *INFO  in complexfreq: if there are problems reading the .eig file this may be due to:\n");
  printf("        1) the nonexistence of the .eig file\n");
  printf("        2) other boundary conditions than in the input deck\n");
  printf("           which created the .eig file\n\n");

  if(fread(&cyclicsymmetry,sizeof(int),1,f1)!=1){
      printf("*ERROR in complexfreq reading the cyclic symmetry flag in the eigenvalue file");
      exit(0);
  }

  if(fread(&nherm,sizeof(int),1,f1)!=1){
      printf("*ERROR in complexfreq reading the Hermitian flag in the eigenvalue file");
      exit(0);
  }

  if(nherm!=1){
      printf("*ERROR in complexfreq: the eigenvectors in the .eig-file result\n");
      printf("       from a non-Hermitian eigenvalue problem. The complex\n");
      printf("       frequency procedure cannot handle that yet\n\n");
      FORTRAN(stop,());
  }

  nsectors=1;

  if(!cyclicsymmetry){

      nkg=*nk;
      neg=*ne;

      if(fread(&nev,sizeof(int),1,f1)!=1){
	  printf("*ERROR in complexfreq reading the number of eigenvalues in the eigenvalue file...");
	  exit(0);
      }
      

      if(nherm==1){
	  d=NNEW(double,nev);
	  if(fread(d,sizeof(double),nev,f1)!=nev){
	      printf("*ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	      exit(0);
	  }
      }else{
	  d=NNEW(double,2*nev);
	  if(fread(d,sizeof(double),2*nev,f1)!=2*nev){
	      printf("*ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	      exit(0);
	  }
      }
      
      ad=NNEW(double,neq[1]);
      adb=NNEW(double,neq[1]);
      au=NNEW(double,nzs[2]);
      aub=NNEW(double,nzs[1]);
      
      /* reading the stiffness matrix */

      if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR in complexfreq reading the diagonal of the stiffness matrix in the eigenvalue file...");
	  exit(0);
      }
      
      if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf("*ERROR in complexfreq reading the off-diagonal terms of the stiffness matrix in the eigenvalue file...");
	  exit(0);
      }
      
      /* reading the mass matrix */

      if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR in complexfreq reading the diagonal of the mass matrix in eigenvalue file...");
	  exit(0);
      }
      
      if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf("*ERROR in complexfreq reading the off-diagonals of the mass matrix in the eigenvalue file...");
	  exit(0);
      }
      
      /* reading the eigenvectors */

      if(nherm==1){
	  z=NNEW(double,neq[1]*nev);
	  if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
	      printf("*ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	      exit(0);
	  }
      }else{
	  z=NNEW(double,2*neq[1]*nev);
	  if(fread(z,sizeof(double),2*neq[1]*nev,f1)!=2*neq[1]*nev){
	      printf("*ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	      exit(0);
	  }
      }

      nm=NNEW(int,nev);
      for(i=0;i<nev;i++){nm[i]=-1;}
  }
  else{

      if(*nmethod==6){
        printf("*ERROR in complexfreq: Coriolis forces cannot\n");
        printf("       be combined with cyclic symmetry\n\n");
        FORTRAN(stop,());
      }

      nev=0;
      do{
	  if(fread(&nmd,sizeof(int),1,f1)!=1){
	      break;
	  }
	  if(fread(&nevd,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq reading the number of eigenvalues in the eigenvalue file...");
	      exit(0);
	      }
	  if(nev==0){
	      if(nherm==1){d=NNEW(double,nevd);}else{d=NNEW(double,2*nevd);}
	      nm=NNEW(int,nevd);
	  }else{
	      printf("*ERROR in complexfreq: flutter forces cannot\n");
	      printf("       be combined with multiple modal diameters\n");
	      printf("       in cyclic symmetry calculations\n\n");
	      FORTRAN(stop,());
	  }
	  
	  if(nherm==1){
	      if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
		  printf("*ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
		  exit(0);
	      }
	  }else{
	      if(fread(&d[nev],sizeof(double),2*nevd,f1)!=2*nevd){
		  printf("*ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
		  exit(0);
	      }
	  }
	  for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	  
	  if(nev==0){
	      adb=NNEW(double,neq[1]);
	      aub=NNEW(double,nzs[1]);

	      if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
		  printf("*ERROR in complexfreq reading the diagonal of the mass matrix in the eigenvalue file...");
		  exit(0);
	      }
	      
	      if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
		  printf("*ERROR in complexfreq reading the off-diagonals of the mass matrix in the eigenvalue file...");
		  exit(0);
	      }
	  }
	  
	  if(nev==0){
	      z=NNEW(double,neq[1]*nevd);
	  }else{
	      RENEW(z,double,(long long)neq[1]*(nev+nevd));
	  }
	  
	  if(fread(&z[neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	      printf("*ERROR in complexfreq reading eigenvectors in the eigenvalue file...");
	      exit(0);
	  }
	  nev+=nevd;
      }while(1);

      /* removing double eigenmodes */

      j=-1;
      for(i=0;i<nev;i++){
	  if((i/2)*2==i){
	      j++;
	      if(nherm==1){
		  d[j]=d[i];
	      }else{
		  d[2*j]=d[2*i];d[2*j+1]=d[2*i+1];
	      }
	      nm[j]=nm[i];
	      for(k=0;k<neq[1];k++){
		  z[j*neq[1]+k]=z[i*neq[1]+k];
	      }
	  }
      }
      nev=j+1;
      if(nherm==1){RENEW(d,double,nev);}else{RENEW(d,double,2*nev);}
      RENEW(nm,int,nev);
      RENEW(z,double,neq[1]*nev);

      /* determining the maximum amount of segments */

      for(i=0;i<*mcs;i++){
	  if(cs[17*i]>nsectors) nsectors=(int)(cs[17*i]+0.5);
      }

        /* determining the maximum number of sectors to be plotted */

      for(j=0;j<*mcs;j++){
	  if(cs[17*j+4]>ngraph) ngraph=(int)cs[17*j+4];
      }
      nkg=*nk*ngraph;
      neg=*ne*ngraph;

  }

  fclose(f1);

  /* assigning nodes and elements to sectors */

  if(cyclicsymmetry){
    inocs=NNEW(int,*nk);
    ielcs=NNEW(int,*ne);
    ielset=cs[12];
    if((*mcs!=1)||(ielset!=0)){
      for(i=0;i<*nk;i++) inocs[i]=-1;
      for(i=0;i<*ne;i++) ielcs[i]=-1;
    }
    
    for(i=0;i<*mcs;i++){
      is=cs[17*i+4];
      if((is==1)&&(*mcs==1)) continue;
      ielset=cs[17*i+12];
      if(ielset==0) continue;
      for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
	if(ialset[i1]>0){
	  iel=ialset[i1]-1;
	  if(ipkon[iel]<0) continue;
	  ielcs[iel]=i;
	  indexe=ipkon[iel];
	  if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
	  else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
	  else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
	  else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
	  else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
	  else {nope=6;}
	  for(i2=0;i2<nope;++i2){
	    node=kon[indexe+i2]-1;
	    inocs[node]=i;
	  }
	}
	else{
	  iel=ialset[i1-2]-1;
	  do{
	    iel=iel-ialset[i1];
	    if(iel>=ialset[i1-1]-1) break;
	    if(ipkon[iel]<0) continue;
	    ielcs[iel]=i;
	    indexe=ipkon[iel];
	    if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
	    else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
	    else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
	    else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
	    else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
	    else {nope=6;}
	    for(i2=0;i2<nope;++i2){
	      node=kon[indexe+i2]-1;
	      inocs[node]=i;
	    }
	  }while(1);
	}
      } 
    }
  }
  
  /* check for rigid body modes 
     if there is a jump of 1.e4 in two subsequent eigenvalues
     all eigenvalues preceding the jump are considered to
     be rigid body modes and their frequency is set to zero */

  if(nherm==1){
      setnull=1.;
      for(i=nev-2;i>-1;i--){
	  if(fabs(d[i])<0.0001*fabs(d[i+1])) setnull=0.;
	  d[i]*=setnull;
      }
  }else{
      setnull=1.;
      for(i=nev-2;i>-1;i--){
	  if(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])<
                 0.0001*sqrt(d[2*i+2]*d[2*i+2]+d[2*i+3]*d[2*i+3])) setnull=0.;
	  d[2*i]*=setnull;
	  d[2*i+1]*=setnull;
      }
  }

  /* determining the frequency ranges corresponding to one
     and the same nodal diameter */

  if(cyclicsymmetry){
      istartnmd=NNEW(int,nev);
      iendnmd=NNEW(int,nev);
      nmd=0;
      inmd=nm[0];
      istartnmd[0]=1;
      for(i=1;i<nev;i++){
	  if(nm[i]==inmd) continue;
	  iendnmd[nmd]=i;
	  nmd++;
	  istartnmd[nmd]=i+1;
	  inmd=nm[i];
      }
      iendnmd[nmd]=nev;
      nmd++;
      RENEW(istartnmd,int,nmd);
      RENEW(iendnmd,int,nmd);
  }

  if(*nmethod==6){

      /* Coriolis */

      neqact=neq[1];

  /* assigning the body forces to the elements */ 

      ifreebody=*ne+1;
      ipobody=NNEW(int,2**ne);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,int,2*(*ne+ifreebody));
      }
      RENEW(ipobody,int,2*(ifreebody-1));
      
      if(cyclicsymmetry){
	  printf("*ERROR in complexfreq: dashpots are not allowed in combination with cyclic symmetry\n");
	  FORTRAN(stop,());
      }

      adc=NNEW(double,neq[1]);
      auc=NNEW(double,nzs[1]);
      FORTRAN(mafillcorio,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
              xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,ncmat_,
	      ttime,&time0,istep,&iinc,ibody));

      /*  zc = damping matrix * eigenmodes */

      zc=NNEW(double,neq[1]*nev);
      for(i=0;i<nev;i++){
	  FORTRAN(op_corio,(&neq[1],aux,&z[i*neq[1]],&zc[i*neq[1]],adc,auc,
	  icol,irow,nzl));
      }
      free(adc);free(auc);

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
	  for(j=i+1;j<nev;j++){
	      cc[i*nev+j]=-cc[j*nev+i];
	  }
      }
      free(zc);

      /* solving for the complex eigenvalues */

      aa=NNEW(double,4*nev*nev);
      bb=NNEW(double,4*nev*nev);
      xx=NNEW(double,4*nev*nev);
      temp=NNEW(double,4*nev*nev);
      eiga=NNEW(double,2*nev);
      eigb=NNEW(double,2*nev);
      eigxx=NNEW(double,2*nev);
      iter=NNEW(int,nev);

      FORTRAN(coriolissolve,(cc,&nev,aa,bb,xx,eiga,eigb,eigxx,
              iter,d,temp));
      
      free(aa);free(bb);free(temp);free(eiga);free(eigb);free(iter);free(cc);

  }else{

      /* flutter */

      /* complex force is being read (e.g. due to fluid flow) */

      if(cyclicsymmetry){
	  neqact=neq[1]/2;
      }else{
	  neqact=neq[1];
      }
      zc=NNEW(double,2*neqact*nev);
      aa=NNEW(double,4*nev*nev);

      FORTRAN(readforce,(zc,&neqact,&nev,nactdof,ikmpc,nmpc,
			 ipompc,nodempc,mi,coefmpc,jobnamef,
                         aa,&igeneralizedforce));

      bb=NNEW(double,4*nev*nev);
      xx=NNEW(double,4*nev*nev);
      eiga=NNEW(double,2*nev);
      eigb=NNEW(double,2*nev);
      eigxx=NNEW(double,2*nev);
      iter=NNEW(int,nev);
      FORTRAN(forcesolve,(zc,&nev,aa,bb,xx,eiga,eigb,eigxx,iter,d,
      		&neq[1],z,istartnmd,iendnmd,&nmd,&cyclicsymmetry,
		&neqact,&igeneralizedforce));
      free(aa);free(bb);free(eiga);free(eigb);free(iter);free(zc);
      
 }

/* sorting the eigenvalues and eigenmodes according to the size of the
   eigenvalues */
      
  ipev=NNEW(int,nev);
  eigxr=NNEW(double,nev);
  aa=NNEW(double,2*nev);
  bb=NNEW(double,4*nev*nev);
  
  FORTRAN(sortev,(&nev,&nmd,eigxx,&cyclicsymmetry,xx,
		  eigxr,ipev,istartnmd,iendnmd,aa,bb));
  
  free(ipev);free(eigxr);free(aa);free(bb);

  /* storing the eigenvalues in the .dat file */

  if(cyclicsymmetry){
      FORTRAN(writeevcscomplex,(eigxx,&nev,nm,&fmin,&fmax));
  }else{
      FORTRAN(writeevcomplex,(eigxx,&nev,&fmin,&fmax));
  }

  /* storing the participation factors */

  eigxr=NNEW(double,nev);
  eigxi=NNEW(double,nev);
  if(nherm==1){
      eig=NNEW(double,nev);
      for(l=0;l<nev;l++){
	  if(d[l]<0.){
	      eig[l]=0.;
	  }else{
	      eig[l]=sqrt(d[l]);
	  }
      }
  }else{

      eig=NNEW(double,2*nev);
      for(l=0;l<nev;l++){
	  eig[2*l]=sqrt(sqrt(d[2*l]*d[2*l]+d[2*l+1]*d[2*l+1])+d[2*l])/sqrt(2.);
	  eig[2*l+1]=sqrt(sqrt(d[2*l]*d[2*l]+d[2*l+1]*d[2*l+1])-d[2*l])/sqrt(2.);
	  if(d[2*l+1]<0.) eig[2*l+1]=-eig[2*l+1];
      }
  }
  for(l=0;l<nev;l++){
      mode=l+1;
      for(k=0;k<nev;k++){
	  eigxr[k]=xx[2*l*nev+2*k];
	  eigxi[k]=xx[2*l*nev+2*k+1];
      }
      FORTRAN(writepf,(eig,eigxr,eigxi,&zero,&nev,&mode,&nherm));
  }
  free(eigxr);free(eigxi);free(eig);free(d);

  if(cyclicsymmetry){
      
       /* assembling the new eigenmodes */

      /* storage in zz: per eigenmode first the complete real part of
         the eigenvector, then the complete imaginary part */
      
      zz=NNEW(double,(long long)2*nev*neqact);
      for(l=0;l<nev;l++){
	  for(i=0;i<neqact;i++){
	      for(k=0;k<nev;k++){

                  /* real part */

		  zz[(long long)2*l*neqact+i]+=
		      (xx[2*l*nev+2*k]*z[(long long)2*k*neqact+i]
		       -xx[2*l*nev+2*k+1]*z[(long long)(2*k+1)*neqact+i]);

		  /* imaginary part */

		  zz[(long long)(2*l+1)*neqact+i]+=
		      (xx[2*l*nev+2*k]*z[(long long)(2*k+1)*neqact+i]
		       +xx[2*l*nev+2*k+1]*z[(long long)2*k*neqact+i]);
	      }
	  }
      }

      /* calculating the scalar product of all old eigenmodes with
         all new eigenmodes => nev x nev matrix */

      xmac=NNEW(double,nev*nev);
      xmaccpx=NNEW(double,4*nev*nev);
      bett=NNEW(double,nev);
      betm=NNEW(double,nev);
      FORTRAN(calcmac,(&neq[1],z,zz,&nev,xmac,xmaccpx,istartnmd,
      		       iendnmd,&nmd,&cyclicsymmetry,&neqact,bett,betm));
      FORTRAN(writemaccs,(xmac,&nev,nm));

      free(xmac);free(bett);free(betm);free(xmaccpx);
      free(z);
      
      /* normalizing the eigenmodes */
      
      z=NNEW(double,neq[1]);
      for(l=0;l<nev;l++){
	  sum=0.;
	  DMEMSET(z,0,neq[1],0.);
	  FORTRAN(op,(&neq[1],aux,&zz[l*neq[1]],z,adb,aub,icol,irow,nzl));
	  for(k=0;k<neq[1];k++){
	      sum+=zz[l*neq[1]+k]*z[k];
	  }
	  
	  sum=sqrt(sum);
	  for(k=0;k<neq[1];k++){
	      zz[l*neq[1]+k]/=sum;
	  }
      }
      free(z);

      /* calculating the mass-weighted internal products (eigenvectors are not 
         necessarily orthogonal, since the matrix of the eigenvalue problem is
         not necessarily Hermitian)
         = orthogonality matrices */

      if(mei[3]==1){
	  
	  xmr=NNEW(double,nev*nev);
	  xmi=NNEW(double,nev*nev);
	  z=NNEW(double,neq[1]);
	  zi=NNEW(double,neq[1]);
	  
	  for(l=0;l<nev;l++){
	      DMEMSET(z,0,neq[1],0.);
	      FORTRAN(op,(&neq[1],aux,&zz[l*neq[1]],z,adb,aub,icol,irow,nzl));
	      for(m=l;m<nev;m++){
		  for(k=0;k<neq[1];k++){
		      xmr[l*nev+m]+=zz[m*neq[1]+k]*z[k];
		  }
	      }
	      
	      memcpy(&zi[0],&zz[(2*l+1)*neqact],sizeof(double)*neqact);
	      for(k=0;k<neqact;k++){zi[neqact+k]=-zz[2*l*neqact+k];}
	      DMEMSET(z,0,neq[1],0.);
	      FORTRAN(op,(&neq[1],aux,zi,z,adb,aub,icol,irow,nzl));
	      for(m=l;m<nev;m++){
		  for(k=0;k<neq[1];k++){
		      xmi[l*nev+m]+=zz[m*neq[1]+k]*z[k];
	      }
	      }
	  }
	  
	  /* Hermitian part of the matrix */
	  
	  for(l=0;l<nev;l++){
	      for(m=0;m<l;m++){
		  xmr[l*nev+m]=xmr[m*nev+l];
		  xmi[l*nev+m]=-xmi[m*nev+l];
	      }
	  }
	  
	  for(l=0;l<nev;l++){
	      for(m=0;m<nev;m++){
		  printf(" %f",xmr[m*nev+l]);
	  }
	      printf("\n");
	  }
	  printf("\n");
	  for(l=0;l<nev;l++){
	      for(m=0;m<nev;m++){
		  printf(" %f",xmi[m*nev+l]);
	      }
	      printf("\n");
	  }
	  free(z);free(zi);
      }

  }else{
      
      /* no cyclic symmmetry */
      
      /* assembling the new eigenmodes */
      
      zz=NNEW(double,2*nev*neq[1]);
      for(l=0;l<nev;l++){
	  for(j=0;j<2;j++){
	      for(i=0;i<neq[1];i++){
		  for(k=0;k<nev;k++){
		      zz[(2*l+j)*neq[1]+i]+=xx[2*l*nev+2*k+j]*
                            z[(long long)k*neq[1]+i];
		  }
	      }
	  }
      }

      /* calculating the scalar product of all old eigenmodes with
         all new eigenmodes => nev x nev matrix */

      xmac=NNEW(double,nev*nev);
      xmaccpx=NNEW(double,4*nev*nev);
      bett=NNEW(double,nev);
      betm=NNEW(double,nev);
      FORTRAN(calcmac,(&neq[1],z,zz,&nev,xmac,xmaccpx,istartnmd,
      		     iendnmd,&nmd,&cyclicsymmetry,&neqact,bett,betm));
      FORTRAN(writemac,(xmac,&nev));
      free(xmac);free(bett);free(betm);free(xmaccpx);

      free(z);
      
      /* normalizing the eigenmodes */
      
      z=NNEW(double,neq[1]);
      for(l=0;l<nev;l++){
	  sum=0.;
	  
	  /* Ureal^T*M*Ureal */
	  
	  DMEMSET(z,0,neq[1],0.);
	  FORTRAN(op,(&neq[1],aux,&zz[2*l*neq[1]],z,adb,aub,icol,irow,nzl));
	  for(k=0;k<neq[1];k++){
	      sum+=zz[2*l*neq[1]+k]*z[k];
	  }
	  
	  /* Uimag^T*M*Uimag */
	  
	  DMEMSET(z,0,neq[1],0.);
	  FORTRAN(op,(&neq[1],aux,&zz[(2*l+1)*neq[1]],z,adb,aub,icol,irow,nzl));
	  for(k=0;k<neq[1];k++){
	      sum+=zz[(2*l+1)*neq[1]+k]*z[k];
	  }
	  
	  sum=sqrt(sum);
	  for(k=0;k<2*neq[1];k++){
	      zz[2*l*neq[1]+k]/=sum;
	  }
      }
      free(z);

      /* calculating the mass-weighted internal products (eigenvectors are not 
         necessarily orthogonal, since the matrix of the eigenvalue problem is
         not necessarily symmetric)
         = orthogonality matrices */

      if(mei[3]==1){
	  
	  xmr=NNEW(double,nev*nev);
	  xmi=NNEW(double,nev*nev);
	  z=NNEW(double,neq[1]);
	  
	  for(l=0;l<nev;l++){
	      sum=0.;
	      
	      /* M*Ureal */
	      
	      DMEMSET(z,0,neq[1],0.);
	      FORTRAN(op,(&neq[1],aux,&zz[2*l*neq[1]],z,adb,aub,icol,irow,nzl));
	      
	      /* Ureal^T*M*Ureal and Uimag^T*M*Ureal */
	      
	      for(m=l;m<nev;m++){
		  for(k=0;k<neq[1];k++){
		      xmr[l*nev+m]+=zz[2*m*neq[1]+k]*z[k];
		  }
		  for(k=0;k<neq[1];k++){
		      xmi[l*nev+m]-=zz[(2*m+1)*neq[1]+k]*z[k];
		  }
	      }
	      
	      /* M*Uimag */
	      
	      DMEMSET(z,0,neq[1],0.);
	      FORTRAN(op,(&neq[1],aux,&zz[(2*l+1)*neq[1]],z,adb,aub,icol,irow,nzl));
	      
	      /* Ureal^T*M*Uimag and Uimag^T*M*Uimag */
	      
	      for(m=l;m<nev;m++){
		  for(k=0;k<neq[1];k++){
		      xmr[l*nev+m]+=zz[(2*m+1)*neq[1]+k]*z[k];
		  }
		  for(k=0;k<neq[1];k++){
		      xmi[l*nev+m]+=zz[2*m*neq[1]+k]*z[k];
		  }
	      }
	  }
	  
	  /* Hermitian part of the matrix */
	  
	  for(l=0;l<nev;l++){
	      for(m=0;m<l;m++){
		  xmr[l*nev+m]=xmr[m*nev+l];
		  xmi[l*nev+m]=-xmi[m*nev+l];
	      }
	  }
	  
	  for(l=0;l<nev;l++){
	      for(m=0;m<nev;m++){
		  printf(" %f",xmr[m*nev+l]);
	      }
	      printf("\n");
	  }
	  printf("\n");
	  for(l=0;l<nev;l++){
	      for(m=0;m<nev;m++){
		  printf(" %f",xmi[m*nev+l]);
	      }
	      printf("\n");
	  }
	  free(z);
      }

  }

 /*storing new eigenmodes and eigenvalues to *.eig-file for later use in 
 steady states dynamic analysis*/

  if(mei[3]==1){

      nherm=0;
      
      if(!cyclicsymmetry){
	  if((f1=fopen(fneig,"wb"))==NULL){
	      printf("*ERROR in complexfreq: cannot open eigenvalue file for writing...");
	      exit(0);
	  }
	  
	  /* storing a zero as indication that this was not a
	     cyclic symmetry calculation */
	  
	  if(fwrite(&cyclicsymmetry,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the cyclic symmetry flag to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* not Hermitian */
	  
	  if(fwrite(&nherm,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the Hermitian flag to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* storing the number of eigenvalues */
	  
	  if(fwrite(&nev,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the number of eigenfrequencies to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* the eigenfrequencies are stored as (radians/time)**2
	     squaring the complexe eigenvalues first */
	  
	  eigx=NNEW(double,2*nev);
	  for(i=0;i<nev;i++){
	      eigx[2*i]=eigxx[2*i]*eigxx[2*i]-eigxx[2*i+1]*eigxx[2*i+1];
	      eigx[2*i+1]=2.*eigxx[2*i]*eigxx[2*i+1];
	  }
	  
	  if(fwrite(eigx,sizeof(double),2*nev,f1)!=2*nev){
	      printf("*ERROR in complexfreq saving the eigenfrequencies to the eigenvalue file...");
	      exit(0);
	  }
	  
	  free(eigx);
	  
	  /* storing the stiffness matrix */
	  
	  if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	      printf("*ERROR in complexfreq saving the diagonal of the stiffness matrix to the eigenvalue file...");
	      exit(0);
	  }
	  if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	      printf("*ERROR in complexfreq saving the off-diagonal entries of the stiffness matrix to the eigenvalue file...");
	      exit(0);
	  }
	  free(ad);free(au);
	  
	  /* storing the mass matrix */
	  
	  if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
	      printf("*ERROR in complexfreq saving the diagonal of the mass matrix to the eigenvalue file...");
	      exit(0);
	  }
	  if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	      printf("*ERROR in complexfreq saving the off-diagonal entries of the mass matrix to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* storing the eigenvectors */
	  
	  lfin=0;
	  lint=0;
	  for(j=0;j<nev;++j){
	      lint=lfin;
	      lfin=lfin+2*neq[1];
	      if(fwrite(&zz[lint],sizeof(double),2*neq[1],f1)!=2*neq[1]){
		  printf("*ERROR in complexfreq saving the eigenvectors to the eigenvalue file...");
		  exit(0);
	      }
	  }
	  
	  /* storing the orthogonality matrices */
	  
	  if(fwrite(xmr,sizeof(double),nev*nev,f1)!=nev*nev){
	      printf("*ERROR in complexfreq saving the real orthogonality matrix to the eigenvalue file...");
	      exit(0);
	  }
	  
	  if(fwrite(xmi,sizeof(double),nev*nev,f1)!=nev*nev){
	      printf("*ERROR in complexfreq saving the imaginary orthogonality matrix to the eigenvalue file...");
	      exit(0);
	  }
	  
      }else{
	  
	  /* opening .eig file for writing */
	  
	  if((f1=fopen(fneig,"wb"))==NULL){
	      printf("*ERROR in complexfreq: cannot open eigenvalue file for writing...");
	      exit(0);
	  }
	  /* storing a one as indication that this was a
	     cyclic symmetry calculation */
	  
	  if(fwrite(&cyclicsymmetry,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the cyclic symmetry flag to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* not Hermitian */
	  
	  if(fwrite(&nherm,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the Hermitian flag to the eigenvalue file...");
	      exit(0);
	  }
	  
	  if(fwrite(&nm[0],sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the nodal diameter to the eigenvalue file...");
	      exit(0);
	  }
	  if(fwrite(&nev,sizeof(int),1,f1)!=1){
	      printf("*ERROR in complexfreq saving the number of eigenvalues to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* the eigenfrequencies are stored as (radians/time)**2
	     squaring the complexe eigenvalues first */
	  
	  eigx=NNEW(double,2*nev);
	  for(i=0;i<nev;i++){
	      eigx[2*i]=eigxx[2*i]*eigxx[2*i]-eigxx[2*i+1]*eigxx[2*i+1];
	      eigx[2*i+1]=2.*eigxx[2*i]*eigxx[2*i+1];
	  }
	  
	  if(fwrite(eigx,sizeof(double),2*nev,f1)!=2*nev){
	      printf("*ERROR in complexfreq saving the eigenfrequencies to the eigenvalue file...");
	      exit(0);
	  }
	  
	  free(eigx);
	  
	  /* storing the mass matrix */
	  
	  if(fwrite(adb,sizeof(double),*neq,f1)!=*neq){
	      printf("*ERROR in complexfreq saving the diagonal of the mass matrix to the eigenvalue file...");
	      exit(0);
	  }
	  if(fwrite(aub,sizeof(double),*nzs,f1)!=*nzs){
	      printf("*ERROR in complexfreq saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* storing the eigenvectors */
	  
	  lfin=0;
	  for(j=0;j<nev;++j){
	      lint=lfin;
	      lfin=lfin+neq[1];
	      if(fwrite(&zz[lint],sizeof(double),neq[1],f1)!=neq[1]){
		  printf("*ERROR in complexfreq saving the eigenvectors to the eigenvalue file...");
		  exit(0);
	      }
	  }
	  
	  /* storing the orthogonality matrices */
	  
	  if(fwrite(xmr,sizeof(double),nev*nev,f1)!=nev*nev){
	      printf("*ERROR in complexfreq saving the real orthogonality matrix to the eigenvalue file...");
	      exit(0);
	  }
	  
	  if(fwrite(xmi,sizeof(double),nev*nev,f1)!=nev*nev){
	      printf("*ERROR in complexfreq saving the imaginary orthogonality matrix to the eigenvalue file...");
	      exit(0);
	  }
      }
      free(adb);free(aub);free(xmr);free(xmi);
      
      fclose(f1);
  }

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,2*mt**nk);
  fn=NNEW(double,2*mt**nk);
  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
     (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
     (strcmp1(&filab[1044],"ERR ")==0)) 
      stn=NNEW(double,12**nk);

  if((strcmp1(&filab[261],"E   ")==0)||(strcmp1(&filab[2523],"MAXE")==0)) 
      een=NNEW(double,12**nk);
  if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,2**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,12**nk);

  inum=NNEW(int,*nk);
  stx=NNEW(double,2*6*mi[0]**ne);
  
  coefmpcnew=NNEW(double,*mpcend);

  cot=NNEW(double,3**nk*ngraph);
  if(*ntrans>0){inotrt=NNEW(int,2**nk*ngraph);}
  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0))

// real and imaginary part of the displacements

    vt=NNEW(double,2*mt**nk*ngraph);
  if(strcmp1(&filab[87],"NT  ")==0)
    t1t=NNEW(double,*nk*ngraph);
  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)||
     (strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0))

// real and imaginary part of the stresses

    stnt=NNEW(double,2*6**nk*ngraph);
  if(strcmp1(&filab[261],"E   ")==0) 
      eent=NNEW(double,2*6**nk*ngraph);
  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0))

// real and imaginary part of the forces

    fnt=NNEW(double,2*mt**nk*ngraph);
  if(strcmp1(&filab[522],"ENER")==0)
    enernt=NNEW(double,*nk*ngraph);
  if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0))
    stxt=NNEW(double,2*6*mi[0]**ne*ngraph);
  if(strcmp1(&filab[2697],"ME  ")==0) 
      emnt=NNEW(double,2*6**nk*ngraph);

  kont=NNEW(int,*nkon*ngraph);
  ipkont=NNEW(int,*ne*ngraph);
  for(l=0;l<*ne*ngraph;l++)ipkont[l]=-1;
  lakont=NNEW(char,8**ne*ngraph);
  inumt=NNEW(int,*nk*ngraph);
  ielmatt=NNEW(int,mi[2]**ne*ngraph);

  nkt=ngraph**nk;
  net=ngraph**ne;

  /* copying the coordinates of the first sector */

  for(l=0;l<3**nk;l++){cot[l]=co[l];}
  if(*ntrans>0){for(l=0;l<*nk;l++){inotrt[2*l]=inotr[2*l];}}
  for(l=0;l<*nkon;l++){kont[l]=kon[l];}
  for(l=0;l<*ne;l++){ipkont[l]=ipkon[l];}
  for(l=0;l<8**ne;l++){lakont[l]=lakon[l];}
  for(l=0;l<*ne;l++){ielmatt[mi[2]*l]=ielmat[mi[2]*l];}

  /* generating the coordinates for the other sectors */

  if(cyclicsymmetry){

    icntrl=1;
    
    FORTRAN(rectcyl,(cot,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));
    
    for(jj=0;jj<*mcs;jj++){
      is=cs[17*jj+4];
      for(i=1;i<is;i++){
	
	theta=i*2.*pi/cs[17*jj];
	
	for(l=0;l<*nk;l++){
	  if(inocs[l]==jj){
	    cot[3*l+i*3**nk]=cot[3*l];
	    cot[1+3*l+i*3**nk]=cot[1+3*l]+theta;
	    cot[2+3*l+i*3**nk]=cot[2+3*l];
	    if(*ntrans>0){inotrt[2*l+i*2**nk]=inotrt[2*l];}
	  }
	}
	for(l=0;l<*nkon;l++){kont[l+i**nkon]=kon[l]+i**nk;}
	for(l=0;l<*ne;l++){
	  if(ielcs[l]==jj){
	    if(ipkon[l]>=0){
	      ipkont[l+i**ne]=ipkon[l]+i**nkon;
	      ielmatt[mi[2]*(l+i**ne)]=ielmat[mi[2]*l];
	      for(l1=0;l1<8;l1++){
		l2=8*l+l1;
		lakont[l2+i*8**ne]=lakon[l2];
	      }
	    }
	  }
	}
      }
    }
    
    icntrl=-1;
    
    FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		     &imag,mi,emnt));
  }

  /* check that the tensor fields which are extrapolated from the
     integration points are requested in global coordinates */

  if(strcmp1(&filab[174],"S   ")==0){
      if((strcmp1(&filab[179],"L")==0)&&(*norien>0)){
      printf("\n*WARNING in complexfreq: element fields in cyclic symmetry calculations\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[179],"G",1);
    }
  }

  if(strcmp1(&filab[261],"E   ")==0){
      if((strcmp1(&filab[266],"L")==0)&&(*norien>0)){
      printf("\n*WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[266],"G",1);
    }
  }

  if(strcmp1(&filab[1479],"PHS ")==0){
      if((strcmp1(&filab[1484],"L")==0)&&(*norien>0)){
      printf("\n*WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[1484],"G",1);
    }
  }

  if(strcmp1(&filab[1653],"MAXS")==0){
      if((strcmp1(&filab[1658],"L")==0)&&(*norien>0)){
      printf("\n*WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[1658],"G",1);
    }
  }   

  if(strcmp1(&filab[2523],"MAXE")==0){
      if((strcmp1(&filab[2528],"L")==0)&&(*norien>0)){
      printf("\n*WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[1658],"G",1);
    }
  }   

  /* allocating fields for magnitude and phase information of
     displacements and stresses */

  if(strcmp1(&filab[870],"PU")==0){
    vr=NNEW(double,mt*nkt);
    vi=NNEW(double,mt*nkt);
  }

  if(strcmp1(&filab[1479],"PHS")==0){
    stnr=NNEW(double,6*nkt);
    stni=NNEW(double,6*nkt);
  }

  if(strcmp1(&filab[1566],"MAXU")==0){
    vmax=NNEW(double,4*nkt);
  }

  if(strcmp1(&filab[1653],"MAXS")==0){
    stnmax=NNEW(double,7*nkt);
  }

  if(strcmp1(&filab[2523],"MAXE")==0){
    eenmax=NNEW(double,7*nkt);
  }

  /* storing the results */

  if(!cyclicsymmetry) (neq[1])*=2;

  lfin=0;
  for(j=0;j<nev;++j){
    lint=lfin;
    lfin=lfin+neq[1];

    /* calculating the cosine and sine */

    for(k=0;k<*mcs;k++){
	theta=nm[j]*2.*pi/cs[17*k];
	cs[17*k+14]=cos(theta);
	cs[17*k+15]=sin(theta);
    }

    if(*nprint>0)FORTRAN(writehe,(&j));

    eei=NNEW(double,6*mi[0]**ne);
    if(*nener==1){
	stiini=NNEW(double,6*mi[0]**ne);
	enerini=NNEW(double,mi[0]**ne);}

    DMEMSET(v,0,2*mt**nk,0.);

    for(k=0;k<neq[1];k+=neq[1]/2){

      for(i=0;i<6*mi[0]**ne;i++){eme[i]=0.;}

      if(k==0) {kk=0;kkv=0;kk6=0;kkx=0;if(*nprint>0)FORTRAN(writere,());}
      else {kk=*nk;kkv=mt**nk;kk6=6**nk;kkx=6*mi[0]**ne;
            if(*nprint>0)FORTRAN(writeim,());}

    /* generating the cyclic MPC's (needed for nodal diameters
       different from 0 */

      for(i=0;i<*nmpc;i++){
	index=ipompc[i]-1;
        /* check whether thermal mpc */
	if(nodempc[3*index+1]==0) continue;
	coefmpcnew[index]=coefmpc[index];
	while(1){
	  index=nodempc[3*index+2];
	  if(index==0) break;
	  index--;

	  icomplex=0;
	  inode=nodempc[3*index];
	  if(strcmp1(&labmpc[20*i],"CYCLIC")==0){
            icomplex=atoi(&labmpc[20*i+6]);}
	  else if(strcmp1(&labmpc[20*i],"SUBCYCLIC")==0){
            for(ij=0;ij<*mcs;ij++){
              lprev=cs[ij*17+13];
              ilength=cs[ij*17+3];
              FORTRAN(nident,(&ics[lprev],&inode,&ilength,&id));
              if(id!=0){
                if(ics[lprev+id-1]==inode){icomplex=ij+1;break;}
              }
            }
	  }

          if(icomplex!=0){
	    idir=nodempc[3*index+1];
	    idof=nactdof[mt*(inode-1)+idir]-1;
            if(idof==-1){xreal=1.;ximag=1.;}
            else{xreal=zz[lint+idof];ximag=zz[lint+idof+neq[1]/2];}
            if(k==0) {
              if(fabs(xreal)<1.e-30)xreal=1.e-30;
	      coefmpcnew[index]=coefmpc[index]*
 		(cs[17*(icomplex-1)+14]+ximag/xreal*cs[17*(icomplex-1)+15]);}
	    else {
              if(fabs(ximag)<1.e-30)ximag=1.e-30;
              coefmpcnew[index]=coefmpc[index]*
		(cs[17*(icomplex-1)+14]-xreal/ximag*cs[17*(icomplex-1)+15]);}
          }
	  else{coefmpcnew[index]=coefmpc[index];}
	}
      }

      if(*iperturb==0){
	results(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,
            &stx[kkx],elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t0,ithermal,
	    prestr,iprestr,filab,eme,&emn[kk6],&een[kk6],iperturb,
	    f,&fn[kkv],nactdof,&iout,qa,vold,&zz[lint+k],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);}
      else{
	results(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,
            &stx[kkx],elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,&emn[kk6],&een[kk6],iperturb,
            f,&fn[kkv],nactdof,&iout,qa,vold,&zz[lint+k],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xload,xloadold,&icfd,inomat);
      }

    }
    free(eei);
    if(*nener==1){free(stiini);free(enerini);}
   
    /* changing the basic results into cylindrical coordinates
       (only for cyclic symmetry structures */

    for(l=0;l<*nk;l++){inumt[l]=inum[l];}

    if(cyclicsymmetry){
      icntrl=2;imag=1;      
      FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));
    }

    /* copying the basis results (real and imaginary) into
       a new field */

    if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
      for(l=0;l<mt**nk;l++){vt[l]=v[l];}
      for(l=0;l<mt**nk;l++){vt[l+mt**nk*ngraph]=v[l+mt**nk];}}
    if(strcmp1(&filab[87],"NT  ")==0)
      for(l=0;l<*nk;l++){t1t[l]=t1[l];};
    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	for(l=0;l<6**nk;l++){stnt[l]=stn[l];}
	for(l=0;l<6**nk;l++){stnt[l+6**nk*ngraph]=stn[l+6**nk];}}
    if(strcmp1(&filab[261],"E   ")==0){
	for(l=0;l<6**nk;l++){eent[l]=een[l];};
	for(l=0;l<6**nk;l++){eent[l+6**nk*ngraph]=een[l+6**nk];}}
    if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
      for(l=0;l<mt**nk;l++){fnt[l]=fn[l];}
      for(l=0;l<mt**nk;l++){fnt[l+mt**nk*ngraph]=fn[l+mt**nk];}}
    if(strcmp1(&filab[522],"ENER")==0)
      for(l=0;l<*nk;l++){enernt[l]=enern[l];};
    if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)){
      for(l=0;l<6*mi[0]**ne;l++){stxt[l]=stx[l];}
      for(l=0;l<6*mi[0]**ne;l++){stxt[l+6*mi[0]**ne*ngraph]=stx[l+6*mi[0]**ne];}}
    if(strcmp1(&filab[2697],"ME  ")==0){
	for(l=0;l<6**nk;l++){emnt[l]=emn[l];};
	for(l=0;l<6**nk;l++){emnt[l+6**nk*ngraph]=emn[l+6**nk];}}
    
    /* mapping the results to the other sectors
       (only for cyclic symmetric structures */

    if(cyclicsymmetry){

      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj+4];
	lprev=cs[17*jj+13];
	for(i=1;i<is;i++){
	  
	  for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
	  
	  theta=i*nm[j]*2.*pi/cs[17*jj];
	  ctl=cos(theta);
	  stl=sin(theta);
	  
	  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      vt[l+mt**nk*i]=v[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  vt[l+mt**nk*i]=ctl*v[l]-stl*v[l+mt**nk];
		}
		
	      }
	    }
	  }
	  
	  /* imaginary part of the displacements in cylindrical
	     coordinates */
	  
	  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      vt[l+mt**nk*(i+ngraph)]=v[l+mt**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  vt[l+mt**nk*(i+ngraph)]=stl*v[l]+ctl*v[l+mt**nk];
		}
	      }
	    }
	  }
	  
	  if(strcmp1(&filab[87],"NT  ")==0){
	    for(l=0;l<*nk;l++){
	      if(inocs[l]==jj) t1t[l+*nk*i]=t1[l];
	    }
	  }
	  
	  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      stnt[l+6**nk*i]=stn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnt[l+6**nk*i]=ctl*stn[l]-stl*stn[l+6**nk];
		}
	      }
	    }
	  }
	  
	  /* imaginary part of the stresses in cylindrical
	     coordinates */
	  
	  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      stnt[l+6**nk*(i+ngraph)]=stn[l+6**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnt[l+6**nk*(i+ngraph)]=stl*stn[l]+ctl*stn[l+6**nk];
		}
	      }
	    }
	  }
	  
	  if(strcmp1(&filab[261],"E   ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      eent[l+6**nk*i]=een[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  eent[l+6**nk*i]=ctl*een[l]-stl*een[l+6**nk];
		}
	      }
	    }
	  }
        
        /* imaginary part of the strains in cylindrical
           coordinates */
        
	  if(strcmp1(&filab[261],"E   ")==0){
	      for(l1=0;l1<*nk;l1++){
		  if(inocs[l1]==jj){
		      
		      /* check whether node lies on axis */
		      
		      ml1=-l1-1;
		      FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		      if(id!=0){
			  if(ics[lprev+id-1]==ml1){
			      for(l2=0;l2<6;l2++){
				  l=6*l1+l2;
				  eent[l+6**nk*(i+ngraph)]=een[l+6**nk];
			      }
			      continue;
			  }
		      }
		      for(l2=0;l2<6;l2++){
			  l=6*l1+l2;
			  eent[l+6**nk*(i+ngraph)]=stl*een[l]+ctl*een[l+6**nk];
		      }
		  }
	      }
	  }
	  
	  if(strcmp1(&filab[2697],"ME  ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      emnt[l+6**nk*i]=emn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  emnt[l+6**nk*i]=ctl*emn[l]-stl*emn[l+6**nk];
		}
	      }
	    }
	  }
        
        /* imaginary part of the mechanical strains in cylindrical
           coordinates */
        
	  if(strcmp1(&filab[2697],"ME  ")==0){
	      for(l1=0;l1<*nk;l1++){
		  if(inocs[l1]==jj){
		      
		      /* check whether node lies on axis */
		      
		      ml1=-l1-1;
		      FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		      if(id!=0){
			  if(ics[lprev+id-1]==ml1){
			      for(l2=0;l2<6;l2++){
				  l=6*l1+l2;
				  emnt[l+6**nk*(i+ngraph)]=emn[l+6**nk];
			      }
			      continue;
			  }
		      }
		      for(l2=0;l2<6;l2++){
			  l=6*l1+l2;
			  emnt[l+6**nk*(i+ngraph)]=stl*emn[l]+ctl*emn[l+6**nk];
		      }
		  }
	      }
	  }
	  
	  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
		
		/* check whether node lies on axis */
		
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      fnt[l+mt**nk*i]=fn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  fnt[l+mt**nk*i]=ctl*fn[l]-stl*fn[l+mt**nk];
		}
	      }
	    }
	  }
        
        /* imaginary part of the forces in cylindrical
           coordinates */

	  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
	      for(l1=0;l1<*nk;l1++){
		  if(inocs[l1]==jj){
		      
		      /* check whether node lies on axis */
		      
		      ml1=-l1-1;
		      FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		      if(id!=0){
			  if(ics[lprev+id-1]==ml1){
			      for(l2=0;l2<4;l2++){
				  l=mt*l1+l2;
				  fnt[l+mt**nk*(i+ngraph)]=fn[l+mt**nk];
			      }
			      continue;
			  }
		      }
		      for(l2=0;l2<4;l2++){
			  l=mt*l1+l2;
			  fnt[l+mt**nk*(i+ngraph)]=stl*fn[l]+ctl*fn[l+mt**nk];
		      }
		  }
	      }
	  }
	  
	  if(strcmp1(&filab[522],"ENER")==0){
	    for(l=0;l<*nk;l++){
	      if(inocs[l]==jj) enernt[l+*nk*i]=0.;
	    }
	  }
	}
      }
      
      icntrl=-2;imag=0;
      
      FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		       &imag,mi,emnt));
      
      FORTRAN(rectcylvi,(cot,&vt[mt**nk*ngraph],&fnt[mt**nk*ngraph],
                         &stnt[6**nk*ngraph],qfnt,&eent[6**nk*ngraph],
                         cs,&nkt,&icntrl,t,filab,&imag,mi,&emnt[6**nk*ngraph]));
      
    }

    /* determining magnitude and phase angle for the displacements */

    if(strcmp1(&filab[870],"PU")==0){
      for(l1=0;l1<nkt;l1++){
	for(l2=0;l2<4;l2++){
	  l=mt*l1+l2;
	  vreal=vt[l];
	  vimag=vt[l+mt**nk*ngraph];
	  vr[l]=sqrt(vreal*vreal+vimag*vimag);
	  if(fabs(vreal)<1.e-10){
	    if(vimag>0){vi[l]=90.;}
	    else{vi[l]=-90.;}
	  }
	  else{
	    vi[l]=atan(vimag/vreal)*constant;
	    if(vreal<0) vi[l]+=180.;
	  }
	}
      }
    }

    /* determining magnitude and phase for the stress */

    if(strcmp1(&filab[1479],"PHS")==0){
      for(l1=0;l1<nkt;l1++){
	for(l2=0;l2<6;l2++){
	  l=6*l1+l2;
	  stnreal=stnt[l];
	  stnimag=stnt[l+6**nk*ngraph];
	  stnr[l]=sqrt(stnreal*stnreal+stnimag*stnimag);
	  if(fabs(stnreal)<1.e-10){
	    if(stnimag>0){stni[l]=90.;}
	    else{stni[l]=-90.;}
	  }
	  else{
	    stni[l]=atan(stnimag/stnreal)*constant;
	    if(stnreal<0) stni[l]+=180.;
	  }
	}
      }
    }

    ++*kode;

    /* storing the real part of the eigenfrequencies in freq */

    freq=eigxx[2*j]/6.283185308;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40*net);ipneigh=NNEW(int,nkt);
    }
    frd(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,
	    kode,filab,eent,t1t,fnt,&freq,epn,ielmatt,matname,enernt,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&j,&nm[j],trab,inotrt,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stxt,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&net,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emnt,
	    thicke,jobnamec,output,qfx);
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
  }

  free(xstiff);if(*nbody>0) free(ipobody);
  free(cstr);free(zz);free(eigxx);free(xx);

  if(cyclicsymmetry){
      free(istartnmd);free(iendnmd);
  }else{
      (neq[1])/=2;
  }

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

