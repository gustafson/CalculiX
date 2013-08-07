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

#ifdef ARPACK

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
#ifdef MATRIXSTORAGE
   #include "matrixstorage.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif

void arpackcs(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int *icol, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,int *mei, double *fei,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double **xstatep, int *npmat_, char *matname, int *mi,
	     int *ics, double *cs, int *mpcend, int *ncmat_,
             int *nstate_, int *mcs, int *nkon, double **enerp,
             char *jobnamec, char *output, char *set, int *nset, 
             int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
             char *cbody, int *ibody, double *xbody, int *nbody,
	     int *nevtot, double *thicke, int *nslavs, double *tietol, 
	      int *mpcinfo,int *ntie,int *istep,
	     int *nnn,char *tieset){

  /* calls the Arnoldi Package (ARPACK) for cyclic symmetry calculations */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A",*lakont=NULL,
      description[13]="            ",fneig[132]="",
      lakonl[2]=" \0";

  int *inum=NULL,k,ido,ldz,iparam[11],ipntr[14],lworkl,idir,nherm=1,
    info,rvec=1,*select=NULL,lfin,j,lint,iout=1,nm,index,inode,id,i,idof,
    ielas=0,icmd=0,kk,l,nkt,icntrl,*kont=NULL,*ipkont=NULL,*inumt=NULL,
    *ielmatt=NULL,net,imag,icomplex,kkv,kk6,iinc=1,nev,ncv,
    mxiter,lprev,ilength,ij,i1,i2,iel,ielset,node,indexe,nope,ml1,
    *inocs=NULL,*ielcs=NULL,jj,l1,l2,ngraph,is,jrow,*ipobody=NULL,
    *inotrt=NULL,symmetryflag=0,inputformat=0,inewton=0,ifreebody,
    mass=1, stiffness=1, buckling=0, rhsi=0, intscheme=0,*ncocon=NULL,
    coriolis=0,iworsttime,l3,iray,mt,kkx,im,ne0,*integerglob=NULL,
    *nshcon=NULL,one=1,irenewxstate,ncont,*itietri=NULL,neq2,*nmc=NULL,
    *koncont=NULL,ismallsliding=0,*itiefac=NULL,*islavsurf=NULL,
    *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,mortar=0,
    *imastop=NULL,*iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,ifacecount,
    mpcfree,memmpc_,icascade,maxlenmpc,nkon0,iit=-1,*irow=NULL,nasym=0,
    kmax1,kmax2,icfd=0,*inomat=NULL;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,*vr=NULL,
    *workl=NULL,*aux=NULL,*d=NULL,sigma,*temp_array=NULL,*vini=NULL,
    *een=NULL,cam[5],*f=NULL,*fn=NULL,qa[3],*fext=NULL,*emn=NULL,
    *epn=NULL,*stiini=NULL,*fnr=NULL,*fni=NULL,fnreal,fnimag,*emeini=NULL,
    *xstateini=NULL,theta,pi,*coefmpcnew=NULL,*xstiff=NULL,*vi=NULL,
    *vt=NULL,*fnt=NULL,*stnt=NULL,*eent=NULL,*cot=NULL,t[3],ctl,stl,
    *t1t=NULL,freq,*stx=NULL,*enern=NULL,*enernt=NULL,*xstaten=NULL,
    *eei=NULL,*enerini=NULL,*cocon=NULL,*qfx=NULL,*qfn=NULL,*qfnt=NULL,
    tol,fmin,fmax,xreal,ximag,*cgr=NULL,*xloadold=NULL,reltime,constant,
    vreal,vimag,*stnr=NULL,*stni=NULL,stnreal,stnimag,*vmax=NULL,
    *stnmax=NULL,vl[4],stnl[6],dd,v1,v2,v3,bb,cc,al[3],cm,cn,tt,
    worstpsmax,vray[3],worstumax,p1[3],p2[3],q[3],tan[3],*springarea=NULL,
    *stxt=NULL,*eenmax=NULL,eenl[6],*emnt=NULL,*xnormastface=NULL,
    *doubleglob=NULL,*shcon=NULL,*cg=NULL,*straight=NULL,
    *xmastnor=NULL,*areaslav=NULL,*dc=NULL,*di=NULL,*xnoels=NULL,
    *workev=NULL,sigmaz[2]={0.,0.},*temp_array2=NULL,*ener=NULL,*xstate=NULL,
    *rwork=NULL,sigmai=0,amp,ampmax,*zstorage=NULL;

  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

  int *ipneigh=NULL,*neigh=NULL;
  
#ifdef SGI
  int token;
#endif
  
  irow=*irowp;ener=*enerp;xstate=*xstatep;
  
  mt=mi[1]+1;
  pi=4.*atan(1.);
  constant=180./pi;
  
  /* copying the frequency parameters */
  
  nev=mei[0];
  ncv=mei[1];
  mxiter=mei[2];
  tol=fei[0];
  fmin=2*pi*fei[1];
  fmax=2*pi*fei[2];
  
  /* assigning the body forces to the elements */ 
  
  if(*nbody>0){
      ifreebody=*ne+1;
      ipobody=NNEW(int,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,int,2*(*ne+ifreebody));
      }
      RENEW(ipobody,int,2*(ifreebody-1));
      if(inewton==1){
	  printf("*ERROR in arpackcs: generalized gravity loading is not allowed in frequency calculations");
	  FORTRAN(stop,());
      }
  }
  
  ne0=*ne;nkon0=*nkon;
  
  /* contact conditions */
  
  if(*iperturb!=0){
      
      memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
      maxlenmpc=mpcinfo[3];
      
      if(*nslavs==0){irenewxstate=1;}else{irenewxstate=0;}
      inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	      lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
	      &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
	      &mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,&ifacecount,
	      nmpc,&mpcfree,&memmpc_,
	      &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	      iperturb,ikboun,nboun,co,istep,&xnoels);
      
      if(ncont!=0){
	  
	  if(*nener==1){
	      RENEW(ener,double,mi[0]*(*ne+*nslavs)*2);
	  }
	  RENEW(ipkon,int,*ne+*nslavs);
	  RENEW(lakon,char,8*(*ne+*nslavs));
	  
	  if(*norien>0){
	      RENEW(ielorien,int,mi[2]*(*ne+*nslavs));
	      for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
	  }
	  RENEW(ielmat,int,mi[2]*(*ne+*nslavs));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;
	  cg=NNEW(double,3*ncont);
	  printf("arpackcs1 cg=%e\n",cg[0]);
	  straight=NNEW(double,16*ncont);
	  
	  /* 11 instead of 10: last position is reserved for the
	     local contact spring element number; needed as
	     pointer into springarea */
	  
	  RENEW(kon,int,*nkon+11**nslavs);
	  if((irenewxstate==1)&&(*nslavs!=0)){
	      RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	      for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
		  xstate[k]=0.;
	      }
	  }
	  xmastnor=NNEW(double,3*nmastnode[*ntie]);
	  xnormastface=NNEW(double,3*9**nslavs);
	  areaslav=NNEW(double,ifacecount);
	  springarea=NNEW(double,2**nslavs);
	  
          /* generating contact spring elements */
	  
	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
		  co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
		  &ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
		  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
		  iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
		  itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
		  imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
		  xnoels);
	  printf("arpackcs2 cg=%e\n",cg[0]);
	  
          /* determining the structure of the stiffness/mass matrix */
	  
	  remastructar(ipompc,&coefmpc,&nodempc,nmpc,
		       &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		       labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		       kon,ipkon,lakon,ne,nnn,nactdof,icol,jq,&irow,isolver,
		       neq,nzs,nmethod,ithermal,iperturb,&mass,mi,ics,cs,
		       mcs);
      }
  }
  
  /* field for initial values of state variables (needed if
     previous static step was viscoplastic and for contact */
  
  if(*nstate_!=0){
      xstateini=NNEW(double,*nstate_*mi[0]*(ne0+*nslavs));
      for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
	  xstateini[k]=xstate[k];
      }
  }
  
  /* determining the internal forces and the stiffness coefficients */
  
  f=NNEW(double,neq[1]);
  
  /* allocating a field for the stiffness matrix */
  
  xstiff=NNEW(double,(long long)27*mi[0]**ne);
  
  iout=-1;
  v=NNEW(double,mt**nk);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  fn=NNEW(double,mt**nk);
  stx=NNEW(double,6*mi[0]**ne);
  inum=NNEW(int,*nk);
  if(*iperturb==0){
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t0,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
	      sideload,xload,xloadold,&icfd,inomat);
  }else{
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
	      sideload,xload,xloadold,&icfd,inomat);
  }
  free(f);free(v);free(fn);free(stx);free(inum);
  iout=1;
  
  /* for the frequency analysis linear strain and elastic properties
     are used */
  
  iperturb[1]=0;ielas=1;
  
  /* determining the maximum number of sectors to be plotted */
  
  ngraph=1;
  for(j=0;j<*mcs;j++){
      if(cs[17*j+4]>ngraph) ngraph=cs[17*j+4];
  }
  
  /* assigning nodes and elements to sectors */
  
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
	      else if (strcmp1(&lakon[8*iel+3],"6")==0)nope=6;
	      else if (strcmp1(&lakon[8*iel],"ES")==0){
		  lakonl[0]=lakon[8*iel+7];
		  nope=atoi(lakonl)+1;}
	      else continue;
	      
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
  
  /* loop over the nodal diameters */
  
  for(nm=cs[1];nm<=cs[2];++nm){
      
      ad=NNEW(double,neq[1]);
      au=NNEW(double,nzs[1]);
      
      adb=NNEW(double,neq[1]);
      aub=NNEW(double,nzs[1]);
      
      fext=NNEW(double,neq[1]);
      if(*iperturb==0){
	  FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		ad,au,fext,nactdof,icol,jq,irow,&neq[1],nzl,nmethod,
		ikmpc,ilmpc,ikboun,ilboun,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,
		t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
		nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		xstiff,npmat_,&dtime,matname,mi,
		ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,&rhsi,
		&intscheme,mcs,&coriolis,ibody,xloadold,&reltime,ielcs,veold,
		springarea,thicke,xnormastface,integerglob,doubleglob,
		tieset,istartset,iendset,ialset,ntie,&nasym));
      }
      else{
	  FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		ad,au,fext,nactdof,icol,jq,irow,&neq[1],nzl,nmethod,
		ikmpc,ilmpc,ikboun,ilboun,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,
		t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
		nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		xstiff,npmat_,&dtime,matname,mi,
		ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,&rhsi,
		&intscheme,mcs,&coriolis,ibody,xloadold,&reltime,ielcs,veold,
		springarea,thicke,xnormastface,integerglob,doubleglob,
		tieset,istartset,iendset,ialset,ntie,&nasym));
	  
	  if(nasym==1){
	      RENEW(au,double,nzs[2]+nzs[1]);
	      RENEW(aub,double,nzs[2]+nzs[1]);
	      symmetryflag=2;
	      inputformat=1;

	      FORTRAN(mafillsmcsas,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
		xboun,nboun,
		ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		ad,au,fext,nactdof,icol,jq,irow,&neq[1],nzl,nmethod,
		ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
		nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		t0,t1,ithermal,prestr,
		iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
		nplicon,plkcon,nplkcon,xstiff,npmat_,&dtime,
		matname,mi,ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,
		&rhsi,&intscheme,mcs,&coriolis,ibody,xloadold,&reltime,ielcs,
		veold,springarea,thicke,xnormastface,integerglob,doubleglob,
		tieset,istartset,iendset,ialset,ntie,&nasym,nstate_,xstateini,
		xstate));
	      
	  }
      }
      
      free(fext);
      
      if(*nmethod==0){
	  
	  /* error occurred in mafill: storing the geometry in frd format */
	  
	  ++*kode;time=0.;
	  inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	  }
	  
	  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	      kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&j,&nm,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);
	  
	  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
	  free(inum);FORTRAN(stop,());
	  
      }
      
      /* LU decomposition of the left hand matrix */
      
      if(nasym==1){sigma=0.;}else{sigma=1.;}
      
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,&symmetryflag,
			 &inputformat,&nzs[2]);
#else
	  printf("*ERROR in arpackcs: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,token);
#else
	  printf("*ERROR in arpackcs: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],nzs);
#else
	  printf("*ERROR in arpackcs: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==6){
#ifdef MATRIXSTORAGE
	  matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
			ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
			lakon,kon,ne,mei,nboun,nmpc,cs,mcs);
#else
	  printf("*ERROR in arpack: the MATRIXSTORAGE library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,
			 &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	  printf("*ERROR in arpackcs: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      
      free(au);free(ad);
      
      /* calculating the ei
	 genvalues and eigenmodes */
      
      printf(" Calculating the eigenvalues and the eigenmodes\n");
      
      ido=0;
      ldz=neq[1];
      for(k=0;k<11;k++)iparam[k]=0;
      iparam[0]=1;
      iparam[2]=mxiter;
      iparam[3]=1;
      iparam[6]=3;
      
      info=0;
      
      resid=NNEW(double,neq[1]);
      z=NNEW(double,(long long)ncv*neq[1]);
      workd=NNEW(double,3*neq[1]);
      
      if(nasym==1){
	  lworkl=3*ncv*(2+ncv);
	  workl=NNEW(double,lworkl);
	  FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
			  workl,&lworkl,&info));
	  /*      neq2=neq[1]/2;
		  lworkl=6*ncv*ncv+10*ncv;
		  workl=NNEW(double,lworkl);
		  rwork=NNEW(double,ncv);
		  FORTRAN(znaupd,(&ido,bmat,&neq2,which,&nev,&tol,resid,&ncv,z,&neq2,
		  iparam,ipntr,workd,workl,&lworkl,rwork,&info));
		  temp_array2=NNEW(double,neq[1]);*/
      }else{
	  lworkl=ncv*(8+ncv);
	  workl=NNEW(double,lworkl);
	  FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
			  iparam,ipntr,workd,workl,&lworkl,&info));
      }
      
      temp_array=NNEW(double,neq[1]);
      
      while((ido==-1)||(ido==1)||(ido==2)){
	  if(ido==-1){
	      if(nasym==1){
		  FORTRAN(opas,(&neq[1],aux,&workd[ipntr[0]-1],temp_array,adb,aub,icol,irow,nzl,nzs));
		  /*	    for(jrow=0;jrow<neq2;jrow++){
			    temp_array2[jrow]=workd[ipntr[0]-1+2*jrow];
			    temp_array2[neq2+jrow]=workd[ipntr[0]+2*jrow];
			    }
			    FORTRAN(opas,(&neq[1],aux,temp_array2,temp_array,adb,aub,
			    icol,irow,nzl,nzs));*/
	      }else{
		  FORTRAN(op,(&neq[1],aux,&workd[ipntr[0]-1],temp_array,adb,aub,
			      icol,irow,nzl));
	      }
	  }
	  
	  /* solve the linear equation system  */
	  
	  if((ido==-1)||(ido==1)){
	      
	      if(ido==-1){
		  if(*isolver==0){
#ifdef SPOOLES
		      spooles_solve(temp_array,&neq[1]);
#endif
		  }
		  else if(*isolver==4){
#ifdef SGI
		      sgi_solve(temp_array,token);
#endif
		  }
		  else if(*isolver==5){
#ifdef TAUCS
		      tau_solve(temp_array,&neq[1]);
#endif
		  }
		  else if(*isolver==7){
#ifdef PARDISO
		      pardiso_solve(temp_array,&neq[1],&symmetryflag);
#endif
		  }
		  /*	    if(nasym==1){
			    for(jrow=0;jrow<neq2;jrow++){
			    workd[ipntr[1]-1+2*jrow]=temp_array[jrow];
			    workd[ipntr[1]+2*jrow]=temp_array[neq2+jrow];
			    }
			    }else{*/
		  for(jrow=0;jrow<neq[1];jrow++){
		      workd[ipntr[1]-1+jrow]=temp_array[jrow];
		      //		}
		  }
	      }
	      else if(ido==1){
		  /*	    if(nasym==1){
			    for(jrow=0;jrow<neq2;jrow++){
			    temp_array2[jrow]=workd[ipntr[2]-1+2*jrow];
			    temp_array2[neq2+jrow]=workd[ipntr[2]+2*jrow];
			    }
			    if(*isolver==0){
			    #ifdef SPOOLES
			    spooles_solve(temp_array2,&neq[1]);
			    #endif
			    }
			    else if(*isolver==4){
			    #ifdef SGI
			    sgi_solve(temp_array2,token);
			    #endif
			    }
			    else if(*isolver==5){
			    #ifdef TAUCS
			    tau_solve(temp_array2,&neq[1]);
			    #endif
			    }
			    else if(*isolver==7){
			    #ifdef PARDISO
			    pardiso_solve(temp_array2,&neq[1],&symmetryflag);
			    #endif
			    }
			    for(jrow=0;jrow<neq2;jrow++){
			    workd[ipntr[1]-1+2*jrow]=temp_array2[jrow];
			    workd[ipntr[1]+2*jrow]=temp_array2[neq2+jrow];
			    }
			    }else{*/
		  if(*isolver==0){
#ifdef SPOOLES
		      spooles_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
		  }
		  else if(*isolver==4){
#ifdef SGI
		      sgi_solve(&workd[ipntr[2]-1],token);
#endif
		  }
		  else if(*isolver==5){
#ifdef TAUCS
		      tau_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
		  }
		  else if(*isolver==7){
#ifdef PARDISO
		      pardiso_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag);
#endif
		  }
		  for(jrow=0;jrow<neq[1];jrow++){
		      workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
		  }
		  //	    }
		  
	      }
	  }
	  
	  if(ido==2){
	      if(nasym==1){
		  FORTRAN(opas,(&neq[1],aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
				adb,aub,icol,irow,nzl,nzs));
		  /*	    for(jrow=0;jrow<neq2;jrow++){
			    temp_array2[jrow]=workd[ipntr[0]-1+2*jrow];
			    temp_array2[neq2+jrow]=workd[ipntr[0]+2*jrow];
			    }
			    FORTRAN(opas,(&neq[1],aux,temp_array2,temp_array,
			    adb,aub,icol,irow,nzl,nzs));
			    for(jrow=0;jrow<neq2;jrow++){
			    workd[ipntr[1]-1+2*jrow]=temp_array[jrow];
			    workd[ipntr[1]+2*jrow]=temp_array[neq2+jrow];
			    }*/
	      }else{
		  FORTRAN(op,(neq,aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
			      adb,aub,icol,irow,nzl));
	      }
	  }
	  
	  if(nasym==1){
	      FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
			      iparam,ipntr,workd,workl,&lworkl,&info));
	      /*	FORTRAN(znaupd,(&ido,bmat,&neq2,which,&nev,&tol,resid,&ncv,z,&neq2,
			iparam,ipntr,workd,workl,&lworkl,rwork,&info));*/
	  }else{
	      FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,
			      z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	  }
      }
      
/*--------------------------------------------------------------------*/
/*
  -----------
  free memory
  -----------
*/
      free(temp_array);free(temp_array2);
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
	  pardiso_cleanup(&neq[1],&symmetryflag);
#endif
      }
      
      if(info!=0){
	  printf("*ERROR in arpackcs: info=%d\n",info);
	  printf("       # of converged eigenvalues=%d\n\n",iparam[4]);
      }         
      
      select=NNEW(int,ncv);
      
      if(nasym==1){
	  d=NNEW(double,nev+1);
	  di=NNEW(double,nev+1);
	  workev=NNEW(double,3*ncv);
	  FORTRAN(dneupd,(&rvec,howmny,select,d,di,z,&ldz,&sigma,&sigmai,
			  workev,bmat,&neq[1],which,&nev,&tol,resid,
			  &ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	  free(workev);
	  dc=NNEW(double,2*nev);
	  
	  /* storing as complex number and taking the square root */
	  
	  for(j=0;j<nev;j++){
	      dc[2*j]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])+d[j])/sqrt(2.);
	      dc[2*j+1]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])-d[j])/sqrt(2.);
	      if(di[j]<0.) dc[2*j+1]=-dc[2*j+1];
	  }
	  FORTRAN(writeevcomplex,(dc,&nev,&fmin,&fmax));
	  free(di);free(dc);
	  /*      d=NNEW(double,2*nev+2);
		  workev=NNEW(double,4*ncv);
		  FORTRAN(zneupd,(&rvec,howmny,select,d,z,&neq2,sigmaz,workev,bmat,&neq2,
		  which,&nev,&tol,resid,&ncv,
		  z,&neq2,iparam,ipntr,workd,workl,&lworkl,rwork,&info));
		  free(workev);free(rwork);*/
	  
	  /* taking the square root */
	  
	  /*      dc=NNEW(double,2*nev);
		  nmc=NNEW(int,nev);
		  for(j=0;j<nev;j++){
		  dc[2*j]=sqrt(sqrt(d[2*j]*d[2*j]+d[2*j+1]*d[2*j+1])+d[2*j])/sqrt(2.);
		  dc[2*j+1]=sqrt(sqrt(d[2*j]*d[2*j]+d[2*j+1]*d[2*j+1])-d[2*j])/sqrt(2.);
		  if(d[2*j+1]<0.) dc[2*j+1]=-dc[2*j+1];
		  nmc[j]=nm;
		  }
		  FORTRAN(writeevcscomplex,(d,&nev,nmc,&fmin,&fmax));
		  FORTRAN(writeevcscomplex,(dc,&nev,nmc,&fmin,&fmax));free(dc);free(nmc);*/
	  
	  /* rearranging the real part of the eigenvalues */
	  
	  /*      for(j=0;j<nev;j++)d[j]=d[2*j];
		  RENEW(d,double,nev);*/
	  
	  /* rearranging the eigenvectors 
	     (first all real parts, then all imaginary parts */
	  
	  /*      v=NNEW(double,neq[1]);
		  for(j=0;j<nev;j++){
		  for(k=0;k<neq[1];k++){
		  v[k]=z[j*neq[1]+k];
		  }
		  for(k=0;k<neq2;k++){
		  z[j*neq[1]+k]=v[2*k];
		  z[j*neq[1]+neq2+k]=v[2*k+1];
		  }
		  }
		  free(v);*/
      }else{
	  d=NNEW(double,nev);
	  FORTRAN(dseupd,(&rvec,howmny,select,d,z,&ldz,&sigma,bmat,&neq[1],which,&nev,
			  &tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	  FORTRAN(writeevcs,(d,&nev,&nm,&fmin,&fmax));
      }
      free(select);free(resid);free(workd);free(workl);
      
      /* for double eigenmodes: the eigenmode for which the largest
	 amplitude has the lowest dof comes first */
      
      neq2=neq[1]/2;
      for (j=0;j<nev;j+=2){
	  
	  ampmax=0.;kmax1=0;
	  for(k=0;k<neq2;k++){
	      amp=z[2*j*neq[1]+k]*z[2*j*neq[1]+k]+
		  z[2*j*neq[1]+neq2+k]*z[2*j*neq[1]+neq2+k];
	      if(amp>ampmax){ampmax=amp;kmax1=k;}
	  }
	  
	  ampmax=0.;kmax2=0;
	  for(k=0;k<neq2;k++){
	      amp=z[(2*j+1)*neq[1]+k]*z[(2*j+1)*neq[1]+k]+
		  z[(2*j+1)*neq[1]+neq2+k]*z[(2*j+1)*neq[1]+neq2+k];
	      if(amp>ampmax){ampmax=amp;kmax2=k;}
	  }
	  
	  if(kmax2<kmax1){
	      printf("exchange!\n");
	      zstorage=NNEW(double,neq[1]);
	      memcpy(zstorage,&z[2*j*neq[1]],sizeof(double)*neq[1]);
	      memcpy(&z[2*j*neq[1]],&z[(2*j+1)*neq[1]],sizeof(double)*neq[1]);
	      memcpy(&z[(2*j+1)*neq[1]],zstorage,sizeof(double)*neq[1]);
	      free(zstorage);
	  }
      }
      
      /* writing the eigenvalues and mass matrix to a binary file */
      
      if(mei[3]==1){
	  
	  strcpy(fneig,jobnamec);
	  strcat(fneig,".eig");
	  
	  /* the first time the file is erased before writing, all subsequent
	     times the data is appended */
	  
	  if(*nevtot==0){
	      if((f1=fopen(fneig,"wb"))==NULL){
		  printf("*ERROR in arpack: cannot open eigenvalue file for writing...");
		  exit(0);
	      }
	      
	      /* storing a one as indication that this was a
		 cyclic symmetry calculation */
	      
	      if(fwrite(&one,sizeof(int),1,f1)!=1){
		  printf("*ERROR saving the cyclic symmetry flag to the eigenvalue file...");
		  exit(0);
	      }
	      
	      /* Hermitian */
	      
	      if(fwrite(&nherm,sizeof(int),1,f1)!=1){
		  printf("*ERROR saving the Hermitian flag to the eigenvalue file...");
		  exit(0);
	      }
	      
	  }else{
	      if((f1=fopen(fneig,"ab"))==NULL){
		  printf("*ERROR in arpack: cannot open eigenvalue file for writing...");
		  exit(0);
	      }
	  }
	  
	  /* nodal diameter */
	  
	  if(fwrite(&nm,sizeof(int),1,f1)!=1){
	      printf("*ERROR saving the nodal diameter to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* number of eigenfrequencies */
	  
	  if(fwrite(&nev,sizeof(int),1,f1)!=1){
	      printf("*ERROR saving the number of eigenfrequencies to the eigenvalue file...");
	      exit(0);
	  }
	  
	  /* the eigenfrequencies are stored as radians/time */
	  
	  if(fwrite(d,sizeof(double),nev,f1)!=nev){
	      printf("*ERROR saving the eigenfrequencies to the eigenvalue file...");
	      exit(0);
	  }
	  if(*nevtot==0){
	      
	      /* mass matrix */
	      
	      if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
		  printf("*ERROR saving the diagonal of the mass matrix to the eigenvalue file...");
		  exit(0);
	      }
	      if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
		  printf("*ERROR saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
		  exit(0);
	      }
	  }
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
      
      temp_array=NNEW(double,neq[1]);
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
      if(strcmp1(&filab[261],"E   ")==0) eent=NNEW(double,2*6**nk*ngraph);
      if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0))
	  
// real and imaginary part of the forces
	  
	  fnt=NNEW(double,2*mt**nk*ngraph);
      if(strcmp1(&filab[2697],"ME  ")==0) emnt=NNEW(double,2*6**nk*ngraph);
      if(strcmp1(&filab[522],"ENER")==0)
	  enernt=NNEW(double,*nk*ngraph);
      if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
	 (strcmp1(&filab[2175],"CONT")==0))
	  stxt=NNEW(double,2*6*mi[0]**ne*ngraph);
      
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
		       &imag,mi,emn));
      
      /* check that the tensor fields which are extrapolated from the
	 integration points are requested in global coordinates */
      
      if(strcmp1(&filab[174],"S   ")==0){
	  if((strcmp1(&filab[179],"L")==0)&&(*norien>0)){
	      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculations\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	      strcpy1(&filab[179],"G",1);
	  }
      }
      
      if(strcmp1(&filab[261],"E   ")==0){
	  if((strcmp1(&filab[266],"L")==0)&&(*norien>0)){
	      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	      strcpy1(&filab[266],"G",1);
	  }
      }
      
      if(strcmp1(&filab[1479],"PHS ")==0){
	  if((strcmp1(&filab[1484],"L")==0)&&(*norien>0)){
	      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	      strcpy1(&filab[1484],"G",1);
	  }
      }
      
      if(strcmp1(&filab[1653],"MAXS")==0){
	  if((strcmp1(&filab[1658],"L")==0)&&(*norien>0)){
	      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	      strcpy1(&filab[1658],"G",1);
	  }
      }   
      
      if(strcmp1(&filab[2523],"MAXE")==0){
	  if((strcmp1(&filab[2528],"L")==0)&&(*norien>0)){
	      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
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
      
      if(strcmp1(&filab[2610],"PRF")==0){
	  fnr=NNEW(double,mt*nkt);
	  fni=NNEW(double,mt*nkt);
      }
      
      /* start of output calculations */
      
      lfin=0;
      for(j=0;j<nev;++j){
	  lint=lfin;
	  lfin=lfin+neq[1];
	  
	  if(mei[3]==1){
	      if(fwrite(&z[lint],sizeof(double),neq[1],f1)!=neq[1]){
		  printf("*ERROR saving data to the eigenvalue file...");
		  exit(0);
	      }
	  }
	  
	  /* check whether the frequency belongs to the requested
	     interval */
	  
	  if(fmin>-0.5){
	      if(fmin*fmin>d[j]) continue;
	  }
	  if(fmax>-0.5){
	      if(fmax*fmax<d[j]) continue;
	  }
	  
	  if(*nprint>0)FORTRAN(writehe,(&j));
	  
	  eei=NNEW(double,6*mi[0]*ne0);
	  if(*nener==1){
	      stiini=NNEW(double,6*mi[0]*ne0);
	      enerini=NNEW(double,mi[0]*ne0);}
	  
//    memset(&v[0],0.,sizeof(double)*2*mt**nk);
	  DMEMSET(v,0,2*mt**nk,0.);
	  
	  for(k=0;k<neq[1];k+=neq[1]/2){
	      
	      for(i=0;i<6*mi[0]*ne0;i++){eme[i]=0.;}
	      
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
			  else{xreal=z[lint+idof];ximag=z[lint+idof+neq[1]/2];}
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
			  f,&fn[kkv],nactdof,&iout,qa,vold,&z[lint+k],
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
			  f,&fn[kkv],nactdof,&iout,qa,vold,&z[lint+k],
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
	  
	  if(strcmp1(&filab[1566],"MAXU")==0){
	      
	      /* determining the ray vector; the components of the
		 ray vector are the coordinates of the node in node set
		 RAY */
	      
	      iray=0;
	      for(i=0;i<*nset;i++){
		  if(strcmp1(&set[81*i],"RAYN")==0){
		      iray=ialset[istartset[i]-1];
		      vray[0]=co[3*iray-3];
		      vray[1]=co[3*iray-2];
		      vray[2]=co[3*iray-1];
		      break;
		  }
	      }
	      if(iray==0){
		  printf("/n*ERROR in arpackcs: no light ray vector/n/n");
		  FORTRAN(stop,());
	      }
	      
	      /* initialization */
	      
	      for(l1=0;l1<4**nk;l1++){vmax[l1]=0.;}
	      
	      /* vector p1 is a point on the rotation axis
		 vector p2 is a unit vector along the axis */
	      
	      for(l2=0;l2<3;l2++){p1[l2]=cs[5+l2];}
	      for(l2=0;l2<3;l2++){p2[l2]=cs[8+l2]-p1[l2];}
	      dd=sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
	      for(l2=0;l2<3;l2++){p2[l2]/=dd;}
	      
	      /* determine the time for the worst displacement
		 orthogonal to a give light ray vector ; */
	      
	      for(l1=0;l1<*nk;l1++){
		  
		  /*  determining a vector through node (l1+1) and
		      orthogonal to the rotation axis */
		  
		  for(l2=0;l2<3;l2++){q[l2]=co[3*l1+l2]-p1[l2];}
		  dd=q[0]*p2[0]+q[1]*p2[1]+q[2]*p2[2];
		  for(l2=0;l2<3;l2++){q[l2]-=dd*p2[l2];}
		  
		  /* determining a vector tan orthogonal to vector q
		     and the ray vector */
		  
		  tan[0]=q[1]*vray[2]-q[2]*vray[1];
		  tan[1]=q[2]*vray[0]-q[0]*vray[2];
		  tan[2]=q[0]*vray[1]-q[1]*vray[0];
		  
		  printf("tangent= %d,%e,%e,%e\n",l1,tan[0],tan[1],tan[2]);
		  
		  worstumax=0.;
		  iworsttime=0;
		  for(l3=0;l3<360;l3++){
		      ctl=cos(l3/constant);
		      stl=sin(l3/constant);
		      for(l2=1;l2<4;l2++){
			  l=mt*l1+l2;
			  vl[l2]=ctl*v[l]-stl*v[l+mt**nk];
		      }
		      
		      /* displacement component along the tangent vector
			 (no absolute value!) */
		      
		      dd=vl[1]*tan[0]+vl[2]*tan[1]+vl[3]*tan[2];
		      if(dd>worstumax){
			  worstumax=dd;
			  iworsttime=l3;
		      }
		  }
		  ctl=cos(iworsttime/constant);
		  stl=sin(iworsttime/constant);
		  for(l2=1;l2<4;l2++){
		      l=mt*l1+l2;
		      vl[l2]=ctl*v[l]-stl*v[l+mt**nk];
		  }
		  vmax[4*l1]=1.*iworsttime;
		  vmax[4*l1+1]=vl[1];
		  vmax[4*l1+2]=vl[2];
		  vmax[4*l1+3]=vl[3];
		  
	      }
	  }
	  
	  /* determine the worst principal stress anywhere
	     in the structure as a function of time; 
	     the worst principal stress is the maximum
	     of the absolute value of the principal stresses */
	  
	  if(strcmp1(&filab[1653],"MAXS")==0){
	      
	      /* determining the set of nodes for the 
		 worst principal stress calculation */
	      
	      ielset=0;
	      for(i=0;i<*nset;i++){
		  if(strcmp1(&set[81*i],"STRESSDOMAINN")==0){
		      ielset=i+1;
		      break;
		  }
	      }
	      if(ielset==0){
		  printf("\n*ERROR in arpackcs: no node set for MAXS\n");
		  printf("       (must have the name STRESSDOMAIN)\n\n");
		  FORTRAN(stop,());
	      }
	      
	      for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
		  if(ialset[i1]>0){
		      l1=ialset[i1]-1;
		      
		      worstpsmax=0.;
		      for(l3=0;l3<360;l3++){
			  ctl=cos(l3/constant);
			  stl=sin(l3/constant);
			  for(l2=0;l2<6;l2++){
			      l=6*l1+l2;
			      stnl[l2]=ctl*stn[l]-stl*stn[l+6**nk];
			  }
			  
			  /* determining the eigenvalues */
			  
			  v1=stnl[0]+stnl[1]+stnl[2];
			  v2=stnl[1]*stnl[2]+stnl[0]*stnl[2]+stnl[0]*stnl[1]-
			      (stnl[5]*stnl[5]+stnl[4]*stnl[4]+stnl[3]*stnl[3]);
			  v3=stnl[0]*(stnl[1]*stnl[2]-stnl[5]*stnl[5])
			      -stnl[3]*(stnl[3]*stnl[2]-stnl[4]*stnl[5])
			      +stnl[4]*(stnl[3]*stnl[5]-stnl[4]*stnl[1]);
			  bb=v2-v1*v1/3.;
			  cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
			  if(fabs(bb)<=1.e-10){
			      if(fabs(cc)>1.e-10){
				  al[0]=-pow(cc,(1./3.));
			      }else{
				  al[0]=0.;
			      }
			      al[1]=al[0];
			      al[2]=al[0];
			  }else{
			      cm=2.*sqrt(-bb/3.);
			      cn=3.*cc/(cm*bb);
			      if(fabs(cn)>1.){
				  if(cn>1.){
				      cn=1.;
				  }else{
				      cn=-1.;
				  }
			      }
			      tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
			      al[0]=cm*cos(tt);
			      al[1]=cm*cos(tt+2.*pi/3.);
			      al[2]=cm*cos(tt+4.*pi/3.);
			  }
			  for(l2=0;l2<3;l2++){
			      al[l2]+=v1/3.;
			  }
			  dd=fabs(al[0]);
			  if(fabs(al[1])>dd) dd=fabs(al[1]);
			  if(fabs(al[2])>dd) dd=fabs(al[2]);
			  if(dd>worstpsmax){
			      worstpsmax=dd;
			      stnmax[7*l1]=dd;
			      for(l2=1;l2<7;l2++){
				  stnmax[7*l1+l2]=stnl[l2-1];
			      }
			  }
		      }
		      
		  }else{
		      l1=ialset[i1-2]-1;
		      do{
			  l1=l1-ialset[i1];
			  if(l1>=ialset[i1-1]-1) break;
			  
			  worstpsmax=0.;
			  for(l3=0;l3<360;l3++){
			      ctl=cos(l3/constant);
			      stl=sin(l3/constant);
			      for(l2=0;l2<6;l2++){
				  l=6*l1+l2;
				  stnl[l2]=ctl*stn[l]-stl*stn[l+6**nk];
			      }
			      
			      /* determining the eigenvalues */
			      
			      v1=stnl[0]+stnl[1]+stnl[2];
			      v2=stnl[1]*stnl[2]+stnl[0]*stnl[2]+stnl[0]*stnl[1]-
				  (stnl[5]*stnl[5]+stnl[4]*stnl[4]+stnl[3]*stnl[3]);
			      v3=stnl[0]*(stnl[1]*stnl[2]-stnl[5]*stnl[5])
				  -stnl[3]*(stnl[3]*stnl[2]-stnl[4]*stnl[5])
				  +stnl[4]*(stnl[3]*stnl[5]-stnl[4]*stnl[1]);
			      bb=v2-v1*v1/3.;
			      cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
			      if(fabs(bb)<=1.e-10){
				  if(fabs(cc)>1.e-10){
				      al[0]=-pow(cc,(1./3.));
				  }else{
				      al[0]=0.;
				  }
				  al[1]=al[0];
				  al[2]=al[0];
			      }else{
				  cm=2.*sqrt(-bb/3.);
				  cn=3.*cc/(cm*bb);
				  if(fabs(cn)>1.){
				      if(cn>1.){
					  cn=1.;
				      }else{
					  cn=-1.;
				      }
				  }
				  tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
				  al[0]=cm*cos(tt);
				  al[1]=cm*cos(tt+2.*pi/3.);
				  al[2]=cm*cos(tt+4.*pi/3.);
			      }
			      for(l2=0;l2<3;l2++){
				  al[l2]+=v1/3.;
			      }
			      dd=fabs(al[0]);
			      if(fabs(al[1])>dd) dd=fabs(al[1]);
			      if(fabs(al[2])>dd) dd=fabs(al[2]);
			      if(dd>worstpsmax){
				  worstpsmax=dd;
				  stnmax[7*l1]=dd;
				  for(l2=1;l2<7;l2++){
				      stnmax[7*l1+l2]=stnl[l2-1];
				  }
			      }
			  }
			  
		      }while(1);
		  }
	      }
	  }
	  
	  /* determine the worst principal strain anywhere
	     in the structure as a function of time; 
	     the worst principal strain is the maximum
	     of the absolute value of the principal strains,
	     times its original sign */
	  
    if(strcmp1(&filab[2523],"MAXE")==0){

      /* determining the set of nodes for the 
         worst principal strain calculation */

      ielset=0;
      for(i=0;i<*nset;i++){
	if(strcmp1(&set[81*i],"STRAINDOMAINN")==0){
	  ielset=i+1;
	  break;
	}
      }
      if(ielset==0){
	printf("\n*ERROR in arpackcs: no node set for MAXE\n");
	printf("       (must have the name STRAINDOMAIN)\n\n");
	FORTRAN(stop,());
      }

      for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
	if(ialset[i1]>0){
	  l1=ialset[i1]-1;

	  worstpsmax=0.;
	  for(l3=0;l3<360;l3++){
	    ctl=cos(l3/constant);
	    stl=sin(l3/constant);
	    for(l2=0;l2<6;l2++){
	      l=6*l1+l2;
	      eenl[l2]=ctl*een[l]-stl*een[l+6**nk];
	    }
	    
	    /* determining the eigenvalues */
	    
	    v1=eenl[0]+eenl[1]+eenl[2];
	    v2=eenl[1]*eenl[2]+eenl[0]*eenl[2]+eenl[0]*eenl[1]-
	      (eenl[5]*eenl[5]+eenl[4]*eenl[4]+eenl[3]*eenl[3]);
	    v3=eenl[0]*(eenl[1]*eenl[2]-eenl[5]*eenl[5])
	      -eenl[3]*(eenl[3]*eenl[2]-eenl[4]*eenl[5])
	      +eenl[4]*(eenl[3]*eenl[5]-eenl[4]*eenl[1]);
	    bb=v2-v1*v1/3.;
	    cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
	    if(fabs(bb)<=1.e-10){
	      if(fabs(cc)>1.e-10){
		al[0]=-pow(cc,(1./3.));
	      }else{
		al[0]=0.;
	      }
	      al[1]=al[0];
	      al[2]=al[0];
	    }else{
	      cm=2.*sqrt(-bb/3.);
	      cn=3.*cc/(cm*bb);
	      if(fabs(cn)>1.){
		if(cn>1.){
		  cn=1.;
		}else{
		  cn=-1.;
		}
	      }
	      tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
	      al[0]=cm*cos(tt);
	      al[1]=cm*cos(tt+2.*pi/3.);
	      al[2]=cm*cos(tt+4.*pi/3.);
	    }
	    for(l2=0;l2<3;l2++){
	      al[l2]+=v1/3.;
	    }
	    dd=fabs(al[0]);
	    if(fabs(al[1])>dd) dd=fabs(al[1]);
	    if(fabs(al[2])>dd) dd=fabs(al[2]);
	    if(dd>worstpsmax){
		worstpsmax=dd;
		eenmax[7*l1]=dd;
		for(l2=1;l2<7;l2++){
		    eenmax[7*l1+l2]=eenl[l2-1];
		}
	    }
	  }
	  
	}else{
	  l1=ialset[i1-2]-1;
	  do{
	    l1=l1-ialset[i1];
	    if(l1>=ialset[i1-1]-1) break;

	    worstpsmax=0.;
	    for(l3=0;l3<360;l3++){
	      ctl=cos(l3/constant);
	      stl=sin(l3/constant);
	      for(l2=0;l2<6;l2++){
		l=6*l1+l2;
		eenl[l2]=ctl*een[l]-stl*een[l+6**nk];
	      }
	      
	      /* determining the eigenvalues */
	      
	      v1=eenl[0]+eenl[1]+eenl[2];
	      v2=eenl[1]*eenl[2]+eenl[0]*eenl[2]+eenl[0]*eenl[1]-
		(eenl[5]*eenl[5]+eenl[4]*eenl[4]+eenl[3]*eenl[3]);
	      v3=eenl[0]*(eenl[1]*eenl[2]-eenl[5]*eenl[5])
		-eenl[3]*(eenl[3]*eenl[2]-eenl[4]*eenl[5])
		+eenl[4]*(eenl[3]*eenl[5]-eenl[4]*eenl[1]);
	      bb=v2-v1*v1/3.;
	      cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
	      if(fabs(bb)<=1.e-10){
		if(fabs(cc)>1.e-10){
		  al[0]=-pow(cc,(1./3.));
		}else{
		  al[0]=0.;
		}
		al[1]=al[0];
		al[2]=al[0];
	      }else{
		cm=2.*sqrt(-bb/3.);
		cn=3.*cc/(cm*bb);
		if(fabs(cn)>1.){
		  if(cn>1.){
		    cn=1.;
		  }else{
		    cn=-1.;
		  }
		}
		tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
		al[0]=cm*cos(tt);
		al[1]=cm*cos(tt+2.*pi/3.);
		al[2]=cm*cos(tt+4.*pi/3.);
	      }
	      for(l2=0;l2<3;l2++){
		al[l2]+=v1/3.;
	      }
	      dd=fabs(al[0]);
	      if(fabs(al[1])>dd) dd=fabs(al[1]);
	      if(fabs(al[2])>dd) dd=fabs(al[2]);
	      if(dd>worstpsmax){
		  worstpsmax=dd;
		  eenmax[7*l1]=dd;
		  for(l2=1;l2<7;l2++){
		      eenmax[7*l1+l2]=eenl[l2-1];
		  }
	      }
	    }
	    
	  }while(1);
	}
      }
    }
    
    /* mapping the results to the other sectors */

    for(l=0;l<*nk;l++){inumt[l]=inum[l];}

    icntrl=2;imag=1;

    FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));

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
    if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
       (strcmp1(&filab[2175],"CONT")==0)){
      for(l=0;l<6*mi[0]**ne;l++){stxt[l]=stx[l];}
      for(l=0;l<6*mi[0]**ne;l++){stxt[l+6*mi[0]**ne*ngraph]=stx[l+6*mi[0]**ne];}}
    if(strcmp1(&filab[2697],"ME  ")==0){
	for(l=0;l<6**nk;l++){emnt[l]=emn[l];};
	for(l=0;l<6**nk;l++){emnt[l+6**nk*ngraph]=emn[l+6**nk];}}

    for(jj=0;jj<*mcs;jj++){
      ilength=cs[17*jj+3];
      is=cs[17*jj+4];
      lprev=cs[17*jj+13];
      for(i=1;i<is;i++){
        
        for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
        
        theta=i*nm*2.*pi/cs[17*jj];
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
        
        /* real part of the mechanical strains */

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
          &stnt[6**nk*ngraph],qfnt,&eent[6**nk*ngraph],cs,&nkt,&icntrl,
          t,filab,&imag,mi,&emnt[6**nk*ngraph]));

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

    /* determining magnitude and phase angle for the forces */

    if(strcmp1(&filab[2610],"PRF")==0){
      for(l1=0;l1<nkt;l1++){
	for(l2=0;l2<4;l2++){
	  l=mt*l1+l2;
	  fnreal=fnt[l];
	  fnimag=fnt[l+mt**nk*ngraph];
	  fnr[l]=sqrt(fnreal*fnreal+fnimag*fnimag);
	  if(fabs(fnreal)<1.e-10){
	    if(fnimag>0){fni[l]=90.;}
	    else{fni[l]=-90.;}
	  }
	  else{
	    fni[l]=atan(fnimag/fnreal)*constant;
	    if(fnreal<0) fni[l]+=180.;
	  }
	}
      }
    }

    ++*kode;
    if(d[j]>=0.){
	freq=sqrt(d[j])/6.283185308;
    }else{
	freq=0.;
    }

    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40*net);ipneigh=NNEW(int,nkt);
    }

    frd(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,
	    kode,filab,eent,t1t,fnt,&freq,epn,ielmatt,matname,enernt,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&j,&nm,trab,inotrt,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stxt,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&net,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emnt,
	    thicke,jobnamec,output,qfx);

    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}

  }

  if((fmax>-0.5)&&(fmax*fmax>d[nev-1])){
    printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
  }

  free(adb);free(aub);free(temp_array);free(coefmpcnew);

  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
     (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
     (strcmp1(&filab[1044],"ERR ")==0)) 
     free(stn);

  free(v);free(fn);free(inum);free(stx);free(z);free(d);

  if((strcmp1(&filab[261],"E   ")==0)||(strcmp1(&filab[2523],"MAXE")==0)) free(een);
  if(strcmp1(&filab[522],"ENER")==0) free(enern);
  if(strcmp1(&filab[2697],"ME  ")==0) free(emn);

  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)) free(vt);
  if(strcmp1(&filab[87],"NT  ")==0) free(t1t);
  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)||
     (strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)) free(stnt);
  if(strcmp1(&filab[261],"E   ")==0) free(eent);
  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)) free(fnt);
  if(strcmp1(&filab[522],"ENER")==0) free(enernt);
  if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
     (strcmp1(&filab[2175],"CONT")==0)) free(stxt);
  if(strcmp1(&filab[2697],"ME  ")==0) free(emnt);

  free(cot);free(kont);free(ipkont);free(lakont);free(inumt);free(ielmatt);
  if(*ntrans>0){free(inotrt);}

  if(mei[3]==1){
      (*nevtot)+=nev;
      fclose(f1);
  }

  /* end loop over the nodal diameters */

  }

  if(*iperturb!=0){
      if(ncont!=0){
	  *ne=ne0;*nkon=nkon0;
	  if(*nener==1){
	      RENEW(ener,double,mi[0]**ne*2);
	  }
	  RENEW(ipkon,int,*ne);
	  RENEW(lakon,char,8**ne);
	  RENEW(kon,int,*nkon);
	  if(*norien>0){
	      RENEW(ielorien,int,mi[2]**ne);
	  }
	  RENEW(ielmat,int,mi[2]**ne);
	  printf("arpackcs3 cg=%e\n",cg[0]);
	  free(cg);
	  free(straight);
	  free(imastop);free(itiefac);free(islavsurf);free(islavnode);
	  free(nslavnode);free(iponoels);free(inoels);free(imastnode);
	  free(nmastnode);free(itietri);free(koncont);
	  free(areaslav);free(springarea);free(xmastnor);free(xnormastface);
      }
  }

  free(inocs);free(ielcs);free(xstiff);
  free(ipobody);

  if(strcmp1(&filab[870],"PU")==0){free(vr);free(vi);}
  if(strcmp1(&filab[1479],"PHS")==0){free(stnr);free(stni);}
  if(strcmp1(&filab[1566],"MAXU")==0){free(vmax);}
  if(strcmp1(&filab[1653],"MAXS")==0){free(stnmax);}
  if(strcmp1(&filab[2523],"MAXE")==0){free(eenmax);}
  if(strcmp1(&filab[2610],"PRF")==0){free(fnr);free(fni);}

  for(i=0;i<6*mi[0]*ne0;i++){eme[i]=0.;}

  if(*iperturb!=0){
      mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
      mpcinfo[3]=maxlenmpc;
  }

  *irowp=irow;*enerp=ener;*xstatep=xstate;

  return;
}

#endif
