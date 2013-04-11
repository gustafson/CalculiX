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

void arpackcs(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int *icol, int *jq, int *irow, int *neq, int *nzl, 
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
             double *xstate, int *npmat_, char *matname, int *mint_,
	     int *ics, double *cs, int *mpcend, int *ncmat_,
             int *nstate_, int *mcs, int *nkon, double *ener,
             char *jobnamec, char *output, char *set, int *nset, 
             int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
             char *cbody, int *ibody, double *xbody, int *nbody,
	     int *nevtot){

  /* calls the Arnoldi Package (ARPACK) for cyclic symmetry calculations */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A",*lakont=NULL,
      description[13]="            ",fneig[132]="";
  int *inum=NULL,k,ido,dz,iparam[11],ipntr[11],lworkl,idir,
    info,rvec=1,*select=NULL,lfin,j,lint,iout=1,nm,index,inode,id,i,idof,
    ielas,icmd,kk,l,nkt,icntrl,*kont=NULL,*ipkont=NULL,*inumt=NULL,
    *ielmatt=NULL,net,imag,icomplex,kk4,kk5,kk6,iinc=1,istep=1,nev,ncv,
    mxiter,lprev,ilength,ij,i1,i2,iel,ielset,node,indexe,nope,ml1,
    *inocs=NULL,*ielcs=NULL,jj,l1,l2,ngraph,is,jrow,*ipobody=NULL,
    *inotrt=NULL,symmetryflag=0,inputformat=0,inewton=0,ifreebody;
  int mass=1, stiffness=1, buckling=0, rhsi=0, intscheme=0,*ncocon=NULL,
      coriolis=0,iworsttime,l3,idummy=1;
  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,*vr=NULL,
    *workl=NULL,*aux=NULL,*d=NULL,sigma=1,*temp_array=NULL,*vini=NULL,
    *een=NULL,cam[3],*f=NULL,*fn=NULL,qa[2],*fext=NULL,*epn=NULL,*stiini=NULL,
    *xstateini=NULL,theta,pi,*coefmpcnew=NULL,*xstiff=NULL,*vi=NULL,
    *vt=NULL,*fnt=NULL,*stnt=NULL,*eent=NULL,*cot=NULL,t[3],ctl,stl,
    *t1t=NULL,freq,*stx=NULL,*enern=NULL,*enernt=NULL,*xstaten=NULL,
    *eei=NULL,*enerini=NULL,*cocon=NULL,*qfx=NULL,*qfn=NULL,*qfnt=NULL,
    tol,fmin,fmax,xreal,ximag,*cgr=NULL,*xloadold=NULL,reltime,constant,
    vreal,vimag,*stnr=NULL,*stni=NULL,stnreal,stnimag,*vmax=NULL,
    *stnmax=NULL,vl[4],stnl[6],dd,v1,v2,v3,bb,cc,al[3],cm,cn,tt,
    worstpsmax;
  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

  int *ipneigh=NULL,*neigh=NULL;

#ifdef SGI
  int token;
#endif

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

  /* determining the internal forces and the stiffness coefficients */

  f=NNEW(double,*neq);

  /* allocating a field for the stiffness matrix */

  xstiff=NNEW(double,27**mint_**ne);

  iout=-1;
  v=NNEW(double,5**nk);
  fn=NNEW(double,4**nk);
  stx=NNEW(double,6**mint_**ne);
  if(*iperturb==0){
     FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t0,ithermal,
	       prestr,iprestr,filab,eme,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               sti,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
               fmpc,nelemload,nload));
  }else{
     FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	       prestr,iprestr,filab,eme,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               sti,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
               fmpc,nelemload,nload));
  }
  free(f);free(v);free(fn);free(stx);
  iout=1;

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

  /* filling in the matrix */

  for(nm=cs[1];nm<=cs[2];++nm){

  ad=NNEW(double,*neq);
  au=NNEW(double,*nzs);

  adb=NNEW(double,*neq);
  aub=NNEW(double,*nzs);

  fext=NNEW(double,*neq);
  if(*iperturb==0){
    FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mint_,
	      ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,&rhsi,
              &intscheme,mcs,&coriolis,ibody,xloadold,&reltime,ielcs));
  }
  else{
    FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mint_,
              ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,&rhsi,
              &intscheme,mcs,&coriolis,ibody,xloadold,&reltime,ielcs));
  }

  free(fext);

  if(*nmethod==0){

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;time=0.;
    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1,
         fn,&time,epn,ielmat,matname,enern,xstaten,nstate_,&istep,&iinc,
	 iperturb,ener,mint_,output,ithermal,qfn,&j,&nm,
         trab,inotr,ntrans,orab,ielorien,norien,description,
	 ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&idummy));
    
    free(inum);free(ipneigh);free(neigh);FORTRAN(stop,());

  }


  /* LU decomposition of the left hand matrix */

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,neq,nzs,&symmetryflag,
                   &inputformat);
#else
    printf("*ERROR in arpackcs: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,neq,nzs,token);
#else
    printf("*ERROR in arpackcs: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,neq,nzs);
#else
    printf("*ERROR in arpackcs: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

  free(au);free(ad);

  /* calculating the eigenvalues and eigenmodes */
  
  printf(" Calculating the eigenvalues and the eigenmodes\n");

  ido=0;
  dz=*neq;
  iparam[0]=1;
  iparam[2]=mxiter;
  iparam[3]=1;
  iparam[6]=3;

  lworkl=ncv*(8+ncv);
  info=0;

  resid=NNEW(double,*neq);
  z=NNEW(double,ncv**neq);
  workd=NNEW(double,3**neq);
  workl=NNEW(double,lworkl);

  FORTRAN(dsaupd,(&ido,bmat,neq,which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
	  workl,&lworkl,&info));

  temp_array=NNEW(double,*neq);

  while((ido==-1)||(ido==1)||(ido==2)){
    if(ido==-1){
      FORTRAN(op,(neq,aux,&workd[ipntr[0]-1],temp_array,adb,aub,icol,irow,nzl));
    }

      /* solve the linear equation system  */

    if((ido==-1)||(ido==1)){

      if(ido==-1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(temp_array,neq);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(temp_array,token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(temp_array,neq);
#endif
        }
        for(jrow=0;jrow<*neq;jrow++){
          workd[ipntr[1]-1+jrow]=temp_array[jrow];
        }
      }
      else if(ido==1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(&workd[ipntr[2]-1],neq);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(&workd[ipntr[2]-1],token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(&workd[ipntr[2]-1],neq);
#endif
        }
        for(jrow=0;jrow<*neq;jrow++){
          workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
        }
      }
      
    }

    if(ido==2){
      FORTRAN(op,(neq,aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],adb,aub,icol,irow,nzl));
    }

    FORTRAN(dsaupd,(&ido,bmat,neq,which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
	    workl,&lworkl,&info));
  }

/*--------------------------------------------------------------------*/
/*
   -----------
   free memory
   -----------
*/
  free(temp_array);
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

  if(info!=0){
    printf("*ERROR in arpackcs: info=%d\n",info);
    printf("       # of converged eigenvalues=%d\n\n",iparam[4]);
  }         

  select=NNEW(int,ncv);
  d=NNEW(double,nev);

  FORTRAN(dseupd,(&rvec,howmny,select,d,z,&dz,&sigma,bmat,neq,which,&nev,&tol,resid,
	  &ncv,z,&dz,iparam,ipntr,workd,workl,&lworkl,&info));

  FORTRAN(writeevcs,(d,&nev,&nm,&fmin,&fmax));

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
      }else{
	  if((f1=fopen(fneig,"ab"))==NULL){
	      printf("*ERROR in arpack: cannot open eigenvalue file for writing...");
	      exit(0);
	  }
      }

      if(fwrite(&nm,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(&nev,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }

      /* the eigenfrequencies are stored as radians/time */

      if(fwrite(d,sizeof(double),nev,f1)!=nev){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(*nevtot==0){
	  if(fwrite(adb,sizeof(double),*neq,f1)!=*neq){
	      printf("*ERROR saving data to the eigenvalue file...");
	      exit(0);
	  }
	  if(fwrite(aub,sizeof(double),*nzs,f1)!=*nzs){
	      printf("*ERROR saving data to the eigenvalue file...");
	      exit(0);
	  }
      }
  }

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,10**nk);
  fn=NNEW(double,8**nk);
  stn=NNEW(double,12**nk);

  if(strcmp1(&filab[18],"E   ")==0) een=NNEW(double,12**nk);
  if(strcmp1(&filab[36],"ENER")==0) enern=NNEW(double,2**nk);

  inum=NNEW(int,*nk);
  stx=NNEW(double,6**mint_**ne);
  
  temp_array=NNEW(double,*neq);
  coefmpcnew=NNEW(double,*mpcend);

  cot=NNEW(double,3**nk*ngraph);
  if(*ntrans>0){inotrt=NNEW(int,2**nk*ngraph);}
  if(strcmp1(&filab[0],"U   ")==0)
    vt=NNEW(double,5**nk*ngraph);
  if(strcmp1(&filab[6],"NT  ")==0)
    t1t=NNEW(double,*nk*ngraph);
  if(strcmp1(&filab[12],"S   ")==0)
    stnt=NNEW(double,6**nk*ngraph);
  if(strcmp1(&filab[18],"E   ")==0)
    eent=NNEW(double,6**nk*ngraph);
  if(strcmp1(&filab[24],"RF  ")==0)
    fnt=NNEW(double,4**nk*ngraph);
  if(strcmp1(&filab[36],"ENER")==0)
    enernt=NNEW(double,*nk*ngraph);

  kont=NNEW(int,*nkon*ngraph);
  ipkont=NNEW(int,*ne*ngraph);
  for(l=0;l<*ne*ngraph;l++)ipkont[l]=-1;
  lakont=NNEW(char,8**ne*ngraph);
  inumt=NNEW(int,*nk*ngraph);
  ielmatt=NNEW(int,*ne*ngraph);

  nkt=ngraph**nk;
  net=ngraph**ne;

  /* copying the coordinates of the first sector */

  for(l=0;l<3**nk;l++){cot[l]=co[l];}
  if(*ntrans>0){for(l=0;l<*nk;l++){inotrt[2*l]=inotr[2*l];}}
  for(l=0;l<*nkon;l++){kont[l]=kon[l];}
  for(l=0;l<*ne;l++){ipkont[l]=ipkon[l];}
  for(l=0;l<8**ne;l++){lakont[l]=lakon[l];}
  for(l=0;l<*ne;l++){ielmatt[l]=ielmat[l];}

  /* generating the coordinates for the other sectors */

  icntrl=1;

  FORTRAN(rectcyl,(cot,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag));

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
            ielmatt[l+i**ne]=ielmat[l];
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

  FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,&imag));

  /* check that the tensor fields which are extrapolated from the
     integration points are requested in global coordinates */

  if(strcmp1(&filab[12],"S   ")==0){
    if(strcmp1(&filab[17],"L")==0){
      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculations\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[17],"G",1);
    }
  }

  if(strcmp1(&filab[18],"E   ")==0){
    if(strcmp1(&filab[23],"L")==0){
      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[23],"G",1);
    }
  }

  if(strcmp1(&filab[102],"PHS ")==0){
    if(strcmp1(&filab[107],"L")==0){
      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[107],"G",1);
    }
  }

  if(strcmp1(&filab[114],"MAXS")==0){
    if(strcmp1(&filab[119],"L")==0){
      printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[119],"G",1);
    }
  }   

  /* allocating fields for magnitude and phase information of
     displacements and stresses */

  if(strcmp1(&filab[60],"PU")==0){
    vr=NNEW(double,5**nk);
    vi=NNEW(double,5**nk);
  }

  if(strcmp1(&filab[102],"PHS")==0){
    stnr=NNEW(double,6**nk);
    stni=NNEW(double,6**nk);
  }

  if(strcmp1(&filab[108],"MAXU")==0){
    vmax=NNEW(double,4**nk);
  }

  if(strcmp1(&filab[114],"MAXS")==0){
    stnmax=NNEW(double,7**nk);
  }

  /* start of output calculations */

  lfin=0;
  for(j=0;j<nev;++j){
    lint=lfin;
    lfin=lfin+*neq;

    if(mei[3]==1){
	if(fwrite(&z[lint],sizeof(double),*neq,f1)!=*neq){
	    printf("*ERROR saving data to the eigenvalue file...");
	    exit(0);
	}
    }

    /* check whether the frequency belongs to the requested
       interval */

    if(fmin>d[j]) continue;
    if(fmax>0.){
      if(fmax<d[j]) break;
    }

    if(*nprint>0)FORTRAN(writehe,(&j));

    /* generating the cyclic MPC's (needed for nodal diameters
       different from 0 */

    eei=NNEW(double,6**mint_**ne);
    if(*nener==1){
	stiini=NNEW(double,6**mint_**ne);
	enerini=NNEW(double,*mint_**ne);}

    for(k=0;k<*neq;k+=*neq/2){

      for(i=0;i<6**mint_**ne;i++){eme[i]=0.;}

      if(k==0) {kk=0;kk4=0;kk5=0;kk6=0;if(*nprint>0)FORTRAN(writere,());}
      else {kk=*nk;kk4=4**nk;kk5=5**nk;kk6=6**nk;if(*nprint>0)FORTRAN(writeim,());}
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
	    idof=nactdof[4*(inode-1)+idir]-1;
            if(idof==-1){xreal=1.;ximag=1.;}
            else{xreal=z[lint+idof];ximag=z[lint+idof+*neq/2];}
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
	FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,&v[kk5],&stn[kk6],inum,
            stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t0,ithermal,
	    prestr,iprestr,filab,eme,&een[kk6],iperturb,
            f,&fn[kk4],nactdof,&iout,qa,vold,&z[lint+k],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,neq,veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],sti,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
            nelemload,nload));}
      else{
	FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,&v[kk5],&stn[kk6],inum,
            stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,&een[kk6],iperturb,
            f,&fn[kk4],nactdof,&iout,qa,vold,&z[lint+k],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,neq,veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,&icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],sti,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
            nelemload,nload));
      }

    }
    free(eei);
    if(*nener==1){free(stiini);free(enerini);}

    /* determining magnitude and phase angle for the displacements */

    if(strcmp1(&filab[60],"PU")==0){
      for(l1=0;l1<*nk;l1++){
	for(l2=0;l2<4;l2++){
	  l=5*l1+l2;
	  vreal=v[l];
	  vimag=v[l+5**nk];
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

    if(strcmp1(&filab[102],"PHS")==0){
      for(l1=0;l1<*nk;l1++){
	for(l2=0;l2<6;l2++){
	  l=6*l1+l2;
	  stnreal=stn[l];
	  stnimag=stn[l+6**nk];
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

    /* determine the time for the worst principal stress anywhere
       in the structure; the worst principal stress is the maximum
       of the absolute value of the principal stresses */

    if((strcmp1(&filab[108],"MAXU")==0)||(strcmp1(&filab[114],"MAXS")==0)){

      worstpsmax=0.;
      iworsttime=0;
      for(l1=0;l1<*nk;l1++){
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
	    iworsttime=l3;
	  }
	}
      }
    }

    /* storing the displacements at the time of the worst
       principal stress anywhere in the structure */

    if(strcmp1(&filab[108],"MAXU")==0){

      for(l1=0;l1<4**nk;l1++){vmax[l1]=0.;}
      for(l1=0;l1<*nk;l1++){
	  ctl=cos(iworsttime/constant);
	  stl=sin(iworsttime/constant);
	  for(l2=1;l2<4;l2++){
	    l=5*l1+l2;
	    vl[l2]=ctl*v[l]-stl*v[l+5**nk];
	  }
	  dd=vl[1]*vl[1]+vl[2]*vl[2]+vl[3]*vl[3];
	  vmax[4*l1]=dd;
	  vmax[4*l1+1]=vl[1];
	  vmax[4*l1+2]=vl[2];
	  vmax[4*l1+3]=vl[3];
      }
      for(l1=0;l1<*nk;l1++){vmax[4*l1]=sqrt(vmax[4*l1]);}
    }

    /* storing the stresses at the time of the worst
       principal stress anywhere in the structure */
        
    if(strcmp1(&filab[114],"MAXS")==0){
      for(l1=0;l1<*nk;l1++){
	ctl=cos(iworsttime/constant);
	stl=sin(iworsttime/constant);
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
	stnmax[7*l1]=dd;
	for(l2=1;l2<7;l2++){
	  stnmax[7*l1+l2]=stnl[l2-1];
	}
      }
    }
    
    /* mapping the results to the other sectors */

    for(l=0;l<*nk;l++){inumt[l]=inum[l];}

    icntrl=2;imag=1;

    FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag));

    if(strcmp1(&filab[0],"U   ")==0)
      for(l=0;l<5**nk;l++){vt[l]=v[l];};
    if(strcmp1(&filab[6],"NT  ")==0)
      for(l=0;l<*nk;l++){t1t[l]=t1[l];};
    if(strcmp1(&filab[12],"S   ")==0)
      for(l=0;l<6**nk;l++){stnt[l]=stn[l];};
    if(strcmp1(&filab[18],"E   ")==0)
      for(l=0;l<6**nk;l++){eent[l]=een[l];};
    if(strcmp1(&filab[24],"RF  ")==0)
      for(l=0;l<4**nk;l++){fnt[l]=fn[l];};
    if(strcmp1(&filab[36],"ENER")==0)
      for(l=0;l<*nk;l++){enernt[l]=enern[l];};

    for(jj=0;jj<*mcs;jj++){
      ilength=cs[17*jj+3];
      is=cs[17*jj+4];
      lprev=cs[17*jj+13];
      for(i=1;i<is;i++){
        
        for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
        
        theta=i*nm*2.*pi/cs[17*jj];
        ctl=cos(theta);
        stl=sin(theta);
        
        if(strcmp1(&filab[0],"U   ")==0){
          for(l1=0;l1<*nk;l1++){
            if(inocs[l1]==jj){

              /* check whether node lies on axis */

	      ml1=-l1-1;
              FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
	      if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		      for(l2=0;l2<4;l2++){
			  l=5*l1+l2;
			  vt[l+5**nk*i]=v[l];
		      }
		      continue;
		  }
	      }
              for(l2=0;l2<4;l2++){
                l=5*l1+l2;
                vt[l+5**nk*i]=ctl*v[l]-stl*v[l+5**nk];
              }
            }
          }
        }
        
        if(strcmp1(&filab[6],"NT  ")==0){
          for(l=0;l<*nk;l++){
	      if(inocs[l]==jj) t1t[l+*nk*i]=t1[l];
          }
        }
        
        if(strcmp1(&filab[12],"S   ")==0){
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
        
        if(strcmp1(&filab[18],"E   ")==0){
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
        
        if(strcmp1(&filab[24],"RF  ")==0){
          for(l1=0;l1<*nk;l1++){
            if(inocs[l1]==jj){

              /* check whether node lies on axis */

	      ml1=-l1-1;
              FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
	      if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		      for(l2=0;l2<4;l2++){
			  l=4*l1+l2;
			  fnt[l+4**nk*i]=fn[l];
		      }
		      continue;
		  }
	      }
              for(l2=0;l2<4;l2++){
                l=4*l1+l2;
                fnt[l+4**nk*i]=ctl*fn[l]-stl*fn[l+4**nk];
              }
            }
          }
        }
        
        if(strcmp1(&filab[36],"ENER")==0){
          for(l=0;l<*nk;l++){
            if(inocs[l]==jj) enernt[l+*nk*i]=0.;
          }
        }
      }
    }

    /* calculating the maxima */

    /*  if(strcmp1(&filab[108],"MAXU")==0){

      for(l1=0;l1<4**nk;l1++){vmax[l1]=0.;}
      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj];
	lprev=cs[17*jj+13];
	for(i=0;i<is;i++){
	  theta=i*nm*2.*pi/cs[17*jj];
	  ctl=cos(theta);
	  stl=sin(theta);
          for(l1=0;l1<*nk;l1++){
	  if(inocs[l1]==jj){*/
	      
              /* check whether node lies on axis */
	      
    /*     ml1=-l1-1;
              FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
	      ionaxis=0;
	      if(id!=0){
		if(ics[lprev+id-1]==ml1){
		  for(l2=1;l2<4;l2++){
		    l=5*l1+l2;
		    vl[l2]=v[l];
		  }
		  ionaxis=1;
		}
	      }
	      if(ionaxis==0){
		for(l2=1;l2<4;l2++){
		  l=5*l1+l2;
		  vl[l2]=ctl*v[l]-stl*v[l+5**nk];
		}
	      }
	      dd=vl[1]*vl[1]+vl[2]*vl[2]+vl[3]*vl[3];
	      if(dd>vmax[4*l1]){
		vmax[4*l1]=dd;
		vmax[4*l1+1]=vl[1];
		vmax[4*l1+2]=vl[2];
		vmax[4*l1+3]=vl[3];
	      }
            }
          }
        }
      }
      for(l1=0;l1<*nk;l1++){vmax[4*l1]=sqrt(vmax[4*l1]);}
      }*/
        
    /*   if(strcmp1(&filab[114],"MAXS")==0){

      for(l1=0;l1<7**nk;l1++){stnmax[l1]=0.;}
      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj];
	lprev=cs[17*jj+13];
	for(i=0;i<is;i++){
	  theta=i*nm*2.*pi/cs[17*jj];
	  ctl=cos(theta);
	  stl=sin(theta);
	  for(l1=0;l1<*nk;l1++){
	  if(inocs[l1]==jj){*/
	      
	      /* check whether node lies on axis */
	      
    /*     ml1=-l1-1;
	      FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
	      ionaxis=0;
	      if(id!=0){
		if(ics[lprev+id-1]==ml1){
		  for(l2=0;l2<6;l2++){
		    l=6*l1+l2;
		    stnl[l2]=stn[l];
		  }
		  ionaxis=1;
		}
	      }
	      if(ionaxis==0){
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnl[l2]=ctl*stn[l]-stl*stn[l+6**nk];
		}
	      }
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
	      dd=al[0];
	      if(fabs(al[1])>fabs(dd)) dd=al[1];
	      if(fabs(al[2])>fabs(dd)) dd=al[2];
	      if(fabs(dd)>fabs(stnmax[7*l1])){
		stnmax[7*l1]=dd;
		for(l2=1;l2<7;l2++){
		  stnmax[7*l1+l2]=stnl[l2-1];
		}
	      }
	    }
	  }
	}
      }
      }*/
	
    icntrl=-2;imag=0;

    FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
      &imag));

    ++*kode;
    freq=d[j]/6.283185308;
    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
    FORTRAN(out,(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,kode,
       filab,eent,t1t,fnt,&freq,epn,ielmatt,matname,enernt,xstaten,nstate_,
       &istep,&iinc,iperturb,ener,mint_,output,ithermal,qfn,&j,&nm,
       trab,inotrt,ntrans,orab,ielorien,norien,description,
       ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&idummy));
    free(ipneigh);free(neigh);

    /*  if(((strcmp1(&filab[60],"PU")==0)||(strcmp1(&filab[102],"PHS")==0)||
        (strcmp1(&filab[108],"MAXU")==0)||(strcmp1(&filab[114],"MAXS")==0))&&
       (2*(j/2)==j)){
      (*kode)++;
      FORTRAN(frdphase,(kode,&freq,nk,inum,vr,vi,stnr,stni,filab,&j,&nm,
			nmethod,vmax,stnmax,&nkcoords));
			}*/

  }

  if((fmax>0.)&&(fmax>d[nev-1])){
    printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
  }

  free(adb);free(aub);free(temp_array);free(coefmpcnew);

  free(v);free(fn);free(stn);free(inum);free(stx);free(resid);
  free(z);free(workd);free(workl);free(select);free(d);

  if(strcmp1(&filab[18],"E   ")==0) free(een);
  if(strcmp1(&filab[36],"ENER")==0) free(enern);

  if(strcmp1(&filab[0],"U   ")==0) free(vt);
  if(strcmp1(&filab[6],"NT  ")==0) free(t1t);
  if(strcmp1(&filab[12],"S   ")==0) free(stnt);
  if(strcmp1(&filab[18],"E   ")==0) free(eent);
  if(strcmp1(&filab[24],"RF  ")==0) free(fnt);
  if(strcmp1(&filab[36],"ENER")==0) free(enernt);

  free(cot);free(kont);free(ipkont);free(lakont);free(inumt);free(ielmatt);
  if(*ntrans>0){free(inotrt);}

  if(mei[3]==1){
      (*nevtot)+=nev;
      fclose(f1);
  }

  }

  free(inocs);free(ielcs);free(xstiff);
  free(ipobody);

  if(strcmp1(&filab[60],"PU")==0){free(vr);free(vi);}
  if(strcmp1(&filab[102],"PHS")==0){free(stnr);free(stni);}
  if(strcmp1(&filab[108],"MAXU")==0){free(vmax);}
  if(strcmp1(&filab[114],"MAXS")==0){free(stnmax);}

  for(i=0;i<6**mint_**ne;i++){eme[i]=0.;}

  return;
}

#endif
