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
#ifdef MATRIXSTORAGE
   #include "matrixstorage.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif

void arpack(double *co, int *nk, int *kon, int *ipkon, char *lakon,
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
	     double *shcon, int *nshcon, double *cocon, int *ncocon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,
	     int *mei, double *fei,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mi,
             int *ncmat_, int *nstate_, double *ener, char *jobnamec,
             char *output, char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, int *isolver, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
             char *cbody, int *ibody,double *xbody, int *nbody){

  /* calls the Arnoldi Package (ARPACK) */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A", fneig[132]="",
      description[13]="            ";

  int *inum=NULL,k,ido,dz,iparam[11],ipntr[11],lworkl,ngraph=1,im,
    info,rvec=1,*select=NULL,lfin,j,lint,iout,ielas=1,icmd=0,mt=mi[1]+1,
    iinc=1,istep=1,nev,ncv,mxiter,jrow,*ipobody=NULL,inewton=0,ifreebody,
    mass[2]={1,1}, stiffness=1, buckling=0, rhsi=0, intscheme=0,noddiam=-1,
    coriolis=0,symmetryflag=0,inputformat=0,*ipneigh=NULL,*neigh=NULL;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,
    *workl=NULL,*aux=NULL,*d=NULL,sigma=1,*temp_array=NULL,
    *een=NULL,sum,cam[5],*f=NULL,*fn=NULL,qa[3],*fext=NULL,*epn=NULL,
    *xstateini=NULL,*xstiff=NULL,*stiini=NULL,*vini=NULL,freq,*stx=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,
    *physcon=NULL,*qfx=NULL,*qfn=NULL,tol,fmin,fmax,pi,*cgr=NULL,
    *xloadold=NULL,reltime,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*cs=NULL,*springarea=NULL;

  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

#ifdef SGI
  int token;
#endif

  if((strcmp1(&filab[870],"PU  ")==0)||
     (strcmp1(&filab[1479],"PHS ")==0)||
     (strcmp1(&filab[1566],"MAXU")==0)||
     (strcmp1(&filab[1653],"MAXS")==0)){
      printf("*ERROR in arpack: PU, PHS, MAXU and MAX was selected in a frequency calculation without cyclic symmetry;\n this is not correct\n");
      FORTRAN(stop,());
  }

  /* copying the frequency parameters */

  pi=4.*atan(1.);

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

  /* field for initial values of state variables (needed if
     previous static step was viscoplastic */

  if(*nstate_!=0){
    xstateini=NNEW(double,*nstate_*mi[0]**ne);
    for(k=0;k<*nstate_*mi[0]**ne;++k){
      xstateini[k]=xstate[k];
    }
  }

  /* determining the internal forces and the stiffness coefficients */

  f=NNEW(double,neq[1]);

  /* allocating a field for the stiffness matrix */

  xstiff=NNEW(double,27*mi[0]**ne);

  iout=-1;
  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stx=NNEW(double,6*mi[0]**ne);
  if(*ithermal>1){
      qfx=NNEW(double,3*mi[0]**ne);
  }
  if(*iperturb==0){
     FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t0,ithermal,
	       prestr,iprestr,filab,eme,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               sti,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
               fmpc,nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea));
  }else{
     FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	       prestr,iprestr,filab,eme,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xboun,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               sti,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
               fmpc,nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea));
  }
  free(f);free(v);free(fn);free(stx);if(*ithermal>1){free(qfx);}
  iout=1;

  /* filling in the matrix */

  ad=NNEW(double,neq[1]);
  au=NNEW(double,nzs[2]);

  adb=NNEW(double,neq[1]);
  aub=NNEW(double,nzs[1]);

  fext=NNEW(double,neq[1]);

  if(*iperturb==0){
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
	      physcon,shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,
		      &coriolis,ibody,xloadold,&reltime,veold,springarea));
  }
  else{
    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
	      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
	      ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	      ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,
	      t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
              physcon,shcon,nshcon,cocon,ncocon,ttime,&time,&istep,&iinc,
		      &coriolis,ibody,xloadold,&reltime,veold,springarea));
  }

  free(fext);

  if(*nmethod==0){

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;time=0.;
    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }
    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1,
         fn,&time,epn,ielmat,matname,enern,xstaten,nstate_,&istep,&iinc,
	 iperturb,ener,mi,output,ithermal,qfn,&j,&noddiam,trab,inotr,ntrans,
         orab,ielorien,norien,description,
	 ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ne,cs,
	 set,nset,istartset,iendset,ialset));
    
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    free(inum);FORTRAN(stop,());

  }

  /* LU decomposition of the left hand matrix */

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
  else if(*isolver==6){
#ifdef MATRIXSTORAGE
    matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
		  ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi);
#else
    printf("*ERROR in arpack: the MATRIXSTORAGE library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1]);
#else
    printf("*ERROR in arpack: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

/*  free(au);free(ad);*/

/* calculating the eigenvalues and eigenmodes */
  
  printf(" Calculating the eigenvalues and the eigenmodes\n\n");

  ido=0;
  dz=neq[1];
  iparam[0]=1;
  iparam[2]=mxiter;
  iparam[3]=1;
  iparam[6]=3;

  lworkl=ncv*(8+ncv);
  info=0;

  resid=NNEW(double,neq[1]);
  long long zsize=ncv*neq[1];
  z=NNEW(double,zsize);
  workd=NNEW(double,3*neq[1]);
  workl=NNEW(double,lworkl);

  FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
	  workl,&lworkl,&info));

  temp_array=NNEW(double,neq[1]);

  while((ido==-1)||(ido==1)||(ido==2)){
    if(ido==-1){
      FORTRAN(op,(&neq[1],aux,&workd[ipntr[0]-1],temp_array,adb,aub,icol,irow,nzl));
    }
    if((ido==-1)||(ido==1)){

      /* solve the linear equation system  */

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
          pardiso_solve(temp_array,&neq[1]);
#endif
        }
        for(jrow=0;jrow<neq[1];jrow++){
          workd[ipntr[1]-1+jrow]=temp_array[jrow];
        }
      }
      else if(ido==1){
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
        if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
        }
        for(jrow=0;jrow<neq[1];jrow++){
          workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
        }
      }

    }

    if(ido==2){
      FORTRAN(op,(&neq[1],aux,&workd[ipntr[0]-1],&workd[ipntr[1]-1],adb,aub,icol,irow,nzl));
    }

    FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&dz,iparam,ipntr,workd,
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
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_cleanup(&neq[1]);
#endif
  }

  if(info!=0){
    printf("*ERROR in arpack: info=%d\n",info);
    printf("       # of converged eigenvalues=%d\n\n",iparam[4]);
  }         

  select=NNEW(int,ncv);
  d=NNEW(double,nev);

  FORTRAN(dseupd,(&rvec,howmny,select,d,z,&dz,&sigma,bmat,&neq[1],which,&nev,&tol,resid,
	  &ncv,z,&dz,iparam,ipntr,workd,workl,&lworkl,&info));

  FORTRAN(writeev,(d,&nev,&fmin,&fmax));

  /* writing the eigenvalues and mass matrix to a binary file */

  if(mei[3]==1){

      strcpy(fneig,jobnamec);
      strcat(fneig,".eig");
      
      if((f1=fopen(fneig,"wb"))==NULL){
	  printf("*ERROR in arpack: cannot open eigenvalue file for writing...");

	  exit(0);
      }
      if(fwrite(&nev,sizeof(int),1,f1)!=1){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }

      /* the eigenfrequencies are stores as radians/time */

      if(fwrite(d,sizeof(double),nev,f1)!=nev){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
      if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf("*ERROR saving data to the eigenvalue file...");
	  exit(0);
      }
  }

  free(au);free(ad);

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stn=NNEW(double,6**nk);
  inum=NNEW(int,*nk);
  stx=NNEW(double,6*mi[0]**ne);
  if(*ithermal>1){
      qfn=NNEW(double,3**nk);
      qfx=NNEW(double,3*mi[0]**ne);
  }
  
  if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);

  temp_array=NNEW(double,neq[1]);

  lfin=0;
  for(j=0;j<nev;++j){
    lint=lfin;
    lfin=lfin+neq[1];

    for(k=0;k<6*mi[0]**ne;k++){eme[k]=0.;}

    /* check whether the frequency belongs to the requested
       interval */

    if(fmin>d[j]) continue;
    if(fmax>0.){
      if(fmax<d[j]) break;
    }

    if(*nprint>0) FORTRAN(writehe,(&j));

    sum=0.;
    for(k=0;k<neq[1];++k)
      temp_array[k]=0.;
    FORTRAN(op,(&neq[1],aux,&z[lint],temp_array,adb,aub,icol,irow,nzl));
    for(k=0;k<neq[1];++k)
      sum+=z[lint+k]*temp_array[k];
    for(k=0;k<neq[1];++k)
      z[lint+k]=z[lint+k]/sqrt(sum);

    if(mei[3]==1){
	if(fwrite(&z[lint],sizeof(double),neq[1],f1)!=neq[1]){
	    printf("*ERROR saving data to the eigenvalue file...");
	    exit(0);
	}
    }

    eei=NNEW(double,6*mi[0]**ne);
    if(*nener==1){
	stiini=NNEW(double,6*mi[0]**ne);
	enerini=NNEW(double,mi[0]**ne);}

//    memset(&v[0],0.,sizeof(double)*mt**nk);
    DMEMSET(v,0,mt**nk,0.);
    if(*iperturb==0){
      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	    stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t0,ithermal,
	    prestr,iprestr,filab,eme,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,&z[lint],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
            sti,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
            nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea));}
    else{
      FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	    stx,elcon,
	    nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	    norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,&z[lint],
	    nodeboun,ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,sti,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
            nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea));
    }
    free(eei);
    if(*nener==1){
	free(stiini);free(enerini);}

    ++*kode;
    freq=d[j]/6.283185308;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }
    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1,
         fn,&freq,epn,ielmat,matname,enern,xstaten,nstate_,&istep,&iinc,
	 iperturb,ener,mi,output,ithermal,qfn,&j,&noddiam,
         trab,inotr,ntrans,orab,ielorien,norien,description,
	 ipneigh,neigh,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ne,cs,
         set,nset,istartset,iendset,ialset));
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
  }

  if((fmax>0.)&&(fmax>d[nev-1])){
    printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
  }

  if(mei[3]==1){
      fclose(f1);
  }

  free(adb);free(aub);free(temp_array);

  free(v);free(fn);free(stn);free(inum);free(stx);free(resid);
  free(z);free(workd);free(workl);free(select);free(d);free(xstiff);
  free(ipobody);

  if(*ithermal>1){free(qfn);free(qfx);}

  if(*nstate_!=0){free(xstateini);}

  if(strcmp1(&filab[261],"E   ")==0) free(een);
  if(strcmp1(&filab[522],"ENER")==0) free(enern);

  for(k=0;k<6*mi[0]**ne;k++){eme[k]=0.;}
    
  return;
}

#endif
