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

void linstatic(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,
	     int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int **icolp, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, 
	     int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old,
	     int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,
	     char *filab, double *eme,
             int *iexpl, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *isolver,
             int *mi, int *ncmat_, int *nstate_, double *cs, int *mcs,
             int *nkon, double *ener, double *xbounold,
	     double *xforcold, double *xloadold,
             char *amname, double *amta, int *namta,
	     int *nam, int *iamforc, int *iamload,
             int *iamt1, int *iamboun, double *ttime, char *output, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, double *trab, 
             int *inotr, int *ntrans, double *fmpc, char *cbody, int *ibody,
	     double *xbody, int *nbody, double *xbodyold, double *tper,
	     double *thicke, char *jobnamec,char *tieset,int *ntie,
             int *istep){
  
  char description[13]="            ";

  int *inum=NULL,k,*icol=NULL,*irow=NULL,ielas,icmd=0,iinc=1,nasym=0,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,*ncocon=NULL,
      *nshcon=NULL,mode=-1,noddiam=-1,*ipobody=NULL,inewton=0,coriolis=0,iout,
      ifreebody,*itg=NULL,ntg=0,symmetryflag=0,inputformat=0,ngraph=1,
      mt=mi[1]+1,ne0,*integerglob=NULL,iglob=0,*ipneigh=NULL,*neigh=NULL,
      icfd=0,*inomat=NULL;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,
         *f=NULL,*fn=NULL,qa[3],*fext=NULL,*epn=NULL,*xstateini=NULL,
         *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
         *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
         *enerini=NULL,*cocon=NULL,*shcon=NULL,*physcon=NULL,*qfx=NULL,
         *qfn=NULL,sigma=0.,*cgr=NULL,*xbodyact=NULL,*vr=NULL,*vi=NULL,
         *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*springarea=NULL,
         *eenmax=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*xnormastface=NULL,
         *emeini=NULL,*doubleglob=NULL;

#ifdef SGI
  int token;
#endif

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time,reltime=1.;

  icol=*icolp;
  irow=*irowp;

  time=*tper;
  dtime=*tper;

  ne0=*ne;

  /* determining the global values to be used as boundary conditions
     for a submodel */

  getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob);

  /* allocating fields for the actual external loading */

  xbounact=NNEW(double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  xforcact=NNEW(double,*nforc);
  xloadact=NNEW(double,2**nload);
  xbodyact=NNEW(double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  if(*ithermal==1){
    t1act=NNEW(double,*nk);
    for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
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
  }

  /* allocating a field for the instantaneous amplitude */

  ampli=NNEW(double,*nam);

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
  *ttime=*ttime+*tper;

  /* determining the internal forces and the stiffness coefficients */

  f=NNEW(double,*neq);

  /* allocating a field for the stiffness matrix */

  xstiff=NNEW(double,(long long)27*mi[0]**ne);

  iout=-1;
  v=NNEW(double,mt**nk);
  fn=NNEW(double,mt**nk);
  stx=NNEW(double,6*mi[0]**ne);
  inum=NNEW(int,*nk);
  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
   	       prestr,iprestr,filab,eme,emn,een,iperturb,
	       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	       ndirboun,xbounact,nboun,ipompc,
	       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
               &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
               emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
               iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	       fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xloadact,xloadold,&icfd,inomat);
  free(v);free(fn);free(stx);free(inum);
  iout=1;
  
  /* determining the system matrix and the external forces */

  ad=NNEW(double,*neq);
  au=NNEW(double,*nzs);
  fext=NNEW(double,*neq);

  FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
	    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	    nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	    nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
	    ikmpc,ilmpc,ikboun,ilboun,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,
	    t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	    nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	    xstiff,npmat_,&dtime,matname,mi,
            ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
            shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,&coriolis,
	    ibody,xloadold,&reltime,veold,springarea,nstate_,
            xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	    tieset,istartset,iendset,ialset,ntie,&nasym));

  /* determining the right hand side */

  b=NNEW(double,*neq);
  for(k=0;k<*neq;++k){
      b[k]=fext[k]-f[k];
  }
  free(fext);free(f);

  if(*nmethod!=0){

    if(*isolver==0){
#ifdef SPOOLES
      spooles(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,&symmetryflag,
              &inputformat,&nzs[2]);
#else
            printf("*ERROR in linstatic: the SPOOLES library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if((*isolver==2)||(*isolver==3)){
      preiter(ad,&au,b,&icol,&irow,neq,nzs,isolver,iperturb);
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,token);
#else
            printf("*ERROR in linstatic: the SGI library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,neq,nzs);
#else
            printf("*ERROR in linstatic: the TAUCS library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
            printf("*ERROR in linstatic: the PARDISO library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }

    free(ad);free(au);

    /* calculating the displacements and the stresses and storing */
    /* the results in frd format for each valid eigenmode */

    v=NNEW(double,mt**nk);
    fn=NNEW(double,mt**nk);
    stn=NNEW(double,6**nk);
    inum=NNEW(int,*nk);
    stx=NNEW(double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
    if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);
    if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);

    eei=NNEW(double,6*mi[0]**ne);
    if(*nener==1){
	stiini=NNEW(double,6*mi[0]**ne);
	enerini=NNEW(double,mi[0]**ne);}

    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,b,nodeboun,ndirboun,xbounact,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,&bet,
            &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
            ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat);

    free(eei);
    if(*nener==1){
	free(stiini);free(enerini);}

    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
    memcpy(&sti[0],&stx[0],sizeof(double)*6*mi[0]**ne);
/*    for(k=0;k<mt**nk;++k){
      vold[k]=v[k];
    }
    for(k=0;k<6*mi[0]**ne;++k){
      sti[k]=stx[k];
      }*/

    ++*kode;

    /* for cyclic symmetric sectors: duplicating the results */

    if(*mcs>0){
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1act,
		   fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
                   nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,
                   qfn,ialset,istartset,iendset,trab,inotr,ntrans,orab,
	           ielorien,norien,sti,veold,&noddiam,set,nset,emn,thicke,
	           jobnamec,&ne0);
    }
    else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	    neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	}
	frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);
	if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    }

    free(v);free(stn);free(inum);
    free(b);free(stx);free(fn);

    if(strcmp1(&filab[261],"E   ")==0) free(een);
    if(strcmp1(&filab[522],"ENER")==0) free(enern);
    if(strcmp1(&filab[2697],"ME  ")==0) free(emn);

  }
  else {

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;
    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
	neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
    }
    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);
    if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    free(inum);FORTRAN(stop,());

  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){xbounold[k]=xbounact[k];}
  for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    for(k=0;k<*nk;++k){vold[mt*k]=t1act[k];}
  }

  free(xbounact);free(xforcact);free(xloadact);free(t1act);free(ampli);
  free(xbodyact);if(*nbody>0) free(ipobody);free(xstiff);

  if(iglob==1){free(integerglob);free(doubleglob);}

  *icolp=icol;
  *irowp=irow;
 
  return;
}
