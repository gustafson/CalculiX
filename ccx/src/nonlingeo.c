/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                          */

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


void nonlingeo(double **cop, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	     ITG *ne, 
	     ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	     ITG **ipompcp, ITG **nodempcp, double **coefmpcp, char **labmpcp,
             ITG *nmpc, 
	     ITG *nodeforc, ITG *ndirforc,double *xforc, ITG *nforc, 
	     ITG *nelemload, char *sideload, double *xload,ITG *nload, 
	     ITG *nactdof, 
	     ITG **icolp, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, 
	     ITG *nmethod, ITG **ikmpcp, ITG **ilmpcp, ITG *ikboun, 
	     ITG *ilboun,
             double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	     double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	     ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     ITG *ithermal,double *prestr, ITG *iprestr, 
	     double **voldp,ITG *iperturb, double *sti, ITG *nzs,  
	     ITG *kode, char *filab, 
             ITG *idrct, ITG *jmax, ITG *jout, double *tinc,
             double *tper, double *tmin, double *tmax, double *eme,
	     double *xbounold, double *xforcold, double *xloadold,
             double *veold, double *accold,
	     char *amname, double *amta, ITG *namta, ITG *nam,
             ITG *iamforc, ITG *iamload,
             ITG *iamt1, double *alpha, ITG *iexpl,
	     ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
             ITG *nplkcon,
             double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
             char *matname, double *qaold, ITG *mi,
             ITG *isolver, ITG *ncmat_, ITG *nstate_, ITG *iumat,
             double *cs, ITG *mcs, ITG *nkon, double **enerp, ITG *mpcinfo,
             char *output,
             double *shcon, ITG *nshcon, double *cocon, ITG *ncocon,
             double *physcon, ITG *nflow, double *ctrl, 
             char *set, ITG *nset, ITG *istartset,
             ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
             char *prset, ITG *nener,ITG *ikforc,ITG *ilforc, double *trab, 
             ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
             ITG *ibody, double *xbody, ITG *nbody, double *xbodyold,
             ITG *ielprop, double *prop, ITG *ntie, char *tieset,
	     ITG *itpamp, ITG *iviewfile, char *jobnamec, double *tietol,
	     ITG *nslavs, double *thicke, ITG *ics, 
	     ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
      *sideface=NULL,*labmpc=NULL,fnffrd[132]=""; 
 
  ITG *inum=NULL,k,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
       icutb=0,istab=0,ifreebody,uncoupled,n1,n2,
       iperturb_sav[2],ilin,*icol=NULL,*irow=NULL,ielas=0,icmd=0,
       memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
       *nodempcref=NULL,memmpcref_,mpcfreeref,*itg=NULL,*ineighe=NULL,
       *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,
       *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,
       ntrit,*inocs=NULL,inewton=0,*ipobody=NULL,*nacteq=NULL,
       *nactdog=NULL,nteq,network,*itietri=NULL,*koncont=NULL,
       ncont,ne0,nkon0,*ipkon=NULL,*kon=NULL,*ielorien=NULL,
       *ielmat=NULL,inext,itp=0,symmetryflag=0,inputformat=0,
       *iruc=NULL,iitterm=0,iturbulent,ngraph=1,ismallsliding=0,
       *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
       *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,*imastnode=NULL,
       *nslavnode=NULL,*nmastnode=NULL,*imastop=NULL,
       *iponoels=NULL,*inoels=NULL,*islavsurfold=NULL,
       *islavact=NULL,mt=mi[1]+1,*nactdofinv=NULL,*ipe=NULL, 
       *ime=NULL,*ikactmech=NULL,nactmech,inode,idir,neold,
       iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
       *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL;

  ITG mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
    coriolis=0,*ipneigh=NULL,*neigh=NULL,maxprevcontel,nslavs_prev_step,
    *nelemface=NULL,*ipoface=NULL,*nodface=NULL,*ifreestream=NULL,iex,
    *isolidsurf=NULL,*neighsolidsurf=NULL,*iponoel=NULL,*inoel=NULL,
    nef=0,nface,nfreestream,nsolidsurf,i,indexe,icfd=0,id,
    node,networknode,iflagact=0,*nodorig=NULL,*ipivr=NULL,iglob=0,
    *inomat=NULL,*ipnei=NULL,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,
    *neifa=NULL,*neiel=NULL,*ielfa=NULL,*ifaext=NULL,nflnei,nfaext;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cg=NULL,
         *cdn=NULL,*vel=NULL,*vfa=NULL,*pslavsurfold=NULL,
         *f=NULL,*fn=NULL,qa[3]={0.,0.,-1.},qam[2]={0.,0.},dtheta,theta,
	 err,ram[8]={0.,0.,0.,0.,0.,0.,0.,0.},*areaslav=NULL,
         *springarea=NULL,ram1[8]={0.,0.,0.,0.,0.,0.,0.,0.},
	 ram2[8]={0.,0.,0.,0.,0.,0.,0.,0.},deltmx,ptime,
         uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,*straight=NULL,
	 *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
	 *xforcact=NULL,*xloadact=NULL,*fext=NULL,*clearini=NULL,
         reltime,time,bet=0.,gam=0.,*aux1=NULL,*aux2=NULL,dtime,*fini=NULL,
         *fextini=NULL,*veini=NULL,*accini=NULL,*xstateini=NULL,
	 *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,
         *xbounini=NULL,dev,*xstiff=NULL,*stx=NULL,*stiini=NULL,
         *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
	 *coefmpcref=NULL,*enerini=NULL,*emn=NULL,
	 *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
	 *adview=NULL,*auview=NULL,*qfx=NULL,
         *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
         *xbodyact=NULL,*cgr=NULL,dthetaref, *vcontu=NULL,*vr=NULL,*vi=NULL,
	 *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
         *f_cm=NULL, *f_cs=NULL,
	 *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
	 *xmastnor=NULL,*emeini=NULL,
	 *doubleglob=NULL,*xnoels=NULL,*au=NULL,
	 *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
	 *x=NULL,*y=NULL,*z=NULL,*xo=NULL,
         *yo=NULL,*zo=NULL,*cdnr=NULL,*cdni=NULL;

#ifdef SGI
  ITG token;
#endif
  
  icol=*icolp;irow=*irowp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  if(*ithermal==4){
      uncoupled=1;
      *ithermal=3;
  }else{
      uncoupled=0;
  }
  
  if(*mortar==0){
      maxprevcontel=*nslavs;
  }else if(*mortar==1){
      maxprevcontel=*nintpoint;
  }
  nslavs_prev_step=*nslavs;

  /* turbulence model 
     iturbulent==0: laminar
     iturbulent==1: k-epsilon
     iturbulent==2: q-omega
     iturbulent==3: SST */
  
  iturbulent=(ITG)physcon[8];
  
  for(k=0;k<3;k++){
      strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];

  if((icascade==2)&&(*iexpl>=2)){
      printf("*ERROR in nonlingeo: linear and nonlinear MPC's depend on each other\n");
      printf("       This is not allowed in a explicit dynamic calculation\n");
      FORTRAN(stop,());
  }
  
  /* check whether the submodel is meant for a fluid-structure
     interaction */
  
  strcpy(fnffrd,jobnamec);
  strcat(fnffrd,"f.frd");
  if((jobnamec[396]!='\0')&&(strcmp1(&jobnamec[396],fnffrd)==0)){
      
      /* fluid-structure interaction: wait till after the compfluid
         call */
      
      integerglob=NNEW(ITG,1);
      doubleglob=NNEW(double,1);
  }else{
      
      /* determining the global values to be used as boundary conditions
	 for a submodel */
      
      getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,xboun,
		       nload,sideload,iamload,&iglob);
  }
  
  /* invert nactdof */
  
  nactdofinv=NNEW(ITG,mt**nk);nodorig=NNEW(ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  free(nodorig);
  
  /* allocating a field for the stiffness matrix */
  
  xstiff=NNEW(double,(long long)27*mi[0]**ne);
  
  /* allocating force fields */
  
  f=NNEW(double,neq[1]);
  fext=NNEW(double,neq[1]);
  
  b=NNEW(double,neq[1]);
  vini=NNEW(double,mt**nk);
  
  aux=NNEW(double,7*maxlenmpc);
  iaux=NNEW(ITG,maxlenmpc);
  
  /* allocating fields for the actual external loading */
  
  xbounact=NNEW(double,*nboun);
  xbounini=NNEW(double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  xforcact=NNEW(double,*nforc);
  xloadact=NNEW(double,2**nload);
  xbodyact=NNEW(double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  
  /* assigning the body forces to the elements */ 
  
  if(*nbody>0){
      ifreebody=*ne+1;
      ipobody=NNEW(ITG,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,ITG,2*(*ne+ifreebody));
      }
      RENEW(ipobody,ITG,2*(ifreebody-1));
      if(inewton==1){cgr=NNEW(double,4**ne);}
  }
  
  /* for mechanical calculations: updating boundary conditions
     calculated in a previous thermal step */
  
  if(*ithermal<2) FORTRAN(gasmechbc,(vold,nload,sideload,
				     nelemload,xload,mi));
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
      itg=NNEW(ITG,*nload+3**nflow);
      ieg=NNEW(ITG,*nflow);
      /* max 6 triangles per face, 4 entries per triangle */
      kontri=NNEW(ITG,24**nload);
      nloadtr=NNEW(ITG,*nload);
      nacteq=NNEW(ITG,4**nk);
      nactdog=NNEW(ITG,4**nk);
      v=NNEW(double,mt**nk);
      FORTRAN(envtemp,(itg,ieg,&ntg,&ntr,sideload,nelemload,
		       ipkon,kon,lakon,ielmat,ne,nload,
                       kontri,&ntri,nloadtr,nflow,ndirboun,nactdog,
                       nodeboun,nacteq,nboun,ielprop,prop,&nteq,
                       v,&network,physcon,shcon,ntmat_,co,
                       vold,set,nshcon,rhcon,nrhcon,mi,nmpc,nodempc,
                       ipompc,labmpc,ikboun,&nasym));
      free(v);
      
      if((*mcs>0)&&(ntr>0)){
	  inocs=NNEW(ITG,*nk);
	  radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
		 iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
      }
      else{ntrit=ntri;}
      
      nzsrad=100*ntr;
      mast1rad=NNEW(ITG,nzsrad);
      irowrad=NNEW(ITG,nzsrad);
      icolrad=NNEW(ITG,ntr);
      jqrad=NNEW(ITG,ntr+1);
      ipointerrad=NNEW(ITG,ntr);
      
      if(ntr>0){
	  mastructrad(&ntr,nloadtr,sideload,ipointerrad,
		      &mast1rad,&irowrad,&nzsrad,
		      jqrad,icolrad);
      }
      
      free(ipointerrad);free(mast1rad);
      RENEW(irowrad,ITG,nzsrad);
      
      RENEW(itg,ITG,ntg);
      ineighe=NNEW(ITG,ntg);
      RENEW(kontri,ITG,4*ntrit);
      RENEW(nloadtr,ITG,ntr);
      
      adview=NNEW(double,ntr);
      auview=NNEW(double,2*nzsrad);
      tarea=NNEW(double,ntr);
      tenv=NNEW(double,ntr);
      fenv=NNEW(double,ntr);
      erad=NNEW(double,ntr);
      
      ac=NNEW(double,nteq*nteq);
      bc=NNEW(double,nteq);
      ipiv=NNEW(ITG,nteq);
      adrad=NNEW(double,ntr);
      aurad=NNEW(double,2*nzsrad);
      bcr=NNEW(double,ntr);
      ipivr=NNEW(ITG,ntr);
  }
  
  /* check for fluid elements */
  
  for(i=0;i<*ne;++i){
      if(ipkon[i]<0) continue;
      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i],"F")==0){icfd=1;nef++;}
  }
  if(icfd==1){
      ipoface=NNEW(ITG,*nk);
      nodface=NNEW(ITG,5*6*nef);
      ipnei=NNEW(ITG,*ne);
      neifa=NNEW(ITG,6**ne);
      neiel=NNEW(ITG,6**ne);
      ielfa=NNEW(ITG,24**ne);
      ifaext=NNEW(ITG,6**ne);
      isolidsurf=NNEW(ITG,6**ne);
      vel=NNEW(double,mi[1]**ne);
      vfa=NNEW(double,mi[1]*6**ne);

      FORTRAN(precfd,(ne,ipkon,kon,lakon,ipnei,neifa,neiel,ipoface,
       nodface,ielfa,&nflnei,&nface,ifaext,&nfaext,
       isolidsurf,&nsolidsurf,set,nset,istartset,iendset,ialset,
       vel,vfa,vold,mi));

      free(ipoface);free(nodface);
      RENEW(neifa,ITG,nflnei);
      RENEW(neiel,ITG,nflnei);
      RENEW(ielfa,ITG,4*nface);
      RENEW(ifaext,ITG,nfaext);
      RENEW(isolidsurf,ITG,nsolidsurf);
      RENEW(vfa,double,mi[1]*nface);
  }
  if(*ithermal>1){qfx=NNEW(double,3*mi[0]**ne);}
  
  /* contact conditions */
  
  inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	  lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
          &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
          mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,ifacecount,
          nmpc,&mpcfree,&memmpc_,
	  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	  iperturb,ikboun,nboun,co,istep,&xnoels);
  
  if(ncont!=0){
      
      cg=NNEW(double,3*ncont);
      straight=NNEW(double,16*ncont);
	  
      /* 11 instead of 10: last position is reserved for the
	 local contact spring element number; needed as
	 pointer into springarea */
      
      if(*mortar==0){
	  RENEW(kon,ITG,*nkon+11**nslavs);
	  springarea=NNEW(double,2**nslavs);
	  if(*nener==1){
	      RENEW(ener,double,mi[0]*(*ne+*nslavs)*2);
	  }
	  RENEW(ipkon,ITG,*ne+*nslavs);
	  RENEW(lakon,char,8*(*ne+*nslavs));
	  
	  if(*norien>0){
	      RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
	      for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
	  }

	  RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;

	  if((maxprevcontel==0)&&(*nslavs!=0)){
	      RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	      for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
		  xstate[k]=0.;
	      }
	  }
	  maxprevcontel=*nslavs;

	  areaslav=NNEW(double,*ifacecount);
      }else if(*mortar==1){
	  islavact=NNEW(ITG,nslavnode[*ntie]);
	  if((*istep==1)||(nslavs_prev_step==0)) clearini=NNEW(double,3*9**nslavs);
      }

      xmastnor=NNEW(double,3*nmastnode[*ntie]);
  }
  
  if((icascade==2)||(ncont!=0)){
      memmpcref_=memmpc_;mpcfreeref=mpcfree;
      nodempcref=NNEW(ITG,3*memmpc_);
      for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
      coefmpcref=NNEW(double,memmpc_);
      for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
  }
  
  if((*ithermal==1)||(*ithermal>=3)){
      t1ini=NNEW(double,*nk);
      t1act=NNEW(double,*nk);
      for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
  /* allocating a field for the instantaneous amplitude */
  
  ampli=NNEW(double,*nam);
  
  fini=NNEW(double,neq[1]);
  
  /* allocating fields for nonlinear dynamics */
  
  if(*nmethod==4){
      mass[0]=1;
      mass[1]=1;
      aux2=NNEW(double,neq[1]);
      fextini=NNEW(double,neq[1]);
      veini=NNEW(double,mt**nk);
      accini=NNEW(double,mt**nk);
      adb=NNEW(double,neq[1]);
      aub=NNEW(double,nzs[1]);
  }

  if((*nstate_!=0)&&(*mortar==0)){
      xstateini=NNEW(double,*nstate_*mi[0]*(*ne+*nslavs));
      for(k=0;k<*nstate_*mi[0]*(*ne+*nslavs);++k){
	  xstateini[k]=xstate[k];
      }
  }
  if((*nstate_!=0)&&(*mortar==1)) xstateini=NNEW(double,1);
  eei=NNEW(double,6*mi[0]**ne);
  stiini=NNEW(double,6*mi[0]**ne);
  emeini=NNEW(double,6*mi[0]**ne);
  if(*nener==1)
      enerini=NNEW(double,mi[0]**ne);
  
  qa[0]=qaold[0];
  qa[1]=qaold[1];
  
  /* normalizing the time */
  
  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,&inext,&itp));
  dtheta=(*tinc)/(*tper);

  /* taking care of a small increment at the end of the step
     for face-to-face penalty contact */

//  if((*mortar==1)&&(dtheta>1.-1.e-6)) dtheta=0.999;
//  if((*mortar==1)&&(dtheta>1.-1.e-6)) dtheta=0.900;
  dthetaref=dtheta;
  if((dtheta<=1.e-6)&&(*iexpl<=1)){
      printf("\n *ERROR in nonlingeo\n");
      printf(" increment size smaller than one millionth of step size\n");
      printf(" increase increment size\n\n");
  }
  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;
  
  /* calculating an initial flux norm */
  
  if(*ithermal!=2){
      if(qau>1.e-10){qam[0]=qau;}
      else if(qa0>1.e-10){qam[0]=qa0;}
      else if(qa[0]>1.e-10){qam[0]=qa[0];}
      else {qam[0]=1.e-2;}
  }
  if(*ithermal>1){
      if(qau>1.e-10){qam[1]=qau;}
      else if(qa0>1.e-10){qam[1]=qa0;}
      else if(qa[1]>1.e-10){qam[1]=qa[1];}
      else {qam[1]=1.e-2;}
  }
  
  /* storing the element and topology information before introducing 
     contact elements */
  
  ne0=*ne;nkon0=*nkon;neold=*ne;
  
  /*********************************************************************/
  
  /* calculating the initial acceleration at the start of the step
     for dynamic calculations */
  
  /*********************************************************************/
  
  if((*nmethod==4)&&(*ithermal!=2)){
      bet=(1.-*alpha)*(1.-*alpha)/4.;
      gam=0.5-*alpha;
      
      /* calculating the stiffness and mass matrix */
      
      reltime=0.;
      time=0.;
      dtime=0.;
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
      
      time=0.;
      dtime=1.;
    
      /*  updating the nonlinear mpc's (also affects the boundary
	  conditions through the nonhomogeneous part of the mpc's)
	  if contact arises the number of MPC's can also change */
      
      cam[0]=0.;cam[1]=0.;
      
      if(icascade==2){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;
	  RENEW(nodempc,ITG,3*memmpcref_);
	  for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	  RENEW(coefmpc,double,memmpcref_);
	  for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
      }

      newstep=0;
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounold,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi));
      if(icascade==2){
	  memmpcref_=memmpc_;mpcfreeref=mpcfree;
	  RENEW(nodempcref,ITG,3*memmpc_);
	  for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	  RENEW(coefmpcref,double,memmpc_);
	  for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
      }
      
      if(icascade>0) remastruct(ipompc,&coefmpc,&nodempc,nmpc,
	      &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	      labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	      kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	      neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	      &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar);
      
      /* invert nactdof */

      free(nactdofinv);nactdofinv=NNEW(ITG,1);
/*      free(nactdofinv);nactdofinv=NNEW(ITG,mt**nk);nodorig=NNEW(ITG,*nk);
      FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			   ipkon,lakon,kon,ne));
			   free(nodorig);*/
      
      iout=-1;
      ielas=1;
      
      fn=NNEW(double,mt**nk);
      stx=NNEW(double,6*mi[0]**ne);
      
      inum=NNEW(ITG,*nk);
      results(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounold,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,&bet,
	      &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	      ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,xstaten,
	      eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	      &ne0,xforc,nforc,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,islavsurf);
      
      free(fn);free(stx);free(inum);
      
      iout=0;
      ielas=0;
      
      reltime=0.;
      time=0.;
      dtime=0.;
      
      if(*iexpl<=1){intscheme=1;}
      
      /* in mafillsm the stiffness and mass matrix are computed;
	 The primary aim is to calculate the mass matrix (not 
	 lumped for an implicit dynamic calculation, lumped for an
	 explicit dynamic calculation). However:
	 - for an implicit calculation the mass matrix is "doped" with
	 a small amount of stiffness matrix, therefore the calculation
	 of the stiffness matrix is needed.
	 - for an explicit calculation the stiffness matrix is not 
	 needed at all. Since the calculation of the mass matrix alone
	 is not possible in mafillsm, the determination of the stiffness
	 matrix is taken as unavoidable "ballast". */
      
      ad=NNEW(double,neq[1]);
      au=NNEW(double,nzs[1]);
      
      FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounold,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	      nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	      nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
	      nmethod,ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	      ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mi,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
	      physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
	      &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
              xstateini,xstate,thicke,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
	      pmastsurf,mortar,clearini));
      
      if(*nmethod==0){
	  
	  /* error occurred in mafill: storing the geometry in frd format */
	  
	  ++*kode;
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(ITG,40**ne);ipneigh=NNEW(ITG,*nk);
	  }
	  
	  ptime=*ttime+time;
	  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	      kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	      thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni);
	  
	  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
	  
	  FORTRAN(stop,());
	  
      }
      
      /* mass x acceleration = f(external)-f(internal) 
	 only for the mechanical loading*/
      
      for(k=0;k<neq[0];++k){
	  b[k]=fext[k]-f[k];
      }
      
      if(*iexpl<=1){
	  
	  /* a small amount of stiffness is added to the mass matrix
	     otherwise the system leads to huge accelerations in 
	     case of discontinuous load changes at the start of the step */
	  
	  dtime=*tinc/10.;
	  scal1=bet*dtime*dtime*(1.+*alpha);
	  for(k=0;k<neq[0];++k){
	      ad[k]=adb[k]+scal1*ad[k];
	  }
	  for(k=0;k<nzs[0];++k){
	      au[k]=aub[k]+scal1*au[k];
	  }
	  if(*isolver==0){
#ifdef SPOOLES
	      spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		      &symmetryflag,&inputformat,&nzs[2]);
#else
	      printf("*ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if((*isolver==2)||(*isolver==3)){
	      preiter(ad,&au,b,&icol,&irow,&neq[0],&nzs[0],isolver,iperturb);
	  }
	  else if(*isolver==4){
#ifdef SGI
	      token=1;
	      sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
#else
	      printf("*ERROR in nonlingeo: the SGI library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[0],&nzs[0]);
#else
	      printf("*ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	      printf("*ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
      }
      
      else{
	  for(k=0;k<neq[0];++k){
	      b[k]=(fext[k]-f[k])/adb[k];
	  }
      }
      
      /* for thermal loading the acceleration is set to zero */
      
      for(k=neq[0];k<neq[1];++k){
	  b[k]=0.;
      }
      
      /* storing the acceleration accold */
      
      for(k=0;k<mt**nk;++k){
	  if(nactdof[k]!=0){accold[k]=b[nactdof[k]-1];}
      }
      
      free(ad);free(au);
      
      /* the mass matrix is kept for subsequent calculations, therefore,
	 no new mass calculation is necessary for the remaining iterations
	 in the present step */
      
      mass[0]=0;intscheme=0;
      
  }
  
  if(*iexpl>1) icmd=3;
  
  /**************************************************************/
  /* starting the loop over the increments                      */
  /**************************************************************/
  
  newstep=1;
  
  while(1.-theta>1.e-6){
      
      if(icutb==0){
	  
	  /* previous increment converged: update the initial values */
	  
	  iinc++;
	  jprint++;
	  
	  /* vold is copied into vini */
	  
	  memcpy(&vini[0],&vold[0],sizeof(double)*mt**nk);
	  
	  for(k=0;k<*nboun;++k){xbounini[k]=xbounact[k];}
	  if((*ithermal==1)||(*ithermal>=3)){
	      for(k=0;k<*nk;++k){t1ini[k]=t1act[k];}
	  }
	  for(k=0;k<neq[1];++k){
	      fini[k]=f[k];
	  }
	  if(*nmethod==4){
	      for(k=0;k<mt**nk;++k){
		  veini[k]=veold[k];
		  accini[k]=accold[k];
	      }
	      for(k=0;k<neq[1];++k){
		  fextini[k]=fext[k];
	      }
	  }
	  if(*ithermal!=2){
	      for(k=0;k<6*mi[0]*ne0;++k){
		  stiini[k]=sti[k];
		  emeini[k]=eme[k];
	      }
	  }
	  if(*nener==1)
	      for(k=0;k<mi[0]*ne0;++k){enerini[k]=ener[k];}

	  if(*mortar==0){
	      if(*nstate_!=0){
		  for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
		      xstateini[k]=xstate[k];
		  }
	      }
	  }	
      }
      
      /* check for max. # of increments */
      
      if(iinc>*jmax){
	  printf(" *ERROR: max. # of increments reached\n\n");
	  FORTRAN(stop,());
      }
      printf(" increment %" ITGFORMAT " attempt %" ITGFORMAT " \n",iinc,icutb+1);
      printf(" increment size= %e\n",dtheta**tper);
      printf(" sum of previous increments=%e\n",theta**tper);
      printf(" actual step time=%e\n",(theta+dtheta)**tper);
      printf(" actual total time=%e\n\n",*ttime+(theta+dtheta)**tper);
      
      printf(" iteration 1\n\n");
      
      qamold[0]=qam[0];
      qamold[1]=qam[1];
      
      /* determining the actual loads at the end of the new increment*/
      
      reltime=theta+dtheta;
      time=reltime**tper;
      dtime=dtheta**tper;
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
              nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
      
      for(i=0;i<3;i++){cam[i]=0.;}for(i=3;i<5;i++){cam[i]=0.5;}
      if(*ithermal>1){radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
       ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,vold,
       shcon,nshcon,ipkon,kon,co,
       kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
       nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
       cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
       ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
       nactdog,nacteq,nodeboun,ndirboun,&network,
       rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,
       ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,ialset,nset,
       ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,iamload,
       jqrad,irowrad,&nzsrad,icolrad,ne);
      }
      
      if(icfd==1){
	  compfluid(&co,nk,&ipkon,&kon,&lakon,ne,&sideface,
            ifreestream,&nfreestream,isolidsurf,neighsolidsurf,&nsolidsurf,
            &iponoel,&inoel,nshcon,shcon,nrhcon,rhcon,&vold,ntmat_,nodeboun,
            ndirboun,nboun,&ipompc,&nodempc,nmpc,&ikmpc,&ilmpc,ithermal,
            ikboun,ilboun,&iturbulent,isolver,iexpl,vcontu,ttime,
            &time,&dtime,nodeforc,ndirforc,xforc,nforc,nelemload,sideload,
            xload,nload,xbody,ipobody,nbody,&ielmat,matname,mi,ncmat_,
            physcon,istep,&iinc,ibody,xloadold,xboun,&coefmpc,
            nmethod,xforcold,xforcact,iamforc,iamload,xbodyold,xbodyact,
            t1old,t1,t1act,iamt1,amta,namta,nam,ampli,xbounold,xbounact,
	    iamboun,itg,&ntg,amname,t0,&nelemface,&nface,cocon,ncocon,xloadact,
	    tper,jmax,jout,set,nset,istartset,iendset,ialset,prset,prlab,
	    nprint,trab,inotr,ntrans,filab,&labmpc,sti,norien,orab,jobnamef,
	    tieset,ntie,mcs,ics,cs,nkon,&mpcfree,&memmpc_,&fmpc,&nef,&inomat,
	    qfx,neifa,neiel,ielfa,ifaext,vfa,vel,ipnei,&nflnei,&nfaext,
            typeboun);

	  /* determining the global values to be used as boundary conditions
	     for a submodel */
	  
	  free(integerglob);free(doubleglob);
	  getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,
                   xboun,nload,sideload,iamload,&iglob);
      }
      
      if((icascade==2)||(ncont!=0)){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;
	  RENEW(nodempc,ITG,3*memmpcref_);
	  for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	  RENEW(coefmpc,double,memmpcref_);
	  for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
      }

      /* generating contact elements */
      
      if((ncont!=0)&&(*mortar<=1)){
	  *ne=ne0;*nkon=nkon0;

	  /* at start of new increment: 
             - copy state variables (node-to-face)
	     - determine slave integration points (face-to-face)
	     - interpolate state variables (face-to-face) */

	  if(icutb==0){
	      if(*mortar==1){

		  if(*nstate_!=0){
		      if(maxprevcontel!=0){
			  islavsurfold=NNEW(ITG,2**ifacecount+2);
			  pslavsurfold=NNEW(double,3**nintpoint);
			  memcpy(&islavsurfold[0],&islavsurf[0],
                                 sizeof(ITG)*(2**ifacecount+2));
			  memcpy(&pslavsurfold[0],&pslavsurf[0],
                                 sizeof(double)*(3**nintpoint));
		      }
		  }

		  *nintpoint=0;

		  precontact(&ncont,ntie,tieset,nset,set,istartset,
                     iendset,ialset,itietri,lakon,ipkon,kon,koncont,ne,
                     cg,straight,co,vold,istep,&iinc,&iit,itiefac,
                     islavsurf,islavnode,imastnode,nslavnode,nmastnode,
                     imastop,mi,ipe,ime,tietol,&iflagact,
		     nintpoint,&pslavsurf,xmastnor,cs,mcs,ics,clearini,
                     nslavs);
		  
		  /* changing the dimension of element-related fields */
		  
		  RENEW(kon,ITG,*nkon+22**nintpoint);
		  RENEW(springarea,double,2**nintpoint);
		  RENEW(pmastsurf,double,6**nintpoint);
		  
		  if(*nener==1){
		      RENEW(ener,double,mi[0]*(*ne+*nintpoint)*2);
		  }
		  RENEW(ipkon,ITG,*ne+*nintpoint);
		  RENEW(lakon,char,8*(*ne+*nintpoint));
		  
		  if(*norien>0){
		      RENEW(ielorien,ITG,mi[2]*(*ne+*nintpoint));
		      for(k=mi[2]**ne;k<mi[2]*(*ne+*nintpoint);k++) ielorien[k]=0;
		  }
		  RENEW(ielmat,ITG,mi[2]*(*ne+*nintpoint));
		  for(k=mi[2]**ne;k<mi[2]*(*ne+*nintpoint);k++) ielmat[k]=1;

	      }
	  }

	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
		  co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
		  &ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
		  &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
		  iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,
                  islavsurf,
		  itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
		  imastnode,nmastnode,xmastnor,filab,mcs,ics,&nasym,
		  xnoels,mortar,pslavsurf,pmastsurf,clearini,&theta);
	  
	  printf("number of contact spring elements=%" ITGFORMAT "\n\n",*ne-ne0);
		  
	  /* interpolating the state variables */

	  if(icutb==0){
	      if(*mortar==1){
		  if(*nstate_!=0){
		      if(maxprevcontel!=0){
			  RENEW(xstateini,double,
                                *nstate_*mi[0]*(ne0+maxprevcontel));
			  for(k=*nstate_*mi[0]*ne0;
                                k<*nstate_*mi[0]*(ne0+maxprevcontel);++k){
			      xstateini[k]=xstate[k];
			  }
		      }
		      
		      RENEW(xstate,double,*nstate_*mi[0]*(ne0+*nintpoint));
		      for(k=*nstate_*mi[0]*ne0;k<*nstate_*mi[0]*(ne0+*nintpoint);k++){
			  xstate[k]=0.;
		      }
		      
		      if((*nintpoint>0)&&(maxprevcontel>0)){
			  iex=2;
			  
			  /* interpolation of xstate */
			  
			  FORTRAN(interpolatestate,(ne,ipkon,kon,lakon,
                               &ne0,mi,xstate,pslavsurf,nstate_,
                               xstateini,islavsurf,islavsurfold,
                               pslavsurfold));
			  
		      }

		      if(maxprevcontel!=0){
			  free(islavsurfold);free(pslavsurfold);
		      }

		      maxprevcontel=*nintpoint;

		      RENEW(xstateini,double,*nstate_*mi[0]*(ne0+*nintpoint));
		      for(k=0;k<*nstate_*mi[0]*(ne0+*nintpoint);++k){
			  xstateini[k]=xstate[k];
		      }
		  }
	      }
	  }
	  
	  if(icascade<1)icascade=1;
      }
      
      /*  updating the nonlinear mpc's (also affects the boundary
	  conditions through the nonhomogeneous part of the mpc's) */
      
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi));
      
      if((icascade==2)||(ncont!=0)){
	  memmpcref_=memmpc_;mpcfreeref=mpcfree;
	  RENEW(nodempcref,ITG,3*memmpc_);
	  for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	  RENEW(coefmpcref,double,memmpc_);
	  for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
      }
      
      if(icascade>0) remastruct(ipompc,&coefmpc,&nodempc,nmpc,
	  &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	  labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	  kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	  neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	  &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar);
      
      /* invert nactdof */
      
      free(nactdofinv);nactdofinv=NNEW(ITG,mt**nk);nodorig=NNEW(ITG,*nk);
      FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			     ipkon,lakon,kon,ne));
      free(nodorig);
      
      /* check whether the forced displacements changed; if so, and
	 if the procedure is static, the first iteration has to be
	 purely linear elastic, in order to get an equilibrium
	 displacement field; otherwise huge (maybe nonelastic)
	 stresses may occur, jeopardizing convergence */
      
      ilin=0;
      
      /* only for iinc=1 a linearized calculation is performed, since
	 for iinc>1 a reasonable displacement field is predicted by using the
	 initial velocity field at the end of the last increment */
      
      if((iinc==1)&&(*ithermal<2)){
	  dev=0.;
	  for(k=0;k<*nboun;++k){
	      err=fabs(xbounact[k]-xbounini[k]);
	      if(err>dev){dev=err;}
	  }
	  if(dev>1.e-5) ilin=1;
      }
      
      /* prediction of the kinematic vectors  */
      
      v=NNEW(double,mt**nk);
      
      prediction(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,accold,v,
		 &iinc,&idiscon,vold,nactdof,mi);
      
      fn=NNEW(double,mt**nk);
      stx=NNEW(double,6*mi[0]**ne);
      
      /* determining the internal forces at the start of the increment
	 
	 for a static calculation with increased forced displacements
	 the linear strains are calculated corresponding to
	 
	 the displacements at the end of the previous increment, extrapolated
	 if appropriate (for nondispersive media) +
	 the forced displacements at the end of the present increment +
	 the temperatures at the end of the present increment (this sum is
	 v) -
	 the displacements at the end of the previous increment (this is vold)
	 
	 these linear strains are converted in stresses by multiplication
	 with the tangent element stiffness matrix and converted into nodal
	 forces. 
	 
	 this boils down to the fact that the effect of forced displacements
	 should be handled in a purely linear way at the
	 start of a new increment, in order to speed up the convergence and
	 (for dissipative media) guarantee smooth loading within the increment.
	 
	 for all other cases the nodal force calculation is based on
	 the true stresses derived from the appropriate strain tensor taking
	 into account the extrapolated displacements at the end of the 
	 previous increment + the forced displacements and the temperatures
	 at the end of the present increment */
      
      iout=-1;
      iperturb_sav[0]=iperturb[0];
      iperturb_sav[1]=iperturb[1];
      
      /* first iteration in first increment: elastic tangent */
      
      if((*nmethod!=4)&&(ilin==1)){
	  
	  ielas=1;
	  
	  iperturb[0]=-1;
	  iperturb[1]=0;
	  
	  for(k=0;k<neq[1];++k){b[k]=f[k];}
	  inum=NNEW(ITG,*nk);
	  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t1ini,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		  &icmd, ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
                  islavsurf);
	  iperturb[0]=0;free(inum);
	  
	  /* check whether any displacements or temperatures are changed
	     in the new increment */
	  
	  for(k=0;k<neq[1];++k){f[k]=f[k]+b[k];}
	  
      }
      else{
	  
	  inum=NNEW(ITG,*nk);
	  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
                  islavsurf);
	  free(inum);
	  
	  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
	  if(*ithermal!=2){
	      for(k=0;k<6*mi[0]*ne0;++k){
		  sti[k]=stx[k];
	      }
	  }
	  
      }
      
      ielas=0;
      iout=0;
      
      free(fn);free(stx);free(v);
      
      /***************************************************************/
      /* iteration counter and start of the loop over the iterations */
      /***************************************************************/

    iit=1;
    icntrl=0;
    ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
    if(uncoupled){
	*ithermal=2;
	iruc=NNEW(ITG,nzs[1]-nzs[0]);
	for(k=0;k<nzs[1]-nzs[0];k++) iruc[k]=irow[k+nzs[0]]-neq[0];
    }

    while(icntrl==0){

    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's) */

      if((iit!=1)||((uncoupled)&&(*ithermal==1))){

	  printf(" iteration %" ITGFORMAT "\n\n",iit);

          FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
              xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
              nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));

	  for(i=0;i<3;i++){cam[i]=0.;}for(i=3;i<5;i++){cam[i]=0.5;}
	  if(*ithermal>1){radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
	     ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,
             ntmat_,vold,shcon,nshcon,ipkon,kon,co,
	     kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
	     nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
             cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
	     ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
	     nactdog,nacteq,nodeboun,ndirboun,&network,
             rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,
	     ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,ialset,
	     nset,ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,
	     iamload,jqrad,irowrad,&nzsrad,icolrad,ne);
	  }

	  if((icascade==2)||
	     ((ncont!=0)&&(ismallsliding==0))){
	      memmpc_=memmpcref_;mpcfree=mpcfreeref;
	      RENEW(nodempc,ITG,3*memmpcref_);
	      for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	      RENEW(coefmpc,double,memmpcref_);
	      for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
	  }

	  if((ncont!=0)&&(*mortar<=1)&&(ismallsliding==0)&&((iit<=8)||(*mortar==1))){
	      neold=*ne;
	      *ne=ne0;*nkon=nkon0;
	      contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
                      straight,nkon,co,vold,ielmat,cs,elcon,istep,
                      &iinc,&iit,ncmat_,ntmat_,&ne0,
                      vini,nmethod,nmpc,&mpcfree,&memmpc_,&ipompc,&labmpc,
                      &ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,iperturb,
                      ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
                      itiefac,areaslav,iponoels,inoels,springarea,tietol,
                      &reltime,imastnode,nmastnode,xmastnor,
                      filab,mcs,ics,&nasym,xnoels,mortar,pslavsurf,pmastsurf,
                      clearini,&theta);

	      if(*mortar==0){
	         if(*ne!=neold){iflagact=1;}
	      }else{
	         if(((*ne-ne0)<(neold-ne0)*0.999)||
                    ((*ne-ne0)>(neold-ne0)*1.001)){iflagact=1;}
              }

	      printf("number of contact spring elements=%" ITGFORMAT "\n\n",*ne-ne0);

	      if(icascade<1)icascade=1;
	  }
	  
	  if(*ithermal==3){
	      for(k=0;k<*nk;++k){t1act[k]=vold[mt*k];}
	  }

	  FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
	        &maxlenmpc,ikmpc,ilmpc,&icascade,
	        kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,&iit,
		&idiscon,&ncont,trab,ntrans,ithermal,mi));

	  if((icascade==2)||
	      ((ncont!=0)&&(ismallsliding==0))){
	      memmpcref_=memmpc_;mpcfreeref=mpcfree;
	      RENEW(nodempcref,ITG,3*memmpc_);
	      for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	      RENEW(coefmpcref,double,memmpc_);
	      for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
	  }

	  if(icascade>0){
	      remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		&mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		&adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar);

	      /* invert nactdof */
	      
	      free(nactdofinv);nactdofinv=NNEW(ITG,mt**nk);nodorig=NNEW(ITG,*nk);
	      FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
				     ipkon,lakon,kon,ne));
	      free(nodorig);
	      
	      v=NNEW(double,mt**nk);
	      stx=NNEW(double,6*mi[0]**ne);
	      fn=NNEW(double,mt**nk);
      
	      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	      iout=-1;
	      
	      inum=NNEW(ITG,*nk);
	      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	        ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
                sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
                islavsurf);

	      /*for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}*/
	      
	      free(v);free(stx);free(fn);free(inum);
	      iout=0;
	      
	  }else{

	      /*for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}*/
	  }
      }
      
      if(*iexpl<=1){

	/* calculating the local stiffness matrix and external loading */

	ad=NNEW(double,neq[1]);
	au=NNEW(double,nzs[1]);

	FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		  nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		  nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		  nmethod,ikmpc,ilmpc,ikboun,ilboun,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		  xstiff,npmat_,&dtime,matname,mi,
                  ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
                  physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		  &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
                  xstateini,xstate,thicke,integerglob,doubleglob,
		  tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		  pmastsurf,mortar,clearini));

	if(nasym==1){
	    RENEW(au,double,2*nzs[1]);
	    if(*nmethod==4) RENEW(aub,double,2*nzs[1]);
	    symmetryflag=2;
	    inputformat=1;

	    FORTRAN(mafillsmas,(co,nk,kon,ipkon,lakon,ne,nodeboun,
                  ndirboun,xbounact,nboun,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		  nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		  nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		  nmethod,ikmpc,ilmpc,ikboun,ilboun,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		  xstiff,npmat_,&dtime,matname,mi,
                  ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
                  physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
                  &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
                  xstateini,xstate,thicke,
                  integerglob,doubleglob,tieset,istartset,iendset,
		  ialset,ntie,&nasym,pslavsurf,pmastsurf,mortar,clearini));
	}
	    

	iperturb[0]=iperturb_sav[0];
	iperturb[1]=iperturb_sav[1];

      }else{

	/* calculating the external loading 

	   This is only done once per increment. In reality, the
           external loading is a function of vold (specifically,
           the body forces and surface loading). This effect is
           neglected, since the increment size in dynamic explicit
           calculations is usually small */

	  FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		  nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		  nbody,cgr,fext,nactdof,&neq[1],
		  nmethod,ikmpc,ilmpc,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1act,ithermal,iprestr,vold,iperturb,
		  iexpl,plicon,nplicon,plkcon,nplkcon,
		  npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		  xbodyold,&reltime,veold,matname,mi,ikactmech,
                  &nactmech));
      }
      
/*      for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
      for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
      for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
      for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}*/

      /* calculating the residual */

      calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	 vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
	 islavnode,nslavnode,mortar,ntie,f_cm,f_cs,mi,
	 nzs,&nasym);
	  
      newstep=0;
      
      if(*nmethod==0){
	  
	  /* error occurred in mafill: storing the geometry in frd format */
	  
	  *nmethod=0;
	  ++*kode;
	  inum=NNEW(ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(ITG,40**ne);ipneigh=NNEW(ITG,*nk);
	  }
	  
	  ptime=*ttime+time;
	  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
		  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni);

	  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
	  free(inum);FORTRAN(stop,());
	  
      }
      
      /* implicit step (static or dynamic */
      
      if(*iexpl<=1){
	  if(*nmethod==4){
	      
	      /* mechanical part */
	      
	      if(*ithermal!=2){
		  scal1=bet*dtime*dtime*(1.+*alpha);
		  for(k=0;k<neq[0];++k){
		      ad[k]=adb[k]+scal1*ad[k];
		  }
		  for(k=0;k<nzs[0];++k){
		      au[k]=aub[k]+scal1*au[k];
		  }
		  
		  /* upper triangle of asymmetric matrix */
		  
		  if(nasym>0){
		      for(k=nzs[2];k<nzs[2]+nzs[0];++k){
			  au[k]=aub[k]+scal1*au[k];
		      }
		  }
	      }
	      
	      /* thermal part */
	      
	      if(*ithermal>1){
		  for(k=neq[0];k<neq[1];++k){
		      ad[k]=adb[k]/dtime+ad[k];
		  }
		  for(k=nzs[0];k<nzs[1];++k){
		      au[k]=aub[k]/dtime+au[k];
		  }
		  
		  /* upper triangle of asymmetric matrix */
		  
		  if(nasym>0){
		      for(k=nzs[2]+nzs[0];k<nzs[2]+nzs[1];++k){
			  au[k]=aub[k]/dtime+au[k];
		      }
		  }
	      }
	  }
	  
	  if(*isolver==0){
#ifdef SPOOLES
	      if(*ithermal<2){
		  spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			  &symmetryflag,&inputformat,&nzs[2]);
	      }else if((*ithermal==2)&&(uncoupled)){
		  n1=neq[1]-neq[0];
		  n2=nzs[1]-nzs[0];
		  spooles(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			  &sigma,&b[neq[0]],&icol[neq[0]],iruc,
			  &n1,&n2,&symmetryflag,&inputformat,&nzs[2]);
	      }else{
		  spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			  &symmetryflag,&inputformat,&nzs[2]);
	      }
#else
	      printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if((*isolver==2)||(*isolver==3)){
	      if(nasym>0){
		  printf(" *ERROR in nonlingeo: the iterative solver cannot be used for asymmetric matrices\n\n");
		  FORTRAN(stop,());
	      }
	      preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
	  }
	  else if(*isolver==4){
#ifdef SGI
	      if(nasym>0){
		  printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
		  FORTRAN(stop,());
	      }
	      token=1;
	      if(*ithermal<2){
		  sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
	      }else if((*ithermal==2)&&(uncoupled)){
		  n1=neq[1]-neq[0];
		  n2=nzs[1]-nzs[0];
		  sgi_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			   &sigma,&b[neq[0]],&icol[neq[0]],iruc,
			   &n1,&n2,token);
	      }else{
		  sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
	      }
#else
	      printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if(*isolver==5){
	      if(nasym>0){
		  printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
		  FORTRAN(stop,());
	      }
#ifdef TAUCS
	      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
	      printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      if(*ithermal<2){
		  pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			       &symmetryflag,&inputformat,jq,&nzs[2]);
	      }else if((*ithermal==2)&&(uncoupled)){
		  n1=neq[1]-neq[0];
		  n2=nzs[1]-nzs[0];
		  pardiso_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			       &sigma,&b[neq[0]],&icol[neq[0]],iruc,
			       &n1,&n2,&symmetryflag,&inputformat,jq,&nzs[2]);
	      }else{
		  pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			       &symmetryflag,&inputformat,jq,&nzs[2]);
	      }
#else
	      printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  
	  if(*mortar<=1){free(ad);free(au);} 
      }
      
      /* explicit dynamic step */
      
      else{
	  if(*ithermal!=2){
	      for(k=0;k<neq[0];++k){
		  b[k]=b[k]/adb[k];
	      }
	  }
	  if(*ithermal>1){
	      for(k=neq[0];k<neq[1];++k){
		  b[k]=b[k]*dtime/adb[k];
	      }
	  }
      }
      
      /* calculating the displacements, stresses and forces */
      
      v=NNEW(double,mt**nk);
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      
      stx=NNEW(double,6*mi[0]**ne);
      fn=NNEW(double,mt**nk);
      
      inum=NNEW(ITG,*nk);
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf);
      free(inum);

      if(*ithermal!=2){
	  if(cam[0]>uam[0]){uam[0]=cam[0];}      
	  if(qau<1.e-10){
	      if(qa[0]>ea*qam[0]){qam[0]=(qamold[0]*jnz+qa[0])/(jnz+1);}
	      else {qam[0]=qamold[0];}
	  }
      }
      if(*ithermal>1){
	  if(cam[1]>uam[1]){uam[1]=cam[1];}      
	  if(qau<1.e-10){
	      if(qa[1]>ea*qam[1]){qam[1]=(qamold[1]*jnz+qa[1])/(jnz+1);}
	      else {qam[1]=qamold[1];}
	  }
      }

      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
      if(*ithermal!=2){
	  for(k=0;k<6*mi[0]*ne0;++k){
	      sti[k]=stx[k];
	  }
      }
      
      free(v);free(stx);free(fn);
      
      /* calculating the residual */
      
      calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	 vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
	 islavnode,nslavnode,mortar,ntie,f_cm,f_cs,mi,
	 nzs,&nasym);
      
      /* calculating the maximum residual */

      for(k=0;k<2;++k){
	  ram2[k]=ram1[k];
	  ram1[k]=ram[k];
	  ram[k]=0.;
      }
      if(*ithermal!=2){
	  for(k=0;k<neq[0];++k){
	      err=fabs(b[k]);
	      if(err>ram[0]){ram[0]=err;ram[2]=k+0.5;}
	  }
      }
      if(*ithermal>1){
	  for(k=neq[0];k<neq[1];++k){
	      err=fabs(b[k]);
	      if(err>ram[1]){ram[1]=err;ram[3]=k+0.5;}
	  }
      }

  /*   Divergence criteria for face-to-face penalty is different */

      if(*mortar==1){
	  for(k=4;k<8;++k){
	      ram2[k]=ram1[k];
	      ram1[k]=ram[k];
	  } 
	  ram[4]=ram[0]+ram1[0];
	  if((iflagact==0)&&(iit>1)){
	      ram[5]=1.5;
	  }else{ram[5]=0.5;}
	  ram[6]=abs((neold-ne0)-(*ne-ne0))+0.5;
	  if(iit>3){
	      if((ram[6]>=ram1[6])&&(ram[6]>=ram2[6])){
		  ram[7]=1.5;
	      }else{ram[7]=0.5;}
	  }
	  
/*	  for(k=2;k<4;++k){
	      ram2[k]=ram1[k];
	      ram1[k]=ram[k+2];
	  } 
	  ram[4]=ram[0]+ram1[0];
	  
	  ram[5]=abs((neold-ne0)-(*ne-ne0));
	  
	  if(iit==2){ram[6]=ram[4];}
	  if((iit>=3)&&(ram[6]>ram[4])){ram[6]=ram[4];}*/
      }
      
      /* next line is inserted to cope with stress-less
	 temperature calculations */
      
      if(*ithermal!=2){
	  if(ram[0]<1.e-6) ram[0]=0.;      
	  printf(" average force= %f\n",qa[0]);
	  printf(" time avg. forc= %f\n",qam[0]);
	  if((ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1==0){
	      printf(" largest residual force= %f\n",
		 ram[0]);
	  }else{
	      inode=(ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1;
	      idir=nactdofinv[(ITG)ram[2]]-mt*(inode-1);
	      printf(" largest residual force= %f in node %" ITGFORMAT " and dof %" ITGFORMAT "\n",
		     ram[0],inode,idir);
	  }
	  printf(" largest increment of disp= %e\n",uam[0]);
	  if((ITG)cam[3]==0){
	      printf(" largest correction to disp= %e\n\n",
                 cam[0]);
	  }else{
	      inode=(ITG)((double)nactdofinv[(ITG)cam[3]]/mt)+1;
	      idir=nactdofinv[(ITG)cam[3]]-mt*(inode-1);
	      printf(" largest correction to disp= %e in node %" ITGFORMAT " and dof %" ITGFORMAT "\n\n",cam[0],inode,idir);
	  }
      }
      if(*ithermal>1){
	  if(ram[1]<1.e-6) ram[1]=0.;      
	  printf(" average flux= %f\n",qa[1]);
	  printf(" time avg. flux= %f\n",qam[1]);
	  if((ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1==0){
	      printf(" largest residual flux= %f\n",
                 ram[1]);
	  }else{
	      inode=(ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1;
	      idir=nactdofinv[(ITG)ram[3]]-mt*(inode-1);
	      printf(" largest residual flux= %f in node %" ITGFORMAT " and dof %" ITGFORMAT "\n",ram[1],inode,idir);
	  }
	  printf(" largest increment of temp= %e\n",uam[1]);
	  if((ITG)cam[4]==0){
	      printf(" largest correction to temp= %e\n\n",
                 cam[1]);
	  }else{
	      inode=(ITG)((double)nactdofinv[(ITG)cam[4]]/mt)+1;
	      idir=nactdofinv[(ITG)cam[4]]-mt*(inode-1);
	      printf(" largest correction to temp= %e in node %" ITGFORMAT " and dof %" ITGFORMAT "\n\n",cam[1],inode,idir);
	  }
      }
      fflush(stdout);

      checkconvergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
	  kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
	  xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
          ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
	  ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
	  vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
	  &theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
	  nactdof,b,tmin,ctrl,amta,namta,itpamp,&inext,&dthetaref,
          &itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
          nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
	  set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,mortar);

    }

    /*********************************************************/
    /*   end of the iteration loop                          */
    /*********************************************************/

    /* icutb=0 means that the iterations in the increment converged,
       icutb!=0 indicates that the increment has to be reiterated with
                another increment size (dtheta) */
   
    if(uncoupled){
	free(iruc);
    }

    if(((qa[0]>ea*qam[0])||(qa[1]>ea*qam[1]))&&(icutb==0)){jnz++;}
    iit=0;

    if(icutb!=0){
      memcpy(&vold[0],&vini[0],sizeof(double)*mt**nk);

      for(k=0;k<*nboun;++k){xbounact[k]=xbounini[k];}
      if((*ithermal==1)||(*ithermal>=3)){
	for(k=0;k<*nk;++k){t1act[k]=t1ini[k];}
      }
	for(k=0;k<neq[1];++k){
	  f[k]=fini[k];
	}
      if(*nmethod==4){
	for(k=0;k<mt**nk;++k){
	  veold[k]=veini[k];
	  accold[k]=accini[k];
	}
	for(k=0;k<neq[1];++k){
	  fext[k]=fextini[k];
	}
      }
      if(*ithermal!=2){
	  for(k=0;k<6*mi[0]*ne0;++k){
	      sti[k]=stiini[k];
	      eme[k]=emeini[k];
	  }
      }
      if(*nener==1)
	  for(k=0;k<mi[0]*ne0;++k){ener[k]=enerini[k];}

      for(k=0;k<*nstate_*mi[0]*(ne0+maxprevcontel);++k){
	  xstate[k]=xstateini[k];
      }	  

      qam[0]=qamold[0];
      qam[1]=qamold[1];
    }
    
    /* face-to-face penalty */

    if((*mortar==1)&&(icutb==0)){
	
	ntrimax=0;
	for(i=0;i<*ntie;i++){	    
	    if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
		ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
	}
	xo=NNEW(double,ntrimax);	    
	yo=NNEW(double,ntrimax);	    
	zo=NNEW(double,ntrimax);	    
	x=NNEW(double,ntrimax);	    
	y=NNEW(double,ntrimax);	    
	z=NNEW(double,ntrimax);	   
	nx=NNEW(ITG,ntrimax);	   
	ny=NNEW(ITG,ntrimax);	    
	nz=NNEW(ITG,ntrimax);
      
	/*  Determination of active nodes (islavact) */
      
	FORTRAN(islavactive,(tieset,ntie,itietri,cg,straight,
		   co,vold,xo,yo,zo,x,y,z,nx,ny,nz,mi,
		   imastop,nslavnode,islavnode,islavact));

	free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);	    
	free(ny);free(nz);
    }

    /* output */

    if((jout[0]==jprint)&&(icutb==0)){

      jprint=0;

      /* calculating the displacements and the stresses and storing */
      /* the results in frd format  */
	
      v=NNEW(double,mt**nk);
      fn=NNEW(double,mt**nk);
      stn=NNEW(double,6**nk);
      if(*ithermal>1) qfn=NNEW(double,3**nk);
      inum=NNEW(ITG,*nk);
      stx=NNEW(double,6*mi[0]**ne);
      
      if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
      if(strcmp1(&filab[435],"PEEQ")==0) epn=NNEW(double,*nk);
      if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);
      if(strcmp1(&filab[609],"SDV ")==0) xstaten=NNEW(double,*nstate_**nk);
      if(strcmp1(&filab[2175],"CONT")==0) cdn=NNEW(double,6**nk);
      if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);

      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);

      iout=2;
      icmd=3;
      
#ifdef COMPANY
      FORTRAN(uinit,());
#endif
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	      ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
              xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
              ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,islavsurf);
      
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);

      iout=0;
      if(*iexpl<=1) icmd=0;
//      FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
//      for(k=0;k<ntg;k++)if(inum[itg[k]-1]>0){inum[itg[k]-1]*=-1;}
      
      ++*kode;
      if(*mcs!=0){
	ptime=*ttime+time;
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
               nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
               ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
               cdn,mortar);
      }
      else{
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(ITG,40**ne);ipneigh=NNEW(ITG,*nk);
	  }

	  ptime=*ttime+time;
	  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni);

	  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
      }
      
      free(v);free(fn);free(stn);free(inum);free(stx);
      if(*ithermal>1){free(qfn);}
      
      if(strcmp1(&filab[261],"E   ")==0) free(een);
      if(strcmp1(&filab[435],"PEEQ")==0) free(epn);
      if(strcmp1(&filab[522],"ENER")==0) free(enern);
      if(strcmp1(&filab[609],"SDV ")==0) free(xstaten);
      if(strcmp1(&filab[2175],"CONT")==0) free(cdn);
      if(strcmp1(&filab[2697],"ME  ")==0) free(emn);
   }
    
  }

  /*********************************************************/
  /*   end of the increment loop                          */
  /*********************************************************/

  if(jprint!=0){

  /* calculating the displacements and the stresses and storing  
     the results in frd format */
  
    v=NNEW(double,mt**nk);
    fn=NNEW(double,mt**nk);
    stn=NNEW(double,6**nk);
    if(*ithermal>1) qfn=NNEW(double,3**nk);
    inum=NNEW(ITG,*nk);
    stx=NNEW(double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
    if(strcmp1(&filab[435],"PEEQ")==0) epn=NNEW(double,*nk);
    if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);
    if(strcmp1(&filab[609],"SDV ")==0) xstaten=NNEW(double,*nstate_**nk);
    if(strcmp1(&filab[2175],"CONT")==0) cdn=NNEW(double,6**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);
    
    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
    iout=2;
    icmd=3;

#ifdef COMPANY
    FORTRAN(uinit,());
#endif
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xbounact,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
            ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
            &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
            mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,islavsurf);
    
    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);

    iout=0;
    if(*iexpl<=1) icmd=0;
//    FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
//    for(k=0;k<ntg;k++)if(inum[itg[k]-1]>0){inum[itg[k]-1]*=-1;}
    
    ++*kode;
    if(*mcs>0){
	ptime=*ttime+time;
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	     t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
             nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
             ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
             cdn,mortar);
    }
    else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	    neigh=NNEW(ITG,40**ne);ipneigh=NNEW(ITG,*nk);
	}

	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni);

	if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    }
    
    free(v);free(fn);free(stn);free(inum);free(stx);
    if(*ithermal>1){free(qfn);}
    
    if(strcmp1(&filab[261],"E   ")==0) free(een);
    if(strcmp1(&filab[435],"PEEQ")==0) free(epn);
    if(strcmp1(&filab[522],"ENER")==0) free(enern);
    if(strcmp1(&filab[609],"SDV ")==0) free(xstaten);
    if(strcmp1(&filab[2175],"CONT")==0) free(cdn);
    if(strcmp1(&filab[2697],"ME  ")==0) free(emn);

  }

  /* setting the velocity to zero at the end of a quasistatic or stationary
     step */

  if(abs(*nmethod)==1){
    for(k=0;k<mt**nk;++k){veold[k]=0.;}
  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){

      /* thermal boundary conditions are updated only if the
         step was thermal or thermomechanical */

      if(ndirboun[k]==0){
	  if(*ithermal<2) continue;

	  /* mechanical boundary conditions are updated only
             if the step was not thermal or the node is a
             network node */

      }else if((ndirboun[k]>0)&&(ndirboun[k]<4)){
	  node=nodeboun[k];
	  FORTRAN(nident,(itg,&node,&ntg,&id));
	  networknode=0;
	  if(id>0){
	      if(itg[id-1]==node) networknode=1;
	  }
	  if((*ithermal==2)&&(networknode==0)) continue;
      }
      xbounold[k]=xbounact[k];
  }
  for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    for(k=0;k<*nk;++k){vold[mt*k]=t1act[k];}
  }
  else if(*ithermal>1){
    for(k=0;k<*nk;++k){t1[k]=vold[mt*k];}
    if(*ithermal>=3){
	for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    }
  }

  qaold[0]=qa[0];
  qaold[1]=qa[1];

  free(f);
  free(b);
  free(xbounact);free(xforcact);free(xloadact);free(xbodyact);
  if(*nbody>0) free(ipobody);if(inewton==1){free(cgr);}
  free(fext);free(ampli);free(xbounini);free(xstiff);
  if((*ithermal==1)||(*ithermal>=3)){free(t1act);free(t1ini);}

  if(*ithermal>1){
      free(itg);free(ieg);free(kontri);free(nloadtr);
      free(nactdog);free(nacteq);free(ineighe);
      free(tarea);free(tenv);free(fenv);free(qfx);
      free(erad);free(ac);free(bc);free(ipiv);
      free(bcr);free(ipivr);free(adview);free(auview);free(adrad);
      free(aurad);free(irowrad);free(jqrad);free(icolrad);
      if((*mcs>0)&&(ntr>0)){free(inocs);}
  }

  if(icfd==1){
      free(neifa);free(neiel);free(ielfa);free(ifaext);
      free(vel);free(vfa);
  }

  free(fini);
  if(*nmethod==4){
    free(aux2);free(fextini);free(veini);free(accini);
    free(adb);free(aub);
  }
  free(eei);free(stiini);free(emeini);
  if(*nener==1)free(enerini);
  if(*nstate_!=0){free(xstateini);}

  free(aux);free(iaux);free(vini);

  if((icascade==2)||(ncont!=0)){
      memmpc_=memmpcref_;mpcfree=mpcfreeref;
      RENEW(nodempc,ITG,3*memmpcref_);
      for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
      RENEW(coefmpc,double,memmpcref_);
      for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
      free(nodempcref);free(coefmpcref);
  }

  if(ncont!=0){
      *ne=ne0;*nkon=nkon0;
      if(*nener==1){
	RENEW(ener,double,mi[0]**ne*2);
      }
      RENEW(ipkon,ITG,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,ITG,*nkon);
      if(*norien>0){
	  RENEW(ielorien,ITG,mi[2]**ne);
      }
      RENEW(ielmat,ITG,mi[2]**ne);

      free(cg);free(straight);
      free(imastop);free(itiefac);free(islavnode);
      free(nslavnode);free(iponoels);free(inoels);free(imastnode);
      free(nmastnode);free(itietri);free(koncont);free(xnoels);
      free(springarea);free(xmastnor);

      if(*mortar==0){
	  free(areaslav);
      }else if(*mortar==1){
	  free(pmastsurf);free(ipe);free(ime);
          free(islavact);
      }
  }

  /* reset icascade */

  if(icascade==1){icascade=0;}

  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;

  if(iglob==1){free(integerglob);free(doubleglob);}

  *icolp=icol;*irowp=irow;*cop=co;*voldp=vold;

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;

  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;

  *islavsurfp=islavsurf;*pslavsurfp=pslavsurf;*clearinip=clearini;

  (*tmin)*=(*tper);
  (*tmax)*=(*tper);

  free(nactdofinv);

  (*ttime)+=(*tper);
  
  return;
}
