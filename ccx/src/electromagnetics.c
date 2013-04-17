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


void electromagnetics(double **cop, ITG *nk, ITG **konp, ITG **ipkonp, 
             char **lakonp,ITG *ne, 
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
	     ITG *kode,char *filab, 
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
             char **setp, ITG *nset, ITG **istartsetp,
             ITG **iendsetp, ITG **ialsetp, ITG *nprint, char *prlab,
             char *prset, ITG *nener,ITG *ikforc,ITG *ilforc, double *trab, 
             ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
             ITG *ibody, double *xbody, ITG *nbody, double *xbodyold,
             ITG *ielprop, double *prop, ITG *ntie, char **tiesetp,
	     ITG *itpamp, ITG *iviewfile, char *jobnamec, double **tietolp,
	     ITG *nslavs, double *thicke, ITG *ics, ITG *nalset, ITG *nmpc_){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
      *labmpc=NULL,kind1[2]="E",kind2[2]="E",*set=NULL,*tieset=NULL,cflag[1]=" "; 
 
  ITG *inum=NULL,k,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
      icutb=0,istab=0,ifreebody,uncoupled=0,maxfaces,
      iperturb_sav[2],*icol=NULL,*irow=NULL,ielas=0,icmd=0,
      memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
      *itg=NULL,*ineighe=NULL,null=0,iactive[3],neqterms,
      *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,index,
      *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,
      ntrit,*inocs=NULL,inewton=0,*ipobody=NULL,*nacteq=NULL,
      *nactdog=NULL,nteq,network,nmastnode,imast,massact[2],
      *ipkon=NULL,*kon=NULL,*ielorien=NULL,nmethodact,
      *ielmat=NULL,inext,itp=0,symmetryflag=0,inputformat=0,
      iitterm=0,ngraph=1,ithermalact=2,*islavact=NULL,
      *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
      *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,mortar=0,
      mt=mi[1]+1,*nactdofinv=NULL, inode,idir,*islavsurf=NULL,
      iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
      *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
      coriolis=0,*ipneigh=NULL,*neigh=NULL,i,icfd=0,id,node,networknode,
      iflagact=0,*nodorig=NULL,*ipivr=NULL,*inomat=NULL,*nodface=NULL,
      *ipoface=NULL,*istartset=NULL,*iendset=NULL,*ialset=NULL;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cdn=NULL,
         *f=NULL,*fn=NULL,qa[3]={0.,0.,-1.},qam[2]={0.,0.},dtheta,theta,
         err,ram[4]={0.,0.,0.,0.},*springarea=NULL,*h0=NULL,
	 ram1[2]={0.,0.},ram2[2]={0.,0.},deltmx,*clearini=NULL,
	 uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,ptime,
	 *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
	 *xforcact=NULL,*xloadact=NULL,*fext=NULL,
         reltime,time,bet=0.,gam=0.,*aux1=NULL,*aux2=NULL,dtime,*fini=NULL,
         *fextini=NULL,*veini=NULL,*xstateini=NULL,
	 *ampli=NULL,*eei=NULL,*t1ini=NULL,
         *xbounini=NULL,*xstiff=NULL,*stx=NULL,
         *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
	 *enerini=NULL,*emn=NULL,*xmastnor=NULL,
	 *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
	 *adview=NULL,*auview=NULL,*qfx=NULL,
         *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
         *xbodyact=NULL,*cgr=NULL,dthetaref,*vr=NULL,*vi=NULL,
	 *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
	 *f_cm=NULL,*f_cs=NULL,*tietol=NULL,
	 *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
         *emeini=NULL,*doubleglob=NULL,*au=NULL,
	 *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
         *cdnr=NULL,*cdni=NULL;

#ifdef SGI
  ITG token;
#endif
 
  /* next line is needed to avoid that elements with negative ipkon
     are taken into account in extrapolate.f */

  strcpy1(&filab[2],"C",1);

  if(*nmethod==8){
      *nmethod=1;
  }else if(*nmethod==9){
      *nmethod=4;
  }else if(*nmethod==10){
      *nmethod=2;
  }
 
  for(k=0;k<3;k++){
      strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];
  
  icol=*icolp;irow=*irowp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

  set=*setp;istartset=*istartsetp;iendset=*iendsetp;ialset=*ialsetp;
  tieset=*tiesetp;tietol=*tietolp;
  
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
  }
  
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
  if(*ithermal>1){qfx=NNEW(double,3*mi[0]**ne);}
  
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
      adb=NNEW(double,neq[1]);
      aub=NNEW(double,nzs[1]);
  }
  
  qa[0]=qaold[0];
  qa[1]=qaold[1];
  
  /* normalizing the time */
  
  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,&inext,&itp));
  dtheta=(*tinc)/(*tper);
  dthetaref=dtheta;
  if(dtheta<=1.e-6){
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
  
  
  /*********************************************************************/
  
  /* calculating the initial quasi-static magnetic intensity due to 
     the coil current */
  
  /*********************************************************************/
  
  /* calculate the current density in the coils

     in this section nload, nforc, nbody and nam are set to zero; the
     electrical potential is supposed to be given (in the form of a
     temperature), the current is calculated (in the form of heat
     flux) by thermal analogy  */
  
  reltime=1.;
  time=0.;
  dtime=0.;
  ithermalact=2;

  nmethodact=1;
  massact[0]=0;
  massact[1]=0;

  if(*ithermal==0){
      qfx=NNEW(double,3*mi[0]**ne);
      t0=NNEW(double,*nk);
  }
  if(strcmp1(&filab[3567],"ECD ")==0){qfn=NNEW(double,3**nk);}
  
  /* the coil current is assumed to be applied at once, i.e. as 
     step loading; the calculation, however, is a quasi-static
     calculation */

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,&null,xloadold,xload,
	      xloadact,iamload,&null,ibody,xbody,&null,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,&null,ampli,&time,&reltime,ttime,&dtime,&ithermalact,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
  
  cam[0]=0.;cam[1]=0.;
  
  /* deactivating all elements except the shells */
  
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[8*i+6],"L")!=0){
	  ipkon[i]=-ipkon[i]-2;
      }
  }
  
  remastruct(ipompc,&coefmpc,&nodempc,nmpc,
	     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	     labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	     kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	     neq,nzs,&nmethodact,&f,&fext,&b,&aux2,&fini,&fextini,
	     &adb,&aub,&ithermalact,iperturb,mass,mi,iexpl,&mortar);
  
  /* invert nactdof */
  
  free(nactdofinv);nactdofinv=NNEW(ITG,mt**nk);nodorig=NNEW(ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  free(nodorig);
  
  iout=-1;
  
  fn=NNEW(double,mt**nk);
  inum=NNEW(ITG,*nk);
  v=NNEW(double,mt**nk);

  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,&ithermalact,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounact,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,&nmethodact,cam,&neq[1],veold,accold,&bet,
          &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
          ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,xstaten,
          eei,enerini,alcon,nalcon,set,nset,istartset,iendset,
          ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	  nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	  ne,xforc,&null,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
          &mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,islavsurf);
  
  free(fn);free(inum);free(v);
  
  iout=1;
  
  ad=NNEW(double,neq[1]);
  au=NNEW(double,nzs[1]);
  
  FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounold,nboun,
		    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		    &null,nelemload,sideload,xloadact,&null,xbodyact,ipobody,
		    &null,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		    &nmethodact,ikmpc,ilmpc,ikboun,ilboun,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		    ielmat,ielorien,norien,orab,ntmat_,
		    t0,t1act,&ithermalact,prestr,iprestr,vold,iperturb,sti,
		    nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		    xstiff,npmat_,&dtime,matname,mi,
		    ncmat_,massact,&stiffness,&buckling,&rhsi,&intscheme,
		    physcon,shcon,nshcon,alcon,nalcon,ttime,&time,istep,&iinc,
		    &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		    xstateini,xstate,thicke,integerglob,doubleglob,
		    tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		    pmastsurf,&mortar,clearini));
  
  if(nmethodact==0){
      
      /* error occurred in mafill: storing the geometry in frd format */
      
      ++*kode;
      
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,&ithermalact,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni);
      
      FORTRAN(stop,());
      
  }
  
  for(k=0;k<neq[1];++k){
      b[k]=fext[k]-f[k];
  }
  
  if(*isolver==0){
#ifdef SPOOLES
      spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
	      &symmetryflag,&inputformat,&nzs[2]);
#else
      printf("*ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
  }
  else if((*isolver==2)||(*isolver==3)){
      preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
  }
  else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
#else
      printf("*ERROR in nonlingeo: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
      printf("*ERROR in nonlingeo: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
      printf("*ERROR in nonlingeo: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
  }

  free(au);free(ad);

  v=NNEW(double,mt**nk);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  
  fn=NNEW(double,mt**nk);
  
  inum=NNEW(ITG,*nk);
  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,&ithermalact,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounact,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,&nmethodact,cam,&neq[1],veold,accold,
	  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	  &icmd,ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,
	  emeini,xstaten,eei,enerini,alcon,nalcon,set,nset,istartset,
	  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	  fmpc,nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,
	  &reltime,ne,xforc,&null,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
          &mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,islavsurf);
  
  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
  
  /* reactivating the deactivated elements */
  
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[8*i+6],"L")!=0){
	  ipkon[i]=-ipkon[i]-2;
      }
  }

  /* createinum is called in order to store the nodes and elements
     of the complete structure, not only of the coil */

  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
	      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi));
  
  ++*kode;
  if(*mcs!=0){
      ptime=*ttime+time;
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,kode,filab,een,
	     t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	     nstate_,istep,&iinc,iperturb,ener,mi,output,&ithermalact,qfn,
	     ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
             cdn,&mortar);
  }else{
      
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,
	  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,&ithermalact,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni);
      
  }
  free(inum);free(v);free(fn);

  inomat=NNEW(ITG,*nk);

  /* calculating the magnetic intensity caused by the current */

  FORTRAN(assigndomtonodes,(ne,lakon,ipkon,kon,ielmat,inomat,
			    elcon,ncmat_,ntmat_,mi));

  h0=NNEW(double,3**nk);

  biosav(ipkon,kon,lakon,ne,co,qfx,h0,mi,inomat,nk);

  if(*ithermal==0){free(qfx);free(t0);}
  if(strcmp1(&filab[3567],"ECD ")==0)free(qfn);
  
  /* deactivating the shell elements */
  
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[8*i+6],"L")==0){
	  ipkon[i]=-ipkon[i]-2;
      }
  }
  
/**************************************************************/
/* creating connecting MPC's between the domains              */
/**************************************************************/

/* creating contact ties between the domains */

  if(*istep==1){
      
      nodface=NNEW(ITG,5*6**ne);
      ipoface=NNEW(ITG,*nk);
      
      RENEW(set,char,81*(*nset+3));
      RENEW(istartset,ITG,*nset+3);
      RENEW(iendset,ITG,*nset+3);
      RENEW(ialset,ITG,*nalset+6**ne);
      RENEW(tieset,char,243*(*ntie+5));
      RENEW(tietol,double,3*(*ntie+5));
      
      FORTRAN(createtiedsurfs,(nodface,ipoface,set,istartset,
	      iendset,ialset,tieset,inomat,ne,ipkon,lakon,kon,ntie,
	      tietol,nalset,nk,nset,iactive));
      
      free(nodface);free(ipoface);
      RENEW(set,char,81**nset);
      RENEW(istartset,ITG,*nset);
      RENEW(iendset,ITG,*nset);
      RENEW(ialset,ITG,*nalset);
      RENEW(tieset,char,243**ntie);
      RENEW(tietol,double,3**ntie);
      
      /* tied contact constraints: generate appropriate MPC's */
      
      tiedcontact(ntie,tieset,nset,set,istartset,iendset,ialset,
       lakon,ipkon,kon,tietol,nmpc, &mpcfree,&memmpc_,
       &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
       ithermal,co,vold,&icfd,nmpc_,mi,nk,istep,ikboun,nboun,
       kind1,kind2);
  }
  
/**************************************************************/
/* creating the A.n MPC                                       */
/**************************************************************/

  FORTRAN(generateeminterfaces,(istartset,iendset,
	  ialset,iactive,ipkon,lakon,kon,ikmpc,nmpc,&maxfaces));
  
  for(i=1;i<3;i++){
      imast=iactive[i];
      if(imast==0) continue;

      /* determining the normals on the face */

      imastnode=NNEW(ITG,8*maxfaces);
      xmastnor=NNEW(double,3*8*maxfaces);

      FORTRAN(normalsoninterface,(istartset,iendset,
       ialset,&imast,ipkon,kon,lakon,imastnode,&nmastnode,
       xmastnor,co));

      /* enlarging the fields for MPC's */

      *nmpc_=*nmpc_+nmastnode;
      RENEW(ipompc,ITG,*nmpc_);
      RENEW(labmpc,char,20**nmpc_+1);
      RENEW(ikmpc,ITG,*nmpc_);
      RENEW(ilmpc,ITG,*nmpc_);
      RENEW(fmpc,double,*nmpc_);
      
      /* determining the maximum number of terms;
	 expanding nodempc and coefmpc to accommodate
	 those terms */
      
      neqterms=3*nmastnode;
      index=memmpc_;
      (memmpc_)+=neqterms;
      RENEW(nodempc,ITG,3*memmpc_);
      RENEW(coefmpc,double,memmpc_);
      for(k=index;k<memmpc_;k++){
	  nodempc[3*k-1]=k+1;
      }
      nodempc[3*memmpc_-1]=0;

      /* creating the A.n MPC's */

      FORTRAN(createinterfacempcs,(imastnode,xmastnor,&nmastnode,
	      ikmpc,ilmpc,nmpc,ipompc,nodempc,coefmpc,labmpc,&mpcfree,
              ikboun,nboun));

      free(imastnode);free(xmastnor);
  }
  
  /* determining the new matrix structure */
      
  remastructem(ipompc,&coefmpc,&nodempc,nmpc,
	     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	     labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	     kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	     neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	     &adb,&aub,ithermal,iperturb,mass,mi,ielmat,elcon,
	     ncmat_,ntmat_,inomat);

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
	  for(k=0;k<neq[1];++k){
	      fini[k]=f[k];
	  }
	  if(*nmethod==4){
	      for(k=0;k<mt**nk;++k){
		  veini[k]=veold[k];
	      }
	      for(k=0;k<neq[1];++k){
		  fextini[k]=fext[k];
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
      
      /* prediction of the next solution (only for temperature) */
      
      v=NNEW(double,mt**nk);
      
      if(*ithermal>2){
	  prediction_em(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,v,
		 &iinc,&idiscon,vold,nactdof,mi);
      }
      
      fn=NNEW(double,mt**nk);
      
      iout=-1;
      iperturb_sav[0]=iperturb[0];
      iperturb_sav[1]=iperturb[1];
      
      /* first iteration in first increment: heat tangent */
      
      inum=NNEW(ITG,*nk);
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
              nslavnode,ntie);
      free(inum);
     
      /* the calculation of the electromagnetic fields is (quasi)linear,
         i.e. the solution of the equations is the fields;
         only the temperature calculation is nonlinear,
         i.e. the solution of the equations is a differential temperature */
 
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
      
      iout=0;
      
      free(fn);free(v);
      
      /***************************************************************/
      /* iteration counter and start of the loop over the iterations */
      /***************************************************************/
      
      iit=1;
      icntrl=0;
      ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
      
      while(icntrl==0){
	  
	  if(iit!=1){
	      
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
	      if(*ithermal>1){radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,
                bcr,ipivr,ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,
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
	      
	  }
	      
	  /* calculating the local stiffness matrix and external loading */
	  
	  ad=NNEW(double,neq[1]);
	  au=NNEW(double,nzs[1]);
	  
	  FORTRAN(mafillem,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
		  xbounact,nboun,
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
		  tieset,istartset,iendset,ialset,ntie,&nasym,iactive,h0,
                  pslavsurf,pmastsurf,&mortar,clearini));
	      
	      iperturb[0]=iperturb_sav[0];
	      iperturb[1]=iperturb_sav[1];

	  /* calculating the residual (f is only for the temperature
             nonzero) */

	  calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	    vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
	    islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
	    nzs,&nasym,ithermal);
	  
	  newstep=0;
	  
	  if(*nmethod==0){
	      
	      /* error occurred in mafill: storing the geometry in frd format */
	      
	      *nmethod=0;
	      ++*kode;
	      inum=NNEW(ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	      
	      ptime=*ttime+time;
	      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
		  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni);
	      
	  }
	  
	  /* implicit step (static or dynamic */
	  
	  if(*nmethod==4){
	      
	      /* mechanical part */
	      
	      if(*ithermal!=2){
		  for(k=0;k<neq[0];++k){
		      ad[k]=adb[k]/dtime+ad[k];
		  }
		  for(k=0;k<nzs[0];++k){
		      au[k]=aub[k]/dtime+au[k];
		  }
		  
		  /* upper triangle of asymmetric matrix */
		  
		  if(nasym>0){
		      for(k=nzs[2];k<nzs[2]+nzs[0];++k){
			  au[k]=aub[k]/dtime+au[k];
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
	      }else{
		  pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			       &symmetryflag,&inputformat,jq,&nzs[2]);
	      }
#else
	      printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  
	  /* calculating the electromagnetic fields and temperatures 
             only the temperature calculation is differential */
	  
	  v=NNEW(double,mt**nk);
	  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	  
	  fn=NNEW(double,mt**nk);
	  
	  inum=NNEW(ITG,*nk);
	  resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		  &icmd,ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,
		  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
                  nslavnode,ntie);
	  free(inum);
	  
	  free(ad);free(au);
	  
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
	  
	  free(v);free(fn);
	  
	  /* calculating the residual */
	  
	  calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	      vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
	      islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
	      nzs,&nasym,ithermal);
	  
	  /* calculating the maximum residual (only thermal part)*/
	  
	  for(k=0;k<2;++k){
	      ram2[k]=ram1[k];
	      ram1[k]=ram[k];
	      ram[k]=0.;
	  }

	  if(*ithermal!=2) ram[0]=0.;

	  if(*ithermal>1){
	      for(k=neq[0];k<neq[1];++k){
		  err=fabs(b[k]);
		  if(err>ram[1]){ram[1]=err;ram[3]=k+0.5;}
	      }
	  }
	  
	  /* printing residuals */
	  
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
	  
	  /* athermal electromagnetic calculations are linear:
             set iit=2 to force convergence */

	  if(*ithermal<=1) iit=2;

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
	    set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
            &mortar);
	  
      }
      
      /*********************************************************/
      /*   end of the iteration loop                          */
      /*********************************************************/
      
      /* icutb=0 means that the iterations in the increment converged,
	 icutb!=0 indicates that the increment has to be reiterated with
	 another increment size (dtheta) */
      
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
	      }
	      for(k=0;k<neq[1];++k){
		  fext[k]=fextini[k];
	      }
	  }
	  
	  qam[0]=qamold[0];
	  qam[1]=qamold[1];
      }
      
      if((jout[0]==jprint)&&(icutb==0)){
	  
	  jprint=0;
	  
	  /* calculating the displacements and the stresses and storing */
	  /* the results in frd format  */
	  
	  v=NNEW(double,mt**nk);
	  fn=NNEW(double,mt**nk);
	  if(*ithermal>1) qfn=NNEW(double,3**nk);
	  if((strcmp1(&filab[3741],"EMFE")==0)||
             (strcmp1(&filab[3828],"EMFB")==0)) stn=NNEW(double,6**nk);
	  inum=NNEW(ITG,*nk);
	  
	  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	  
	  iout=2;
	  icmd=3;
	  
	  resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,
		  xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		  nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
                  nslavnode,ntie);
	  
	  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
	  iout=0;
	  icmd=0;
	  FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
	  
	  ++*kode;
	  if(*mcs!=0){
	      ptime=*ttime+time;
	      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
		 t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
		 nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
		 ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
		  norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
                 cdn,&mortar);
	  }else{
	      
	      ptime=*ttime+time;
	      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
		  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,
                  enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni);
	      
	  }
	  
	  free(v);free(fn);free(inum);
	  if(*ithermal>1){free(qfn);}
	  if(strcmp1(&filab[3741],"EMF ")==0) free(stn);
	  
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
      if(*ithermal>1) qfn=NNEW(double,3**nk);
      if(strcmp1(&filab[3741],"EMF ")==0) stn=NNEW(double,6**nk);
      inum=NNEW(ITG,*nk);
      
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      iout=2;
      icmd=3;
      
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	      ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,
	      xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
              nslavnode,ntie);
      
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
      
      iout=0;
      icmd=0;
      FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
      
      ++*kode;
      if(*mcs>0){
	  ptime=*ttime+time;
	  frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
		 t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
		 nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
		 ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
		 norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
                 cdn,&mortar);
      }else{
	  
	  ptime=*ttime+time;
	  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	      kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	      thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni);
	  
      }
      
      free(v);free(fn);free(inum);
      if(*ithermal>1){free(qfn);}
      if(strcmp1(&filab[3741],"EMF ")==0) free(stn);
      
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
  if(*nbody>0) free(ipobody);
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
  
  free(fini);
  if(*nmethod==4){
      free(aux2);free(fextini);free(veini);
      free(adb);free(aub);
  }
  
  free(aux);free(iaux);free(vini);free(h0);free(inomat);
  
  /* reset icascade */
  
  if(icascade==1){icascade=0;}
  
  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;
  
  *icolp=icol;*irowp=irow;*cop=co;*voldp=vold;
  
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
  
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;

  *setp=set;*istartsetp=istartset;*iendsetp=iendset;*ialsetp=ialset;
  *tiesetp=tieset;*tietolp=tietol;
  
  (*tmin)*=(*tper);
  (*tmax)*=(*tper);
  
  free(nactdofinv);
 
  if(*nmethod==1){
      *nmethod=8;
  }else if(*nmethod==4){
      *nmethod=9;
  }else if(*nmethod==2){
      *nmethod=10;
  }

  (*ttime)+=(*tper);
  
  return;
}