/*     CalculiX - A 3-dimensional finite element program                 */
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


void electromagnetics(double **cop, int *nk, int **konp, int **ipkonp, 
             char **lakonp,int *ne, 
	     int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int **ipompcp, int **nodempcp, double **coefmpcp, char **labmpcp,
             int *nmpc, 
	     int *nodeforc, int *ndirforc,double *xforc, int *nforc, 
	     int *nelemload, char *sideload, double *xload,int *nload, 
	     double *ad, double *au, double *b, int *nactdof, 
	     int **icolp, int *jq, int **irowp, int *neq, int *nzl, 
	     int *nmethod, int **ikmpcp, int **ilmpcp, int *ikboun, 
	     int *ilboun,
             double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int **ielmatp,
	     int **ielorienp, int *norien, double *orab, int *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     int *ithermal,double *prestr, int *iprestr, 
	     double **voldp,int *iperturb, double *sti, int *nzs,  
	     int *kode, double *adb, double *aub,char *filab, 
             int *idrct, int *jmax, int *jout, double *tinc,
             double *tper, double *tmin, double *tmax, double *eme,
	     double *xbounold, double *xforcold, double *xloadold,
             double *veold, double *accold,
	     char *amname, double *amta, int *namta, int *nam,
             int *iamforc, int *iamload,
             int *iamt1, double *alpha, int *iexpl,
	     int *iamboun, double *plicon, int *nplicon, double *plkcon,
             int *nplkcon,
             double **xstatep, int *npmat_, int *istep, double *ttime,
             char *matname, double *qaold, int *mi,
             int *isolver, int *ncmat_, int *nstate_, int *iumat,
             double *cs, int *mcs, int *nkon, double **enerp, int *mpcinfo,
             int *nnn, char *output,
             double *shcon, int *nshcon, double *cocon, int *ncocon,
             double *physcon, int *nflow, double *ctrl, 
             char *set, int *nset, int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener,int *ikforc,int *ilforc, double *trab, 
             int *inotr, int *ntrans, double **fmpcp, char *cbody,
             int *ibody, double *xbody, int *nbody, double *xbodyold,
             int *ielprop, double *prop, int *ntie, char *tieset,
	     int *itpamp, int *iviewfile, char *jobnamec, double *tietol,
	     int *nslavs, double *thicke, int *ics, int *nalset, int *nmpc_){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
      *labmpc=NULL,kind1[2]="E",kind2[2]="E"; 
 
  int *inum=NULL,k,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
      icutb=0,istab=0,ifreebody,uncoupled,n1,n2,maxfaces,
      iperturb_sav[2],*icol=NULL,*irow=NULL,ielas=0,icmd=0,
      memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
      *itg=NULL,*ineighe=NULL,null=0,iactive[3],neqterms,
      *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,index,
      *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,
      ntrit,*inocs=NULL,inewton=0,*ipobody=NULL,*nacteq=NULL,
      *nactdog=NULL,nteq,network,nmastnode,imast,
      *ipkon=NULL,*kon=NULL,*ielorien=NULL,
      *ielmat=NULL,inext,itp=0,symmetryflag=0,inputformat=0,
      *iruc=NULL,iitterm=0,ngraph=1,
      *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
      *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,mortar=0,
      mt=mi[1]+1,*nactdofinv=NULL, inode,idir,
      iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
      *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
      coriolis=0,*ipneigh=NULL,*neigh=NULL,i,icfd=0,id,node,networknode,
      iflagact=0,*nodorig=NULL,*ipivr=NULL,*inomat=NULL,*nodface=NULL,
      *ipoface=NULL;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,
         *f=NULL,*fn=NULL,qa[3]={0.,0.,-1.},qam[2]={0.,0.},dtheta,theta,
         err,ram[4]={0.,0.,0.,0.},*springarea=NULL,*h0=NULL,
	 ram1[2]={0.,0.},ram2[2]={0.,0.},deltmx,
         uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,
	 *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
	 *xforcact=NULL,*xloadact=NULL,*fext=NULL,
         reltime,time,bet=0.,gam=0.,*aux1=NULL,*aux2=NULL,dtime,*fini=NULL,
         *fextini=NULL,*veini=NULL,*xstateini=NULL,
	 *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,
         *xbounini=NULL,*xstiff=NULL,*stx=NULL,
         *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
	 *enerini=NULL,*emn=NULL,*xmastnor=NULL,
	 *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
	 *adview=NULL,*auview=NULL,*qfx=NULL,
         *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
         *xbodyact=NULL,*cgr=NULL,dthetaref,*vr=NULL,*vi=NULL,
	 *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
         *f_cm=NULL,*f_cs=NULL,
	 *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
         *xnormastface=NULL,*emeini=NULL,*doubleglob=NULL;

#ifdef SGI
  int token;
#endif
 
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
  
  /* invert nactdof */
  
  nactdofinv=NNEW(int,mt**nk);nodorig=NNEW(int,*nk);
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
  iaux=NNEW(int,maxlenmpc);
  
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
      ipobody=NNEW(int,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,int,2*(*ne+ifreebody));
      }
      RENEW(ipobody,int,2*(ifreebody-1));
  }
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
      itg=NNEW(int,*nload+3**nflow);
      ieg=NNEW(int,*nflow);
      /* max 6 triangles per face, 4 entries per triangle */
      kontri=NNEW(int,24**nload);
      nloadtr=NNEW(int,*nload);
      nacteq=NNEW(int,4**nk);
      nactdog=NNEW(int,4**nk);
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
	  inocs=NNEW(int,*nk);
	  radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
               iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
      }
      else{ntrit=ntri;}
      
      nzsrad=100*ntr;
      mast1rad=NNEW(int,nzsrad);
      irowrad=NNEW(int,nzsrad);
      icolrad=NNEW(int,ntr);
      jqrad=NNEW(int,ntr+1);
      ipointerrad=NNEW(int,ntr);
      
      if(ntr>0){
	  mastructrad(&ntr,nloadtr,sideload,ipointerrad,
		      &mast1rad,&irowrad,&nzsrad,
		      jqrad,icolrad);
      }
      
      free(ipointerrad);free(mast1rad);
      RENEW(irowrad,int,nzsrad);
      
      RENEW(itg,int,ntg);
      ineighe=NNEW(int,ntg);
      RENEW(kontri,int,4*ntrit);
      RENEW(nloadtr,int,ntr);
      
      adview=NNEW(double,ntr);
      auview=NNEW(double,2*nzsrad);
      tarea=NNEW(double,ntr);
      tenv=NNEW(double,ntr);
      fenv=NNEW(double,ntr);
      erad=NNEW(double,ntr);
      
      ac=NNEW(double,nteq*nteq);
      bc=NNEW(double,nteq);
      ipiv=NNEW(int,nteq);
      adrad=NNEW(double,ntr);
      aurad=NNEW(double,2*nzsrad);
      bcr=NNEW(double,ntr);
      ipivr=NNEW(int,ntr);
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
  
  /* calculating the initial magnetic intensity due to the coil current */
  
  /*********************************************************************/
  
  /* calculate the current density in the coils

     in this section nload, nforc, nbody and nam are set to zero; the
     electrical potential is supposed to be given (in the form of a
     temperature), the current is calculated (in the form of heat
     flux) by thermal analogy  */
  
  reltime=1.;
  time=0.;
  dtime=0.;
  
  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,&null,xloadold,xload,
	      xloadact,iamload,&null,ibody,xbody,&null,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,&null,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
  
  cam[0]=0.;cam[1]=0.;
  
  /* deactivating all elements except the shells */
  
  for(i=0;i<*ne;i++){
      if(strcmp1(&lakon[8*i+6],"L")!=0){
	  ipkon[i]=-ipkon[i]-1;
      }
  }
  
  remastruct(ipompc,&coefmpc,&nodempc,nmpc,
	     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	     labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	     kon,ipkon,lakon,ne,nnn,nactdof,icol,jq,&irow,isolver,
	     neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	     &adb,&aub,ithermal,iperturb,mass,mi);
  
  /* invert nactdof */
  
  free(nactdofinv);nactdofinv=NNEW(int,mt**nk);nodorig=NNEW(int,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  free(nodorig);
  
  iout=-1;
  
  fn=NNEW(double,mt**nk);
  
  inum=NNEW(int,*nk);
  resultsinduction(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,stx,
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
	  nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	  ne,xforc,&null,thicke,xnormastface,shcon,nshcon,
	     sideload,xloadact,xloadold,&icfd,inomat,h0);
  
  free(fn);free(inum);
  
  iout=0;
  
  ad=NNEW(double,neq[1]);
  au=NNEW(double,nzs[1]);
  
  FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounold,nboun,
		    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		    &null,nelemload,sideload,xloadact,&null,xbodyact,ipobody,
		    &null,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		    nmethod,ikmpc,ilmpc,ikboun,ilboun,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		    ielmat,ielorien,norien,orab,ntmat_,
		    t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		    nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		    xstiff,npmat_,&dtime,matname,mi,
		    ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		    physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		    &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		    xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
		    tieset,istartset,iendset,ialset,ntie,&nasym));
  
  if(nmethod==0){
      
      /* error occurred in mafill: storing the geometry in frd format */
      
      ++*kode;
      
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx);
      
      FORTRAN(stop,());
      
  }
  
  for(k=0;k<neq[0];++k){
      b[k]=fext[k]-f[k];
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

  v=NNEW(double,mt**nk);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  
  fn=NNEW(double,mt**nk);
  
  inum=NNEW(int,*nk);
  resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
	  fmpc,nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,
	  &reltime,ne,xforc,&null,thicke,xnormastface,shcon,nshcon,
	     sideload,xloadact,xloadold,&icfd,inomat,h0);
  free(inum);
  
  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
  
  ++*kode;
  if(*mcs!=0){
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	     t1act,fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	     nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
	     ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne);
  }else{
      
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx);
      
  }

  /* calculating the magnetic intensity caused by the current */

  FORTRAN(assigndomtonodes,(ne,lakon,ipkon,kon,ielmat,inomat,
			    elcon,ncmat_,ntmat_,mi));

  h0=NNEW(double,3**nk);

  biosav(ipkon,kon,lakon,ne,co,qfx,h0,mi,inomat,nk);
  
/* reverting the active elements */
  
  for(i=0;i<*ne;i++){
      ipkon[i]=-ipkon[i]-2;
  }
  
  free(ad);free(au);
  
/**************************************************************/
/* creating connecting MPC's between the domains              */
/**************************************************************/

/* creating contact ties between the domains */

  if(*istep==1){
      
      nodface=NNEW(int,5*6**ne);
      ipoface=NNEW(int,*nk);
      
      RENEW(set,char,81*(*nset+3));
      RENEW(istartset,int,*nset+3);
      RENEW(iendset,int,*nset+3);
      RENEW(ialset,int,*nalset+6**ne);
      RENEW(tieset,char,243*(*ntie+5));
      RENEW(tietol,double,2*(*ntie+5));
      
      FORTRAN(createtiedsurfs,(nodface,ipoface,set,istartset,
	      iendset,ialset,tieset,inomat,ne,ipkon,lakon,kon,ntie,
	      tietol,nalset,nk,nset,iactive));
      
      free(nodface);free(ipoface);
      RENEW(set,char,81**nset);
      RENEW(istartset,int,*nset);
      RENEW(iendset,int,*nset);
      RENEW(ialset,int,*nalset);
      RENEW(tieset,char,243**ntie);
      RENEW(tietol,double,2**ntie);
      
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

      /* determining the normals on the face */

      imastnode=NNEW(int,8*maxfaces);
      xmastnor=NNEW(double,3*8*maxfaces);

      FORTRAN(normalsoninterface,(istartset,iendset,
       ialset,&imast,ipkon,kon,lakon,imastnode,&nmastnode,
       xmastnor,co));

      /* enlarging the fields for MPC's */

      nmpc_=nmpc_+nmastnode;
      RENEW(ipompc,int,*nmpc_);
      RENEW(labmpc,char,20**nmpc_+1);
      RENEW(ikmpc,int,*nmpc_);
      RENEW(ilmpc,int,*nmpc_);
      RENEW(fmpc,double,*nmpc_);
      
      /* determining the maximum number of terms;
	 expanding nodempc and coefmpc to accommodate
	 those terms */
      
      neqterms=3*nmastnode;
      index=memmpc_;
      (memmpc_)+=neqterms;
      RENEW(nodempc,int,3*memmpc_);
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
      printf(" increment %d attempt %d \n",iinc,icutb+1);
      printf(" increment size= %e\n",dtheta**tper);
      printf(" sum of previous increments=%e\n",theta**tper);
      printf(" actual step time=%e\n",(theta+dtheta)**tper);
      printf(" actual total time=%e\n\n",*ttime+dtheta**tper);
      
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
      
      /* prediction of the kinematic vectors  */
      
      v=NNEW(double,mt**nk);
      
      prediction_em(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,v,
		 &iinc,&idiscon,vold,nactdof,mi);
      
      fn=NNEW(double,mt**nk);
      
      iout=-1;
      iperturb_sav[0]=iperturb[0];
      iperturb_sav[1]=iperturb[1];
      
      /* first iteration in first increment: elastic tangent */
      
      inum=NNEW(int,*nk);
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
	      &reltime,ne,xforc,nforc,thicke,xnormastface,shcon,nshcon,
		 sideload,xloadact,xloadold,&icfd,inomat,h0);
      free(inum);
      
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
	  
	  if((iit!=1)||((uncoupled)&&(*ithermal==1))){
	      
	      printf(" iteration %d\n\n",iit);
	      
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
	  
	  FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
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
                  xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
		  tieset,istartset,iendset,ialset,ntie,&nasym));
	      
	      iperturb[0]=iperturb_sav[0];
	      iperturb[1]=iperturb_sav[1];

	  /* calculating the residual */

	  calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	    vini,&dtime,accold,nk,adb,aub,icol,irow,nzl,alpha,fextini,fini,
	    islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
	    nzs,&nasym,ad,au);
	  
	  newstep=0;
	  
	  if(*nmethod==0){
	      
	      /* error occurred in mafill: storing the geometry in frd format */
	      
	      *nmethod=0;
	      ++*kode;
	      inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	      
	      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
		  kode,filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx);
	      
	  }
	  
	  /* implicit step (static or dynamic */
	  
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
	  
	  /* calculating the displacements, stresses and forces */
	  
	  v=NNEW(double,mt**nk);
	  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	  
	  fn=NNEW(double,mt**nk);
	  
	  inum=NNEW(int,*nk);
	  resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
		  &reltime,ne,xforc,nforc,thicke,xnormastface,shcon,nshcon,
		     sideload,xloadact,xloadold,&icfd,inomat,h0);
	  free(inum);
	  
	  /* calculating the residual */
	  
	  calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	      vini,&dtime,accold,nk,adb,aub,icol,irow,nzl,alpha,fextini,fini,
	      islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
	      nzs,&nasym,ad,au);
	  
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
	  
	  /* next line is inserted to cope with stress-less
	     temperature calculations */
	  
	  if(*ithermal!=2){
	      if(ram[0]<1.e-6) ram[0]=0.;      
	      printf(" average force= %f\n",qa[0]);
	      printf(" time avg. forc= %f\n",qam[0]);
	      if((int)((double)nactdofinv[(int)ram[2]]/mt)+1==0){
		  printf(" largest residual force= %f\n",
			 ram[0]);
	      }else{
		  inode=(int)((double)nactdofinv[(int)ram[2]]/mt)+1;
		  idir=nactdofinv[(int)ram[2]]-mt*(inode-1);
		  printf(" largest residual force= %f in node %d and dof %d\n",
			 ram[0],inode,idir);
	      }
	      printf(" largest increment of disp= %e\n",uam[0]);
	      if((int)cam[3]==0){
		  printf(" largest correction to disp= %e\n\n",
			 cam[0]);
	      }else{
		  inode=(int)((double)nactdofinv[(int)cam[3]]/mt)+1;
		  idir=nactdofinv[(int)cam[3]]-mt*(inode-1);
		  printf(" largest correction to disp= %e in node %d and dof %d\n\n",cam[0],inode,idir);
	      }
	  }
	  if(*ithermal>1){
	      if(ram[1]<1.e-6) ram[1]=0.;      
	      printf(" average flux= %f\n",qa[1]);
	      printf(" time avg. flux= %f\n",qam[1]);
	      if((int)((double)nactdofinv[(int)ram[3]]/mt)+1==0){
		  printf(" largest residual flux= %f\n",
			 ram[1]);
	      }else{
		  inode=(int)((double)nactdofinv[(int)ram[3]]/mt)+1;
		  idir=nactdofinv[(int)ram[3]]-mt*(inode-1);
		  printf(" largest residual flux= %f in node %d and dof %d\n",ram[1],inode,idir);
	      }
	      printf(" largest increment of temp= %e\n",uam[1]);
	      if((int)cam[4]==0){
		  printf(" largest correction to temp= %e\n\n",
			 cam[1]);
	      }else{
		  inode=(int)((double)nactdofinv[(int)cam[4]]/mt)+1;
		  idir=nactdofinv[(int)cam[4]]-mt*(inode-1);
		  printf(" largest correction to temp= %e in node %d and dof %d\n\n",cam[1],inode,idir);
	      }
	  }
	  
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
	    set,nset,istartset,iendset,ialset,emn,thicke,jobnamec);
	  
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
	  inum=NNEW(int,*nk);
	  
	  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	  
	  iout=2;
	  icmd=3;
	  
	  resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
		  &reltime,ne,xforc,nforc,thicke,xnormastface,shcon,nshcon,
		     sideload,xloadact,xloadold,&icfd,inomat,h0);
	  
	  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
	  iout=0;
	  icmd=0;
	  FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
	  
	  ++*kode;
	  if(*mcs!=0){
	      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
		 t1act,fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
		 nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
		 ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
		  norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne);
	  }else{
	      
	      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
		  kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,
                  enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx);
	      
	  }
	  
	  free(v);free(fn);free(inum);
	  if(*ithermal>1){free(qfn);}
	  
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
      inum=NNEW(int,*nk);
      
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      iout=2;
      icmd=3;
      
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
	      &reltime,ne,xforc,nforc,thicke,xnormastface,shcon,nshcon,
		 sideload,xloadact,xloadold,&icfd,inomat,h0);
      
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
      
      iout=0;
      icmd=0;
      FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
      
      ++*kode;
      if(*mcs>0){
	  frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
		 t1act,fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
		 nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
		 ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
		 norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne);
      }else{
	  
	  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	      kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	      thicke,jobnamec,output,qfx);
	  
      }
      
      free(v);free(fn);free(inum);
      if(*ithermal>1){free(qfn);}
      
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
  
  free(aux);free(iaux);free(vini);free(h0);
  
  /* reset icascade */
  
  if(icascade==1){icascade=0;}
  
  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;
  
  *icolp=icol;*irowp=irow;*cop=co;*voldp=vold;
  
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
  
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;
  
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
  
  return;
}
