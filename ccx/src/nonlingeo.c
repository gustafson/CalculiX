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


void nonlingeo(double **cop, int *nk, int **konp, int **ipkonp, char **lakonp,
	     int *ne, 
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
	     int *nslavs, double *thicke, int *ics){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
      *sideface=NULL,*labmpc=NULL,fnffrd[132]=""; 
 
  int *inum=NULL,k,l,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
       icutb=0,istab=0,ifreebody,uncoupled,n1,n2,nzlc,im,regmode,
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
       *nslavnode=NULL,*nmastnode=NULL,mortar=0,*imastop=NULL,
       *iponoels=NULL,*inoels=NULL,nzsc,*irowc=NULL,*jqc=NULL,
       *islavact=NULL,*irowqdt=NULL,*jqqdt=NULL,nzsqdt,*icolc=NULL,
       *irowbd=NULL,*jqbd=NULL,mt=mi[1]+1,*nactdofinv=NULL,*ipe=NULL, 
       *ime=NULL,*ikactmech=NULL,nactmech,ifacecount,inode,idir,neold,
       iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
       *jqrad=NULL,*ipointerrad=NULL,*iwatchactiv=NULL,*integerglob=NULL;

  int mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
    coriolis=0,*ipneigh=NULL,*neigh=NULL,irenewxstate,j,
    *nelemface=NULL,*ipoface=NULL,*nodface=NULL,*ifreestream=NULL,
    *isolidsurf=NULL,*neighsolidsurf=NULL,*iponoel=NULL,*inoel=NULL,
    nef=0,nface,nfreestream,nsolidsurf,inoelfree,i,indexe,icfd=0,id,
    node,networknode,*jqtemp=NULL,*icoltemp=NULL,*irowtemp=NULL,
    nzstemp[3],iflagact=0,*nodorig=NULL,*ipivr=NULL,iflag_fric,iglob=0,
    *nslavspc=NULL,*islavspc=NULL,nsspc,*nslavmpc=NULL,*islavmpc=NULL,nsmpc,
    *nmastspc=NULL,*imastspc=NULL,nmspc,*nmastmpc=NULL,*imastmpc=NULL,nmmpc,
    *islavborder=NULL,*jqb=NULL,*irowb=NULL,nzsbd2, *islavactdof=NULL,
    *islavactini=NULL,iflagact_old,*inomat=NULL,*islavactdoftie=NULL;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cg=NULL,mu,fkinv,
         *f=NULL,*fn=NULL,qa[3]={0.,0.,-1.},qam[2]={0.,0.},dtheta,theta,
	 err,ram[4]={0.,0.,0.,0.},*areaslav=NULL,*springarea=NULL,p0,beta,
	 ram1[2]={0.,0.},ram2[2]={0.,0.},deltmx,*auc=NULL,*adc=NULL,
         uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,*straight=NULL,
	 *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,*bdd=NULL,
	 *xforcact=NULL,*xloadact=NULL,*fext=NULL,*gap=NULL,
         reltime,time,bet=0.,gam=0.,*aux1=NULL,*aux2=NULL,dtime,*fini=NULL,
         *fextini=NULL,*veini=NULL,*accini=NULL,*xstateini=NULL,
	 *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,*auqdt=NULL,
         *xbounini=NULL,dev,*xstiff=NULL,*stx=NULL,*stiini=NULL,
         *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
	 *coefmpcref=NULL,*enerini=NULL,*slavnor=NULL,*emn=NULL,
	 *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
	 *adview=NULL,*auview=NULL,*qfx=NULL,*bhat=NULL,
         *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
         *xbodyact=NULL,*cgr=NULL,dthetaref, *vcontu=NULL,*vr=NULL,*vi=NULL,
	 *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
         *cstress=NULL,*cdisp=NULL,*aubd=NULL, *f_cm=NULL, *f_cs=NULL,
	 *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
         *xmastnor=NULL,*xnormastface=NULL,*emeini=NULL,*slavtan=NULL,
	 *bp=NULL,*doubleglob=NULL,*pslavdual=NULL,*Bd=NULL,*Dd=NULL,
	 *dhinv=NULL,*bpini=NULL,*cstressini=NULL,*xnoels=NULL;

#ifdef SGI
  int token;
#endif

  if(*ithermal==4){
      uncoupled=1;
      *ithermal=3;
  }else{
      uncoupled=0;
  }
  
  /* turbulence model 
     iturbulent==0: laminar
     iturbulent==1: k-epsilon
     iturbulent==2: q-omega
     iturbulent==3: SST */
  
  iturbulent=(int)physcon[8];
  
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
  
  /* check whether the submodel is meant for a fluid-structure
     interaction */
  
  strcpy(fnffrd,jobnamec);
  strcat(fnffrd,"f.frd");
  if(strcmp1(&jobnamec[396],fnffrd)==0){
      
      /* fluid-structure interaction: wait till after the compfluid
         call */
      
      integerglob=NNEW(int,1);
      doubleglob=NNEW(double,1);
  }else{
      
      /* determining the global values to be used as boundary conditions
	 for a submodel */
      
      getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,xboun,
		       nload,sideload,iamload,&iglob);
  }
  
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
      if(inewton==1){cgr=NNEW(double,4**ne);}
  }
  
  /* for mechanical calculations: updating boundary conditions
     calculated in a previous thermal step */
  
  if(*ithermal<2) FORTRAN(gasmechbc,(vold,nload,sideload,
				     nelemload,xload,mi));
  
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
  
  /* check for fluid elements */
  
  for(i=0;i<*ne;++i){
      if(ipkon[i]<0) continue;
      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i],"F")==0){icfd=1;nef++;}
  }
  if(icfd==1){
      sideface=NNEW(char,6*nef);
      nelemface=NNEW(int,6*nef);
      ipoface=NNEW(int,*nk);
      nodface=NNEW(int,5*6*nef);
      ifreestream=NNEW(int,*nk);
      isolidsurf=NNEW(int,*nk);
      neighsolidsurf=NNEW(int,*nk);
      iponoel=NNEW(int,*nk);
      inoel=NNEW(int,3*20*nef);
      inomat=NNEW(int,*nk);
      FORTRAN(precfd,(nelemface,sideface,&nface,ipoface,nodface,
	      ne,ipkon,kon,lakon,ikboun,ilboun,xboun,nboun,nk,isolidsurf,
	      &nsolidsurf,ifreestream,&nfreestream,neighsolidsurf,
	      iponoel,inoel,&inoelfree,&nef,co,ipompc,nodempc,ikmpc,ilmpc,
	      nmpc,set,istartset,iendset,ialset,nset,&iturbulent,inomat,
              ielmat));
      RENEW(sideface,char,nface);
      RENEW(nelemface,int,nface);
      free(ipoface);free(nodface);
      RENEW(ifreestream,int,nfreestream);
      RENEW(isolidsurf,int,nsolidsurf);
      RENEW(neighsolidsurf,int,nsolidsurf);
      RENEW(inoel,int,3*inoelfree);
      vcontu=NNEW(double,2**nk);
  }
  if(*ithermal>1){qfx=NNEW(double,3*mi[0]**ne);}
  
  /* contact conditions */
  
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
      straight=NNEW(double,16*ncont);
      
      if(mortar==1){
	  
	  /* adding one element per slave node, similar to 
	     spring elements;
	     needed for output in frd-format of CDISP and CSTRES */
	  
	  RENEW(kon,int,*nkon+*nslavs);
	  if((irenewxstate==1)&&(*nslavs!=0)){
	      RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	      for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
		  xstate[k]=0.;
	      }
	  }      
	  for(k=0;k<*nslavs;k++){
	      ipkon[*ne+k]=*nkon+k;
	      kon[*nkon+k]=islavnode[k];
	      strcpy1(&lakon[8*(*ne+k)]," S    C1",8);
	  }
	  islavactdoftie=NNEW(int,nslavnode[*ntie]);
	  bp=NNEW(double,nslavnode[*ntie]);
	  islavact=NNEW(int,nslavnode[*ntie]);
	  gap=NNEW(double,nslavnode[*ntie]);
	  slavnor=NNEW(double,3*nslavnode[*ntie]);
	  slavtan=NNEW(double,6*nslavnode[*ntie]);
	  cdisp=NNEW(double,6*nslavnode[*ntie]);
	  
	  irowb=NNEW(int,*nk);
	  Bd=NNEW(double,*nk);
	  nzsbd2=*nk;
	  Dd=NNEW(double,*nk);
	  jqb=NNEW(int, *nk+1); 
	  
	  /* allocation of temperary fields: stores the structure
	     of the stiffness matrix without mortar contact */
	  
	  cstress=NNEW(double,3*nslavnode[*ntie]);
	  
	  /* fields for cutback  */
	  bpini=NNEW(double,nslavnode[*ntie]);
	  islavactini=NNEW(int,nslavnode[*ntie]);  
	  cstressini=NNEW(double,3*nslavnode[*ntie]);
	  
	  /* check for friction */
	  iflag_fric=0;
	  printf("\t***mortar***\n");
	  for (i=0;i<*ntie;i++){
	      if(tieset[i*(81*3)+80]=='C'){
		  FORTRAN(getcontactparams,(&mu,&regmode,&fkinv,&p0,&beta,
					    tietol,elcon,&i,ncmat_,ntmat_));
		  printf("\ttie %d mu %e regmode %d kinv %e p0 %e beta %e \n",
			 i+1,mu,regmode,fkinv, p0, beta );
		  for(j=nslavnode[i];j<nslavnode[i+1];j++){
		      islavactdoftie[j]=i;
		      if(mu >1.e-10){iflag_fric=1;}
		  }
	      }
	  }	
	  printf("\tiflag_firc %d \n \t************\n",iflag_fric);
	  /* coeffs for dual basis functions */
	  pslavdual=NNEW(double,16*itiefac[2**ntie-1]);
	  
	  /* checking for SPC's and MPC's on slave and master surface */
	  nslavspc=NNEW(int,2*nslavnode[*ntie]);
	  islavspc=NNEW(int,2**nboun);
	  nslavmpc=NNEW(int,2*nslavnode[*ntie]);
	  islavmpc=NNEW(int,2**nmpc);
	  nmastspc=NNEW(int,2*nmastnode[*ntie]);
	  imastspc=NNEW(int,2**nboun);
	  nmastmpc=NNEW(int,2*nmastnode[*ntie]);
	  imastmpc=NNEW(int,2**nmpc);
	  islavborder=NNEW(int,nslavnode[*ntie]);      
	  FORTRAN(conttiemortar,(lakon,ipkon,kon,ntie,tieset,nset,set,
		itiefac,islavsurf,islavnode,imastnode,nslavnode,nmastnode,
		nslavs,iponoels,inoels,ipoface,nodface,nk,nboun,ndirboun,
		nodeboun,xboun,nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,
                ikmpc,ilmpc,nslavspc,islavspc,&nsspc,nslavmpc,islavmpc,&nsmpc,
		nmastspc,imastspc,&nmspc,nmastmpc,imastmpc,&nmmpc,
	        islavborder));       
	  RENEW(islavspc,int,2*nsspc+1);
	  RENEW(islavmpc,int,2*nsmpc+1);
	  RENEW(imastspc,int,2*nmspc+1);
	  RENEW(imastmpc,int,2*nmmpc+1);       
	  dhinv=NNEW(double,9*nslavnode[*ntie]);      
      }else{
	  
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
      }
  }
  
  if((icascade==2)||(ncont!=0)){
      memmpcref_=memmpc_;mpcfreeref=mpcfree;
      nodempcref=NNEW(int,3*memmpc_);
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
  
  if(*nstate_!=0){
      xstateini=NNEW(double,*nstate_*mi[0]*(*ne+*nslavs));
      for(k=0;k<*nstate_*mi[0]*(*ne+*nslavs);++k){
	  xstateini[k]=xstate[k];
      }
  }
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
  
  ne0=*ne;nkon0=*nkon;
  
  /*********************************************************************/
  
  /* calculating the initial acceleration at the start of the step
     for dynamic calculations */
  
  /*********************************************************************/
  
  if((*nmethod==4)&&(*ithermal!=2)){
      bet=(1.-*alpha)*(1.-*alpha)/4.;
      gam=0.5-*alpha;
      
      /* calculating the stiffness and mass matrix 
	 the stress must be determined to calculate the 
	 stiffness matrix*/
      
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
      if(*ithermal>1){radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
       ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,vold,
       shcon,nshcon,ipkon,kon,co,
       kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
       nflow,ikboun,xboun,nboun,ithermal,&iinc,&iit,
       cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
       ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
       nactdog,nacteq,nodeboun,ndirboun,&network,rhcon,
       nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,ctrl,
       xloadold,&reltime,nmethod,set,mi,istartset,iendset,ialset,nset,
       ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,iamload,
       jqrad,irowrad,&nzsrad,icolrad,ne);}
      
      if((icascade==2)||(ncont!=0)){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;
	  RENEW(nodempc,int,3*memmpcref_);
	  for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	  RENEW(coefmpc,double,memmpcref_);
	  for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
      }
      
      if((ncont!=0)&&(mortar==0)){
	  *ne=ne0;*nkon=nkon0;
	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
	       ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	       co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
	       &ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
	       &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	       iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
	       itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
	       imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
	       xnoels);
	  if(icascade<1)icascade=1;
      }
      newstep=0;
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounold,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi));
      if((icascade==2)||(ncont!=0)){
	  memmpcref_=memmpc_;mpcfreeref=mpcfree;
	  RENEW(nodempcref,int,3*memmpc_);
	  for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	  RENEW(coefmpcref,double,memmpc_);
	  for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
      }
      
      if(icascade>0) remastruct(ipompc,&coefmpc,&nodempc,nmpc,
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
      ielas=1;
      
      fn=NNEW(double,mt**nk);
      stx=NNEW(double,6*mi[0]**ne);
      
      inum=NNEW(int,*nk);
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
	  &ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
          sideload,xloadact,xloadold,&icfd,inomat);
      
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
              xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
	      tieset,istartset,iendset,ialset,ntie,&nasym));
      
      if(nmethod==0){
	  
	  /* error occurred in mafill: storing the geometry in frd format */
	  
	  ++*kode;
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	  }
	  
	  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	      kode,filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	      thicke,jobnamec,output,qfx);
	  
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
	  if(*nstate_!=0){
	      for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
		  xstateini[k]=xstate[k];
	      }
	  }	
	  if(mortar==1){
	      for (i=0;i<*ntie;i++){
		  for(j=nslavnode[i];j<nslavnode[i+1];j++){
		      islavactini[j]=islavact[j];
		      bpini[j]=bp[j];
		      for(k=0;k<3;k++){
			  cstressini[3*j+k]=cstress[3*j+k];
		      }
		  }    
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
            qfx);
	  
	  /* determining the global values to be used as boundary conditions
	     for a submodel */
	  
	  free(integerglob);free(doubleglob);
	  getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,
                   xboun,nload,sideload,iamload,&iglob);
      }
      
      if((icascade==2)||
	 ((ncont!=0)&&((iinc==1)||(ismallsliding<2)))){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;
	  RENEW(nodempc,int,3*memmpcref_);
	  for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	  RENEW(coefmpc,double,memmpcref_);
	  for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
      }
      
      if((ncont!=0)&&(mortar==0)&&((iinc==1)||(ismallsliding<2))){
	  *ne=ne0;*nkon=nkon0;
	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
	     ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	     co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
	     &ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
	     &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
	     iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
	     itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
	     imastnode,nmastnode,xmastnor,xnormastface,filab,mcs,ics,&nasym,
             xnoels);
	  
	  printf("number of contact spring elements=%d\n\n",*ne-ne0);
	  
	  if(icascade<1)icascade=1;
      }
      
      /*  updating the nonlinear mpc's (also affects the boundary
	  conditions through the nonhomogeneous part of the mpc's) */
      
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi));
      
      if((icascade==2)||
	 ((ncont!=0)&&((iinc==1)||(ismallsliding<2)))){
	  memmpcref_=memmpc_;mpcfreeref=mpcfree;
	  RENEW(nodempcref,int,3*memmpc_);
	  for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	  RENEW(coefmpcref,double,memmpc_);
	  for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
      }
      
      if(icascade>0) remastruct(ipompc,&coefmpc,&nodempc,nmpc,
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
	  inum=NNEW(int,*nk);
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
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xloadact,xloadold,&icfd,inomat);
	  iperturb[0]=0;free(inum);
	  
	  /* check whether any displacements or temperatures are changed
	     in the new increment */
	  
	  for(k=0;k<neq[1];++k){f[k]=f[k]+b[k];}
	  
      }
      else{
	  
	  inum=NNEW(int,*nk);
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
	       &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
               sideload,xloadact,xloadold,&icfd,inomat);
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
      if(mortar==1){
	  iwatchactiv=NNEW(int,2*30);
	  bdd=NNEW(double,neq[1]);
	  bhat=NNEW(double,neq[1]);
	  jqqdt=NNEW(int,neq[1]+1);
	  islavactdof=NNEW(int,neq[1]);
      }
      
      /***************************************************************/
      /* iteration counter and start of the loop over the iterations */
      /***************************************************************/

    if(mortar==1){iflagact_old=1;}
    iit=1;
    icntrl=0;
    ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
    if(uncoupled){
	*ithermal=2;
	iruc=NNEW(int,nzs[1]-nzs[0]);
	for(k=0;k<nzs[1]-nzs[0];k++) iruc[k]=irow[k+nzs[0]]-neq[0];
    }

    while(icntrl==0){

    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's) */

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
	      RENEW(nodempc,int,3*memmpcref_);
	      for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	      RENEW(coefmpc,double,memmpcref_);
	      for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
	  }

	  if((ncont!=0)&&(mortar==0)&&(ismallsliding==0)&&(iit<=8)){
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
                      &reltime,imastnode,nmastnode,xmastnor,xnormastface,
                      filab,mcs,ics,&nasym,xnoels);
	      if(*ne!=neold){iflagact=1;}

	      printf("number of contact spring elements=%d\n\n",*ne-ne0);

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
	      RENEW(nodempcref,int,3*memmpc_);
	      for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	      RENEW(coefmpcref,double,memmpc_);
	      for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
	  }

	  if(icascade>0){
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
	      
	      v=NNEW(double,mt**nk);
	      stx=NNEW(double,6*mi[0]**ne);
	      fn=NNEW(double,mt**nk);
      
	      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	      iout=-1;
	      
	      inum=NNEW(int,*nk);
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
		&reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
                sideload,xloadact,xloadold,&icfd,inomat);

	      /*for(k=0;k<neq[1];++k){printf("f=%d,%f\n",k,f[k]);}*/
	      
	      free(v);free(stx);free(fn);free(inum);
	      iout=0;
	      
	  }else{

	      /*for(k=0;k<neq[1];++k){printf("f=%d,%f\n",k,f[k]);}*/
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
                  xstateini,xstate,thicke,xnormastface,integerglob,doubleglob,
		  tieset,istartset,iendset,ialset,ntie,&nasym));

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
                  xnormastface,integerglob,doubleglob,tieset,istartset,iendset,
		  ialset,ntie,&nasym));
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
      
/*      for(k=0;k<neq[1];++k){printf("f=%d,%f\n",k,f[k]);}
      for(k=0;k<neq[1];++k){printf("fext=%d,%f\n",k,fext[k]);}
      for(k=0;k<neq[1];++k){printf("ad=%d,%f\n",k,ad[k]);}
      for(k=0;k<nzs[1];++k){printf("au=%d,%f\n",k,au[k]);}*/

      /* calculating the residual */
      if(mortar==1){
//       free(f_cm);free(f_cs);
       f_cm=NNEW(double,neq[1]);f_cs=NNEW(double,neq[1]);
      }
      calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
	 vini,&dtime,accold,nk,adb,aub,icol,irow,nzl,alpha,fextini,fini,
	 islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
	 nzs,&nasym);
    
      /*****************************************************************/
      /* if(mortar==1) start loop to find correct active set, else run */ 
      /* do loop only once */
      /*****************************************************************/
     
      if(mortar==1){
        iflagact=0;
	ismallsliding=1;
	nzsc=nzs[1];
	auc=NNEW(double,nzsc);
	adc=NNEW(double,neq[1]);
	irowc=NNEW(int,nzsc);
	icolc=NNEW(int,neq[1]);
	jqc=NNEW(int,neq[1]+1);
        if(iit==1 || ismallsliding==0){
	    aubd=NNEW(double, nzs[1]);
	    irowbd=NNEW(int, nzs[1]);
	    jqbd=NNEW(int, neq[1]+1);
	}
        if(iit==1 || ismallsliding==0){	
	free(irowb);free(Bd);free(Dd);free(jqb);
        irowb=NNEW(int,*nk);
        Bd=NNEW(double,*nk);
        nzsbd2=*nk;
        Dd=NNEW(double,*nk);
        jqb=NNEW(int, *nk+1);
	}
        /* storing the original stiffness matrix */
        jqtemp=NNEW(int,neq[1]+1);
        irowtemp=NNEW(int,nzs[1]);
        icoltemp=NNEW(int,neq[1]);
	for(i=0;i<3;i++){nzstemp[i]=nzs[i];}
	for (i=0;i<neq[1];i++){jqtemp[i]=jq[i];icoltemp[i]=icol[i];}
	jqtemp[neq[1]]=jq[neq[1]];
	for (i=0;i<nzs[1];i++){irowtemp[i]=irow[i];}

        if(iit==1 || ismallsliding==0){
	    DMEMSET(slavnor,0,3*nslavnode[*ntie],0.);
	    DMEMSET(slavtan,0,6*nslavnode[*ntie],0.);
	}
      }

      do{      
	  if(mortar==1){
	    if(iit>1 && ismallsliding==1){iflagact=1;}
	    contactmortar(&ncont,ntie,tieset,nset,set,istartset,iendset,
	       ialset,itietri,lakon,
	       ipkon,kon,koncont,ne,cg,straight,co,vold,ielmat,cs,
	       elcon,istep,&iinc,&iit,ncmat_,ntmat_,&ne0,vini,nmethod,neq,
	       nzs,nactdof,itiefac,islavsurf,islavnode,imastnode,nslavnode,
               nmastnode,&ifacecount,ad,
	       &au,b,&irow,icol,jq,imastop,iponoels,inoels,&nzsc,&auc,
	       adc,&irowc,jqc,islavact,gap,bdd,&auqdt,&irowqdt,jqqdt,&nzsqdt,
	       &nzlc,slavnor,slavtan,bhat,icolc,&aubd,&irowbd,jqbd,mi,ipe,
	       ime,tietol,&iflagact,cstress,bp,&iflag_fric,nk,nboun,
	       ndirboun,nodeboun,xbounact,nmpc,ipompc,nodempc,coefmpc,
               ikboun,ilboun,ikmpc,
	       ilmpc,nslavspc,islavspc,&nsspc,nslavmpc,islavmpc,&nsmpc,
               nmastspc,imastspc,&nmspc,
	       nmastmpc,imastmpc,&nmmpc,islavborder,pslavdual,&Bd,Dd,jqb,
               &irowb,&nzsbd2,
	       islavactdof,dhinv,islavactdoftie,plicon,nplicon,npmat_,nelcon);

	    nzs[0]=nzs[1];
	    symmetryflag=2;
	    inputformat=3;
	  }
	  
	  newstep=0;
	  
	  if(*nmethod==0){
	      
	      /* error occurred in mafill: storing the geometry in frd format */
	      
	      *nmethod=0;
	      ++*kode;
	      inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	      if(strcmp1(&filab[1044],"ZZS")==0){
		  neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	      }

	      frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
		  kode,filab,een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,
		  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
		  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
		  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
		  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
		  thicke,jobnamec,output,qfx);

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
	      
	      if(mortar==0){free(ad);free(au);} 
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
	  
	  if(mortar==1){
	      
	      /* restoring the structure of the original stiffness
		 matrix */
	      
	      for(i=0;i<3;i++){nzs[i]=nzstemp[i];}
	      for (i=0;i<neq[1];i++){jq[i]=jqtemp[i];icol[i]=icoltemp[i];}
	      jq[neq[1]]=jqtemp[neq[1]];
	      for (i=0;i<nzs[1];i++){irow[i]=irowtemp[i];}
	      free(jqtemp);free(irowtemp);free(icoltemp);
	      iflagact=iflagact_old;

	      contactstress_fric3(bhat,adc,auc,jqc,irowc,neq,gap,bdd,aubd,
		  jqbd,irowbd,b,islavact,auqdt,irowqdt,jqqdt,ntie,nslavnode,
		  islavnode,nmastnode,imastnode,slavnor,slavtan,icolc,&nzlc,
                  nactdof,&iflagact,cstress,mi,cdisp,f_cs,f_cm,&iit,
                  iwatchactiv,vold,bp,nk,nboun,ndirboun,nodeboun,xbounact,
                  nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
                  nslavspc,islavspc,&nsspc,nslavmpc,islavmpc,&nsmpc,
                  nmastspc,imastspc,&nmspc,nmastmpc,imastmpc,&nmmpc,
		  pslavdual,ipkon,kon,lakon,islavsurf,itiefac,iponoels,
                  inoels,islavactdof,dhinv,Dd,Bd,jqb,irowb,&nzsbd2,tieset,
		  elcon,tietol,ncmat_,ntmat_,plicon,nplicon,npmat_,nelcon);

	     iflagact_old=iflagact;
	    
	  }
      }while(mortar==1 && iflagact==2 && iit==1);
      
      /*****************************************************************/
      /* if(mortar==1) correct active set founded and  */ 
      /* system of equations solved, only once for penalty*/
      /*****************************************************************/
      
	/* calculating the slave contact forces 
           Ph.D. Thesis Stefan Hartmann eqn. (6.26) */

      if(mortar==1){
	  
	  free(auc);free(adc);free(irowc);free(icolc);free(jqc);
	  free(au);free(ad);free(auqdt);free(irowqdt);
	  
	  if(ismallsliding==0){
	      free(aubd);free(jqbd);free(irowbd);
	  }
      }
      
      /* calculating the displacements, stresses and forces */
      
      v=NNEW(double,mt**nk);
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      
      stx=NNEW(double,6*mi[0]**ne);
      fn=NNEW(double,mt**nk);
      
      inum=NNEW(int,*nk);
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
	 &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
         sideload,xloadact,xloadold,&icfd,inomat);
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
	 vini,&dtime,accold,nk,adb,aub,icol,irow,nzl,alpha,fextini,fini,
	 islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
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

      if(mortar==1){
	  free(f_cs);free(f_cm);
      }  
    }

    /*********************************************************/
    /*   end of the iteration loop                          */
    /*********************************************************/

    /* icutb=0 means that the iterations in the increment converged,
       icutb!=0 indicates that the increment has to be reiterated with
                another increment size (dtheta) */

    if(mortar==1){
      if(ismallsliding==1){      
	free(aubd);free(jqbd);free(irowbd);
      }
      free(bdd);free(bhat);free(jqqdt);	  
      free(iwatchactiv);free(islavactdof);
   }
   
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
      if(*nstate_!=0){
	for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
	  xstate[k]=xstateini[k];
	}	
      }

      qam[0]=qamold[0];
      qam[1]=qamold[1];
      if(mortar==1){
       for (i=0;i<*ntie;i++){
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	  islavact[j]=islavactini[j];
	  bp[j]=bpini[j];
	  for(k=0;k<3;k++){
	    cstress[3*j+k]=cstressini[3*j+k];
          }
	}    
       }        
      }      
    }
    
    if((jout[0]==jprint)&&(icutb==0)){

      jprint=0;

      /* calculating the displacements and the stresses and storing */
      /* the results in frd format  */
	
      v=NNEW(double,mt**nk);
      fn=NNEW(double,mt**nk);
      stn=NNEW(double,6**nk);
      if(*ithermal>1) qfn=NNEW(double,3**nk);
      inum=NNEW(int,*nk);
      stx=NNEW(double,6*mi[0]**ne);
      
      if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
      if(strcmp1(&filab[435],"PEEQ")==0) epn=NNEW(double,*nk);
      if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);
      if(strcmp1(&filab[609],"SDV ")==0) xstaten=NNEW(double,*nstate_**nk);
      if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);

      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);

      iout=2;
      icmd=3;
      
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
              &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat);
      
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);

      iout=0;
      if(*iexpl<=1) icmd=0;
      FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
//      for(k=0;k<ntg;k++)if(inum[itg[k]-1]>0){inum[itg[k]-1]*=-1;}
      
      ++*kode;
      if(*mcs!=0){
         if(mortar==1){
	      RENEW(stx,double,6*mi[0]*(*ne+*nslavs));
	      for(k=0;k<*nslavs;k++){
		  for(l=0;l<6;l++){stx[6*mi[0]*(*ne+k)+l]=cdisp[6*k+l];}
	      }
	      *ne+=*nslavs;
	      *nkon+=*nslavs;
	  }	
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
               nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
               ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0);
	  if(mortar==1){*ne-=*nslavs;*nkon-=*nslavs;}	       
      }
      else{
	  if(strcmp1(&filab[1044],"ZZS")==0){
	      neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	  }
	  if(mortar==1){
	      RENEW(stx,double,6*mi[0]*(*ne+*nslavs));
	      for(k=0;k<*nslavs;k++){
		  for(l=0;l<6;l++){stx[6*mi[0]*(*ne+k)+l]=cdisp[6*k+l];}
	      }
	      *ne+=*nslavs;
	  }

	  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);

	  if(mortar==1){*ne-=*nslavs;}
	  if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
      }
      
      free(v);free(fn);free(stn);free(inum);free(stx);
      if(*ithermal>1){free(qfn);}
      
      if(strcmp1(&filab[261],"E   ")==0) free(een);
      if(strcmp1(&filab[435],"PEEQ")==0) free(epn);
      if(strcmp1(&filab[522],"ENER")==0) free(enern);
      if(strcmp1(&filab[609],"SDV ")==0) free(xstaten);
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
    inum=NNEW(int,*nk);
    stx=NNEW(double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) een=NNEW(double,6**nk);
    if(strcmp1(&filab[435],"PEEQ")==0) epn=NNEW(double,*nk);
    if(strcmp1(&filab[522],"ENER")==0) enern=NNEW(double,*nk);
    if(strcmp1(&filab[609],"SDV ")==0) xstaten=NNEW(double,*nstate_**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) emn=NNEW(double,6**nk);
    
    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
    iout=2;
    icmd=3;

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
            &reltime,&ne0,xforc,nforc,thicke,xnormastface,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat);
    
    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);

    iout=0;
    if(*iexpl<=1) icmd=0;
    FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
//    for(k=0;k<ntg;k++)if(inum[itg[k]-1]>0){inum[itg[k]-1]*=-1;}
    
    ++*kode;
    if(*mcs>0){
      	if(mortar==1){
	    RENEW(stx,double,6*mi[0]*(*ne+*nslavs));
	    for(k=0;k<*nslavs;k++){
		for(l=0;l<6;l++){stx[6*mi[0]*(*ne+k)+l]=cdisp[6*k+l];}
	    }
	    *ne+=*nslavs;
	    *nkon+=*nslavs;
	}       
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	     t1act,fn,ttime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
             nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
             ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0);
	if(mortar==1){*ne-=*nslavs;*nkon-=*nslavs;} 	     
    }
    else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	    neigh=NNEW(int,40**ne);ipneigh=NNEW(int,*nk);
	}
	if(mortar==1){
	    RENEW(stx,double,6*mi[0]*(*ne+*nslavs));
	    for(k=0;k<*nslavs;k++){
		for(l=0;l<6;l++){stx[6*mi[0]*(*ne+k)+l]=cdisp[6*k+l];}
	    }
	    *ne+=*nslavs;
	}

	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx);

	if(mortar==1){*ne-=*nslavs;}
	if(strcmp1(&filab[1044],"ZZS")==0){free(ipneigh);free(neigh);}
    }
    
    free(v);free(fn);free(stn);free(inum);free(stx);
    if(*ithermal>1){free(qfn);}
    
    if(strcmp1(&filab[261],"E   ")==0) free(een);
    if(strcmp1(&filab[435],"PEEQ")==0) free(epn);
    if(strcmp1(&filab[522],"ENER")==0) free(enern);
    if(strcmp1(&filab[609],"SDV ")==0) free(xstaten);
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
      free(sideface);free(nelemface);free(ifreestream);
      free(isolidsurf);free(neighsolidsurf);free(iponoel);free(inoel);
      free(vcontu);free(inomat);
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
      RENEW(nodempc,int,3*memmpcref_);
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
      RENEW(ipkon,int,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,int,*nkon);
      if(*norien>0){
	  RENEW(ielorien,int,mi[2]**ne);
      }
      RENEW(ielmat,int,mi[2]**ne);
      free(cg);free(straight);
      free(imastop);free(itiefac);free(islavsurf);free(islavnode);
      free(nslavnode);free(iponoels);free(inoels);free(imastnode);
      free(nmastnode);free(itietri);free(koncont);free(xnoels);

    /* deleting contact MPC's (not for modal dynamics calculations) */

      remcontmpc(nmpc,labmpc,&mpcfree,nodempc,ikmpc,ilmpc,coefmpc,ipompc);

      if(mortar==1){
	free(islavact);free(gap);free(slavnor);free(slavtan);
        free(cstress);free(ipe);free(ime);
        free(cdisp);free(bp);free(islavactdoftie);
        free(nslavspc);free(islavspc);free(nslavmpc);free(islavmpc);
        free(nmastspc);free(imastspc);free(nmastmpc);free(imastmpc);
        free(islavborder); free(pslavdual);
	free(Dd); free(Bd);free(irowb);free(jqb);free(dhinv);
	free(cstressini);free(bpini);free(islavactini);
      }else{
	  free(areaslav);free(springarea);free(xmastnor);free(xnormastface);
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

  (*tmin)*=(*tper);
  (*tmax)*=(*tper);

  free(nactdofinv);
  
  return;
}
