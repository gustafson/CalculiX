/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef CALCULIX_MPI
#include <spoolesMPI.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

#ifdef CALCULIX_MPI
int myid = 0, nproc = 0;
#endif

int main(int argc,char *argv[])
{
  
int *kon=NULL, *nodeboun=NULL, *ndirboun=NULL, *ipompc=NULL,
	*nodempc=NULL, *nodeforc=NULL, *ndirforc=NULL,
	*nelemload=NULL,
	*nnn=NULL, *nactdof=NULL, *icol=NULL,*ics=NULL,
	*jq=NULL, *mast1=NULL, *irow=NULL, *rig=NULL,
	*ikmpc=NULL, *ilmpc=NULL, *ikboun=NULL, *ilboun=NULL,
	*npn=NULL, *adj=NULL, *xadj=NULL, *iw=NULL, *nreorder=NULL,
	*mmm=NULL, *xnpn=NULL, *ipointer=NULL,
	*istartset=NULL, *iendset=NULL, *ialset=NULL, *ielmat=NULL,
	*ielorien=NULL, *nrhcon=NULL, *nodebounold=NULL, *ndirbounold=NULL,
	*nelcon=NULL, *nalcon=NULL, *iamforc=NULL,  *iamload=NULL,
	*iamt1=NULL, *namta=NULL, *ipkon=NULL, *iamboun=NULL,
	*nplicon=NULL, *nplkcon=NULL, *inotr=NULL, *iponor=NULL, *knor=NULL,
	*ikforc=NULL, *ilforc=NULL, *iponoel=NULL, *inoel=NULL, *nshcon=NULL,
	*ncocon=NULL,*ibody=NULL, *inum1=NULL,*ielprop=NULL,
    *inum2=NULL,*ipoinpc=NULL,cfd=0,mt;
    
double *co=NULL, *xboun=NULL, *coefmpc=NULL, *xforc=NULL,
	*xload=NULL, *ad=NULL, *au=NULL, *xbounold=NULL, *xforcold=NULL,
	*b=NULL, *vold=NULL, *sti=NULL, *xloadold=NULL, *xnor=NULL,
	*reorder=NULL,*dcs=NULL, *thickn=NULL, *thicke=NULL, *offset=NULL,
	*elcon=NULL, *rhcon=NULL, *alcon=NULL, *alzero=NULL, *t0=NULL, *t1=NULL,
	*prestr=NULL, *orab=NULL, *amta=NULL, *veold=NULL, *accold=NULL,
	*adb=NULL, *aub=NULL, *t1old=NULL, *eme=NULL, *plicon=NULL, *plkcon=NULL,
	*xstate=NULL, *trab=NULL, *ener=NULL, *shcon=NULL, *cocon=NULL,
	*cs=NULL,*tietol=NULL,*fmpc=NULL,*prop=NULL,
	*xbody=NULL,*xbodyold=NULL;
    
double ctrl[27]={4.5,8.5,9.5,16.5,10.5,4.5,0.,5.5,0.,0.,0.25,0.5,0.75,0.85,0.,0.,1.5,0.,0.005,0.01,0.,0.,0.02,1.e-5,1.e-3,1.e-8,1.e30};
    
char *sideload=NULL, *set=NULL, *matname=NULL, *orname=NULL, *amname=NULL,
	*filab=NULL, *lakon=NULL, *labmpc=NULL, *prlab=NULL, *prset=NULL, 
	jobnamec[396]="",jobnamef[132]="",output[4]="frd", *typeboun=NULL,
	*inpc=NULL,*tieset=NULL,*cbody=NULL;
    
int nk,ne,nboun,nmpc,nforc,nload,nprint,nset,nalset,nentries=14,
  nmethod,neq[3]={0,0,0},i,mpcfree=1,mei[4],j,nzl,nam,nbounold=0,
  nforcold=0,nloadold=0,nbody,nbody_=0,nbodyold=0,
  k,nzs[3],nmpc_=0,nload_=0,nforc_=0,istep,istat,nboun_=0,
  iperturb[2]={0,0},nmat,ntmat_=0,norien,ithermal[2]={0,0},
  iprestr,kode,isolver=0,
  jout[2]={1,1},nlabel,nkon=0,idrct,jmax[2],iexpl,nevtot=0,
//  iplas=0,npmat_=0,mi[2]={0,0},ntrans,mpcend=-1,namtot_=0,iumat=0,mpcmult,
  iplas=0,npmat_=0,mi[2]={0,3},ntrans,mpcend=-1,namtot_=0,iumat=0,mpcmult,
  icascade=0,maxlenmpc,mpcinfo[4],ne1d=0,ne2d=0,infree[4]={0,0,0,0},
  callfrommain,nflow=0,jin=0,irstrt=0,nener=0,jrstrt=0,nenerold,
  nline,ipoinp[2*nentries],*inp=NULL,ntie,ntie_=0,mcs=0,kflag=2,nprop_=0,
  nprop=0,itpamp=0,iviewfile,nkold,nevdamp_=0;

int *meminset=NULL,*rmeminset=NULL;

int nzs_,nk_=0,ne_=0,nset_=0,nalset_=0,nmat_=0,norien_=0,nam_=0,
    ntrans_=0,ncs_=0,nstate_=0,ncmat_=0,memmpc_=0,nprint_=0;

double fei[3],tinc,tper,tmin,tmax,*xmodal=NULL,
    alpha,ttime=0.,qaold[2]={0.,0.},physcon[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};


#ifdef CALCULIX_MPI
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
#endif


if(argc==1){printf("Usage: CalculiX.exe -i jobname\n");FORTRAN(stop,());}
else{
  for(i=1;i<argc;i++){
    if(strcmp1(argv[i],"-i")==0) {
    strcpy(jobnamec,argv[i+1]);strcpy1(jobnamef,argv[i+1],132);jin++;break;}
    if(strcmp1(argv[i],"-v")==0) {
	printf("\nThis is version version 2.1\n\n");
	FORTRAN(stop,());
    }
  }
  if(jin==0){strcpy(jobnamec,argv[1]);strcpy1(jobnamef,argv[1],132);}

  for(i=1;i<argc;i++){
    if(strcmp1(argv[i],"-o")==0) {
    strcpy(output,argv[i+1]);break;}
  }
}

FORTRAN(openfile,(jobnamef,output));

printf("\n************************************************************\n\n");
printf("CalculiX version 2.1, Copyright(C) 1998-2007 Guido Dhondt\n");
printf("CalculiX comes with ABSOLUTELY NO WARRANTY. This is free\n");
printf("software, and you are welcome to redistribute it under\n");
printf("certain conditions, see gpl.htm\n\n");
printf("************************************************************\n\n");
printf("You are using an executable made on Do 4. Mär 20:53:32 CET 2010\n");
fflush(stdout);

istep=0;
istat=0;
iprestr=0;
kode=0;

/* default solver */

#if defined(SGI)
 isolver=4;
#elif defined(PARDISO)
 isolver=7;
#elif defined(SPOOLES)
 isolver=0;
#elif defined(TAUCS)
 isolver=5;
#else
 isolver=3;
#endif

/* conservative estimate of the fields to be allocated */

readinput(jobnamec,&inpc,&nline,&nset_,ipoinp,&inp,&ipoinpc,ithermal); 

set=NNEW(char,81*nset_);
meminset=NNEW(int,nset_);
rmeminset=NNEW(int,nset_);

FORTRAN(allocation,(&nload_,&nforc_,&nboun_,&nk_,&ne_,&nmpc_,&nset_,&nalset_,
   &nmat_,&ntmat_,&npmat_,&norien_,&nam_,&nprint_,mi,&ntrans_,
   set,meminset,rmeminset,&ncs_,&namtot_,&ncmat_,&memmpc_,&ne1d,
   &ne2d,&nflow,jobnamec,&irstrt,ithermal,&nener,&nstate_,&istep,
   inpc,ipoinp,inp,&ntie_,&nbody_,&nprop_,ipoinpc,&nevdamp_));

free(set);free(meminset);free(rmeminset);mt=mi[1]+1;

nzs_=20000000;

nload=0;nbody=0;nforc=0;nboun=0;nk=0;nmpc=0;nam=0;

/* caveat: change nlabel in radmatrix.f and expand.c as well 
   if changing next line, as well as in storeresidual.f */
nlabel=27;

while(istat>=0) {

  fflush(stdout);

  /* in order to reduce the number of variables to be transferred to
     the subroutines, the max. field sizes are (for most fields) copied
     into the real sizes */

  nzs[1]=nzs_;
  nprint=nprint_;

  if((istep == 0)||(irstrt<0)) {
    ne=ne_;
    nset=nset_;
    nalset=nalset_;
    nmat=nmat_;
    norien=norien_;
    ntrans=ntrans_;
    ntie=ntie_;

    /* allocating space before the first step */

    /* coordinates and topology */

    co=NNEW(double,3*nk_);
    kon=NNEW(int,28*ne_);
    ipkon=NNEW(int,ne_);
    lakon=NNEW(char,8*ne_);

    /* property cards */

    ielprop=NNEW(int,ne_);
    for(i=0;i<ne_;i++) ielprop[i]=-1;
    prop=NNEW(double,nprop_);

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0)){
	iponor=NNEW(int,56*ne_);
	for(i=0;i<56*ne_;i++) iponor[i]=-1;
	xnor=NNEW(double,36*ne1d+24*ne2d);
	knor=NNEW(int,24*(ne1d+ne2d));
	thickn=NNEW(double,2*nk_);
	thicke=NNEW(double,56*ne_);
	offset=NNEW(double,2*ne_);
	iponoel=NNEW(int,nk_);
	inoel=NNEW(int,9*ne1d+24*ne2d);
	rig=NNEW(int,nk_);
	infree[2]=1;
    }

    /* SPC's */

    nodeboun=NNEW(int,nboun_);
    ndirboun=NNEW(int,nboun_);
    typeboun=NNEW(char,nboun_);
    iamboun=NNEW(int,nboun_);
    xboun=NNEW(double,nboun_);
    ikboun=NNEW(int,nboun_);
    ilboun=NNEW(int,nboun_);

    /* MPC's */

    ipompc=NNEW(int,nmpc_);
    nodempc=NNEW(int,3*memmpc_);
    for(i=0;i<3*memmpc_;i+=3){nodempc[i+2]=i/3+2;}
    nodempc[3*memmpc_-1]=0;
    coefmpc=NNEW(double,memmpc_);
    labmpc=NNEW(char,20*nmpc_+1);
    ikmpc=NNEW(int,nmpc_);
    ilmpc=NNEW(int,nmpc_);
    fmpc=NNEW(double,nmpc_);

    /* nodal loads */

    nodeforc=NNEW(int,2*nforc_);
    ndirforc=NNEW(int,nforc_);
    iamforc=NNEW(int,nforc_);
    xforc=NNEW(double,nforc_);
    ikforc=NNEW(int,nforc_);
    ilforc=NNEW(int,nforc_);

    /* distributed facial loads */

    nelemload=NNEW(int,2*nload_);
    iamload=NNEW(int,2*nload_);
    sideload=NNEW(char,20*nload_);
    xload=NNEW(double,2*nload_);

    /* distributed volumetric loads */

    cbody=NNEW(char,81*nbody_);
    ibody=NNEW(int,3*nbody_);
    xbody=NNEW(double,7*nbody_);
    xbodyold=NNEW(double,7*nbody_);

    /* printing output */

    prlab=NNEW(char,6*nprint_);
    prset=NNEW(char,81*nprint_);

    /* set definitions */

    set=NNEW(char,81*nset);
    istartset=NNEW(int,nset);
    iendset=NNEW(int,nset);
    ialset=NNEW(int,nalset);

    /* (hyper)elastic constants */

    elcon=NNEW(double,(ncmat_+1)*ntmat_*nmat);
    nelcon=NNEW(int,2*nmat);

    /* density */

    rhcon=NNEW(double,2*ntmat_*nmat);
    nrhcon=NNEW(int,nmat);

    /* specific heat */

    shcon=NNEW(double,4*ntmat_*nmat);
    nshcon=NNEW(int,nmat);

    /* thermal expansion coefficients */

    alcon=NNEW(double,7*ntmat_*nmat);
    nalcon=NNEW(int,2*nmat);
    alzero=NNEW(double,nmat);

    /* conductivity */

    cocon=NNEW(double,7*ntmat_*nmat);
    ncocon=NNEW(int,2*nmat);

    /* isotropic and kinematic hardening coefficients*/

    plicon=NNEW(double,(2*npmat_+1)*ntmat_*nmat);
    nplicon=NNEW(int,(ntmat_+1)*nmat);
    plkcon=NNEW(double,(2*npmat_+1)*ntmat_*nmat);
    nplkcon=NNEW(int,(ntmat_+1)*nmat);

    /* linear dynamic properties */
    
    xmodal=NNEW(double,11+nevdamp_);
    xmodal[10]=nevdamp_+0.5;

    /* internal state variables */

    if(nstate_>0){xstate=NNEW(double,nstate_*mi[0]*ne);}

    /* material orientation */

    orname=NNEW(char,80*norien);
    orab=NNEW(double,7*norien);
    ielorien=NNEW(int,ne_);

    /* transformations */

    trab=NNEW(double,7*ntrans);
    inotr=NNEW(int,2*nk_);

    /* amplitude definitions */

    amname=NNEW(char,80*nam_);
    amta=NNEW(double,2*namtot_);
    namta=NNEW(int,3*nam_);

    /* temperatures */

    if((ne1d==0)&&(ne2d==0)){
	t0=NNEW(double,nk_);
	t1=NNEW(double,nk_);}
    else{
	t0=NNEW(double,3*nk_);
	t1=NNEW(double,3*nk_);}
    iamt1=NNEW(int,nk_);

    prestr=NNEW(double,6*mi[0]*ne_);
    vold=NNEW(double,mt*nk_);
    veold=NNEW(double,mt*nk_);

    ielmat=NNEW(int,ne_);

    matname=NNEW(char,80*nmat);

    filab=NNEW(char,87*nlabel);

    /* tied constraints */

    if(ntie_>0){
      tieset=NNEW(char,243*ntie_);
      tietol=NNEW(double,ntie_);
      cs=NNEW(double,17*ntie_);
    }

    /* temporary fields for cyclic symmetry calculations */

    if(ncs_>0){
      ics=NNEW(int,24*ncs_);
      dcs=NNEW(double,30*ncs_);
    }

  }
  else {

    /* allocating and reallocating space for subsequent steps */

    if((nmethod != 4) && ((nmethod != 1) || (iperturb[0] < 2))){
      veold=NNEW(double,mt*nk_);
    }
    else{
      RENEW(veold,double,mt*nk_);
      memset(&veold[mt*nk],0,sizeof(double)*mt*(nk_-nk));
    }
    RENEW(vold,double,mt*nk_);
    memset(&vold[mt*nk],0,sizeof(double)*mt*(nk_-nk));

 /*   if(nmethod != 4){free(accold);}*/

    RENEW(nodeboun,int,nboun_);
    RENEW(ndirboun,int,nboun_);
    RENEW(typeboun,char,nboun_);
    RENEW(xboun,double,nboun_);
    RENEW(ikboun,int,nboun_);
    RENEW(ilboun,int,nboun_);

    RENEW(nodeforc,int,2*nforc_);
    RENEW(ndirforc,int,nforc_);
    RENEW(xforc,double,nforc_);
    RENEW(ikforc,int,nforc_);
    RENEW(ilforc,int,nforc_);

    RENEW(nelemload,int,2*nload_);
    RENEW(sideload,char,20*nload_);
    RENEW(xload,double,2*nload_);

    RENEW(cbody,char,81*nbody_);
    RENEW(ibody,int,3*nbody_);
    RENEW(xbody,double,7*nbody_);
    RENEW(xbodyold,double,7*nbody_);
    for(i=7*nbodyold;i<7*nbody_;i++) xbodyold[i]=0;

    if(nam > 0) {
      RENEW(iamforc,int,nforc_);
      RENEW(iamload,int,2*nload_);
      RENEW(iamboun,int,nboun_);
      RENEW(amname,char,80*nam_);
      RENEW(amta,double,2*namtot_);
      RENEW(namta,int,3*nam_);
    }

    RENEW(ipompc,int,nmpc_);

    RENEW(labmpc,char,20*nmpc_+1);
    RENEW(ikmpc,int,nmpc_);
    RENEW(ilmpc,int,nmpc_);
    RENEW(fmpc,double,nmpc_);

    if(ntrans > 0){
      RENEW(inotr,int,2*nk_);
    }

    RENEW(co,double,3*nk_);

    if(ithermal[0] != 0){
	if((ne1d==0)&&(ne2d==0)){
	    RENEW(t0,double,nk_);
	    RENEW(t1,double,nk_);
	}
      if(nam > 0) {RENEW(iamt1,int,nk_);}
    }

  }

  /* allocation of fields in the restart file */

  if(irstrt<0){
    nodebounold=NNEW(int,nboun_);
    ndirbounold=NNEW(int,nboun_);
    xbounold=NNEW(double,nboun_);
    xforcold=NNEW(double,nforc_);
    xloadold=NNEW(double,2*nload_);
    if(ithermal[0]!=0) t1old=NNEW(double,nk_); 
    sti=NNEW(double,6*mi[0]*ne);
    eme=NNEW(double,6*mi[0]*ne);
    if(nener==1)ener=NNEW(double,mi[0]*ne*2);
    nnn=NNEW(int,nk_);
  }

  nenerold=nener;
  nkold=nk;

  /* reading the input file */

  FORTRAN(calinput,(co,&nk,kon,ipkon,lakon,&nkon,&ne,
            nodeboun,ndirboun,xboun,&nboun,
	    ipompc,nodempc,coefmpc,&nmpc,&nmpc_,nodeforc,ndirforc,xforc,&nforc,
	    &nforc_,nelemload,sideload,xload,&nload,&nload_,
	    &nprint,prlab,prset,&mpcfree,&nboun_,mei,set,istartset,iendset,
	    ialset,&nset,&nalset,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
	    alzero,t0,t1,matname,ielmat,orname,orab,ielorien,amname,
            amta,namta,&nam,&nmethod,iamforc,iamload,iamt1,
	    ithermal,iperturb,&istat,&istep,&nmat,&ntmat_,&norien,prestr,
	    &iprestr,&isolver,fei,veold,&tinc,&tper,
	    xmodal,filab,jout,&nlabel,&idrct,
	    jmax,&tmin,&tmax,&iexpl,&alpha,iamboun,plicon,nplicon,
	    plkcon,nplkcon,&iplas,&npmat_,mi,&nk_,trab,inotr,&ntrans,
	    ikboun,ilboun,ikmpc,ilmpc,ics,dcs,&ncs_,&namtot_,cs,&nstate_,
	    &ncmat_,&iumat,&mcs,labmpc,iponor,xnor,knor,thickn,thicke,
	    ikforc,ilforc,offset,iponoel,inoel,rig,infree,nshcon,shcon,
            cocon,ncocon,physcon,&nflow,
            ctrl,&memmpc_,&maxlenmpc,&ne1d,&ne2d,&nener,vold,nodebounold,
            ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
            sti,ener,xstate,jobnamec,nnn,&irstrt,&ttime,
            qaold,output,typeboun,inpc,&nline,ipoinp,inp,tieset,tietol,
            &ntie,fmpc,cbody,ibody,xbody,&nbody,&nbody_,xbodyold,&nam_,
	    ielprop,&nprop,&nprop_,prop,&itpamp,&iviewfile,ipoinpc,&cfd));

/*	FORTRAN(writeboun,(nodeboun,ndirboun,xboun,typeboun,&nboun));*/

  if(istat<0) break;

  /*RENEW(inpc,char,(long long)132*nline);*/
  /* RENEW(inp,int,3*ipoinp[23]); */

  if(istep == 1) {

    /* reallocating space in the first step */

    /* allocating and initializing fields pointing to the previous step */

    RENEW(vold,double,mt*nk);
    sti=NNEW(double,6*mi[0]*ne);

    /* strains */

    eme=NNEW(double,6*mi[0]*ne);

    /* residual stresses/strains */

    if(iprestr==1) {
	RENEW(prestr,double,6*mi[0]*ne);
	for(i=0;i<ne;i++){
	    for(j=0;j<mi[0];j++){
		for(k=0;k<6;k++){
		    sti[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		}
	    }
	}
    }
    else if(iprestr==2){
	RENEW(prestr,double,6*mi[0]*ne);
	for(i=0;i<ne;i++){
	    for(j=0;j<mi[0];j++){
		for(k=0;k<6;k++){
		    eme[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		}
	    }
	}
    }
    else {
	free(prestr);
    }

    nodebounold=NNEW(int,nboun);
    ndirbounold=NNEW(int,nboun);
    xbounold=NNEW(double,nboun);
    xforcold=NNEW(double,nforc);
    xloadold=NNEW(double,2*nload);

    /* initial temperatures: store in the "old" boundary conditions */

    if(ithermal[0]>1){
      for(i=0;i<nboun;i++){
	if(ndirboun[i]==0){
	    xbounold[i]=vold[mt*(nodeboun[i]-1)];
	}
      }
    }

    /* initial temperatures: store in the "old" temperature field */

    if(ithermal[0]!=0){
      t1old=NNEW(double,nk);
      for(i=0;i<nk;i++) t1old[i]=t0[i];
    }

    /* element definition */

    RENEW(kon,int,nkon);
    RENEW(ipkon,int,ne);
    RENEW(lakon,char,8*ne);

    /* property cards */

    if(nprop>0){
	RENEW(ielprop,int,ne);
	RENEW(prop,double,nprop);
    }else{
	free(ielprop);free(prop);
    }

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0)){
	RENEW(iponor,int,2*nkon);
	RENEW(xnor,double,infree[0]);
	RENEW(knor,int,infree[1]);
	free(thickn);
	RENEW(thicke,double,2*nkon);
	RENEW(offset,double,2*ne);
	RENEW(inoel,int,3*(infree[2]-1));
	RENEW(iponoel,int,infree[3]);
	RENEW(rig,int,infree[3]);
    }

    /* set definitions */ 

    RENEW(set,char,81*nset);
    RENEW(istartset,int,nset);
    RENEW(iendset,int,nset);
    RENEW(ialset,int,nalset);

    /* material properties */

    RENEW(elcon,double,(ncmat_+1)*ntmat_*nmat);
    RENEW(nelcon,int,2*nmat);

    RENEW(rhcon,double,2*ntmat_*nmat);
    RENEW(nrhcon,int,nmat);

    RENEW(shcon,double,4*ntmat_*nmat);
    RENEW(nshcon,int,nmat);

    RENEW(cocon,double,7*ntmat_*nmat);
    RENEW(ncocon,int,2*nmat);

    RENEW(alcon,double,7*ntmat_*nmat);
    RENEW(nalcon,int,2*nmat);
    RENEW(alzero,double,nmat);

    RENEW(matname,char,80*nmat);
    RENEW(ielmat,int,ne);

    /* allocating space for the state variables */

    /*    if(nstate_>0){
      xstate=NNEW(double,nstate_*mi[0]*ne);
      }*/

    /* next statements for plastic materials and nonlinear springs */

    if(npmat_>0){
	RENEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	RENEW(nplicon,int,(ntmat_+1)*nmat);
    }else{
	free(plicon);free(nplicon);
    }
    /* next statements only for plastic materials */

    if(iplas!=0){
      RENEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
      RENEW(nplkcon,int,(ntmat_+1)*nmat);
    }
    else{
      free(plkcon);free(nplkcon);
    }

    /* material orientation */

    if(norien > 0) {
      RENEW(orname,char,80*norien);
      RENEW(ielorien,int,ne);
      RENEW(orab,double,7*norien);
    }
    else {
      free(orname);
      free(ielorien);
      free(orab);
    }

    /* amplitude definitions */

    if(nam > 0) {
      RENEW(amname,char,80*nam);
      RENEW(namta,int,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }
    else {
      free(amname);
      free(amta);
      free(namta);
      free(iamforc);
      free(iamload);
      free(iamboun);
    }

    if(ntrans > 0){
      RENEW(trab,double,7*ntrans);
    }
    else{free(trab);free(inotr);}

    if(ithermal[0] == 0){free(t0);free(t1);}
    if((ithermal[0] == 0)||(nam<=0)){free(iamt1);}

    if(ncs_>0){
      RENEW(ics,int,ncs_);
      free(dcs);}

 
  /* tied contact constraints: generate appropriate MPC's */

  tiedcontact(&ntie, tieset, &nset, set,istartset, iendset, ialset,
       lakon, ipkon, kon,tietol,&nmpc, &mpcfree, &memmpc_,
       &ipompc, &labmpc, &ikmpc, &ilmpc,&fmpc, &nodempc, &coefmpc,
       ithermal, co, vold,&cfd,&nmpc_,mi);

 }else{

    /* reallocating space in all but the first step (>1) */

    RENEW(vold,double,mt*nk);

    /* if the SPC boundary conditions were changed in the present step,
       they have to be rematched with those in the last step. Removed SPC 
       boundary conditions do not appear any more (this is different from
       forces and loads, where removed forces or loads are reset to zero;
       a removed SPC constraint does not have a numerical value any more) */
       
    reorder=NNEW(double,nboun);
    nreorder=NNEW(int,nboun);
    if(nbounold<nboun){
      RENEW(xbounold,double,nboun);
      RENEW(nodebounold,int,nboun);
      RENEW(ndirbounold,int,nboun);
    }
    FORTRAN(spcmatch,(xboun,nodeboun,ndirboun,&nboun,xbounold,nodebounold,
		      ndirbounold,&nbounold,ikboun,ilboun,vold,reorder,nreorder,
                      mi));
    RENEW(xbounold,double,nboun);
    RENEW(nodebounold,int,nboun);
    RENEW(ndirbounold,int,nboun);
    free(reorder); free(nreorder);

    /* for additional forces or loads in the present step, the
       corresponding slots in the force and load fields of the
       previous steps are initialized */

    RENEW(xforcold,double,nforc);
    for(i=nforcold;i<nforc;i++) xforcold[i]=0;

    RENEW(xloadold,double,2*nload);
    for(i=2*nloadold;i<2*nload;i++) xloadold[i]=0;

    if(ithermal[0]!=0){
      RENEW(t1old,double,nk);
    }

    if(nam > 0) {
      RENEW(amname,char,80*nam);
      RENEW(namta,int,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }
    
  }

  /* reallocating fields for all steps (>=1) */

  RENEW(co,double,3*nk);

  RENEW(nodeboun,int,nboun);
  RENEW(ndirboun,int,nboun);
  RENEW(typeboun,char,nboun);
  RENEW(xboun,double,nboun);
  RENEW(ikboun,int,nboun);
  RENEW(ilboun,int,nboun);
    
  RENEW(nodeforc,int,2*nforc);
  RENEW(ndirforc,int,nforc);
  RENEW(xforc,double,nforc);
  RENEW(ikforc,int,nforc);
  RENEW(ilforc,int,nforc);

  RENEW(nelemload,int,2*nload);
  RENEW(sideload,char,20*nload);
  RENEW(xload,double,2*nload);

  RENEW(cbody,char,81*nbody);
  RENEW(ibody,int,3*nbody);
  RENEW(xbody,double,7*nbody);
  RENEW(xbodyold,double,7*nbody);

  RENEW(ipompc,int,nmpc);
  RENEW(labmpc,char,20*nmpc+1);
  RENEW(ikmpc,int,nmpc);
  RENEW(ilmpc,int,nmpc);
  RENEW(fmpc,double,nmpc);

  /* energy */

  if((nener==1)&&(nenerold==0)){
    ener=NNEW(double,mi[0]*ne*2);
    if((istep>1)&&(iperturb[0]>1)){
      printf("*ERROR in CalculiX: in nonlinear calculations");
      printf("       energy output must be selected in the first step");
      FORTRAN(stop,());
    }
  }

  /* initial velocities and accelerations */

  if((nmethod == 4) || ((nmethod == 1) && (iperturb[0] >= 2))) {
    RENEW(veold,double,mt*nk);
  }
  else {free(veold);}

  if((nmethod == 4)&&(iperturb[0]>1)) {
    accold=NNEW(double,mt*nk);
    }

  if(nam > 0) {
    RENEW(iamforc,int,nforc);
    RENEW(iamload,int,2*nload);
    RENEW(iamboun,int,nboun);
  }

  /* temperature loading */
  
  if(ithermal[0] != 0){
      if((ne1d==0)&&(ne2d==0)){
	  RENEW(t0,double,nk);
	  RENEW(t1,double,nk);
      }
    if(nam > 0) {RENEW(iamt1,int,nk);}
  }

  if(ntrans > 0){
    RENEW(inotr,int,2*nk);
  }

  /*  sorting the elements with distributed loads */

  if(nload>0){
      if(nam>0){
	  FORTRAN(isortiddc2,(nelemload,iamload,xload,xloadold,sideload,&nload,&kflag));
      }else{
	  FORTRAN(isortiddc1,(nelemload,xload,xloadold,sideload,&nload,&kflag));
      }
  }
  
  /*   calling the user routine ufaceload (can be empty) */

  FORTRAN(ufaceload,(co,ipkon,kon,lakon,nelemload,sideload,&nload));
  
  /* decascading MPC's and renumbering the equations: only necessary
  if MPC's changed */

  if(((istep == 1)||(ntrans>0)||(mpcend<0)||(nk!=nkold))&&(icascade==0)) {

    /* decascading the MPC's */

    printf(" Decascading the MPC's\n\n");

    callfrommain=1;
    cascade(ipompc,&coefmpc,&nodempc,&nmpc,
	    &mpcfree,nodeboun,ndirboun,&nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,&mpcmult,
	    labmpc,&nk,&memmpc_,&icascade,&maxlenmpc,
            &callfrommain,iperturb,ithermal);

    if(istep==1) nnn=NNEW(int,nk);
    else RENEW(nnn,int,nk);
    for(i=1;i<=nk;++i)
	nnn[i-1]=i;
	
//    if((icascade==0)&&(isolver!=6)){
    if((icascade==10)&&(isolver!=6)){

	/* renumbering the nodes */
	
	printf(" Renumbering the nodes to decrease the profile:\n");
	fflush(stdout);

	npn=NNEW(int,20*ne+mpcend);
	adj=NNEW(int,380*ne+mpcmult);
	xadj=NNEW(int,nk+1);
	iw=NNEW(int,3*nk+1);
	mmm=NNEW(int,nk);
	xnpn=NNEW(int,ne+nmpc+1);
	inum1=NNEW(int,nk);
	inum2=NNEW(int,nk);
	
	FORTRAN(renumber,(&nk,kon,ipkon,lakon,&ne,ipompc,nodempc,&nmpc,nnn,
	npn,adj,xadj,iw,mmm,xnpn,inum1,inum2));
	
	free(npn);free(adj);free(xadj);free(iw);free(mmm);free(xnpn);
	free(inum1);free(inum2);
    }

  }

  /* determining the matrix structure: changes if SPC's have changed */
  
  if(icascade==0) printf(" Determining the structure of the matrix:\n");
  
  nactdof=NNEW(int,mt*nk);  
  mast1=NNEW(int,nzs[1]);
  irow=NNEW(int,nzs[1]);
  
  if((mcs==0)||(cs[1]<0)){
      
      icol=NNEW(int,4*nk);
      jq=NNEW(int,4*nk+1);
      ipointer=NNEW(int,4*nk);
      
      if(icascade==0){
	  mastruct(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,ipompc,
		   nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,neq,nnn,
		   ikmpc,ilmpc,ipointer,nzs,&nmethod,ithermal,
                   ikboun,ilboun,iperturb,mi);
      }
      else{neq[0]=1;neq[1]=1;neq[2]=1;}
  }
  else{
      
      icol=NNEW(int,8*nk);
      jq=NNEW(int,8*nk+1);
      ipointer=NNEW(int,8*nk);
      
      mastructcs(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,
		 ipompc,nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,
		 neq,nnn,ikmpc,ilmpc,ipointer,nzs,&nmethod,
		 ics,cs,labmpc,&mcs,mi);
  }
  
  free(ipointer);free(mast1);
  if(icascade==0)RENEW(irow,int,nzs[2]);

  /* nmethod=1: static analysis   */
  /* nmethod=2: frequency analysis  */
  /* nmethod=3: buckling analysis */
  /* nmethod=4: linear dynamic analysis */

  if((nmethod<=1)||(iperturb[0]>1))
    {
	if(iperturb[0]<2){
	
	prespooles(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     ad,au,b,nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,adb,aub,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,&tper);

      }

      else{

	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	nonlingeo(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     ad,au,b,nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,&ikmpc, 
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     &vold,iperturb,sti,nzs,&kode,adb,aub,filab,&idrct,jmax,
	     jout,&tinc,&tper,&tmin,&tmax,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,iamload,iamt1,&alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,nnn,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,tieset,&itpamp,&iviewfile,jobnamec,tietol);

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];


      }
    }
  else if(nmethod==2)
    {
      /* FREQUENCY ANALYSIS */

      if((mcs==0)||(cs[1]<0)){
#ifdef ARPACK
	arpack(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     ad,au,b,nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     shcon,nshcon,cocon,ncocon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,vold,iperturb,sti,nzs,
	     &kode,adb,aub,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,&ncmat_,&nstate_,ener,jobnamec,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody);}
#else
            printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
            FORTRAN(stop,());}
#endif

      else{
#ifdef ARPACK
	arpackcs(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     ad,au,b,nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,adb,aub,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,ics,cs,&mpcend,&ncmat_,
             &nstate_,&mcs,&nkon,ener,jobnamec,output,set,&nset,istartset,
             iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody,&nevtot);}
#else
            printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
            FORTRAN(stop,());}
#endif

    }
  else if(nmethod==3)
    {
#ifdef ARPACK
	arpackbu(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, 
	     nelemload,sideload,xload,&nload, 
	     ad,au,b,nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,adb,aub,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,&ncmat_,&nstate_,ener,output,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody);
#else
            printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
  else if(nmethod==4)
    {
	if((ne1d!=0)||(ne2d!=0)){
	    printf(" *WARNING: 1-D or 2-D elements may cause problems in modal dynamic calculations\n");
	    printf("           ensure that point loads defined in a *MODAL DYNAMIC step\n");
	    printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	    printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	    printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");}

      printf(" Composing the dynamic response from the eigenmodes\n\n");

      dyna(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,&sti,nzs,
	    &tinc,&tper,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	   &ntie,tieset,&idrct,jmax,&tmin,&tmax,ctrl,&itpamp,tietol,&nalset,&nnn);
    }
  else if(nmethod==5)
    {
	  if((ne1d!=0)||(ne2d!=0)){
	      printf(" *WARNING: 1-D or 2-D elements may cause problems in steady state calculations\n");
	      printf("           ensure that point loads defined in a *STEADY STATE DYNAMICS step\n");
	      printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	      printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	      printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");}

      printf(" Composing the steady state response from the eigenmodes\n\n");

      steadystate(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,sti,nzs,
	    &tinc,&tper,xmodal,veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,ics,cs,&mpcend,&nnn);
    }

  free(nactdof);
  free(icol);
  free(jq);
  free(irow);

  /* deleting the perturbation loads and temperatures */

  if((iperturb[0] == 1)&&(nmethod==3)) {
    nforc=0;
    nload=0;
    nbody=0;
    if(ithermal[0] == 1) {
      for(k=0;k<nk;++k){
	t1[k]=t0[k];
      }
    }
  }

  else{
    nbounold=nboun;
    for (i=0;i<nboun;i++) {
	nodebounold[i]=nodeboun[i];
	ndirbounold[i]=ndirboun[i];
    }
    nforcold=nforc;
    nloadold=nload;
    nbodyold=nbody;

    /* resetting the amplitude to none except for time=total time amplitudes */

    if(nam > 0) {
	for (i=0;i<nboun;i++) {
	    if(iamboun[i]>0){
		if(namta[3*iamboun[i]-1]>0){
		    iamboun[i]=0;
		    xboun[i]=xbounold[i];}
	    }
	}
	for (i=0;i<nforc;i++){
	    if(iamforc[i]>0){
		if(namta[3*iamforc[i]-1]>0){
		    iamforc[i]=0;
		    xforc[i]=xforcold[i];}
	    }
	}
	for (i=0;i<2*nload;i++){
	    if(iamload[i]>0){
		if(namta[3*iamload[i]-1]>0){
		    iamload[i]=0;
		    xload[i]=xloadold[i];}
	    }
	}
	for (i=1;i<3*nbody;i=i+3){
	    if(ibody[i]>0){
		if(namta[3*ibody[i]-1]>0){
		    ibody[i]=0;
		    xbody[7*(i-1)/3]=xbodyold[7*(i-1)/3];}
	    }
	}
	if(ithermal[0]==1) {
	    if(iamt1[i]>0){
		if(namta[3*iamt1[i]-1]>0){
		    iamt1[i]=0;
		    t1[i]=t1old[i];}
	    }
	}
    }
  }


  if((nmethod == 4)&&(iperturb[0]>1)) free(accold);

  if(irstrt>0){
    jrstrt++;
    if(jrstrt==irstrt){
      jrstrt=0;
      FORTRAN(restartwrite,(&istep, &nset, &nload, &nforc, &nboun, &nk, &ne, 
        &nmpc, &nalset, &nmat, &ntmat_, &npmat_, &norien, &nam, &nprint,  
        mi, &ntrans, &ncs_, &namtot_, &ncmat_, &mpcend,&maxlenmpc, &ne1d, 
        &ne2d, &nflow, &nlabel, &iplas, &nkon,ithermal,&nmethod,iperturb, 
        &nstate_,&nener, set, istartset, iendset, ialset, co, kon, ipkon, 
        lakon, nodeboun, ndirboun, iamboun, xboun, ikboun, ilboun, ipompc, 
        nodempc, coefmpc, labmpc, ikmpc, ilmpc, nodeforc, ndirforc, iamforc, 
        xforc, ikforc, ilforc, nelemload, iamload, sideload, xload, 
        elcon, nelcon, rhcon, nrhcon, alcon, nalcon, 
        alzero, plicon, nplicon, plkcon, nplkcon, orname, orab, ielorien, 
        trab, inotr, amname, amta, namta, t0, t1, iamt1, veold, 
        ielmat,matname, prlab,prset,filab, vold,nodebounold, 
        ndirbounold, xbounold, xforcold, xloadold, t1old, eme, 
        iponor, xnor, knor, thickn, thicke, offset, iponoel, inoel, rig, 
        shcon, nshcon, cocon, ncocon, ics, 
	sti, ener, xstate, jobnamec,infree,nnn,prestr,&iprestr,cbody, 
	ibody,xbody,&nbody,xbodyold,&ttime,qaold,cs,&mcs,output,
        physcon,ctrl,typeboun,fmpc,tieset,&ntie));
    }
  } 
	  
}

 FORTRAN(closefile,());

#ifdef CALCULIX_MPI
MPI_Finalize();
#endif

 return 0;
      
}




