/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                          */

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

#ifdef __WIN32
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

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
  
FILE *f1;
  
int *kon=NULL, *nodeboun=NULL, *ndirboun=NULL, *ipompc=NULL,
	*nodempc=NULL, *nodeforc=NULL, *ndirforc=NULL,
        *nelemload=NULL,im,*inodesd=NULL,nload1,*idefforc=NULL,
	*nactdof=NULL, *icol=NULL,*ics=NULL,
        *jq=NULL, *mast1=NULL, *irow=NULL, *rig=NULL,*idefbody=NULL,
	*ikmpc=NULL, *ilmpc=NULL, *ikboun=NULL, *ilboun=NULL,
	*nreorder=NULL,*ipointer=NULL,*idefload=NULL,
	*istartset=NULL, *iendset=NULL, *ialset=NULL, *ielmat=NULL,
	*ielorien=NULL, *nrhcon=NULL, *nodebounold=NULL, *ndirbounold=NULL,
	*nelcon=NULL, *nalcon=NULL, *iamforc=NULL,  *iamload=NULL,
	*iamt1=NULL, *namta=NULL, *ipkon=NULL, *iamboun=NULL,
	*nplicon=NULL, *nplkcon=NULL, *inotr=NULL, *iponor=NULL, *knor=NULL,
	*ikforc=NULL, *ilforc=NULL, *iponoel=NULL, *inoel=NULL, *nshcon=NULL,
        *ncocon=NULL,*ibody=NULL,*ielprop=NULL,*islavsurf=NULL,
        *ipoinpc=NULL,icfd=0,mt,nxstate,nload0;
    
double *co=NULL, *xboun=NULL, *coefmpc=NULL, *xforc=NULL,*clearini=NULL,
	*xload=NULL, *xbounold=NULL, *xforcold=NULL,
	*vold=NULL, *sti=NULL, *xloadold=NULL, *xnor=NULL,
	*reorder=NULL,*dcs=NULL, *thickn=NULL, *thicke=NULL, *offset=NULL,
	*elcon=NULL, *rhcon=NULL, *alcon=NULL, *alzero=NULL, *t0=NULL, *t1=NULL,
	*prestr=NULL, *orab=NULL, *amta=NULL, *veold=NULL, *accold=NULL,
        *t1old=NULL, *eme=NULL, *plicon=NULL, *pslavsurf=NULL, *plkcon=NULL,
	*xstate=NULL, *trab=NULL, *ener=NULL, *shcon=NULL, *cocon=NULL,
        *cs=NULL,*tietol=NULL,*fmpc=NULL,*prop=NULL,*t0g=NULL,*t1g=NULL,
	*xbody=NULL,*xbodyold=NULL;
    
double ctrl[27]={4.5,8.5,9.5,16.5,10.5,4.5,0.,5.5,0.,0.,0.25,0.5,0.75,0.85,0.,0.,1.5,0.,0.005,0.01,0.,0.,0.02,1.e-5,1.e-3,1.e-8,1.e30};
    
char *sideload=NULL, *set=NULL, *matname=NULL, *orname=NULL, *amname=NULL,
     *filab=NULL, *lakon=NULL, *labmpc=NULL, *prlab=NULL, *prset=NULL, 
     jobnamec[528]="",jobnamef[132]="",output[4]="asc", *typeboun=NULL,
     *inpc=NULL,*tieset=NULL,*cbody=NULL,fneig[132]="",*sideloadtemp=NULL,
     kind1[2]="T",kind2[2]="T";
    
int nk,ne,nboun,nmpc,nforc,nload,nprint,nset,nalset,nentries=15,
  nmethod,neq[3]={0,0,0},i,mpcfree=1,mei[4],j,nzl,nam,nbounold=0,
  nforcold=0,nloadold=0,nbody,nbody_=0,nbodyold=0,network=0,
  k,nzs[3],nmpc_=0,nload_=0,nforc_=0,istep,istat,nboun_=0,nintpoint=0,
  iperturb[2]={0,0},nmat,ntmat_=0,norien,ithermal[2]={0,0},
  iprestr,kode,isolver=0,nslavs=0,nkon_=0,ne0,nkon0,mortar=0,
  jout[2]={1,1},nlabel,nkon=0,idrct,jmax[2],iexpl,nevtot=0,ifacecount=0,
  iplas=0,npmat_=0,mi[3]={0,3,1},ntrans,mpcend=-1,namtot_=0,iumat=0,mpcmult,
  icascade=0,maxlenmpc,mpcinfo[4],ne1d=0,ne2d=0,infree[4]={0,0,0,0},
  callfrommain,nflow=0,jin=0,irstrt=0,nener=0,jrstrt=0,nenerold,
  nline,ipoinp[2*nentries],*inp=NULL,ntie,ntie_=0,mcs=0,nprop_=0,
  nprop=0,itpamp=0,iviewfile,nkold,nevdamp_=0,npt_=0,cyclicsymmetry;

int *meminset=NULL,*rmeminset=NULL;

int nzs_,nk_=0,ne_=0,nset_=0,nalset_=0,nmat_=0,norien_=0,nam_=0,
    ntrans_=0,ncs_=0,nstate_=0,ncmat_=0,memmpc_=0,nprint_=0;

double fei[3],tinc,tper,tmin,tmax,*xmodal=NULL,
    alpha,ttime=0.,qaold[2]={0.,0.},physcon[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

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
	printf("\nThis is version version 2.7\n\n");
	FORTRAN(stop,());
    }
  }
  if(jin==0){strcpy(jobnamec,argv[1]);strcpy1(jobnamef,argv[1],132);}

  for(i=1;i<argc;i++){
    if(strcmp1(argv[i],"-o")==0) {
    strcpy(output,argv[i+1]);break;}
  }
}

#ifdef BAM
int lop=0,lrestart=0,kstep=1,kinc=1;
double time[2],dtime;
FORTRAN(uexternaldb,(&lop,&lrestart,time,&dtime,&kstep,&kinc));
#endif

FORTRAN(openfile,(jobnamef,output));

printf("\n************************************************************\n\n");
printf("CalculiX version 2.7, Copyright(C) 1998-2014 Guido Dhondt\n");
printf("CalculiX comes with ABSOLUTELY NO WARRANTY. This is free\n");
printf("software, and you are welcome to redistribute it under\n");
printf("certain conditions, see gpl.htm\n\n");
printf("************************************************************\n\n");
printf("You are using an executable made on So 23. Feb 13:14:40 CET 2014\n");
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
   inpc,ipoinp,inp,&ntie_,&nbody_,&nprop_,ipoinpc,&nevdamp_,&npt_,&nslavs,
   &nkon_,&mcs,&mortar,&ifacecount,&nintpoint));

free(set);free(meminset);free(rmeminset);mt=mi[1]+1;

nzs_=20000000;

nload=0;nbody=0;nforc=0;nboun=0;nk=0;nmpc=0;nam=0;

/* caveat: if changing next line:
   - change noelfiles appropriately
   - change nlabel in geomview.f, expand.c, storeresidual.f
     and createmddof.f
   - change the dimension of label in geomview.f */

nlabel=46;

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
    kon=NNEW(int,nkon_);
    ipkon=NNEW(int,ne_);
    lakon=NNEW(char,8*ne_);

    /* property cards */

    if(nprop_>0){
	ielprop=NNEW(int,ne_);
	for(i=0;i<ne_;i++) ielprop[i]=-1;
	prop=NNEW(double,nprop_);
    }

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0)){
	iponor=NNEW(int,2*nkon_);
	for(i=0;i<2*nkon_;i++) iponor[i]=-1;
	xnor=NNEW(double,36*ne1d+24*ne2d);
	knor=NNEW(int,24*(ne1d+ne2d)*(mi[2]+1));
	thickn=NNEW(double,2*nk_);
	thicke=NNEW(double,mi[2]*nkon_);
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
    if((istep == 0)||((irstrt<0)&&(nam_>0)))iamboun=NNEW(int,nboun_);
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
    if((istep == 0)||((irstrt<0)&&(nam_>0)))iamforc=NNEW(int,nforc_);
    idefforc=NNEW(int,nforc_);
    xforc=NNEW(double,nforc_);
    ikforc=NNEW(int,nforc_);
    ilforc=NNEW(int,nforc_);

    /* distributed facial loads */

    nelemload=NNEW(int,2*nload_);
    if((istep == 0)||((irstrt<0)&&(nam_>0)))iamload=NNEW(int,2*nload_);
    idefload=NNEW(int,nload_);
    sideload=NNEW(char,20*nload_);
    xload=NNEW(double,2*nload_);

    /* distributed volumetric loads */

    cbody=NNEW(char,81*nbody_);
    idefbody=NNEW(int,nbody_);
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

    if(npmat_>0){
	plicon=NNEW(double,(2*npmat_+1)*ntmat_*nmat);
	nplicon=NNEW(int,(ntmat_+1)*nmat);
	plkcon=NNEW(double,(2*npmat_+1)*ntmat_*nmat);
	nplkcon=NNEW(int,(ntmat_+1)*nmat);
    }

    /* linear dynamic properties */
    
    xmodal=NNEW(double,11+nevdamp_);
    xmodal[10]=nevdamp_+0.5;

    /* internal state variables (nslavs is needed for restart
       calculations) */

    if(mortar==0){
	xstate=NNEW(double,nstate_*mi[0]*(ne+nslavs));
	nxstate=nstate_*mi[0]*(ne+nslavs);
    }else if(mortar==1){
	xstate=NNEW(double,nstate_*mi[0]*(ne+nintpoint));
	nxstate=nstate_*mi[0]*(ne+nintpoint);
    }

    /* material orientation */

    if((istep == 0)||((irstrt<0)&&(norien>0))) {
	orname=NNEW(char,80*norien);
	orab=NNEW(double,7*norien);
	ielorien=NNEW(int,mi[2]*ne_);
    }

    /* transformations */

    if((istep == 0)||((irstrt<0)&&(ntrans>0))) {
	trab=NNEW(double,7*ntrans);
	inotr=NNEW(int,2*nk_);
    }

    /* amplitude definitions */

    if((istep == 0)||((irstrt<0)&&(nam_>0))) {
	amname=NNEW(char,80*nam_);
	amta=NNEW(double,2*namtot_);
	namta=NNEW(int,3*nam_);
    }

    if((istep == 0)||((irstrt<0)&&(ithermal[0]>0))) {
	t0=NNEW(double,nk_);
	t1=NNEW(double,nk_);
	if((ne1d!=0)||(ne2d!=0)){
	    t0g=NNEW(double,2*nk_);
	    t1g=NNEW(double,2*nk_);
	}
    }

    /* the number in next line is NOT 1.2357111317 -> points
       to user input; instead it is a generic nonzero
       initialization */

    if(istep==0){
	DMEMSET(t0,0,nk_,1.2357111319);
	DMEMSET(t1,0,nk_,1.2357111319);
    }
    
    if((istep == 0)||((irstrt<0)&&(ithermal[0]>0)&&(nam_>0)))iamt1=NNEW(int,nk_);

    if((istep==0)||((irstrt<0)&&(iprestr>0)))prestr=NNEW(double,6*mi[0]*ne_);

    vold=NNEW(double,mt*nk_);
    veold=NNEW(double,mt*nk_);

    ielmat=NNEW(int,mi[2]*ne_);

    matname=NNEW(char,80*nmat);

    filab=NNEW(char,87*nlabel);

    /* tied constraints */

    if(ntie_>0){
      tieset=NNEW(char,243*ntie_);
      tietol=NNEW(double,3*ntie_);
      cs=NNEW(double,17*ntie_);
    }

    /* temporary fields for cyclic symmetry calculations */

    if((ncs_>0)||(npt_>0)){
      if(2*npt_>24*ncs_){
	ics=NNEW(int,2*npt_);
      }else{
	ics=NNEW(int,24*ncs_);
      }
      if(npt_>30*ncs_){
	  dcs=NNEW(double,npt_);
      }else{
	  dcs=NNEW(double,30*ncs_);
      }
    }

  }
  else {

      /* allocating and reallocating space for subsequent steps */
      
//    if(ncont==1)nload=nload0;
      
      if((nmethod!=4)&&(nmethod!=5)&&(nmethod!=8)&&(nmethod!=9)&& 
       ((abs(nmethod)!=1)||(iperturb[0]<2))){
        veold=NNEW(double,mt*nk_);
    }
    else{
      RENEW(veold,double,mt*nk_);
      DMEMSET(veold,mt*nk,mt*nk_,0.);
    }
    RENEW(vold,double,mt*nk_);
    DMEMSET(vold,mt*nk,mt*nk_,0.);

 /*   if(nmethod != 4){free(accold);}*/

    RENEW(nodeboun,int,nboun_);
    RENEW(ndirboun,int,nboun_);
    RENEW(typeboun,char,nboun_);
    RENEW(xboun,double,nboun_);
    RENEW(ikboun,int,nboun_);
    RENEW(ilboun,int,nboun_);

    RENEW(nodeforc,int,2*nforc_);
    RENEW(ndirforc,int,nforc_);
    idefforc=NNEW(int,nforc_);
    RENEW(xforc,double,nforc_);
    RENEW(ikforc,int,nforc_);
    RENEW(ilforc,int,nforc_);

    RENEW(nelemload,int,2*nload_);
    idefload=NNEW(int,nload_);
    RENEW(sideload,char,20*nload_);
    RENEW(xload,double,2*nload_);

    RENEW(cbody,char,81*nbody_);
    idefbody=NNEW(int,nbody_);
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
      RENEW(inotr,int,2*nk_);DMEMSET(inotr,2*nk,2*nk_,0);
    }

    RENEW(co,double,3*nk_);DMEMSET(co,3*nk,3*nk_,0.);

    if(ithermal[0] != 0){
	RENEW(t0,double,nk_);DMEMSET(t0,nk,nk_,0.);
	RENEW(t1,double,nk_);DMEMSET(t1,nk,nk_,0.);
	if((ne1d!=0)||(ne2d!=0)){
	    RENEW(t0g,double,2*nk_);DMEMSET(t0g,2*nk,2*nk_,0.);
	    RENEW(t1g,double,2*nk_);DMEMSET(t1g,2*nk,2*nk_,0.);
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
    if(mcs>ntie_) RENEW(cs,double,17*mcs);
    if(mortar==1){
	islavsurf=NNEW(int,2*ifacecount+2);
	pslavsurf=NNEW(double,3*nintpoint);
	clearini=NNEW(double,3*9*nslavs);
    }
			
  }

  nenerold=nener;
  nkold=nk;

  /* opening the eigenvalue file and checking for cyclic symmetry */

  strcpy(fneig,jobnamec);
  strcat(fneig,".eig");
  cyclicsymmetry=0;
  if((f1=fopen(fneig,"rb"))!=NULL){
      if(fread(&cyclicsymmetry,sizeof(int),1,f1)!=1){
	  printf("*ERROR reading the information whether cyclic symmetry is involved in the eigenvalue file");
	  exit(0);
      }
      fclose(f1);
  }

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
            ctrl,&maxlenmpc,&ne1d,&ne2d,&nener,vold,nodebounold,
            ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
            sti,ener,xstate,jobnamec,&irstrt,&ttime,
            qaold,output,typeboun,inpc,ipoinp,inp,tieset,tietol,
            &ntie,fmpc,cbody,ibody,xbody,&nbody,&nbody_,xbodyold,&nam_,
	    ielprop,&nprop,&nprop_,prop,&itpamp,&iviewfile,ipoinpc,&icfd,
	    &nslavs,t0g,t1g,&network,&cyclicsymmetry,idefforc,idefload,
	    idefbody,&mortar,&ifacecount,islavsurf,pslavsurf,clearini));

  nload0=nload;free(idefforc);free(idefload);free(idefbody);

  if((abs(nmethod)!=1)||(iperturb[0]<2))icascade=0;

/*	FORTRAN(writeboun,(nodeboun,ndirboun,xboun,typeboun,&nboun));*/

  if(istat<0) break;

  if(istep == 1) {

  /* tied contact constraints: generate appropriate MPC's */

    tiedcontact(&ntie, tieset, &nset, set,istartset, iendset, ialset,
       lakon, ipkon, kon,tietol,&nmpc, &mpcfree, &memmpc_,
       &ipompc, &labmpc, &ikmpc, &ilmpc,&fmpc, &nodempc, &coefmpc,
       ithermal, co, vold,&icfd,&nmpc_,mi,&nk,&istep,ikboun,&nboun,
       kind1,kind2);

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
	RENEW(thicke,double,mi[2]*nkon);
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
    RENEW(ielmat,int,mi[2]*ne);

    /* allocating space for the state variables */

    if(mortar==0){
	RENEW(xstate,double,nstate_*mi[0]*(ne+nslavs));
	for(i=nxstate;i<nstate_*mi[0]*(ne+nslavs);i++){xstate[i]=0.;}
    }else if(mortar==1){
	RENEW(xstate,double,nstate_*mi[0]*(ne+nintpoint));
	for(i=nxstate;i<nstate_*mi[0]*(ne+nintpoint);i++){xstate[i]=0.;}
    }

    /* next statements for plastic materials and nonlinear springs */

    if(npmat_>0){
	RENEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	RENEW(nplicon,int,(ntmat_+1)*nmat);
	RENEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
	RENEW(nplkcon,int,(ntmat_+1)*nmat);
    }

    /* material orientation */

    if(norien > 0) {
      RENEW(orname,char,80*norien);
      RENEW(ielorien,int,mi[2]*ne);
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

    if(ithermal[0] == 0){free(t0);free(t1);free(t0g);free(t1g);}
    if((ithermal[0] == 0)||(nam<=0)){free(iamt1);}

    if(ncs_>0){
      RENEW(ics,int,ncs_);
      free(dcs);
    }else if(npt_>0){free(ics);}

    if(mcs>0){
	RENEW(cs,double,17*mcs);
    }else{
	free(cs);
    }

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

  /* temperature loading */
  
  if(ithermal[0] != 0){
      RENEW(t0,double,nk);
      RENEW(t1,double,nk);
      if((ne1d!=0)||(ne2d!=0)){
	  RENEW(t0g,double,2*nk);
	  RENEW(t1g,double,2*nk);
      }
      if(nam > 0) {RENEW(iamt1,int,nk);}
  }

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

  if((nmethod==4)||(nmethod==5)||(nmethod==8)||(nmethod==9)||
     ((abs(nmethod)==1)&&(iperturb[0]>=2))){
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

  /* generate force convection elements */

  if(network==1){
      ne0=ne;nkon0=nkon;nload1=nload;
      RENEW(ipkon,int,ne+nload);
      RENEW(lakon,char,8*(ne+nload));
      RENEW(kon,int,nkon+9*nload);
      inodesd=NNEW(int,nk);
      RENEW(nelemload,int,4*nload);
      RENEW(sideload,char,40*nload);
      
      FORTRAN(genadvecelem,(inodesd,ipkon,&ne,lakon,kon,&nload,
			    sideload,nelemload,&nkon));
      
      free(inodesd);
      RENEW(ipkon,int,ne);
      RENEW(lakon,char,8*ne);
      RENEW(kon,int,nkon);
      RENEW(sti,double,6*mi[0]*ne);
      RENEW(eme,double,6*mi[0]*ne);
      if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
      if(nprop>0) RENEW(ielprop,int,ne);
      if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);
      RENEW(nelemload,int,2*nload);
      RENEW(sideload,char,20*nload);
      RENEW(xload,double,2*nload);
      RENEW(xloadold,double,2*nload);
      if(nam>0){
	  RENEW(iamload,int,2*nload);
	  for(i=2*nload1;i<2*nload;i++)iamload[i]=0;
      }
      if(nener==1)RENEW(ener,double,mi[0]*ne*2);
      if(norien>0)RENEW(ielorien,int,mi[2]*ne);
      RENEW(ielmat,int,mi[2]*ne);
      for(i=mi[2]*ne0;i<mi[2]*ne;i++)ielmat[i]=1;
  }

  if(ntrans > 0){
    RENEW(inotr,int,2*nk);
  }
  
  /*   calling the user routine ufaceload (can be empty) */

  if(ithermal[1]>=2){
      sideloadtemp=NNEW(char,20*nload);
      for(i=0;i<nload;i++){
	  strcpy1(&sideloadtemp[20*i],&sideload[20*i],20);
	  if((strcmp1(&sideload[20*i]," ")==0)&&
	     (strcmp1(&sideload[20*i+1]," ")!=0)){
	      strcpy1(&sideloadtemp[20*i],"F",1);
	  }
      }
      FORTRAN(ufaceload,(co,ipkon,kon,lakon,&nboun,nodeboun,
              nelemload,sideloadtemp,&nload,&ne,&nk));
      free(sideloadtemp);
  }
  
  /* decascading MPC's only necessary if MPC's changed */

  if(((istep == 1)||(ntrans>0)||(mpcend<0)||(nk!=nkold))&&(icascade==0)) {

    /* decascading the MPC's */

    printf(" Decascading the MPC's\n\n");

    callfrommain=1;
    cascade(ipompc,&coefmpc,&nodempc,&nmpc,
	    &mpcfree,nodeboun,ndirboun,&nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,&mpcmult,
	    labmpc,&nk,&memmpc_,&icascade,&maxlenmpc,
            &callfrommain,iperturb,ithermal);
  }

  /* determining the matrix structure: changes if SPC's have changed */
  
  if((icascade==0)&&(nmethod<8)) printf(" Determining the structure of the matrix:\n");
  
  nactdof=NNEW(int,mt*nk);  
  mast1=NNEW(int,nzs[1]);
  irow=NNEW(int,nzs[1]);
  
  if((mcs==0)||(cs[1]<0)){
      
      icol=NNEW(int,mt*nk);
      jq=NNEW(int,mt*nk+1);
      ipointer=NNEW(int,mt*nk);
      
      if((icascade==0)&&(nmethod<8)){
	  mastruct(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,ipompc,
		   nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,neq,
		   ikmpc,ilmpc,ipointer,nzs,&nmethod,ithermal,
                   ikboun,ilboun,iperturb,mi,&mortar);
      }
      else{neq[0]=1;neq[1]=1;neq[2]=1;}
  }
  else{
      
      icol=NNEW(int,8*nk);
      jq=NNEW(int,8*nk+1);
      ipointer=NNEW(int,8*nk);
      
      mastructcs(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,
		 ipompc,nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,
		 neq,ikmpc,ilmpc,ipointer,nzs,&nmethod,
		 ics,cs,labmpc,&mcs,mi,&mortar);
  }
  
  free(ipointer);free(mast1);
  if((icascade==0)&&(nmethod<8))RENEW(irow,int,nzs[2]);

  /* nmethod=1: static analysis   */
  /* nmethod=2: frequency analysis  */
  /* nmethod=3: buckling analysis */
  /* nmethod=4: (linear or nonlinear) dynamic analysis */
  /* nmethod=5: steady state dynamics analysis */
  /* nmethod=6: Coriolis frequency calculation */
  /* nmethod=7: flutter frequency calculation */
     

  if((nmethod<=1)||((iperturb[0]>1)&&(nmethod<8)))
    {
	if(iperturb[0]<2){
	
	linstatic(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,&tper,thicke,jobnamec,tieset,&ntie,&istep);

      }

      else{

	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	nonlingeo(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,&ikmpc, 
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,&tinc,&tper,&tmin,&tmax,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,iamload,iamt1,&alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,tieset,&itpamp,&iviewfile,jobnamec,tietol,&nslavs,thicke,
	     ics,&nintpoint,&mortar,
	     &ifacecount,typeboun,&islavsurf,&pslavsurf,&clearini);

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];

      }
    }else if(nmethod==2){
      
      /* FREQUENCY ANALYSIS */
      
      if((mcs==0)||(cs[1]<0)){
#ifdef ARPACK
	  
	  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	  mpcinfo[3]=maxlenmpc;
	  
	  arpack(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     shcon,nshcon,cocon,ncocon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,vold,iperturb,sti,nzs,
	     &kode,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,&ncmat_,&nstate_,&ener,jobnamec,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,&nslavs,tietol,&nkon,mpcinfo,
	     &ntie,&istep,&mcs,ics,tieset,cs,&nintpoint,&mortar,&ifacecount,
	     &islavsurf,&pslavsurf,&clearini);

	  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
	  maxlenmpc=mpcinfo[3];

#else
	  printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	  FORTRAN(stop,());
#endif

      }else{

#ifdef ARPACK

	  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	  mpcinfo[3]=maxlenmpc;
	  
	  arpackcs(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,
             xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,ics,cs,&mpcend,&ncmat_,
             &nstate_,&mcs,&nkon,&ener,jobnamec,output,set,&nset,istartset,
             iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody,&nevtot,thicke,&nslavs,tietol,mpcinfo,
	     &ntie,&istep,tieset,&nintpoint,&mortar,&ifacecount,&islavsurf,
             &pslavsurf,&clearini);

	  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
	  maxlenmpc=mpcinfo[3];

#else
	  printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	  FORTRAN(stop,());
#endif

      }
  }else if(nmethod==3){
    
#ifdef ARPACK
	arpackbu(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, 
	     nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,&ncmat_,&nstate_,ener,output,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,jobnamec);
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
            nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	    &ntie,tieset,&idrct,jmax,&tmin,&tmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,&nslavs);
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
	    &tinc,&tper,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
	    xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,ics,cs,&mpcend,
	    ctrl,ikforc,ilforc,thicke);
    }
  else if((nmethod==6)||(nmethod==7))
    {

      printf(" Composing the complex eigenmodes from the real eigenmodes\n\n");

      complexfreq(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
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
	    &ntie,tieset,&idrct,jmax,&tmin,&tmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,jobnamef,mei);
    }
  else if(nmethod>7){

	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	electromagnetics(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,
             ndirboun,xboun,&nboun, 
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,&ikmpc, 
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,&tinc,&tper,&tmin,&tmax,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,iamload,iamt1,&alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             &set,&nset,&istartset,&iendset,&ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,&tieset,&itpamp,&iviewfile,jobnamec,&tietol,&nslavs,thicke,
	     ics,&nalset,&nmpc_);

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];
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
  }else{
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

  /* removing the advective elements, if any */

  if(network==1){
      ne=ne0;nkon=nkon0;
      RENEW(ipkon,int,ne);
      RENEW(lakon,char,8*ne);
      RENEW(kon,int,nkon);
      RENEW(sti,double,6*mi[0]*ne);
      RENEW(eme,double,6*mi[0]*ne);
      if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
      if(nprop>0) RENEW(ielprop,int,ne);
      if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);
      if(nener==1)RENEW(ener,double,mi[0]*ne*2);
      if(norien>0)RENEW(ielorien,int,mi[2]*ne);
      RENEW(ielmat,int,mi[2]*ne);
  }

  nload=nload0;

  if((nmethod == 4)&&(iperturb[0]>1)) free(accold);

  if(irstrt>0){
    jrstrt++;
    if(jrstrt==irstrt){
      jrstrt=0;
      FORTRAN(restartwrite,(&istep,&nset,&nload,&nforc,&nboun,&nk,&ne,
        &nmpc,&nalset,&nmat,&ntmat_,&npmat_,&norien,&nam,&nprint, 
        mi,&ntrans,&ncs_,&namtot_,&ncmat_,&mpcend,&maxlenmpc,&ne1d,
        &ne2d,&nflow,&nlabel,&iplas,&nkon,ithermal,&nmethod,iperturb,
        &nstate_,&nener,set,istartset,iendset,ialset,co,kon,ipkon,
        lakon,nodeboun,ndirboun,iamboun,xboun,ikboun,ilboun,ipompc,
        nodempc,coefmpc,labmpc,ikmpc,ilmpc,nodeforc,ndirforc,iamforc,
        xforc,ikforc,ilforc,nelemload,iamload,sideload,xload,
        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
        alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,ielorien,
        trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
        ielmat,matname,prlab,prset,filab,vold,nodebounold,
        ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
        iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
        shcon,nshcon,cocon,ncocon,ics,
	sti,ener,xstate,jobnamec,infree,prestr,&iprestr,cbody,
	ibody,xbody,&nbody,xbodyold,&ttime,qaold,cs,&mcs,output,
	physcon,ctrl,typeboun,fmpc,tieset,&ntie,tietol,&nslavs,t0g,t1g,
	&nprop,ielprop,prop,&mortar,&nintpoint,&ifacecount,islavsurf,
	pslavsurf,clearini));
    }
  } 
	  
}

FORTRAN(closefile,());

strcpy(fneig,jobnamec);
strcat(fneig,".frd");
if((f1=fopen(fneig,"ab"))==NULL){
    printf("*ERROR in frd: cannot open frd file for writing...");
    exit(0);
}
fprintf(f1," 9999\n");
fclose(f1);

/* deallocating the fields
   this section is addressed immediately after leaving calinput */

free(ipoinpc);free(inpc);free(inp);

if(ncs_>0) free(ics);
if((ncs_<=0)&&(npt_>0)) free(dcs);
if(mcs>0) free(cs);
free(tieset);free(tietol);

free(co);free(kon);free(ipkon);free(lakon);

free(nodeboun);free(ndirboun);free(typeboun);free(xboun);free(ikboun);
free(ilboun);free(nodebounold);free(ndirbounold);free(xbounold);

free(ipompc);free(labmpc);free(ikmpc);free(ilmpc);free(fmpc);
free(nodempc);free(coefmpc);

free(nodeforc);free(ndirforc);free(xforc);free(ikforc);free(ilforc);
free(xforcold);

free(nelemload);free(sideload);free(xload);free(xloadold);

free(cbody);free(ibody);free(xbody);free(xbodyold);

if(nam>0){free(iamboun);free(iamforc);free(iamload);free(amname);
    free(amta);free(namta);}

free(set);free(istartset);free(iendset);free(ialset);

free(elcon);free(nelcon);free(rhcon);free(nrhcon);free(shcon);free(nshcon);
free(cocon);free(ncocon);free(alcon);free(nalcon);free(alzero);
if(nprop>0){free(ielprop);free(prop);}
if(npmat_>0){free(plicon);free(nplicon);free(plkcon);free(nplkcon);}

if(norien>0){free(orname);free(orab);free(ielorien);}
if(ntrans>0){free(trab);free(inotr);}
if(iprestr>0){free(prestr);}

if(ithermal[0]!=0){
    free(t0);free(t1);free(t1old);
    if(nam>0) free(iamt1);
    if((ne1d!=0)||(ne2d!=0)){free(t0g);free(t1g);}
}

free(prlab);free(prset);free(filab);free(xmodal);

free(ielmat);free(matname);

free(sti);free(eme);free(ener);free(xstate);

free(vold);free(veold);

if((ne1d!=0)||(ne2d!=0)){
    free(iponor);free(xnor);free(knor);free(thicke);free(offset);
    free(iponoel);free(inoel);free(rig);
}

if((mortar==1)&&(nstate_!=0)){
    free(islavsurf);free(pslavsurf);
}

#ifdef CALCULIX_MPI
MPI_Finalize();
#endif

 return 0;
      
}




