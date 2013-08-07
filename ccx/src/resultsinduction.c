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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static char *lakon1,*matname1,*sideload1;

static int *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iperturb1,*iout1,*nmethod1,
    *nplkcon1,*npmat1_,*mi1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_fn1,calcul_qa1,
    *nal=NULL,*ipompc1,*nodempc1,*nmpc1,*ncocon1,*ikmpc1,*ilmpc1,
    num_cpus,mt1,*nk1,*nshcon1,*nelemload1,*nload1;

static double *co1,*v1,*elcon1,*rhcon1,*alcon1,*orab1,*t01,
    *fn1=NULL,*qa1=NULL,*vold1,*dtime1,*time1,
    *ttime1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*sti1,
    *eei1,*springarea1,*reltime1,*coefmpc1,
    *cocon1,*qfx1,*shcon1,*xload1,
    *xloadold1,*h01;

void resultsinduction(double *co,int *nk,int *kon,int *ipkon,char *lakon,
       int *ne,
       double *v,double *stn,int *inum,double *stx,double *elcon,int *nelcon,
       double *rhcon,int *nrhcon,double *alcon,int *nalcon,double *alzero,
       int *ielmat,int *ielorien,int *norien,double *orab,int *ntmat_,
       double *t0,
       double *t1,int *ithermal,double *prestr,int *iprestr,char *filab,
       double *eme,double *emn,
       double *een,int *iperturb,double *f,double *fn,int *nactdof,int *iout,
       double *qa,double *vold,double *b,int *nodeboun,int *ndirboun,
       double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
       char *labmpc,int *nmpc,int *nmethod,double *cam,int *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,int *nplicon,double *plkcon,
       int *nplkcon,double *xstateini,double *xstiff,double *xstate,int *npmat_,
       double *epn,char *matname,int *mi,int *ielas,int *icmd,int *ncmat_,
       int *nstate_,
       double *sti,double *vini,int *ikboun,int *ilboun,double *ener,
       double *enern,double *emeini,double *xstaten,double *eei,double *enerini,
       double *cocon,int *ncocon,char *set,int *nset,int *istartset,
       int *iendset,
       int *ialset,int *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       int *inotr,int *ntrans,double *fmpc,int *nelemload,int *nload,
       int *ikmpc,int *ilmpc,
       int *istep,int *iinc,double *springarea,double *reltime, int *ne0,
       double *xforc, int *nforc, double *thicke,double *xnormastface,
       double *shcon,int *nshcon,char *sideload,double *xload,
       double *xloadold,int *icfd,int *inomat,double *h0){

    int intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,iener,ikin,
        intpointvart,mt=mi[1]+1,i,j,*ithread=NULL;

    /*

     calculating integration point values (strains, stresses,
     heat fluxes, material tangent matrices and nodal forces)

     storing the nodal and integration point results in the
     .dat file

     iout=-2: v is assumed to be known and is used to
              calculate strains, stresses..., no result output
              corresponds to iout=-1 with in addition the
              calculation of the internal energy density
     iout=-1: v is assumed to be known and is used to
              calculate strains, stresses..., no result output;
              is used to take changes in SPC's and MPC's at the
              start of a new increment or iteration into account
     iout=0: v is calculated from the system solution
             and strains, stresses.. are calculated, no result output
     iout=1:  v is calculated from the system solution and strains,
              stresses.. are calculated, requested results output
     iout=2: v is assumed to be known and is used to 
             calculate strains, stresses..., requested results output */
      
    /* variables for multithreading procedure */
    
    int sys_cpus;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_RESULTS");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
	
    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }

// next line is to be inserted in a similar way for all other paralell parts

    if(*ne<num_cpus) num_cpus=*ne;
    
    pthread_t tid[num_cpus];
    
    /* 1. nodewise storage of the primary variables
       2. determination which derived variables have to be calculated */

    FORTRAN(resultsini,(nk,v,ithermal,filab,iperturb,f,fn,
       nactdof,iout,qa,vold,b,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
       veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,&iener,
       &ikin,&intpointvart,xforc,nforc));

   /* next statement allows for storing the displacements in each
      iteration: for debugging purposes */

    if((strcmp1(&filab[3],"I")==0)&&(*iout==0)){
	FORTRAN(frditeration,(co,nk,kon,ipkon,lakon,ne,v,
			      ttime,ielmat,matname,mi,istep,iinc));
    }

    /* calculating the stresses and material tangent at the 
       integration points; calculating the internal forces */

    if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){

	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;v1=v;elcon1=elcon;
        nelcon1=nelcon;ielmat1=ielmat;ntmat1_=ntmat_;vold1=vold;dtime1=dtime;
        matname1=matname;mi1=mi;ncmat1_=ncmat_;sti1=sti;eei1=eei;alcon1=alcon;
	nalcon1=nalcon;h01=h0;

	/* calculating the stresses */
	
	if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
		printf(" Using up to %d cpu(s) for the stress calculation.\n\n", num_cpus);
	}
	
	/* create threads and wait */
	
	ithread=NNEW(int,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsemmt, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	free(ithread);
    }

    /* calculating the thermal flux and material tangent at the 
       integration points; calculating the internal point flux */

    if((ithermal[0]>=2)&&(intpointvart==1)){

	fn1=NNEW(double,num_cpus*mt**nk);
	qa1=NNEW(double,num_cpus*3);
	nal=NNEW(int,num_cpus);

	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;v1=v;
        elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;nrhcon1=nrhcon;
	ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
        ntmat1_=ntmat_;t01=t0;iperturb1=iperturb;iout1=iout;vold1=vold;
        ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
        dtime1=dtime;time1=time;ttime1=ttime;plkcon1=plkcon;
        nplkcon1=nplkcon;xstateini1=xstateini;xstiff1=xstiff;
        xstate1=xstate;npmat1_=npmat_;matname1=matname;mi1=mi;
        ncmat1_=ncmat_;nstate1_=nstate_;cocon1=cocon;ncocon1=ncocon;
        qfx1=qfx;ikmpc1=ikmpc;ilmpc1=ilmpc;istep1=istep;iinc1=iinc;
        springarea1=springarea;calcul_fn1=calcul_fn;calcul_qa1=calcul_qa;
        mt1=mt;nk1=nk;shcon1=shcon;nshcon1=nshcon;ithermal1=ithermal;
        nelemload1=nelemload;nload1=nload;nmethod1=nmethod;reltime1=reltime;
        sideload1=sideload;xload1=xload;xloadold1=xloadold;

	/* calculating the heat flux */
	
	printf(" Using up to %d cpu(s) for the heat flux calculation.\n\n", num_cpus);
	
	/* create threads and wait */
	
	ithread=NNEW(int,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsthermemmt, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	
	for(i=0;i<*nk;i++){
		fn[mt*i]=fn1[mt*i];
	}
	for(i=0;i<*nk;i++){
	    for(j=1;j<num_cpus;j++){
		fn[mt*i]+=fn1[mt*i+j*mt**nk];
	    }
	}
	free(fn1);free(ithread);
	
        /* determine the internal concentrated heat flux */

	qa[1]=qa1[1];
	for(j=1;j<num_cpus;j++){
	    qa[1]+=qa1[1+j*3];
	}
	
	free(qa1);
	
	for(j=1;j<num_cpus;j++){
	    nal[0]+=nal[j];
	}

	if(calcul_qa==1){
	    if(nal[0]>0){
		qa[1]/=nal[0];
	    }
	}
	free(nal);
    }

    /* storing results in the .dat file
       extrapolation of integration point values to the nodes
       interpolation of 3d results for 1d/2d elements */

    FORTRAN(resultsprint,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
       stx,ielorien,norien,orab,t1,ithermal,filab,een,iperturb,fn,
       nactdof,iout,vold,nodeboun,ndirboun,nboun,nmethod,ttime,xstate,
       epn,mi,
       nstate_,ener,enern,xstaten,eei,set,nset,istartset,iendset,
       ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
       nelemload,nload,&ikin,ielmat,thicke,eme,emn,rhcon,nrhcon,shcon,
       nshcon,cocon,ncocon,ntmat_,sideload,icfd,inomat));
  
  return;

}

/* subroutine for multithreading of resultsem */

void *resultsemmt(int *i){

    int indexfn,indexqa,indexnal,nea,neb,nedelta;

    indexfn=*i*mt1**nk1;
    indexqa=*i*3;
    indexnal=*i;

    nedelta=(int)floor(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;

    FORTRAN(resultsem,(co1,kon1,ipkon1,lakon1,v1,elcon1,nelcon1,ielmat1,
       ntmat1_,vold1,dtime1,matname1,mi1,ncmat1_,&nea,&neb,sti1,eei1,alcon1,
       nalcon1,h01));

    return NULL;
}

/* subroutine for multithreading of resultsmech */

void *resultsthermemmt(int *i){

    int indexfn,indexqa,indexnal,nea,neb,nedelta;

    indexfn=*i*mt1**nk1;
    indexqa=*i*3;
    indexnal=*i;
    
    nedelta=(int)floor(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;
//    if(neb>*ne1) neb=*ne1;

//    printf("therm i=%d,nea=%d,neb=%d\n",i,nea,neb);
//    printf("indexfn=%d,indexqa=%d,indexnal=%d\n",indexfn,indexqa,indexnal);

    FORTRAN(resultstherm,(co1,kon1,ipkon1,lakon1,ne1,v1,
	   elcon1,nelcon1,rhcon1,nrhcon1,ielmat1,ielorien1,norien1,orab1,
	   ntmat1_,t01,iperturb1,&fn1[indexfn],shcon1,nshcon1,
	   iout1,&qa1[indexqa],vold1,ipompc1,nodempc1,coefmpc1,nmpc1,
           dtime1,time1,ttime1,plkcon1,nplkcon1,xstateini1,xstiff1,xstate1,
           npmat1_,
           matname1,mi1,ncmat1_,nstate1_,cocon1,ncocon1,
           qfx1,ikmpc1,ilmpc1,istep1,iinc1,springarea1,
	   &calcul_fn1,&calcul_qa1,&nal[indexnal],&nea,&neb,ithermal1,
           nelemload1,nload1,nmethod1,reltime1,sideload1,xload1,xloadold1));

    return NULL;
}
