/*     CalculiX - A 3-dimensional finite element program                 */
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


void checkconvergence(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	  int *ne, double *stn, int *nmethod, 
	  int *kode, char *filab, double *een, double *t1act,
          double *time, double *epn,int *ielmat,char *matname,
          double *enern, double *xstaten, int *nstate_, int *istep,
          int *iinc, int *iperturb, double *ener, int *mint_, char *output,
          int *ithermal, double *qfn, int *mode, int *noddiam, double *trab,
          int *inotr, int *ntrans, double *orab, int *ielorien, int *norien,
          char *description,double *sti,
	  int *icutb, int *iit, double *dtime, double *qa, double *vold,
          double *qam, double *ram1, double *ram2, double *ram,
          double *cam, double *uam, int *ntg, double *ttime,
          int *icntrl, double *theta, double *dtheta, double *veold,
          double *vini, int *idrct, double *tper,int *istab, double *tmax, 
          int *nactdof, double *b, double *tmin, double *ctrl, double *amta,
          int *namta, int *itpamp, int *inext, double *dthetaref, int *itp,
          int *jprint, int *jout, int *uncoupled, double *t1, int *iitterm,
          int *nelemload, int *nload, int *nodeboun, int *nboun, int *itg,
          int *ndirboun, double *deltmx){

    int k,*ipneigh=NULL,*neigh=NULL,*inum=NULL,id,istart,iend,inew,i,j;
    double c1[2],c2[2],reftime,*fn=NULL;

    int i0,ir,ip,ic,il,ig,ia,iest,iest1=0,iest2=0,iconvergence,idivergence,
      ngraph=1;
    double df,dc,db,dd,ran,can,rap,ea,cae,ral,da,*vr=NULL,*vi=NULL,*stnr=NULL,
    *stni=NULL,*vmax=NULL,*stnmax=NULL;

    i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
    df=ctrl[10];dc=ctrl[11];db=ctrl[12];da=ctrl[13];dd=ctrl[16];
    ran=ctrl[18];can=ctrl[19];rap=ctrl[22];
    ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];

    if(*ithermal!=2){
	if(qa[0]>ea*qam[0]){
	    if(*iit<=ip){c1[0]=ran;}
	    else{c1[0]=rap;}
	    c2[0]=can;
	}
	else{
	    c1[0]=ea;
	    c2[0]=cae;
	}
	if(ram1[0]<ram2[0]){ram2[0]=ram1[0];}
    }
    if(*ithermal>1){
	if(qa[1]>ea*qam[1]){
	    if(*iit<=ip){c1[1]=ran;}
	    else{c1[1]=rap;}
	    c2[1]=can;
	}
	else{
	    c1[1]=ea;
	    c2[1]=cae;
	}
	if(ram1[1]<ram2[1]){ram2[1]=ram1[1];}
    }

    iconvergence=0;

    /* mechanical */

    if(*ithermal<2){
	if((ram[0]<=c1[0]*qam[0])&&
	   ((cam[0]<=c2[0]*uam[0])||
	    (((ram[0]*cam[0]<c2[0]*uam[0]*ram2[0])||(ram[0]<=ral*qam[0])||(qa[0]<=ea*qam[0]))&&(*ntg==0))||
	    (cam[0]<1.e-8))) iconvergence=1;
    }

    /* thermal */

    if(*ithermal==2){
	if((ram[1]<=c1[1]*qam[1])&&
           (cam[2]<*deltmx)&&
	   ((cam[1]<=c2[1]*uam[1])||
	    (((ram[1]*cam[1]<c2[1]*uam[1]*ram2[1])||(ram[1]<=ral*qam[1])||(qa[1]<=ea*qam[1]))&&(*ntg==0))||
	    (cam[1]<1.e-8)))iconvergence=1;
    }

    /* thermomechanical */

    if(*ithermal==3){
	if(((ram[0]<=c1[0]*qam[0])&&
	    ((cam[0]<=c2[0]*uam[0])||
	     (((ram[0]*cam[0]<c2[0]*uam[0]*ram2[0])||(ram[0]<=ral*qam[0])||(qa[0]<=ea*qam[0]))&&(*ntg==0))||
	     (cam[0]<1.e-8)))&&
	   ((ram[1]<=c1[1]*qam[1])&&
            (cam[2]<*deltmx)&&
	    ((cam[1]<=c2[1]*uam[1])||
	     (((ram[1]*cam[1]<c2[1]*uam[1]*ram2[1])||(ram[1]<=ral*qam[1])||(qa[1]<=ea*qam[1]))&&(*ntg==0))||
	     (cam[1]<1.e-8))))iconvergence=1;
    }

	
	/* increment convergence reached */
	
    if(iconvergence==1){
	*ttime=*ttime+*dtime;
	/*	if((*uncoupled)&&(*ithermal==1))(*iit)-=(*iitterm);*/
	FORTRAN(writesummary,(istep,iinc,icutb,iit,ttime,time,dtime));
	if(*uncoupled){
	    if(*ithermal==2){
	        *iitterm=*iit;
		*ithermal=1;
		for(k=0;k<*nk;++k){t1[k]=vold[5*k];}
		*ttime=*ttime-*dtime;
		*iit=1;
		(ctrl[0])*=4;
		printf(" thermal convergence\n\n");
		return;
	    }else{
		*ithermal=3;
		*iit=*iitterm;
		(ctrl[0])/=4;
	    }
	}
	
	*icntrl=1;
	*icutb=0;
	*theta=*theta+*dtheta;
	
	/* defining a mean "velocity" for static calculations: is used to
	   extrapolate the present results for next increment */
	
	if(*nmethod != 4){
	    for(i=0;i<*nk;i++){
		for(j=1;j<4;j++){
		    veold[4*i+j]=(vold[5*i+j]-vini[5*i+j])/(*dtime);
		}
	    }
	    /* for(k=0;k<4**nk;++k){
		veold[k]=(vold[k]-vini[k])/(*dtime);
		}*/
	}
	
	/* check whether next increment size must be decreased */
	
	if((*iit>il)&&(*idrct==0)){
	    if(*idrct==0){
		*dtheta=*dthetaref*db;
		*dthetaref=*dtheta;
		printf(" convergence; the increment size is decreased to %e\n\n",*dtheta**tper);
		if(*dtheta<*tmin){
		    printf("\n *ERROR: increment size smaller than minimum\n");
		    printf(" best solution and residuals are in the frd file\n\n");
		    fn=NNEW(double,4**nk);
		    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
		    FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
                      nk,sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
		      ielorien,co,nelemload,nload,nodeboun,nboun,itg,ntg,
                      vold,ndirboun));
		    ++*kode;
		    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
		    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,
                         nmethod,kode,
			 filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			 xstaten,nstate_,istep,iinc,iperturb,ener,mint_,output,
			 ithermal,qfn,mode,noddiam,
			 trab,inotr,ntrans,orab,ielorien,norien,description,
			 ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
		    free(ipneigh);free(neigh);
		    FORTRAN(uout,(vold));
		    FORTRAN(stop,());
		}
	    }
	    else{
		printf("convergence\n\n");}
	}
	
	/* check whether next increment size can be increased */
	
	else if(*iit<=ig){
	    if((*istab==1)&&(*idrct==0)){
		*dtheta=*dthetaref*dd;
		*dthetaref=*dtheta;
		printf(" convergence; the increment size is increased to %e\n\n",*dtheta**tper);
	    }
	    else{
		*istab=1;
		printf(" convergence\n\n");
		*dtheta=*dthetaref;
	    }
	}
	else{
	    *istab=0;
	    printf(" convergence\n\n");
	    *dtheta=*dthetaref;
	}
	
	if((*dtheta>*tmax)&&(*idrct==0)){
	    *dtheta=*tmax;
	    *dthetaref=*dtheta;
	    printf(" the increment size exceeds thetamax and is decreased to %e\n\n",*dtheta**tper);
	}

        /* if itp=1 the increment just finished ends at a time point */

	if((*itpamp>0)&&(*idrct==0)){
	    if(*itp==1){
		*jprint=*jout;
	    }else{
		*jprint=*jout+1;
	    }
	    if(namta[3**itpamp-1]<0){
		reftime=*ttime+*dtheta**tper+1.01e-6;
	    }else{
		reftime=*time+*dtheta**tper+1.01e-6;
	    }
	    istart=namta[3**itpamp-3];
	    iend=namta[3**itpamp-2];
	    FORTRAN(identamta,(amta,&reftime,&istart,&iend,&id));
	    if(id<istart){
		inew=istart;
	    }else{
		inew=id+1;
	    }
	    if(*inext<inew){
		if(namta[3**itpamp-1]<0){
		    *dtheta=(amta[2**inext-2]-*ttime)/(*tper);
		}else{
		    *dtheta=(amta[2**inext-2]-*time)/(*tper);
		}
		(*inext)++;
		*itp=1;
		printf(" the increment size exceeds a time point and is decreased to %e\n\n",*dtheta**tper);
	    }else{*itp=0;}
	}

	if(*dtheta>1.-*theta){
	    *dtheta=1.-*theta;
	    *dthetaref=*dtheta;
	    printf(" the increment size exceeds the remainder of the step and is decreased to %e\n\n",*dtheta**tper);
	    if(*dtheta<=1.e-6){(*ttime)+=(*dtheta**tper);}
	}
    }
    else{
	
	/* check for the amount of iterations */
	
	if(*iit>ic){
	    printf("\n *ERROR: too many iterations needed\n");
	    printf(" best solution and residuals are in the frd file\n\n");
	    fn=NNEW(double,4**nk);
	    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	    FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,nk,sti,stn,
		ipkon,inum,kon,lakon,ne,mint_,orab,ielorien,co,
                nelemload,nload,nodeboun,nboun,itg,ntg,vold,ndirboun));
	    ++*kode;
	    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
	    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,nmethod,kode,
			 filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			 xstaten,nstate_,istep,iinc,iperturb,ener,mint_,output,
			 ithermal,qfn,mode,noddiam,
			 trab,inotr,ntrans,orab,ielorien,norien,description,
			 ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
	    free(ipneigh);free(neigh);
	    FORTRAN(uout,(vold));
	    FORTRAN(stop,());
	}	
	
	/* check for diverging residuals */
	
	if((*iit>=i0)||(fabs(ram[0])>1.e20)||(fabs(cam[0])>1.e20)||
	               (fabs(ram[1])>1.e20)||(fabs(cam[1])>1.e20)||
                       (cam[2]>*deltmx)){
	    idivergence=0;
	    if(*ithermal!=2){
		if((ram1[0]>ram2[0])&&(ram[0]>ram2[0])&&(ram[0]>c1[0]*qam[0]))
		    idivergence=1;
	    }

            /* for thermal calculations the maximum temperature change
               is checked as well */

            if(*ithermal>1){
	        if((ram1[1]>ram2[1])&&(ram[1]>ram2[1])&&(ram[1]>c1[1]*qam[1]))
		    idivergence=1;
		if(cam[2]>*deltmx) idivergence=2;
	    }
	    if(idivergence>0){
		if(*idrct==1) {
		    printf("\n *ERROR: solution seems to diverge; please try \n");
		    printf(" automatic incrementation; program stops\n");
		    printf(" best solution and residuals are in the frd file\n\n");
		    fn=NNEW(double,4**nk);
		    inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
		    FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,nk,
                       sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
		       ielorien,co,nelemload,nload,nodeboun,nboun,itg,ntg,
                       vold,ndirboun));
		    ++*kode;
		    ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
		    FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,
                       inum,nmethod,kode,
		       filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
		       xstaten,nstate_,istep,iinc,iperturb,ener,mint_,output,
		       ithermal,qfn,mode,noddiam,
		       trab,inotr,ntrans,orab,ielorien,norien,description,
		       ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
		    free(ipneigh);free(neigh);
		    FORTRAN(uout,(vold));
		    FORTRAN(stop,());
		}
		else {
		    if(idivergence==1){
		      *dtheta=*dtheta*df;
		    }else{
		      *dtheta=*dtheta**deltmx/cam[2]*da;
		    }
		    *dthetaref=*dtheta;
		    printf(" divergence; the increment size is decreased to %e\n",*dtheta**tper);
		    printf(" the increment is reattempted\n\n");
		    *istab=0;
		    if(*itp==1){
		      *itp=0;
		      (*inext)--;
		    }
		    if(*dtheta<*tmin){
			printf("\n *ERROR: increment size smaller than minimum\n");
			printf(" best solution and residuals are in the frd file\n\n");
			fn=NNEW(double,4**nk);
			inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
			FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
                           nk,sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
			   ielorien,co,nelemload,nload,nodeboun,nboun,
                           itg,ntg,vold,ndirboun));
			++*kode;
			ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
			FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,
                           inum,nmethod,kode,
			   filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			   xstaten,nstate_,istep,iinc,iperturb,ener,mint_,
                           output,ithermal,qfn,mode,noddiam,
			   trab,inotr,ntrans,orab,ielorien,norien,description,
			   ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
			free(ipneigh);free(neigh);
			FORTRAN(uout,(vold));
			FORTRAN(stop,());
		    }
		    *icntrl=1;
		    (*icutb)++;
		    if(*icutb>ia){
			printf("\n *ERROR: too many cutbacks\n");
			printf(" best solution and residuals are in the frd file\n\n");
			fn=NNEW(double,4**nk);
			inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
			FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
                           nk,sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
			   ielorien,co,nelemload,nload,nodeboun,nboun,
                           itg,ntg,vold,ndirboun));
			++*kode;
			ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
			FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,
                           inum,nmethod,kode,
			   filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			   xstaten,nstate_,istep,iinc,iperturb,ener,mint_,
                           output,ithermal,qfn,mode,noddiam,
			   trab,inotr,ntrans,orab,ielorien,norien,description,
			   ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
			free(ipneigh);free(neigh);
			FORTRAN(uout,(vold));
			FORTRAN(stop,());
		    }
		    if(*uncoupled){
		      if(*ithermal==1){
			(ctrl[0])/=4;
		      }
		      *ithermal=3;
		    }
		    return;
		}
	    }
	}
	
	/* check for too slow convergence */
	
	if(*iit>=ir){
	    if(*ithermal!=2){
		iest1=(int)ceil(*iit+log(ran*qam[0]/(ram[0]))/log(ram[0]/(ram1[0])));
	    }
	    if(*ithermal>1){
		iest2=(int)ceil(*iit+log(ran*qam[1]/(ram[1]))/log(ram[1]/(ram1[1])));
	    }
	    if(iest1>iest2){iest=iest1;}else{iest=iest2;}
	    if(iest>0){
	    printf(" estimated number of iterations till convergence = %d\n",
		   iest);
	    }
	    if((iest>ic)||(*iit==ic)){
		
		if(*idrct!=1){
		    *dtheta=*dtheta*dc;
		    *dthetaref=*dtheta;
		    printf(" too slow convergence; the increment size is decreased to %e\n",*dtheta**tper);
		    printf(" the increment is reattempted\n\n");
		    *istab=0;
		    if(*dtheta<*tmin){
			printf("\n *ERROR: increment size smaller than minimum\n");
			printf(" best solution and residuals are in the frd file\n\n");
			fn=NNEW(double,4**nk);
			inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
			FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
                           nk,sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
			   ielorien,co,nelemload,nload,nodeboun,nboun,
                           itg,ntg,vold,ndirboun));
			++*kode;
			ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
			FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,
                           inum,nmethod,kode,
			   filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			   xstaten,nstate_,istep,iinc,iperturb,ener,
                           mint_,output,ithermal,qfn,mode,noddiam,
			   trab,inotr,ntrans,orab,ielorien,norien,description,
			   ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
			free(ipneigh);free(neigh);
			FORTRAN(uout,(vold));
			FORTRAN(stop,());
		    }
		    *icntrl=1;
		    (*icutb)++;
		    if(*icutb>ia){
			printf("\n *ERROR: too many cutbacks\n");
			printf(" best solution and residuals are in the frd file\n\n");
			fn=NNEW(double,4**nk);
			inum=NNEW(int,*nk);for(k=0;k<*nk;k++) inum[k]=1;
			FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
                           nk,sti,stn,ipkon,inum,kon,lakon,ne,mint_,orab,
			   ielorien,co,nelemload,nload,nodeboun,nboun,
                           itg,ntg,vold,ndirboun));
			++*kode;
			ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
			FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,vold,stn,
                           inum,nmethod,kode,
			   filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
			   xstaten,nstate_,istep,iinc,iperturb,ener,mint_,
                           output,ithermal,qfn,mode,noddiam,
			   trab,inotr,ntrans,orab,ielorien,norien,description,
			   ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph));
			free(ipneigh);free(neigh);
			FORTRAN(uout,(vold));
			FORTRAN(stop,());
		    }
		    if(*uncoupled){
		      if(*ithermal==1){
			(ctrl[0])/=4;
		      }
		      *ithermal=3;
		    }
		    return;
		}
	    }
	}
	
	printf(" no convergence\n\n");
	
	(*iit)++;
	
    }
    
    return;
}
