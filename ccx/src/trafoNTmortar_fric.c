/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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
#include <time.h>
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/** changing au due to N and T (normal and tangential
   *    direction at the slave surface) 
   * 	changing b due to N and T (normal and tangential
   *	direction at the slave surface) 
 * @param [out] au
 * @param [out] b
*/
void trafoNTmortar_fric(int *neq,int *nzs, int *islavactdof,int *islavact,int *nslavnode, int *nmastnode, int *ncone, 
        double *ad, double *au, double *b, int *irow, int *jq,
        int *nzsc, double *auc,
        double *adc, int *irowc, int *jqc,
        double *gap, double *bdd, double *auqdt, int *irowqdt,
        int *jqqdt, int *nzsqdt, int *nzlc,double *slavnor, double *slavtan,double *bhat,
	double *aubd, int *irowbd, int *jqbd, double *vold, double *cstress,double *bp_old,int *nactdof,
	int *islavnode, int *ntie, int *mi,int *nk,double *friccoeff){

    int i,j,k,l,m,numb,ntrimax,debug,idof1,idof2,idof3,imodification,
        nzsbd,jrow,jcol,islavnodeentry,jslavnodeentry,mt=mi[1]+1;

    double t1,t2,t3,e1,e2,e3,help,s,dd1,dd2,upt1,upt2,n1,n2,n3,*rstick=NULL,
         *bhat2=NULL, *auqdt2=NULL, *unitmatrix=NULL,bp, *t=NULL,*n=NULL, *up_t=NULL,
         *tslip=NULL,*vp=NULL,*mp=NULL,*fp=NULL,ep,constant=1.E10,*lmhat=NULL,nlmhat,
         alpha,beta,delta,lambda_n,*lambda_t=NULL,nlambda_t,det,*rp=NULL,*hp=NULL,*rphat=NULL,
         *lp=NULL,*ltslip=NULL,aux,*mp2=NULL,*mp3=NULL,*u_old=NULL,scal,*ltu=NULL,up_n;
    
    debug=0;
    imodification=1;
    

    u_old=NNEW(double,neq[1]);
    lmhat=NNEW(double,2);
    lambda_t=NNEW(double,2);
    fp=NNEW(double,4);
    mp=NNEW(double,4);
    mp2=NNEW(double,4);
    mp3=NNEW(double,4);    
    lp=NNEW(double,4);
    ltslip=NNEW(double,6);
    rp=NNEW(double,2);
    rphat=NNEW(double,2);
    hp=NNEW(double,2);
    n=NNEW(double,3);
    t=NNEW(double,6);
    up_t=NNEW(double,2);
    vp=NNEW(double,2);
    rstick=NNEW(double, 2*3);
    tslip=NNEW(double, 2*3);
    ltu=NNEW(double,2);
    
    
    bhat2=NNEW(double,neq[1]);
    auqdt2=NNEW(double,jqqdt[neq[1]]);
    unitmatrix=NNEW(double,neq[1]);
    for(j=1;j<*nk+1;j++){
	    idof1=nactdof[mt*j-3]-1;
	    idof2=nactdof[mt*j-2]-1;
	    idof3=nactdof[mt*j-1]-1; 
	    if(idof1>-1) u_old[idof1]=vold[mt*j-3];
	    if(idof2>-1) u_old[idof2]=vold[mt*j-2];
	    if(idof3>-1) u_old[idof3]=vold[mt*j-1];
//	    if(idof1>-1 && idof2>-1 && idof3>-1){
//	    printf("u(%d)=%e %e %e \n",j,u_old[idof1],u_old[idof2],u_old[idof3]);
//	    }
    }  
    
    for(j=0;j<neq[1];j++){
	    unitmatrix[j]=1.;
    }
    for(j=0;j<jqqdt[neq[1]];j++){ 
           auqdt2[j]=-auqdt[j];
    }      
	FORTRAN(opnonsym, (&neq[1], &aux, u_old, bhat2, unitmatrix, auqdt2, jqqdt, irowqdt));
	
	free(unitmatrix);free(auqdt2);
	
//    if(debug==1){  
    for (i=0;i<*ntie;i++){
//            FORTRAN(getmu,(&mu,tietol,elcon,i,ncmat_,ntmat_));
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    islavnodeentry=islavnode[j];
//	    friccoeff[j]=mu;
	    /* mechanical degrees of freedom */
	    
	    idof1=nactdof[mt*islavnodeentry-3]-1;
	    idof2=nactdof[mt*islavnodeentry-2]-1;
	    idof3=nactdof[mt*islavnodeentry-1]-1;
	    
//	    printf("idof %d %d %d\n",idof1,idof2,idof3);
//	    printf("u(%d)=%e %e %e \n",islavnodeentry,u_old[idof1],u_old[idof2],u_old[idof3]); 	    
//	    printf("uhat(%d)=%e %e %e \n",islavnodeentry,bhat2[idof1],bhat2[idof2],bhat2[idof3]);    
//	    printf("bp= %e\n",bp_old[j]);
	}
    }
//    }	


     for(j=0;j<neq[1];j++){ //loop over columns 
        for(i=jq[j]-1;i<jq[j+1]-1;i++){ //loop over rows
	  jslavnodeentry = floor(islavactdof[j]/10.);
	  jcol= islavactdof[j]-10*jslavnodeentry;
	  k=irow[i]-1;
          islavnodeentry = floor(islavactdof[k]/10.);
	  jrow= islavactdof[k]-10*islavnodeentry;
	  
	  if(islavactdof[k]>0 && islavact[islavnodeentry-1]==1){ // row stick
	    scal= 1.0/bdd[k];
//	    scal=1.0;
            for(l=0;l<3;l++){
	     n[l]=slavnor[3*(islavnodeentry-1)+l];
	    }
	    for(l=0;l<6;l++){
	     t[l]=slavtan[6*(islavnodeentry-1)+l];
	    }
	    if(jrow==1){
            up_t[0]=bhat2[k]*t[0]+bhat2[k+1]*t[1]+bhat2[k+2]*t[2];
	    up_t[1]=bhat2[k]*t[3]+bhat2[k+1]*t[4]+bhat2[k+2]*t[5];
	    lambda_n=n[0]*cstress[k]+n[1]*cstress[k+1]+n[2]*cstress[k+2];
	    lambda_t[0]=t[0]*cstress[k]+t[1]*cstress[k+1]+t[2]*cstress[k+2];
	    lambda_t[1]=t[3]*cstress[k]+t[4]*cstress[k+1]+t[5]*cstress[k+2];
	    }else if(jrow==2){
            up_t[0]=bhat2[k-1]*t[0]+bhat2[k]*t[1]+bhat2[k+1]*t[2];
	    up_t[1]=bhat2[k-1]*t[3]+bhat2[k]*t[4]+bhat2[k+1]*t[5];
	    lambda_n=n[0]*cstress[k-1]+n[1]*cstress[k]+n[2]*cstress[k+1];
	    lambda_t[0]=t[0]*cstress[k-1]+t[1]*cstress[k]+t[2]*cstress[k+1];
	    lambda_t[1]=t[3]*cstress[k-1]+t[4]*cstress[k]+t[5]*cstress[k+1];
	    } 
	    nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
	    bp=bp_old[islavnodeentry-1];

	    for(l=0;l<3;l++){
	     for(m=0;m<2;m++){
	      rstick[m*3+l]=friccoeff[islavnodeentry-1]*up_t[m]*n[l];
	     }
	    }
	    if(debug==1){
	    printf("j=%d k=%d jslavnode=%d kslavnode=%d jcol=%d jrow=%d kact=%d\n",j,k,jslavnodeentry,islavnodeentry,jcol,jrow,islavact[islavnodeentry-1]);
	    printf("\t n= %e %e %e \n",n[0],n[1],n[2]);
	    printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);
	    printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);
	    printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);
	    if(jrow==1){
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k],bhat2[k+1],bhat2[k+2],up_t[0],up_t[1],neq[1]);
	    }else{
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k-1],bhat2[k],bhat2[k+1],up_t[0],up_t[1],neq[1]);
	    }
	    printf("\t rstick1: %e %e %e\n",rstick[0],rstick[1],rstick[2]);
	    printf("\t rstick2: %e %e %e\n",rstick[3],rstick[4],rstick[5]);
	    }
	    if(islavnodeentry==jslavnodeentry){ // row stick and column stick
//	        printf("\t -bp*t[jcol]=%e -bp*t[3+jcol]=%e\n",-bp*t[jcol-1],-bp*t[3+jcol-1]);
	      if(jrow==2 && jcol==1){ 
	        e1=adc[j];
	        e2=auc[i];
	        e3=auc[i+1];
	        ad[j]=n[0];
		
	        au[i]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3 )-bp*t[jcol-1];
	        au[i+1]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3) -bp*t[3+jcol-1];
		if(debug==1){
	        printf("\t au= %e %e %e\n",e1,e2,e3);
	        printf("\t au= %e %e %e\n",ad[j],au[i],au[i+1]);
		}
	      }else if(jrow==1 && jcol==2){
	        e1=auc[i];
	        e2=ad[j];
	        e3=auc[i+1];
	        au[i]=n[1];
	        ad[j]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3) -bp*t[jcol-1];
	        au[i+1]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3) -bp*t[3+jcol-1];
		if(debug==1){
	        printf("\t au= %e %e %e\n",e1,e2,e3);
	        printf("\t au= %e %e %e\n",au[i],ad[j],au[i+1]);
		}
	      }else if(jrow==1 && jcol==3){
	        e1=auc[i];
	        e2=auc[i+1];
	        e3=ad[j];
	        au[i]=n[2];
	        au[i+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3) -bp*t[jcol-1];
	        ad[j]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3) -bp*t[3+jcol-1];
		if(debug==1){
	        printf("\t au= %e %e %e\n",e1,e2,e3);
	        printf("\t au= %e %e %e\n",au[i],au[i+1],ad[j]);
		}
	      }else{
	        printf("trafoNTmortar_fric: somethings wrong with diagonal stop!\n");
	        FORTRAN(stop,());		
	      }	
	      i=i+1;
	    }else if(jrow==1){ // row stick and column normal or slip
	      e1=auc[i];
	      e2=auc[i+1];
	      e3=auc[i+2];
	      au[i]=0.0;
	      au[i+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3);
	      au[i+2]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3);	      
              if(debug==1){
	      printf("\t au= %e %e %e\n",e1,e2,e3);
	      printf("\t au= %e %e %e\n",au[i],au[i+1],au[i+2]);
	      }
	      i=i+2;	      
	    }else{
	        printf("trafoNTmortar_fric: somethings wrong with idofs stop!\n");
	        FORTRAN(stop,());		
	    }
	    
	  }else if(islavactdof[k]>0 && islavact[islavnodeentry-1]==2){ //row slip 
	    scal= 1.0/bdd[k];
//	    scal=1.0;
            for(l=0;l<3;l++){
	     n[l]=slavnor[3*(islavnodeentry-1)+l];
	    }
	    for(l=0;l<6;l++){
	     t[l]=slavtan[6*(islavnodeentry-1)+l];
	    }
	    if(jrow==1){
            up_t[0]=bhat2[k]*t[0]+bhat2[k+1]*t[1]+bhat2[k+2]*t[2];
	    up_t[1]=bhat2[k]*t[3]+bhat2[k+1]*t[4]+bhat2[k+2]*t[5];
	    lambda_n=n[0]*cstress[k]+n[1]*cstress[k+1]+n[2]*cstress[k+2];
	    lambda_t[0]=t[0]*cstress[k]+t[1]*cstress[k+1]+t[2]*cstress[k+2];
	    lambda_t[1]=t[3]*cstress[k]+t[4]*cstress[k+1]+t[5]*cstress[k+2];
	    }else if(jrow==2){
            up_t[0]=bhat2[k-1]*t[0]+bhat2[k]*t[1]+bhat2[k+1]*t[2];
	    up_t[1]=bhat2[k-1]*t[3]+bhat2[k]*t[4]+bhat2[k+1]*t[5];
	    lambda_n=n[0]*cstress[k-1]+n[1]*cstress[k]+n[2]*cstress[k+1];
	    lambda_t[0]=t[0]*cstress[k-1]+t[1]*cstress[k]+t[2]*cstress[k+1];
	    lambda_t[1]=t[3]*cstress[k-1]+t[4]*cstress[k]+t[5]*cstress[k+1];
	    } 
	    

	    nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
//	    if(lambda_n>1E-10){
//	      bp=friccoeff*(lambda_n);
//	    }else{
	      bp=bp_old[islavnodeentry-1];
//	    }	    
	    lmhat[0]=lambda_t[0]+constant*up_t[0];
	    lmhat[1]=lambda_t[1]+constant*up_t[1];
	    nlmhat=sqrt(lmhat[0]*lmhat[0]+lmhat[1]*lmhat[1]);
	  
	    ep=bp/nlmhat;
	    if(imodification==1){
	    for(l=0;l<2;l++){
	      for(m=0;m<2;m++){
	        fp[2*m+l]=(1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[m]*lmhat[l];
	      }
	    }
	    }else if(imodification==2){
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[0]*lmhat[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	     fp[3]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[1]*lmhat[1];
	    }else if(imodification==3){
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[1]*lambda_t[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[1];
	     fp[3]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[1]*lambda_t[1];	      
	    }  
	    for(l=0;l<4;l++){mp[l]=0.0;mp2[l]=0.0;mp3[l]=0.0;}
	    for(l=0;l<2;l++){
	        mp[2*l+l]=1.0;
	      for(m=0;m<2;m++){
	        mp[2*m+l]=(mp[2*m+l]-fp[2*m+l])*ep;
	      }
	    }
	  
	    alpha=(lambda_t[0]*lmhat[0]+lambda_t[1]*lmhat[1])/(nlambda_t*nlmhat);
	    delta=min(1.0,nlambda_t/bp);
	    if(imodification==1){
	    if(alpha<0.0){
	     beta=1.0/(1.0-alpha*delta);
	    }else{
	     beta=1.0;
	    }
	    }else if(imodification==2){
	      beta=2/(2-(alpha-1)*delta);
	    }else if(imodification==3){
	      beta=1.0;
	    }  
	    mp2[0]=1.0;
	    mp2[3]=1.0;
	    for(l=0;l<4;l++){
	     mp2[l]=mp2[l]-beta*mp[l]; 
	    }  
	    det=mp2[0]*mp2[3]-mp2[1]*mp2[2];
	    mp3[0]=(1.0/det)*mp2[3];
	    mp3[1]=-(1.0/det)*mp2[1];
	    mp3[2]=-(1.0/det)*mp2[2];
	    mp3[3]=(1.0/det)*mp2[0];
	    
	    vp[0]=(1.0/(nlmhat))*(mp3[0]*lmhat[0]+mp3[1]*lmhat[1]);
	    vp[1]=(1.0/(nlmhat))*(mp3[2]*lmhat[0]+mp3[3]*lmhat[1]);
	    
	    if(nlmhat<1.e-10){
	     vp[0]=0.0;
	     vp[1]=0.0;
	    }  
	  
	    for(l=0;l<3;l++){
	      for(m=0;m<2;m++){
	       tslip[m*3+l]=t[m*3+l]-friccoeff[islavnodeentry-1]*vp[m]*n[l];
	      }
	    } 	    
	    if(debug==1){
	    printf("j=%d k=%d jslavnode=%d kslavnode=%d jcol=%d jrow=%d kact=%d\n",j,k,jslavnodeentry,islavnodeentry,jcol,jrow,islavact[islavnodeentry-1]);
	    printf("\t n= %e %e %e \n",n[0],n[1],n[2]);
	    printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);
	    printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);	   
	    printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);
	    if(jrow==1){
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k],bhat2[k+1],bhat2[k+2],up_t[0],up_t[1],neq[1]);
	    }else{
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k-1],bhat2[k],bhat2[k+1],up_t[0],up_t[1],neq[1]);
	    }
	    printf("\t ep= %e alpha=%e delta=%e beta=%e\n",ep,alpha,delta,beta);
	    printf("\t fp= %e %e %e %e\n",fp[0],fp[1],fp[2],fp[3]);
	    printf("\t mp= %e %e %e %e\n",mp[0],mp[1],mp[2],mp[3]);	
	    printf("\t mp2= %e %e %e %e det=%e\n",mp2[0],mp2[1],mp2[2],mp2[3],det);
	    printf("\t mp3= %e %e %e %e\n",mp3[0],mp3[1],mp3[2],mp3[3]);	
	    printf("\t vp= %e %e\n",vp[0],vp[1]);
	    printf("\t tslip1: %e %e %e\n",tslip[0],tslip[1],tslip[2]);
	    printf("\t tslip2: %e %e %e\n",tslip[3],tslip[4],tslip[5]);
	    }
	    if(islavnodeentry==jslavnodeentry){ //row slip and column slip
	      lp[0]=constant*(mp3[0]-1.0);
	      lp[1]=constant*(mp3[1]);
	      lp[2]=constant*(mp3[2]);	
	      lp[3]=constant*(mp3[3]-1.0);
	      for(l=0;l<3;l++){
	       for(m=0;m<2;m++){
	         ltslip[m*3+l]=lp[2*m]*t[l]+lp[m*2+1]*t[3+l];
	       }
	      } 
	      if(nlmhat<1.E-10){
	      for(l=0;l<3;l++){
	       for(m=0;m<2;m++){
	         ltslip[m*3+l]=0.0;
	       }
	      }		
	      }	
	      if(debug==1){
	      printf("\t lp= %e %e %e %e\n",lp[0],lp[1],lp[2],lp[3]);
	      printf("\t ltslip1:%e %e %e \n",ltslip[0],ltslip[1],ltslip[2]);
	      printf("\t ltslip2:%e %e %e \n",ltslip[3],ltslip[4],ltslip[5]);
	      }
	      if(jrow==2 && jcol==1){ 
	        e1=ad[j];
	        e2=auc[i];
	        e3=auc[i+1];
	        ad[j]=n[0];
	        au[i]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3 )+ltslip[jcol-1];
	        au[i+1]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3) +ltslip[3+jcol-1];
		if(debug==1){
	        printf("\t au= %e %e %e\n",e1,e2,e3);
	        printf("\t au2= %e %e %e\n",ad[j],au[i],au[i+1]);
		}
	      }else if(jrow==1 && jcol==2){
	        e1=auc[i];
	        e2=ad[j];
	        e3=auc[i+1];
	        au[i]=n[1];
	        ad[j]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3) +ltslip[jcol-1];
	        au[i+1]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3) +ltslip[3+jcol-1];
		if(debug==1){
	        printf("\t au= %e %e %e\n",e1,e2,e3);
	        printf("\t au2= %e %e %e\n",au[i],ad[j],au[i+1]);
		}
	      }else if(jrow==1 && jcol==3){
	        e1=auc[i];
	        e2=auc[i+1];
	        e3=ad[j];
	        au[i]=n[2];
	        au[i+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3) +ltslip[jcol-1];
	        ad[j]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3) +ltslip[3+jcol-1];
		if(debug==1){
			      printf("\t au= %e %e %e\n",e1,e2,e3);
	      printf("\t au2= %e %e %e\n",au[i],au[i+1],ad[j]);
		}
	      }else{
	        printf("trafoNTmortar_fric: somethings wrong with diagonal stop!\n");
	        FORTRAN(stop,());		
	      }	      
	      i=i+1;
	    }else if(jrow==1){ //row slip and column stick or normal
	      e1=auc[i];
	      e2=auc[i+1];
	      e3=auc[i+2];
	      au[i]=0.0;
	      au[i+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3);
	      au[i+2]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3);	
	      if(debug==1){
	      printf("\t au= %e %e %e\n",e1,e2,e3);
	      printf("\t au2= %e %e %e\n",au[i],au[i+1],au[i+2]);	      
	      }
	      i=i+2;	      
	    }else{
	        printf("trafoNTmortar_fric: somethings wrong with idofs stop!\n");
	        FORTRAN(stop,());		
	    }		          
	  }   
      }
     }  

    /* changing b due to N and T (normal and tangential
       direction at the slave surface */
    
    for(k=0;k<neq[1];k++){
      if(islavactdof[k]>0){
	islavnodeentry = floor(islavactdof[k]/10.);
	jrow= islavactdof[k]-10*islavnodeentry;
	scal=1.0/bdd[k];
//	scal=1.0;
        for(l=0;l<3;l++){
	   n[l]=slavnor[3*(islavnodeentry-1)+l];
	}
	for(l=0;l<6;l++){
	   t[l]=slavtan[6*(islavnodeentry-1)+l];
	}

        up_t[0]=bhat2[k]*t[0]+bhat2[k+1]*t[1]+bhat2[k+2]*t[2];
	up_t[1]=bhat2[k]*t[3]+bhat2[k+1]*t[4]+bhat2[k+2]*t[5];
	up_n=bhat2[k]*n[0]+bhat2[k+1]*n[1]+bhat2[k+2]*n[2];
	lambda_n=n[0]*cstress[k]+n[1]*cstress[k+1]+n[2]*cstress[k+2];
	lambda_t[0]=t[0]*cstress[k]+t[1]*cstress[k+1]+t[2]*cstress[k+2];
	lambda_t[1]=t[3]*cstress[k]+t[4]*cstress[k+1]+t[5]*cstress[k+2];
	nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
//	if(lambda_n>1E-10){
//	  bp=friccoeff*(lambda_n);
//	}else{
	  bp=bp_old[islavnodeentry-1];
//	}


	if(islavact[islavnodeentry-1]==1 && jrow==1){ //stick
	  for(l=0;l<3;l++){
	    for(m=0;m<2;m++){
	      rstick[m*3+l]=friccoeff[islavnodeentry-1]*up_t[m]*n[l];
	    }
	  }  
	  e1=bhat[k];
	  e2=bhat[k+1];
	  e3=bhat[k+2];
         /* right side if solving Ku=f,see hüeber */	  
//	  b[k]=gap[islavnodeentry-1];
//	  b[k+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3)-bp*up_t[0];
//	  b[k+2]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3)-bp*up_t[1];
         /* right side if solving K du=f */
	  b[k]=gap[islavnodeentry-1]-up_n;
	  b[k+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3);
	  b[k+2]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3);
          if(debug==1){
	    printf("k %d %d %d jrow %d activ %d",k,k+1,k+2,jrow,islavact[islavnodeentry-1]);
	    printf("\t bp=%e bp_alt=%e \n",friccoeff[islavnodeentry-1]*(lambda_n),bp);
	    printf("\t n= %e %e %e \n",n[0],n[1],n[2]);
	    printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);
	    printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);
	    printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k],bhat2[k+1],bhat2[k+2],up_t[0],up_t[1],neq[1]);
	    printf("\t rstick1: %e %e %e\n",rstick[0],rstick[1],rstick[2]);
	    printf("\t rstick2: %e %e %e\n",rstick[3],rstick[4],rstick[5]);
	  printf("k=%d activ=%d\n",k,islavact[islavnodeentry-1]);
	  printf("\t b_old= %e %e %e \n",e1,e2,e3);
	  printf("\t b_new= %e %e %e \n",b[k],b[k+1],b[k+2]);
	  }
	  k=k+2;
	  
	}else if(islavact[islavnodeentry-1]==2 && jrow==1){ //slip
	  lmhat[0]=lambda_t[0]+constant*up_t[0];
	  lmhat[1]=lambda_t[1]+constant*up_t[1];
	  nlmhat=sqrt(lmhat[0]*lmhat[0]+lmhat[1]*lmhat[1]);
	  
	  ep=bp/nlmhat;
	  if(imodification==1){
	    for(l=0;l<2;l++){
	      for(m=0;m<2;m++){
	        fp[2*m+l]=(1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[m]*lmhat[l];
	      }
	    }
	    }else if(imodification==2){
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[0]*lmhat[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	     fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	     fp[3]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[1]*lmhat[1];
	    }else if(imodification==3){
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[1]*lambda_t[0];
	     fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[1];
	     fp[3]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[1]*lambda_t[1];	      
	    } 
	  for(l=0;l<4;l++){mp[l]=0.0;mp2[l]=0.0;}
	  for(l=0;l<2;l++){
	      mp[2*l+l]=1.0;
	    for(m=0;m<2;m++){
	      mp[2*m+l]=(mp[2*m+l]-fp[2*m+l])*ep;
	    }
	  }
	  
	  alpha=(lambda_t[0]*lmhat[0]+lambda_t[1]*lmhat[1])/(nlambda_t*nlmhat);
	  delta=min(1.0,nlambda_t/bp);
	  if(imodification==1){
	    if(alpha<0.0){
	     beta=1.0/(1.0-alpha*delta);
	    }else{
	     beta=1.0;
	    }
	  }else if(imodification==2){
	      beta=2/(2-(alpha-1)*delta);
	  }else if(imodification==3){
	      beta=1.0;
	  }  

	  mp2[0]=1.0;
	  mp2[1]=0.0;
	  mp2[3]=1.0;
	  mp2[2]=0.0;
	  for(l=0;l<4;l++){
	     mp2[l]=mp2[l]-beta*mp[l]; 
	  } 
	  det=mp2[0]*mp2[3]-mp2[1]*mp2[2];
	  mp3[0]=(1.0/det)*mp2[3];
	  mp3[1]=-(1.0/det)*mp2[1];
	  mp3[2]=-(1.0/det)*mp2[2];
	  mp3[3]=(1.0/det)*mp2[0];
	    
	  vp[0]=(1.0/(nlmhat))*(mp3[0]*lmhat[0]+mp3[1]*lmhat[1]);
	  vp[1]=(1.0/(nlmhat))*(mp3[2]*lmhat[0]+mp3[3]*lmhat[1]);
	  
	  for(l=0;l<3;l++){
	    for(m=0;m<2;m++){
	     tslip[m*3+l]=t[m*3+l]-friccoeff[islavnodeentry-1]*vp[m]*n[l];
	    }
	  } 
//	  det=fp[0]*fp[3]-fp[1]*fp[2];
	  hp[0]=ep*(fp[0]*lmhat[0]+fp[1]*lmhat[1]);
	  hp[1]=ep*(fp[2]*lmhat[0]+fp[3]*lmhat[1]);
	  
	  rp[0]=-(mp3[0]*hp[0]+mp3[1]*hp[1]);
	  rp[1]=-(mp3[2]*hp[0]+mp3[3]*hp[1]);
	  
	  rphat[0]=rp[0]+bp*vp[0];
	  rphat[1]=rp[1]+bp*vp[1];
	  
	  if(nlmhat<1e-10){
	     rphat[0]=0.0;
	     rphat[1]=0.0;
	  }
	  /* solving for du -> need more arrays*/    
	  lp[0]=constant*(mp3[0]-1-0);
	  lp[1]=constant*(mp3[1]);
	  lp[2]=constant*(mp3[2]);	
	  lp[3]=constant*(mp3[3]-1.0);

	  ltu[0]=lp[0]*up_t[0]+lp[1]*up_t[1];
	  ltu[1]=lp[2]*up_t[0]+lp[3]*up_t[1];
	  
	  e1=bhat[k];
	  e2=bhat[k+1];
	  e3=bhat[k+2];
	  
	  /* right side if solving Ku=f, see hüeber*/
//	  b[k]=gap[islavnodeentry-1];
//	  b[k+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3)+rphat[0];
//	  b[k+2]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3)+rphat[1];

	  /* right side if solving K du=f*/ 
	  b[k]=gap[islavnodeentry-1]-up_n;
	  b[k+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3)+rphat[0]-ltu[0];
	  b[k+2]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3)+rphat[1]-ltu[1];	  
         if(debug==1){
	  printf("k=%d activ=%d\n",k,islavact[islavnodeentry-1]);
	    printf("\t n= %e %e %e \n",n[0],n[1],n[2]);
	    printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);
	    printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);	   
	    printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);
	    printf("\t lmhat= %e %e  nlmhat=%e \n",lmhat[0],lmhat[1],nlmhat);	    
	    printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,bhat2[k],bhat2[k+1],bhat2[k+2],up_t[0],up_t[1],neq[1]);
	    printf("\t ep= %e alpha=%e delta=%e beta=%e\n",ep,alpha,delta,beta);
	    printf("\t fp= %e %e %e %e\n",fp[0],fp[1],fp[2],fp[3]);
	    printf("\t mp= %e %e %e %e\n",mp[0],mp[1],mp[2],mp[3]);	
	    printf("\t mp2= %e %e %e %e det=%e\n",mp2[0],mp2[1],mp2[2],mp2[3],det);
	    printf("\t mp3= %e %e %e %e\n",mp3[0],mp3[1],mp3[2],mp3[3]);	    
	    printf("\t vp= %e %e\n",vp[0],vp[1]);
	    printf("\t lp= %e * (%e %e %e %e)\n",constant,lp[0]/constant,lp[1]/constant,lp[2]/constant,lp[3]/constant);
	    printf("\t tslip1: %e %e %e\n",tslip[0],tslip[1],tslip[2]);
	    printf("\t tslip2: %e %e %e\n",tslip[3],tslip[4],tslip[5]);	
	    
	    printf("\t hp= %e %e\n",hp[0],hp[1]);
	    printf("\t rp= %e %e \n",rp[0],rp[1]);
	    printf("\t rphat= %e %e \n",rphat[0],rphat[1]);
	    printf("\t ltu= %e %e ehat1=%e ehat2=%e\n",ltu[0],ltu[1],scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3),scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3));
	  printf("\t b_old= %e %e %e \n",e1,e2,e3);
	  printf("\t b_new= %e %e %e \n",b[k],b[k+1],b[k+2]);
	  }
	  k=k+2;	  
	}else{
	  printf("trafoNTmortar_fric: somethings wrong stop!\n");
	  FORTRAN(stop,());
	}  

      }
    }
/*    for (i=0;i<neq[1];i++){
  	printf("TNT:bhat(%d)= %e \t b(%d)0 %e \n",i,bhat[i],i,b[i]);
   }*/
    free(lmhat);
    free(lambda_t);
    free(fp);
    free(mp);
    free(lp);
    free(ltslip);
    free(rp);
    free(rphat);
    free(hp);
    free(n);
    free(t);
    free(up_t);
    free(vp);
    free(rstick);
    free(tslip);
    free(ltu);
    
    free(bhat2);

}