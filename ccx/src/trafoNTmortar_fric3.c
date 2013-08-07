/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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
   * 
   * author: Saskia Sitzmann
 * @param [out] au
 * @param [out] b
*/
void trafoNTmortar_fric3(int *neq,int *nzs, int *islavactdof,int *islavact,int *nslavnode, int *nmastnode, int *ncone, 
        double *ad, double *au, double *b, int *irow, int *jq,
        int *nzsc, double *auc,
        double *adc, int *irowc, int *jqc,
        double *gap, double *bdd, double *auqdt, int *irowqdt,
        int *jqqdt, int *nzsqdt, int *nzlc,double *slavnor, double *slavtan,double *bhat,
	double *aubd, int *irowbd, int *jqbd, double *vold, double *cstress,double *bp_old,int *nactdof,
	int *islavnode, int *ntie, int *mi,int *nk,
        int *nboun,int *ndirboun,int *nodeboun,double *xboun,
        int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
        int *ikboun,int *ilboun,int *ikmpc,int *ilmpc,
        int *nslavspc,int *islavspc,int *nsspc,int *nslavmpc,int *islavmpc,int *nsmpc,
        int *nmastspc,int *imastspc,int *nmspc,int *nmastmpc,int *imastmpc,int *nmmpc,
        double *Bd,double *Dd,int *jqb,int *irowb, int *nzsbd2,char *tieset,
	int *islavactdoftie,int *nelcon, double  *elcon, double *tietol,int *ncmat_,int *ntmat_,
	double *plicon,int *nplicon, int *npmat_){

    int i,j,jj,k,kk,l,m,debug,idof1,idof2,idof3,imodification,
        jrow,jcol,islavnodeentry,jslavnodeentry,mt=mi[1]+1,nodes,dim,node,islavk2,
        dirblock,dirind,dirdep,ist,node1,id,dir,index,perturbed=1,derivmode,regmode;

    double t1,t2,e1,e2,e3,*rstick=NULL,
         bp, *t=NULL,*n=NULL, *up_t=NULL,
         *tslip=NULL,*vp=NULL,*mp=NULL,*fp=NULL,ep,constant=1.E10,*lmhat=NULL,nlmhat,
         alpha,beta,delta,lambda_n,*lambda_t=NULL,nlambda_t,det,*rp=NULL,*hp=NULL,*rphat=NULL,
         *lp=NULL,*ltslip=NULL,*mp2=NULL,*mp3=NULL,*u_old=NULL,scal,*ltu=NULL,up_n,*n2=NULL,
         coefdep,nt2,nt1,nn,that[6],dist,n11,n22,aninvloc,gnc,dgnc,mu,p0,beta_e;
    
    debug=0;
    imodification=1;
    
    u_old=NNEW(double,3*nslavnode[*ntie]);
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
    n2=NNEW(double,3);
    t=NNEW(double,6);
    up_t=NNEW(double,2);
    vp=NNEW(double,2);
    rstick=NNEW(double, 2*3);
    tslip=NNEW(double, 2*3);
    ltu=NNEW(double,2);
    
    /** get uhat_k-1 **/    
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){	
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	  nodes=islavnode[j];	    
	  u_old[j*3]=vold[mt*(nodes)-3];	    
	  u_old[j*3+1]=vold[mt*(nodes)-2];	    	    
	  u_old[j*3+2]=vold[mt*(nodes)-1];	    
          if(debug==1) printf("j %d node %d uold %e %e %e\n",j,nodes, u_old[j*3], u_old[j*3+1], u_old[j*3+2] );		
	}	
      }
    }
    for(i=0;i<*nk; i++){     
      for (j=jqb[i]-1;j<jqb[i+1]-1;j++){	
	nodes=irowb[j];	
	islavk2=0;	
	for(jj=0;jj<*ntie;jj++){  
	  if(tieset[jj*(81*3)+80]=='C'){
	  dim=nslavnode[jj+1]-nslavnode[jj];	          
	  FORTRAN(nident,(&islavnode[nslavnode[jj]], &nodes,&dim, &id));	 
	  if(id>0 && islavnode[nslavnode[jj]+id-1]==nodes){islavk2=nslavnode[jj]+id-1;}	
	  }
	}        
        for(k=0;k<3;k++){	       
	  idof1=nactdof[mt*(i+1)-3+k]-1;	       
	  idof2=nactdof[mt*nodes-3+k]-1;	       
	  jj=-1;	       
	  if(idof1>-1){	       
	    for(kk=jqqdt[idof1]-1;kk<jqqdt[idof1+1]-1;kk++){		
	      if(irowqdt[kk]-1==idof2){jj=kk;} 	       
	    }	       
	  }	       	                           
          if(jj>-1) u_old[(islavk2)*3+k]=u_old[(islavk2)*3+k]-auqdt[jj]*vold[mt*(i+1)-3+k];	       
        }
      }
    }
    /** loop over columns **/     
    for(j=0;j<neq[1];j++){ ///loop over columns         
      for(i=jq[j]-1;i<jq[j+1]-1;i++){ ///loop over rows	  
	jslavnodeentry = floor(islavactdof[j]/10.);	  
	jcol= islavactdof[j]-10*jslavnodeentry;	  
	k=irow[i]-1;          
	islavnodeentry = floor(islavactdof[k]/10.);	  
	jrow= islavactdof[k]-10*islavnodeentry;	  
	if(islavactdof[k]>0 && (islavact[islavnodeentry-1]==1||islavact[islavnodeentry-1]==2)){	   
	  node=islavnode[islavnodeentry-1];	   
	  idof1=nactdof[mt*node-3]-1;	   
	  idof2=nactdof[mt*node-2]-1;	   
	  idof3=nactdof[mt*node-1]-1;	    	   
	  scal= 1.0/Dd[node-1];
          /** get normal and tangetials **/           
	  for(l=0;l<3;l++){	     
	    n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	    n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
	  }
	  for(l=0;l<6;l++){	     
	    t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	    that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
	  }
	  /** modify then due to SPC's/MPC's **/
	  /** check whether a direction is blocked with an SPC/MPC **/	   
	 dirblock=0; 	   
	 for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){	      
	   if(islavmpc[2*jj+1]==1){		 
	     dirblock=1;	      
	   }	   
	 } 
//	 debug=0;
	   
	 if(debug==1){           
	   printf("trafoNT: j:%d, k:%d, slavnodeenty,j,i: %d, %d jcol:%d, jrow:%d\n ",j,k, jslavnodeentry,islavnodeentry,jcol,jrow);	   
	   printf("node %d idof %d %d %d\n",node,idof1,idof2,idof3);	  	   
	 }
	 /** handle SPCs and MPCs -> modify tangentials and normal **/ 	   
	 for(jj=nslavspc[2*(islavnodeentry-1)];jj<nslavspc[2*(islavnodeentry-1)+1];jj++){	    
	   t[0]=0;t[1]=0;t[2]=0;	    
	   ist=islavspc[2*jj];            
	   dir=ndirboun[ist-1];	    
	   node1=nodeboun[ist-1];	    
	   if(debug==1){printf("jj %d ist %d dir %d node %d\n",jj,ist,dir,node1);}
           t[dir-1]=1.0;	    
	   dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];	    
	   n[0]=n[0]-dist*t[0];	    
	   n[1]=n[1]-dist*t[1];	    
	   n[2]=n[2]-dist*t[2];	    
	   nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	    
	   n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;	   
	 }
	 for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){	    	   
	   t[0]=0;t[1]=0;t[2]=0;             
	   ist=islavmpc[2*jj];            
	   dirdep=nodempc[3*(ist-1)+1];            
	   coefdep=coefmpc[ist-1];            
	   index=nodempc[3*(ist-1)+2];	    
	   node1=nodempc[3*(ist-1)]; 	    
	   t[dirdep-1]=coefdep;	    
	   if(debug==1){printf("jj %d ist %d dir %d node %d coef %e\n",jj,ist,dirdep,node1,coefdep);}	                
            while(index!=0){
               dirind=nodempc[3*(index-1)+1];
	       node1=nodempc[3*(index-1)]; 
	       if(node1==node)t[dirind-1]=coefmpc[index-1];
               index=nodempc[3*(index-1)+2];
	       if(debug==1){printf("index %d \n",index);}
	    }
            nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
	    t[0]=t[0]/nt1;t[1]=t[1]/nt1;t[2]=t[2]/nt1;
	    dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];
	    n[0]=n[0]-dist*t[0];
	    n[1]=n[1]-dist*t[1];
	    n[2]=n[2]-dist*t[2];
	    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	    n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;	    
	   
	 }
	 /** get second tangential and modify it due to SPC's/MPC's**/
	if(nslavspc[2*(islavnodeentry-1)+1]-nslavspc[2*(islavnodeentry-1)]>0 ||nslavmpc[2*(islavnodeentry-1)+1]-nslavmpc[2*(islavnodeentry-1)]>0){	       	    
	  /** caution: t[3-5]= hat(t2)=n x hat(t1) **/	    
	  t[3]=n[1]*t[2]-n[2]*t[1];            
	  t[4]=n[2]*t[0]-n[0]*t[2];            
	  t[5]=n[0]*t[1]-n[1]*t[0];	    
	  for(kk=0;kk<6;kk++){that[kk]=t[kk];}	    
	    for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){             
	      ist=islavmpc[2*jj];             
	      dirdep=nodempc[3*(ist-1)+1];             
	      coefdep=coefmpc[ist-1];             
	      index=nodempc[3*(ist-1)+2];	     
	      node1=nodempc[3*(ist-1)];                
	      while(index!=0){               
		dirind=nodempc[3*(index-1)+1];	       
		node1=nodempc[3*(index-1)]; 	       
		if(node1==node){		 
		  t[dirind-1+3]=t[dirind-1+3]-coefmpc[index-1]*t[dirdep-1+3]/coefdep;		 
		  n[dirind-1]=n[dirind-1]-coefmpc[index-1]*n[dirdep-1]/coefdep;		 	       
		}
		index=nodempc[3*(index-1)+2];	    
	      }
	      t[dirdep-1+3]=0.0;	    
	      n[dirdep-1]=0.0;	    	   
	    }	     
	    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	
	    nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
	    nt2=sqrt(t[3]*t[3]+t[4]*t[4]+t[5]*t[5]);
	    if(debug==1){
	    printf("n= %e %e %e \n",n[0],n[1],n[2]);
	    printf("t1= %e %e %e \n",t[0],t[1],t[2]);
	    printf("t2= %e %e %e \n",t[3],t[4],t[5]);
	    }
	    /** caution: t[3-5]= tilde(t2), hat(t2) with modifications due to SPC/MPCs
	        *        that[0-5] =hat(t1) ,hat(t2)
	        *	n[0-2]= n*hat(D)/D in case of directinal blocking
		*       n[0-2]= n otherwise 
		*       n2[0-2]= original n
		**/
	   }	    	    
	    /** calculate fields needed for Coulomb friction **/ 
	   up_n=   u_old[(islavnodeentry-1)*3]*n2[0]+u_old[(islavnodeentry-1)*3+1]*n2[1]+u_old[(islavnodeentry-1)*3+2]*n2[2];
           up_t[0]=u_old[(islavnodeentry-1)*3]*that[0]+u_old[(islavnodeentry-1)*3+1]*that[1]+u_old[(islavnodeentry-1)*3+2]*that[2];
	   up_t[1]=u_old[(islavnodeentry-1)*3]*that[3]+u_old[(islavnodeentry-1)*3+1]*that[4]+u_old[(islavnodeentry-1)*3+2]*that[5];
	   lambda_n=n2[0]*cstress[(islavnodeentry-1)*3]+n2[1]*cstress[(islavnodeentry-1)*3+1]+n2[2]*cstress[(islavnodeentry-1)*3+2];
	   lambda_t[0]=that[0]*cstress[(islavnodeentry-1)*3]+that[1]*cstress[(islavnodeentry-1)*3+1]+that[2]*cstress[(islavnodeentry-1)*3+2];
	   lambda_t[1]=that[3]*cstress[(islavnodeentry-1)*3]+that[4]*cstress[(islavnodeentry-1)*3+1]+that[5]*cstress[(islavnodeentry-1)*3+2];
	   nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
	   bp=bp_old[islavnodeentry-1];
	   /// perturebed lagrange
           FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&p0,&beta_e,tietol,elcon,&islavactdoftie[islavnodeentry-1],ncmat_,ntmat_));
	   derivmode=1;
	   FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,&aninvloc,&p0,&beta_e,elcon,nelcon,&islavactdoftie[islavnodeentry-1],ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
	   
	   if(islavact[islavnodeentry-1]==1){ /// case: row stick	    
	     for(l=0;l<3;l++){	     
	       for(m=0;m<2;m++){	      
		 rstick[m*3+l]=mu*up_t[m]*n2[l];	     
	       }	    
	     }
	     if(debug==1){	    
	       printf("j=%d k=%d jslavnode=%d kslavnode=%d jcol=%d jrow=%d kact=%d\n",j,k,jslavnodeentry,islavnodeentry,jcol,jrow,islavact[islavnodeentry-1]);	    
//	       printf("\t n= %e %e %e \n",n[0],n[1],n[2]);	    
//	       printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);	    
//	       printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);	    
	       printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);	   
	       printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,u_old[(islavnodeentry-1)*3],u_old[(islavnodeentry-1)*3+1],u_old[(islavnodeentry-1)*3+2],up_t[0],up_t[1],neq[1]);	    
	       printf("\t rstick1: %e %e %e\n",rstick[0],rstick[1],rstick[2]);	    
	       printf("\t rstick2: %e %e %e\n",rstick[3],rstick[4],rstick[5]);	    
	     }
	     if(islavnodeentry==jslavnodeentry){ 
	       /// row stick and column stick (row and colum belong to same slave node)
	       if(idof1>-1 && idof2>-1 && idof3>-1){ /// 3D	      
		 if(jrow==2 && jcol==1){ 	        
		   e1=adc[j];	        
		   e2=auc[i];	        
		   e3=auc[i+1];
		   if(perturbed==1){
		     ad[j]=n[0]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     ad[j]=n[0];
		   }
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
		   if(perturbed==1){
		     au[i]=n[1]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     au[i]=n[1];
		   }		   
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
		   if(perturbed==1){
		     au[i]=n[2]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     au[i]=n[2];
		   }			   
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
	       }else if(idof2>-1 && idof3>-1){ ///2D auf 3D
		 if(jrow==3){                     	        
		   e1=adc[j];                     	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                     	        
		   e2=adc[j];		   		              
		 }
		 if(nn<0.5){                                
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 }  
		 t1=rstick[4];			
		 t2=rstick[5];
		 n11=n2[1];
		 n22=n2[2];
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}							       
	         if(jrow==3){					
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }
		   au[i]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{				
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }		   
		   ad[j]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 } 	       	     
	       }else if(idof1>-1 && idof3>-1){ ///2D auf 3D		              
		 if(jrow==3){                     	        
		   e1=adc[j];                     	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                     	        
		   e2=adc[j];		   		              
		 } 
		 if(nn<0.3){                                
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 } 
		 t1=rstick[3];				
		 t2=rstick[5];
		 n11=n2[0];
		 n22=n2[2];
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}				
	         if(jrow==3){					
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }
		   au[i]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{					
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }		   
		   ad[j]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 }	       	     
	       }else if(idof1>-1 && idof2>-1){ ///2D auf 3D		              
		 if(jrow==2){                    	        
		   e1=adc[j];                     	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                     	        
		   e2=adc[j];		   		              
		 } 
		 if(nn<0.3){                               
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 } 		        
		 t1=rstick[3];				
		 t2=rstick[4];
		 n11=n2[0];
		 n22=n2[1];
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}
	         if(jrow==2){						
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }		   
		   au[i]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{				
		   au[i]=n[jcol-1];	
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }		   
		   ad[j]=scal*(t1*e1+t2*e2)-bp*t[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 } 			       	       	     
	       }else {                    
		 printf("trafoNT: Idofs ordered the wrong way...stop");				
		 FORTRAN(stop,());	     
	       }   	    
	     }else{ /// row stick and column normal or slip	     
	       if(idof1>-1 && idof2>-1 && idof3>-1){	      
		 e1=auc[i];	      
		 e2=auc[i+1];	      
		 e3=auc[i+2];	      	
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);		   
		 }else{		     
		   au[i]=0.0;		   
		 }
		 au[i+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3);	      
		 au[i+2]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3);	                   
		 if(debug==1){	      
		   printf("\t au= %e %e %e\n",e1,e2,e3);	      
		   printf("\t au= %e %e %e\n",au[i],au[i+1],au[i+2]);	      
		 }
		 i=i+2;	     
	       }else if(idof2>-1 && idof3>-1){ ///2D auf 3D			                        	
		 e1=auc[i];                     		
		 e2=auc[i+1];                     		
		 t1=rstick[4];				
		 t2=rstick[5];
		 n11=n2[1];
		 n22=n2[2];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 if(debug==1){				
		   printf("\t nn %e nt %e %e\n",nn,nt1,nt2);				
		   printf("\t t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],au[i+1]);				
		 }
		 i=i+1;	       	      	     	       
	       }else if(idof1>-1 && idof3>-1){ ///2D auf 3D			                            		
		 e1=auc[i];                     		
		 e2=auc[i+1];                     		
		 t1=rstick[3];				
		 t2=rstick[5];
		 n11=n2[0];
		 n22=n2[2];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 i=i+1;	       	      	       	    
	       }else if(idof1>-1 && idof2>-1){ ///2D auf 3D			                             		
		 e1=auc[i];                     		
		 e2=auc[i+1];                     		
		 t1=rstick[3];				
		 t2=rstick[4];	
		 n11=n2[1];
		 n22=n2[2];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 i=i+1;      	     
	       }else{ ///1D auf 3D
		 e1=auc[i];
		 if(idof1>-1){		 
		   n11=n2[0];
		 }else if(idof2>-1){
		   n11=n2[1];
		 }else{
		   n11=n2[2];
		 }
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1);		   
		 }else{		     
		   au[i]=0.0;		   
		 }	     
	       } 	      	    
	     }	    	   	   
	   }else if(islavact[islavnodeentry-1]==2){ ///row slip 	     
	     /** fields needed for slip case **/
	     lmhat[0]=lambda_t[0]+constant*up_t[0];
	     lmhat[1]=lambda_t[1]+constant*up_t[1];
	     nlmhat=sqrt(lmhat[0]*lmhat[0]+lmhat[1]*lmhat[1]);	  
	     ep=bp/nlmhat;
	    
	     if(mu>1.E-10){ /// contact tie with friction	     	     
	       if(imodification==1){	      
		 for(l=0;l<2;l++){	       
		   for(m=0;m<2;m++){	        
		     fp[2*m+l]=(1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[m]*lmhat[l];	       
		   }	      
		 }	     
	       }else if(imodification==2){	      
		 fp[0]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[0]*lmhat[0];	      
		 fp[1]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);	      
		 fp[2]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);	      
		 fp[3]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[1]*lmhat[1];	     
	       }else if(imodification==3){	      
		 fp[0]= (1.0/(max(bp,nlambda_t)*nlambda_t))*lambda_t[0]*lambda_t[0];	      
		 fp[1]= (1.0/(max(bp,nlambda_t)*nlambda_t))*lambda_t[1]*lambda_t[0];	      
		 fp[2]= (1.0/(max(bp,nlambda_t)*nlambda_t))*lambda_t[0]*lambda_t[1];	      
		 fp[3]= (1.0/(max(bp,nlambda_t)*nlambda_t))*lambda_t[1]*lambda_t[1];	      	     
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
		   tslip[m*3+l]=that[m*3+l]-mu*vp[m]*n2[l];	       
		 }	     
	       } 	    
	     }else{ /// contact tie without friction	     
	       for(l=0;l<3;l++){	       
		 for(m=0;m<2;m++){	        
		   tslip[m*3+l]=t[m*3+l];	       
		 }	     
	       }
//	       scal=1.0;
	    
	     }  
	     if(debug==1){	    
	       printf("j=%d k=%d jslavnode=%d kslavnode=%d jcol=%d jrow=%d kact=%d\n",j,k,jslavnodeentry,islavnodeentry,jcol,jrow,islavact[islavnodeentry-1]);	    
//	       printf("\t n= %e %e %e \n",n[0],n[1],n[2]);	    	   
//	       printf("\t t1= %e %e %e \n",t[0],t[1],t[2]);	    
//	       printf("\t t2= %e %e %e \n",t[3],t[4],t[5]);	   	    
	       printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);	    
	       printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,u_old[(islavnodeentry-1)*3],u_old[(islavnodeentry-1)*3+1],u_old[(islavnodeentry-1)*3+2],up_t[0],up_t[1],neq[1]);	    
	       if(mu>1.E-10){	    
		 printf("\t ep= %e alpha=%e delta=%e beta=%e\n",ep,alpha,delta,beta);	    
		 printf("\t fp= %e %e %e %e\n",fp[0],fp[1],fp[2],fp[3]);	    
		 printf("\t mp= %e %e %e %e\n",mp[0],mp[1],mp[2],mp[3]);		    
		 printf("\t mp2= %e %e %e %e det=%e\n",mp2[0],mp2[1],mp2[2],mp2[3],det);	    
		 printf("\t mp3= %e %e %e %e\n",mp3[0],mp3[1],mp3[2],mp3[3]);		    
		 printf("\t vp= %e %e\n",vp[0],vp[1]);	    
	       }
	       printf("\t tslip1: %e %e %e\n",tslip[0],tslip[1],tslip[2]);	    
	       printf("\t tslip2: %e %e %e\n",tslip[3],tslip[4],tslip[5]);	    
	     }
	     if(islavnodeentry==jslavnodeentry){ ///row slip and column slip	      
	       lp[0]=constant*(mp3[0]-1.0);	      
	       lp[1]=constant*(mp3[1]);	      
	       lp[2]=constant*(mp3[2]);		      
	       lp[3]=constant*(mp3[3]-1.0);	      
	       for(l=0;l<3;l++){	       
		 for(m=0;m<2;m++){	         
		   ltslip[m*3+l]=lp[2*m]*t[l]+lp[m*2+1]*t[3+l];	       
		 }	      
	       } 
	       if(mu<1.E-10){ ///contact tie without friction	      
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
	       if(idof1>-1 && idof2>-1 && idof3>-1){	     
		 if(jrow==2 && jcol==1){ 	        
		   e1=ad[j];	        
		   e2=auc[i];	        
		   e3=auc[i+1];	        
		   if(perturbed==1){
		     ad[j]=n[0]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     ad[j]=n[0];
		   }		   
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
		   if(perturbed==1){
		     au[i]=n[1]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     au[i]=n[1];
		   }		   
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
		   if(perturbed==1){
		     au[i]=n[2]+scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
		   }else{
		     au[i]=n[2];
		   }		   
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
	       }else if(idof2>-1 && idof3>-1){ ///2D auf 3D		              
		 if(jrow==3){                     	        
		   e1=adc[j];                     	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                    	        
		   e2=adc[j];		   		              
		 } 
		 if(nn<0.5){                                
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 }  
                 t1=tslip[4];				
		 t2=tslip[5];
		 n11=n2[1];
		 n22=n2[2];
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}				
	         if(jrow==3){					
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }		   
		   au[i]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{				
		   au[i]=n[jcol-1];
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }
		   ad[j]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==19){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 } 	       	     
	       }else if(idof1>-1 && idof3>-1){ ///2D auf 3D		              
		 if(jrow==3){                     	        
		   e1=adc[j];                    	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                     	        
		   e2=adc[j];		   		              
		 } 
		 if(nn<0.3){                                
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 } 		              
                 t1=tslip[3];				
		 t2=tslip[5];	
		 n11=n2[0];
		 n22=n2[2];		 
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}				
	         if(jrow==3){					
		   ad[j]=n[jcol-1];
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }		   
		   au[i]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{				
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }		   
		   ad[j]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 }	       	     	       
	       }else if(idof1>-1 && idof2>-1){ ///2D auf 3D		              
		 if(jrow==2){                     	        
		   e1=adc[j];                     	        
		   e2=auc[i];		              
		 }else{                     	        
		   e1=auc[i];                     	        
		   e2=adc[j];		   		              
		 } 
		 if(nn<0.3){                                
		   printf("trafoNT: normal direction can't be blocked...stop");				
		   FORTRAN(stop,());				  				
		 } 		               
                 t1=tslip[3];				
		 t2=tslip[4];
		 n11=n2[0];
		 n22=n2[1];		 
		 if(debug==1){printf("nn %e nt %e %e\n",nn,nt1,nt2);}			       
		 if(jrow==2){					
		   if(perturbed==1){
		     ad[j]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     ad[j]=n[jcol-1];
		   }		   
		   au[i]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e %e \n",t1,t2,e1,e2,ad[j],au[i]);}			       
		 }else{				
		   if(perturbed==1){
		     au[i]=n[jcol-1]+scal*(dgnc)*(n11*e1+n22*e2);
		   }else{
		     au[i]=n[jcol-1];
		   }		   
		   ad[j]=scal*(t1*e1+t2*e2)+ltslip[3+jcol-1];				
		   if(debug==1){printf("t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],ad[j]);}				 			       
		 } 			       	       	     
	       }else {                   
		 printf("trafoNT: Idofs ordered the wrong way...stop");				
		 FORTRAN(stop,());	     
	       }   	    
	     }else{ ///row slip and column stick or normal	     
	       if(idof1>-1 && idof2>-1 && idof3>-1){	      	      
		 e1=auc[i];	      
		 e2=auc[i+1];	      
		 e3=auc[i+2];	      
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3);	      
		 au[i+2]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3);		      
		 if(debug==1){	      
		   printf("\t au= %e %e %e\n",e1,e2,e3);	      
		   printf("\t au2= %e %e %e\n",au[i],au[i+1],au[i+2]);	      	      
		 }
	         i=i+2;	     	       
	       }else if(idof2>-1 && idof3>-1){ ///2D auf 3D			                             		
		 e1=auc[i];                     		
		 e2=auc[i+1];                    		
		 t1=tslip[4];				
		 t2=tslip[5];
		 n11=n2[1];
		 n22=n2[2];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 if(debug==1){				
		   printf("\t nn %e nt %e %e\n",nn,nt1,nt2);				
		   printf("\t t %e %e  e %e %e a(i+1)= %e  %e \n",t1,t2,e1,e2,au[i],au[i+1]);				
		 }
		 i=i+1;	       	      	     	       
	       }else if(idof1>-1 && idof3>-1){ ///2D auf 3D			                            		
		 e1=auc[i];                     		
		 e2=auc[i+1];                     		
		 t1=tslip[3];				
		 t2=tslip[5];
		 n11=n2[0];
		 n22=n2[2];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 i=i+1;	       	      	       	     
	       }else if(idof1>-1 && idof2>-1){ ///2D auf 3D			                             		
		 e1=auc[i];                     		
		 e2=auc[i+1];                     		
		 t1=tslip[3];				
		 t2=tslip[4];
		 n11=n2[0];
		 n22=n2[1];
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1+n22*e2);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 
		 au[i+1]=scal*(t1*e1+t2*e2);				
		 i=i+1;      	     
	       }else {
		 e1=auc[i];
		 if(idof1>-1){		 
		   n11=n2[0];
		 }else if(idof2>-1){
		   n11=n2[1];
		 }else{
		   n11=n2[2];
		 }
		 if(perturbed==1){		     
		   au[i]=scal*(dgnc)*(n11*e1);		   
		 }else{		     
		   au[i]=0.0;		   
		 }		 	     
	       }	      	    
	     }	   	   
	   }  	  	
	}      
      }     
    }  

    /** changing b due to N and T (normal and tangential
       direction at the slave surface **/
    
    for(k=0;k<neq[1];k++){      
      if(islavactdof[k]>0){	
	islavnodeentry = floor(islavactdof[k]/10.);
	jrow= islavactdof[k]-10*islavnodeentry;
	node=islavnode[islavnodeentry-1];
	idof1=nactdof[mt*node-3]-1;
	idof2=nactdof[mt*node-2]-1;
	idof3=nactdof[mt*node-1]-1;		   
	scal= 1.0/Dd[node-1];	
//        debug=0;
        /** get normal and tangetials **/
        for(l=0;l<3;l++){	     
	  n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	  n2[l]=slavnor[3*(islavnodeentry-1)+l];	     
	}
	for(l=0;l<6;l++){	     
	  t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	  that[l]=slavtan[6*(islavnodeentry-1)+l];	     
	}	    
	/** modify then due to SPC's/MPC's **/
	/** check whether a direction is blocked with an SPC/MPC **/
	dirblock=0;	     
	for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){	       
	  if(islavmpc[2*jj+1]==1){		 
	    dirblock=1;	       
	  }	     
	}  
	if(debug==1){	   	   
	  printf("node %d idof %d %d %d act %d\n",node,idof1,idof2,idof3,islavact[islavnodeentry-1]);	  	   
	}
	/** handle SPCs and MPCs -> modify tangentials and normal **/
	   
	if(nslavspc[2*(islavnodeentry-1)+1]-nslavspc[2*(islavnodeentry-1)]>0 ||nslavmpc[2*(islavnodeentry-1)+1]-nslavmpc[2*(islavnodeentry-1)]>0){	    
	  t[0]=0;t[1]=0;t[2]=0; 	   
	}  
	for(jj=nslavspc[2*(islavnodeentry-1)];jj<nslavspc[2*(islavnodeentry-1)+1];jj++){	    
	  t[0]=0;t[1]=0;t[2]=0; 	    
	  ist=islavspc[2*jj];            
	  dir=ndirboun[ist-1];	    
	  node1=nodeboun[ist-1];	    
	  if(debug==1){	      
	    printf("jj %d ist %d dir %d node %d\n",jj,ist,dir,node1);	    
	  }
          t[dir-1]=1.0;	    
	  dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];	    
	  n[0]=n[0]-dist*t[0];	    
	  n[1]=n[1]-dist*t[1];	    
	  n[2]=n[2]-dist*t[2];	    
	  nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	    
	  n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;	   
	}  	  
	for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){	    
	  t[0]=0;t[1]=0;t[2]=0;             
	  ist=islavmpc[2*jj];           
	  dirdep=nodempc[3*(ist-1)+1];            
	  coefdep=coefmpc[ist-1];            
	  index=nodempc[3*(ist-1)+2];	    
	  node1=nodempc[3*(ist-1)]; 	   
	  t[dirdep-1]=coefdep;    	    
	  if(debug==1){printf("jj %d ist %d dir %d node %d coef %e\n",jj,ist,dirdep,node1,coefdep);}	    
          while(index!=0){               
	    dirind=nodempc[3*(index-1)+1];	       
	    node1=nodempc[3*(index-1)]; 	       
	    if(node1==node) t[dirind-1]=coefmpc[index-1];               
               index=nodempc[3*(index-1)+2];
	       if(debug==1){printf("index %d \n",index);}
	    }
            nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
	    t[0]=t[0]/nt1;t[1]=t[1]/nt1;t[2]=t[2]/nt1;
	    dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];
	    n[0]=n[0]-dist*t[0];
	    n[1]=n[1]-dist*t[1];
	    n[2]=n[2]-dist*t[2];
	    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	    n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;
	   
	}
	/** get second tangential and modify it due to SPC's/MPC's**/	   
	if(nslavspc[2*(islavnodeentry-1)+1]-nslavspc[2*(islavnodeentry-1)]>0 ||nslavmpc[2*(islavnodeentry-1)+1]-nslavmpc[2*(islavnodeentry-1)]>0){	       	    
	  /** caution: t[3-5]= hat(t2)=n x hat(t1) **/	     
	  t[3]=n[1]*t[2]-n[2]*t[1];            
	  t[4]=n[2]*t[0]-n[0]*t[2];            
	  t[5]=n[0]*t[1]-n[1]*t[0];	   
	  dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];	  		    
	  for(kk=0;kk<6;kk++){that[kk]=t[kk];}	    
	    for(jj=nslavmpc[2*(islavnodeentry-1)];jj<nslavmpc[2*(islavnodeentry-1)+1];jj++){             
	      ist=islavmpc[2*jj];             
	      dirdep=nodempc[3*(ist-1)+1];             
	      coefdep=coefmpc[ist-1];             
	      index=nodempc[3*(ist-1)+2];	     
	      node1=nodempc[3*(ist-1)];                 
	      while(index!=0){              
		dirind=nodempc[3*(index-1)+1];	       
		node1=nodempc[3*(index-1)];	       
		if(node==node1){	        
		  t[dirind-1+3]=t[dirind-1+3]-coefmpc[index-1]*t[dirdep-1+3]/coefdep;		
		  n[dirind-1]=n[dirind-1]-coefmpc[index-1]*n[dirdep-1]/coefdep;			       
		}
                index=nodempc[3*(index-1)+2];	    	      
	      }
	      t[dirdep-1+3]=0.0;	    
	      n[dirdep-1]=0.0;	    	   
	    }	     
	    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	
	    nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
	    nt2=sqrt(t[3]*t[3]+t[4]*t[4]+t[5]*t[5]);
	    if(debug==1){	     
	      printf("nmod= %e %e %e \n",n[0],n[1],n[2]);	    
	      printf("t1= %e %e %e \n",t[0],t[1],t[2]);	    
	      printf("t2= %e %e %e \n",t[3],t[4],t[5]);	    
	      printf("nhat= %e %e %e \n",n2[0],n2[1],n2[2]);	    
	      printf("that1= %e %e %e \n",that[0],that[1],that[2]);	    
	      printf("that2= %e %e %e \n",that[3],that[4],that[5]);	    	    
	    }
	    /** caution: t[3-5]= tilde(t2), hat(t2) with modifications due to SPC/MPCs
	        *        that[0-5] =hat(t1) ,hat(t2)
	        *	n[0-2]= n*hat(D)/D in case of directinal blocking
		*       n[0-2]= n otherwise 
		*       n2[0-2]= original n
		**/	   
	}		    
	/** calculate needed fields for coulombfriction **/ 
	up_n=   u_old[(islavnodeentry-1)*3]*n2[0]+u_old[(islavnodeentry-1)*3+1]*n2[1]+u_old[(islavnodeentry-1)*3+2]*n2[2];
        up_t[0]=u_old[(islavnodeentry-1)*3]*that[0]+u_old[(islavnodeentry-1)*3+1]*that[1]+u_old[(islavnodeentry-1)*3+2]*that[2];
	up_t[1]=u_old[(islavnodeentry-1)*3]*that[3]+u_old[(islavnodeentry-1)*3+1]*that[4]+u_old[(islavnodeentry-1)*3+2]*that[5];
	lambda_n=n2[0]*cstress[(islavnodeentry-1)*3]+n2[1]*cstress[(islavnodeentry-1)*3+1]+n2[2]*cstress[(islavnodeentry-1)*3+2];
	lambda_t[0]=that[0]*cstress[(islavnodeentry-1)*3]+that[1]*cstress[(islavnodeentry-1)*3+1]+that[2]*cstress[(islavnodeentry-1)*3+2];
	lambda_t[1]=that[3]*cstress[(islavnodeentry-1)*3]+that[4]*cstress[(islavnodeentry-1)*3+1]+that[5]*cstress[(islavnodeentry-1)*3+2];
	nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
	bp=bp_old[islavnodeentry-1];
	/// perturebed lagrange
        FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&p0,&beta_e,tietol,elcon,&islavactdoftie[islavnodeentry-1],ncmat_,ntmat_));	
	derivmode=1;
	FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,&aninvloc,&p0,&beta_e,elcon,nelcon,&islavactdoftie[islavnodeentry-1],ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
        derivmode=0;
	FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&gnc,&aninvloc,&p0,&beta_e,elcon,nelcon,&islavactdoftie[islavnodeentry-1],ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
	if(islavact[islavnodeentry-1]==1){ ///stick	  
	  for(l=0;l<3;l++){
	    for(m=0;m<2;m++){
	      rstick[m*3+l]=mu*up_t[m]*n2[l];
	    }
	  }  
	  if(idof1>-1 && idof2>-1 && idof3>-1 &&jrow==1){
	    e1=bhat[k];
	    e2=bhat[k+1];
	    e3=bhat[k+2];
            /**right side if solving K du=f **/	     
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n[0]*e1+n[1]*e2+n[2]*e3);	     
	    b[k+1]=scal*(rstick[0]*e1+rstick[1]*e2+rstick[2]*e3);	     
	    b[k+2]=scal*(rstick[3]*e1+rstick[4]*e2+rstick[5]*e3);         
	    if(debug==1){	    
	      printf("\t bp=%e bp_alt=%e \n",mu*(lambda_n),bp);	    
	      printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);	   
	      printf("\t uold(%d)= %e %e %e u_t=%e %e max=%d \n",k,u_old[(islavnodeentry-1)*3],u_old[(islavnodeentry-1)*3+1],u_old[(islavnodeentry-1)*3+2],up_t[0],up_t[1],neq[1]);	    
	      printf("\t rstick1: %e %e %e\n",rstick[0],rstick[1],rstick[2]);	    
	      printf("\t rstick2: %e %e %e\n",rstick[3],rstick[4],rstick[5]);	  
	      printf("k=%d activ=%d\n",k,islavact[islavnodeentry-1]);	  
	      printf("\t b_old= %e %e %e \n",e1,e2,e3);	  
	      printf("\t b_new= %e %e %e \n",b[k],b[k+1],b[k+2]);	  
	    }
	    k=k+2;	  	  
	  }else if(idof2>-1 && idof3>-1){ ///2D auf 3D	    	    
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=rstick[4];
	    t2=rstick[5];
            b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[1]*e1+n2[2]*e2);	    
	    b[k+1]=scal*(t1*e1+t2*e2);
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }
	    k=k+1;  	  
	  }else if(idof1>-1 && idof3>-1){ ///2D auf 3D	    
	    b[k]=gap[islavnodeentry-1];
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=rstick[3];
	    t2=rstick[5];
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[0]*e1+n2[2]*e2);
	    b[k+1]=scal*(t1*e1+t2*e2);
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }	  
	    k=k+1; 	  
	  }else if(idof1>-1 && idof2>-1){ ///2D auf 3D	    
	    b[k]=gap[islavnodeentry-1];
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=rstick[3];
	    t2=rstick[4];
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[0]*e1+n2[1]*e2);
	    b[k+1]=scal*(t1*e1+t2*e2);
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }	 
	    k=k+1;  					  
	  }else{	    
	    e1=adc[k];
	    e2=bhat[k];
	    if(idof1>-1){
	      n11=n2[0];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[0]+scal*(dgnc)*(n11*e1);
	    }else if(idof2>-1) {
	      n11=n2[1];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[1]+scal*(dgnc)*(n11*e1);	    
	    }else{
	      n11=n2[2];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[2]+scal*(dgnc)*(n11*e1);	      
	    }	    
	  }  		
	}else if(islavact[islavnodeentry-1]==2){ ///slip	  
	  if(mu>1.E-10){ /// contact tie with friction	    
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
	      fp[1]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	      fp[2]= (1.0/(max(bp,nlambda_t)*nlmhat)*2)*(lambda_t[1]*lmhat[0]+lambda_t[0]*lmhat[1]);
	      fp[3]= (1.0/(max(bp,nlambda_t)*nlmhat))*lambda_t[1]*lmhat[1];
	    }else if(imodification==3){
	      fp[0]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[0];
	      fp[1]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[1]*lambda_t[0];
	      fp[2]= (1.0/(max(bp,nlambda_t)*max(bp,nlambda_t)))*lambda_t[0]*lambda_t[1];
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
	       tslip[m*3+l]=that[m*3+l]-mu*vp[m]*n2[l];
	      }
	    } 
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
	  }else{///contact tie without friction	      
	    for(l=0;l<3;l++){	    
	      for(m=0;m<2;m++){	      
		tslip[m*3+l]=t[m*3+l];	       
	      }	      	      	      
	    }
//	    scal=1.0;
	    ltu[0]=0.0;
	    ltu[1]=0.0;
	    rphat[0]=0.0;
	    rphat[1]=0.0;
	 }  
         if(debug==1){
	    printf("k=%d activ=%d\n",k,islavact[islavnodeentry-1]);
//	    printf("\t n= %e %e %e \n",n2[0],n2[1],n2[2]);
//	    printf("\t t1= %e %e %e \n",that[0],that[1],that[2]);
//	    printf("\t t2= %e %e %e \n",that[3],that[4],that[5]);
	    printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],lambda_t[1],nlambda_t,bp);
	    printf("\t lmhat= %e %e  nlmhat=%e \n",lmhat[0],lmhat[1],nlmhat);	    
	    printf("\t uold(%d)= %e %e %e u_t= %e %e %e \n",k,u_old[(islavnodeentry-1)*3],u_old[(islavnodeentry-1)*3+1],u_old[(islavnodeentry-1)*3+2],up_n,up_t[0],up_t[1]);
	    if(mu>1.E-10){	    
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
	    printf("\t ltu= %e %e \n",ltu[0],ltu[1]);
	    }else{
	      printf("\t no friction\n");
	    printf("\t tslip1: %e %e %e\n",tslip[0],tslip[1],tslip[2]);
	    printf("\t tslip2: %e %e %e\n",tslip[3],tslip[4],tslip[5]);	
	    printf("\t ltu= %e %e gap* %e  gap %e\n",ltu[0],ltu[1],gap[islavnodeentry-1]-up_n,gap[islavnodeentry-1]);	    
	    } 	
	 }
	 if(idof1>-1 && idof2>-1 && idof3>-1 &&jrow==1){
	   e1=bhat[k];
	   e2=bhat[k+1];
	   e3=bhat[k+2];
	  /* right side if solving K du=f*/ 
	   b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n[0]*e1+n[1]*e2+n[2]*e3);
	   b[k+1]=scal*(tslip[0]*e1+tslip[1]*e2+tslip[2]*e3)+rphat[0]-ltu[0];
	   b[k+2]=scal*(tslip[3]*e1+tslip[4]*e2+tslip[5]*e3)+rphat[1]-ltu[1];
	   if(debug==1){
	   printf("\t b_old= %e %e %e \n",e1,e2,e3);
	   printf("\t b_new= %e %e %e \n",b[k],b[k+1],b[k+2]);
	   }
	   k=k+2; 
	 }else if(idof2>-1 && idof3>-1){ ///2D auf 3D
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=tslip[4];
	    t2=tslip[5];
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[1]*e1+n2[2]*e2);
	    b[k+1]=scal*(t1*e1+t2*e2)-ltu[1];
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }
	    k=k+1;   
	 }else if(idof1>-1 && idof3>-1){ ///2D auf 3D
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=tslip[3];
	    t2=tslip[5];
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[0]*e1+n2[2]*e2);
	    b[k+1]=scal*(t1*e1+t2*e2)-ltu[1];
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }	  
	    k=k+1;  
	 }else if(idof1>-1 && idof2>-1){ ///2D auf 3D
	    e1=bhat[k];
	    e2=bhat[k+1];
	    t1=tslip[3];
	    t2=tslip[4];
	    b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+(dgnc)*scal*(n2[0]*e1+n2[1]*e2);
	    b[k+1]=scal*(t1*e1+t2*e2)-ltu[1];
	    if(debug==1){
	      printf("\t node %d bold %e %e b %e %e \n",node,bhat[k],bhat[k+1], b[k],b[k+1]);
	    }	 
	    k=k+1;  				 
	 }else{ ///1D auf 3D 
	    b[k]=gap[islavnodeentry-1];
	    e1=adc[k];
	    e2=bhat[k];
	    if(idof1>-1){
	      n11=n2[0];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[0]+scal*(dgnc)*(n11*e1);
	    }else if(idof2>-1) {
	      n11=n2[1];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[1]+scal*(dgnc)*(n11*e1);	    
	    }else{
	      n11=n2[2];
	      b[k]=gap[islavnodeentry-1]+gnc-dgnc*lambda_n+scal*(dgnc)*(n11*e2);
	      ad[k]=n[2]+scal*(dgnc)*(n11*e1);	      
	    }	    
	 }
	}

      }
    }
    free(lmhat);
    free(lambda_t);
    free(fp);
    free(mp);
    free(mp2);
    free(mp3);    
    free(lp);
    free(ltslip);
    free(rp);
    free(rphat);
    free(hp);
    free(n);
    free(n2);    
    free(t);
    free(up_t);
    free(vp);
    free(rstick);
    free(tslip);
    free(ltu);
    free(u_old);
//    FORTRAN(stop,());
}
