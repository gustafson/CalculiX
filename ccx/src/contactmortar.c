/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

void contactmortar(int *ncont, int *ntie, char *tieset, int *nset, char *set,
        int *istartset, int *iendset, int *ialset, int *itietri,
        char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
        double *cg, double *straight, double *co,
        double *vold, int *ielmat, double *cs, double *elcon,
        int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
        int *ifcont1, int *ifcont2, int *ne0, double *vini,
        int *nmethod,int *neq, int *nzs, int *nactdof, int *itiefac,
        int *islavsurf, int *islavnode, int *imastnode,
        int *nslavnode, int *nmastnode, int *ncone, double *ad,
        double **aup, double *b, int **irowp, int *icol, int *jq, int *imastop,
        int *iponoels, int *inoels, int *nzsc, double **aucp,
        double *adc, int **irowcp, int *jqc, int *islavact,
        double *gap, double *bdd, double **auqdtp, int **irowqdtp,
        int *jqqdt, int *nzsqdt, int *nzlc,double *slavnor,double *bhat,
        int *icolc){
    
    int i,j,k,l,m,numb,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,ifree=0,kflag, 
        nzsbd,*mast1=NULL,*ipointer=NULL,*irowbd=NULL,*jqbd=NULL, 
        *irowc=NULL,*imastsurf=NULL,jrow,jcol,islavnodeentry,
        *islavactdof=NULL,*irow=NULL,*irowqdt=NULL,* jqctemp=NULL, 
        * irowctemp=NULL,jslavnodeentry;

    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*aubd=NULL, 
      *auc=NULL, *pmastsurf=NULL,*auqdt=NULL,
      *slavtan=NULL,t1,t2,t3,e1,e2,e3,*au=NULL,* auctemp=NULL;
	  

    irow = *irowp; au=*aup; auc=*aucp; irowc=*irowcp; auqdt=*auqdtp;
    irowqdt=*irowqdtp;
	
    FORTRAN(updatecont,(koncont,ncont,co,vold,
            cg,straight));

    /* determining the size of the auxiliary fields */

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
    nx=NNEW(int,ntrimax);
    ny=NNEW(int,ntrimax);
    nz=NNEW(int,ntrimax);

    imastsurf=NNEW(int,9**ncone);
    pmastsurf=NNEW(double,18**ncone);
    slavtan=NNEW(double,6*nslavnode[*ntie]);
    
    FORTRAN(gencontrel,(tieset,ntie,itietri,ipkon,kon,
        lakon,set,cg,straight,&ifree,
        koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,cs,
        elcon,istep,iinc,iit,ncmat_,ntmat_,
        vini,nmethod,islavsurf,imastsurf,pmastsurf,itiefac,
        islavnode,nslavnode,slavnor,slavtan,imastop,gap,
        islavact));
	
    free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);
    free(ny);free(nz);

    RENEW(imastsurf,int,ifree);
    RENEW(pmastsurf,double,2*ifree);
	
    /* coupling the active slave degrees of freedom with the corresponding slave
       node */
    
    islavactdof=NNEW(int,neq[1]);
    
    FORTRAN(genislavactdof,(ntie,neq,nactdof,nslavnode,islavact,islavactdof,islavnode));
    
    /* calculating the coupling matrices bd (described by aubd, irowbd, jqbd and nzsbd; 
       nonsymmetric, zero diagonal) and dd (in bdd; diagonal matrix */

    nzsbd = nzs[1];

    aubd = NNEW(double, nzsbd);
    irowbd = NNEW(int, nzsbd);
    jqbd = NNEW(int, neq[1]+1);

    bdfill(&irowbd, jqbd, &aubd, bdd, &nzsbd, ntie,
        ipkon, kon, lakon, nslavnode, nmastnode, imastnode, islavnode, 
        islavsurf, imastsurf, pmastsurf, itiefac, neq, nactdof,co,vold,
        iponoels, inoels); 

    free(imastsurf);free(pmastsurf);
    
    /* modifying the stiffnes matrix with the coupling matrices; the
       modified (symmetric) matrix is described in asymmetric form by
       the fields auc, adc, irowc, jqc and nzsc */ 

    *nzsc = nzs[1];
    RENEW(auc,double, *nzsc);
    RENEW(irowc,int, *nzsc);
    
    *nzsqdt=*nzsc;
    auqdt=NNEW(double,*nzsqdt);
    irowqdt=NNEW(int,*nzsqdt);

    multimortar(au, ad, irow, jq, nzs,
	   aubd, bdd, irowbd, jqbd, &nzsbd,
	   &auc, adc, &irowc, jqc, nzsc,
	   auqdt,irowqdt,jqqdt,nzsqdt,
	   neq,b,bhat);

    /* sorting the row numbers within each column in irowc */

    for (j=0;j<neq[1];j++){
      kflag = 2;
      if(jqc[j+1]-jqc[j]>0){
	numb=(jqc[j+1]-jqc[j]);
	FORTRAN(isortid,(&irowc[jqc[j]-1],&auc[jqc[j]-1],&numb,&kflag));
      }
    }
	
    /* copying auc,adc,irowc, jqc and bhat into
       au,ad,irow,jq and b */

    RENEW(au,double,*nzsc);
    RENEW(irow,int,*nzsc);
    for(i=0;i<neq[1];i++){
      ad[i]=adc[i];
      jq[i]=jqc[i];
      b[i]=bhat[i];
    }
    jq[neq[1]]=jqc[neq[1]];
    for(i=0;i<*nzsc;i++){
      au[i]=auc[i];
      irow[i]=irowc[i];
    }

	nzs[1]=jq[neq[1]]-1;
	
	for(i=0;i<neq[1];i++){
	  printf("islavactdof=%d %d\n",i+1,islavactdof[i]);
	}


    /* changing au due to N and T (normal and tangential
       direction at the slave surface */
    
    for(j=0;j<neq[1];j++){
      for(i=jq[j]-1;i<jq[j+1]-1;i++){
	k=irow[i]-1;

	/*  k is the row number, j is the column number */

	if(islavactdof[k]>0){
	  islavnodeentry = floor(islavactdof[k]/10.);
	  jrow= islavactdof[k]-10*islavnodeentry;

	  /* islavnodeentry is the slave node entry in field
             islavnode corresponding to row k,
             jslavnodeentry is the slave node entry in field
             islavnode corresponding to column j */

	  if (jrow==1){
	    e1=auc[i];
	    jslavnodeentry=floor(islavactdof[j]/10.);
	    if (islavnodeentry!=jslavnodeentry){
	      au[i]=0;
	    }else{
	      jcol=islavactdof[j]-10*jslavnodeentry;
	      	      au[i]=slavnor[3*(islavnodeentry-1)+jcol-1];
	      //au[i]=slavnor[3*(islavnodeentry-1)+jcol-1]*bdd[k]; //normer avec d le N..
	    }
	    printf("au,%d,%d,%d,%f\n",k+1,j+1,jrow,au[i]);
	  }
	  else if (jrow==2){
	    t1=slavtan[6*(islavnodeentry-1)];
	    t2=slavtan[6*(islavnodeentry-1)+1];
	    t3=slavtan[6*(islavnodeentry-1)+2];
	    e2=auc[i];
	    if((k+1)==j){e3=adc[j];}else{e3=auc[i+1];}
	    au[i]=t1*e1+t2*e2+t3*e3;
	    printf("au,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,j+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	  }
	  else{
	    t1=slavtan[6*(islavnodeentry-1)+3];
	    t2=slavtan[6*(islavnodeentry-1)+4];
	    t3=slavtan[6*(islavnodeentry-1)+5];
	    //	    e3=au[i];
	    au[i]=t1*e1+t2*e2+t3*e3;
	    printf("au,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,j+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	  }
	}

        /* diagonal terms */

	if(k==(j-1)){
	  if(islavactdof[j]>0){
	    jslavnodeentry=floor(islavactdof[j]/10.);
	    jrow=islavactdof[j]-10*jslavnodeentry;
	    if(jrow==1){
	      e1=adc[j];
	      	      ad[j]=slavnor[3*(jslavnodeentry-1)+jrow-1];
								    //ad[j]=slavnor[3*(jslavnodeentry-1)+jrow-1]*bdd[j];
	    printf("ad,%d,%d,%f\n",j+1,jrow,ad[j]);
	    }else if(jrow==2){
	      t1=slavtan[6*(jslavnodeentry-1)];
	      t2=slavtan[6*(jslavnodeentry-1)+1];
	      t3=slavtan[6*(jslavnodeentry-1)+2];
	      e2=adc[j];
	      e3=auc[i+1];
	      ad[j]=t1*e1+t2*e2+t3*e3;
	    printf("ad,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",j+1,jrow,t1,t2,t3,e1,e2,e3,ad[j]);
	    }else{
	      t1=slavtan[6*(jslavnodeentry-1)+3];
	      t2=slavtan[6*(jslavnodeentry-1)+4];
	      t3=slavtan[6*(jslavnodeentry-1)+5];
	      //	      e3=ad[j];
	      ad[j]=t1*e1+t2*e2+t3*e3;
	    printf("ad,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",j+1,jrow,t1,t2,t3,e1,e2,e3,ad[j]);
	    }
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
	  if (jrow==1){
	    e1=bhat[k];
	    b[k]=gap[islavnodeentry-1];
	    printf("b,%d,%d,%f\n",k+1,jrow,ad[j]);
	  }
	  else if (jrow==2){
	    t1=slavtan[6*(islavnodeentry-1)];
	    t2=slavtan[6*(islavnodeentry-1)+1];
	    t3=slavtan[6*(islavnodeentry-1)+2];
	    e2=bhat[k];
	    e3=bhat[k+1];
	    b[k]=t1*e1+t2*e2+t3*e3;
	    printf("b,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	  }
	  else{
	    t1=slavtan[6*(islavnodeentry-1)+3];
	    t2=slavtan[6*(islavnodeentry-1)+4];
	    t3=slavtan[6*(islavnodeentry-1)+5];
	    //	    e3=b[k];
	    b[k]=t1*e1+t2*e2+t3*e3;
	    printf("b,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	  }
	}
      }

  int number=10;

  FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));

   number=7;

  FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));


  // number=8;

//  FORTRAN(writematrix,(au,bhat,irow,jq,&neq[1],&number));


 // number=9;

 // FORTRAN(writematrix,(au,b,irow,jq,&neq[1],&number));

    
    free(islavactdof);
	
    free(slavtan);
    free(aubd);free(jqbd);free(irowbd);
    
    /* So far every nonzero in Auc was stored; however,
       Auc is symmetric. To reduce the computational effort
       in subroutine contactstress Auc is now stored in
       a symmetrical form, i.e. only half of the matrix
       is stored */
    
    auctemp=NNEW(double, *nzsc);
    jqctemp=NNEW(int,neq[1]+1);
    irowctemp=NNEW(int,*nzsc);
    
    k=0;
    jqctemp[0]=1;
    for (i=0;i<neq[1];i++){ 
      
      for (j=jqc[i]-1;j<jqc[i+1]-1;j++){ 
	if(irowc[j]>i+1){ 
	  auctemp[k]=auc[j]; 
	  irowctemp[k++]=irowc[j];	
	  
	}			
      }
      jqctemp[i+1]=k+1; 
      
    }
    jqctemp[neq[1]]=k+1;
    
    *nzsc=k;
    for(i=0;i<*nzsc;i++){
      auc[i]=auctemp[i];
      irowc[i]=irowctemp[i];
    }
    
    for(i=0;i<neq[1]+1;i++){
      jqc[i]=jqctemp[i];
    }
    
    
    free(auctemp);
    free(irowctemp);
    
    RENEW(auc,double,k);
    RENEW(irowc,int,k);
    
    /* calculating icol and icolc (needed for SPOOLES) */
    
    for(i=0; i<neq[1]; i++){
      icol[i] = jq[i+1]-jq[i];
    }
    
    for(i=0; i<neq[1]; i++){
      icolc[i] = jqc[i+1]-jqc[i];
    }	
    
    free(jqctemp);

    /* nzlc is the number of the rightmost column with 
       nonzero off-diagonal terms */

    //  number=10;

//  FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));

    *nzlc=0;
    for(i=neq[1]-1;i>-1;i--){
      if(icolc[i]>0){
	*nzlc=i+1;
	break;
      }
    }
    
    *irowp = irow; *aup=au; *aucp=auc; *irowcp=irowc; *auqdtp=auqdt;
    *irowqdtp=irowqdt;

    return;
}
