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

/** changing au due to N and T (normal and tangential
   *    direction at the slave surface) 
   * 	changing b due to N and T (normal and tangential
   *	direction at the slave surface) 
 * @param [out] au
 * @param [out] b
*/
void trafoNTmortar(int *neq,int *nzs, int *islavactdof,int *nslavnode, int *nmastnode, int *ncone, 
        double *ad, double *au, double *b, int *irow, int *jq,
        int *nzsc, double *auc,
        double *adc, int *irowc, int *jqc,
        double *gap, double *bdd, double *auqdt, int *irowqdt,
        int *jqqdt, int *nzsqdt, int *nzlc,double *slavnor, double *slavtan,double *bhat,
	double *aubd, int *irowbd, int *jqbd){

    int i,j,k,numb,ntrimax,debug,
        nzsbd,l,jrow,jcol,islavnodeentry,jslavnodeentry;

    double t1,t2,t3,e1,e2,e3,help,s,dd1,dd2;
    
    debug=0;


/*  	for(j=nslavnode[0]; j<nslavnode[1]; j++){
           printf("trafoNT: slavnor(%d): %e, %e, %e \n", j,slavnor[3*(j-1)],slavnor[3*(j-1)+1],slavnor[3*(j-1)+2]);
        }
*/ 
     for(j=0;j<neq[1];j++){
//        printf("trafoNT: j: %d, jq[j]: %d, jq[j+1]:%d \n", j,jq[j],jq[j+1] );
        for(i=jq[j]-1;i<jq[j+1]-1;i++){
	  jslavnodeentry = floor(islavactdof[j]/10.);
	  jcol= islavactdof[j]-10*jslavnodeentry;
	  k=irow[i]-1;
          islavnodeentry = floor(islavactdof[k]/10.);
	  jrow= islavactdof[k]-10*islavnodeentry;
          s=1.0;
 //         if(k==5742 || k==5743 || k==5744){debug=1;}else{debug=0;}
//         if(k==13627 || k==13628 || k==13629){debug=1;}else{debug=0;}
          if(islavactdof[j]>0){// jcol activ, jrow activ
            if(islavactdof[k]>0){
            if(debug==1){
		printf("slavnor: %e %e %e \n",slavnor[3*(islavnodeentry-1)+1-1],slavnor[3*(islavnodeentry-1)+2-1],slavnor[3*(islavnodeentry-1)+3-1]);
		help=slavnor[3*(islavnodeentry-1)+1-1]*slavtan[6*(islavnodeentry-1)]+slavnor[3*(islavnodeentry-1)+2-1]*slavtan[6*(islavnodeentry-1)+1]+slavnor[3*(islavnodeentry-1)+3-1]*slavtan[6*(islavnodeentry-1)+2];
		dd1=slavnor[3*(islavnodeentry-1)]*slavtan[6*(islavnodeentry-1)]+slavnor[3*(islavnodeentry-1)+1]*slavtan[6*(islavnodeentry-1)+1]+slavnor[3*(islavnodeentry-1)+2]*slavtan[6*(islavnodeentry-1)+2];
		dd2=slavtan[6*(islavnodeentry-1)]*slavtan[6*(islavnodeentry-1)+3]+slavtan[6*(islavnodeentry-1)+1]*slavtan[6*(islavnodeentry-1)+4]+slavtan[6*(islavnodeentry-1)+2]*slavtan[6*(islavnodeentry-1)+5];
		printf("scalarproducts: %e %e %e \n",help,dd1,dd2);
         printf("trafoNT: j:%d, k:%d, slavnodeenty,j,i: %d, %d \n ",j,k, jslavnodeentry,islavnodeentry);}
               switch(jcol){
                 case 1:
                   if(k==j+1){//1. obere Nebendiagonale-> speichere gleich auch Werte für Diagonale  und 2. Nebendiagonale
                     	e1=adc[j];
                     	e2=auc[i];
                     	e3=auc[i+1];
                        help=slavnor[3*(islavnodeentry-1)+jcol-1];
                     	ad[j]=s*slavnor[3*(islavnodeentry-1)+jcol-1];
                     	t1=slavtan[6*(islavnodeentry-1)];
			t2=slavtan[6*(islavnodeentry-1)+1];
			t3=slavtan[6*(islavnodeentry-1)+2];
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
			au[i]=t1*e1+t2*e2+t3*e3;
                     	t1=slavtan[6*(islavnodeentry-1)+3];
			t2=slavtan[6*(islavnodeentry-1)+4];
			t3=slavtan[6*(islavnodeentry-1)+5];
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
			au[i+1]=t1*e1+t2*e2+t3*e3; //Achtung: i wird hochgezählt! 
        if(debug==1){
	printf("trafoNT: K+d(1-3):%e,\t %e,\t %e \n ",adc[j],auc[i],auc[i+1]);
          printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d K(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, ad[j],au[i],au[i+1]);
	printf("\n");}
			i=i+1; 

                   }else{
			if(jrow==1){
				if(islavnodeentry!=jslavnodeentry){
				au[i]=0;
				}else{
				au[i]=s*slavnor[3*(islavnodeentry-1)+jcol-1];
				}
                     		e1=auc[i];
                     		e2=auc[i+1];
                     		e3=auc[i+2];
                     		t1=slavtan[6*(islavnodeentry-1)];
				t2=slavtan[6*(islavnodeentry-1)+1];
				t3=slavtan[6*(islavnodeentry-1)+2];
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
				au[i+1]=t1*e1+t2*e2+t3*e3;
                     		t1=slavtan[6*(islavnodeentry-1)+3];
				t2=slavtan[6*(islavnodeentry-1)+4];
				t3=slavtan[6*(islavnodeentry-1)+5];
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
				au[i+2]=t1*e1+t2*e2+t3*e3;
        if(debug==1){
	printf("trafoNT: K(1-3):%e,\t %e,\t %e \n ",auc[i],auc[i+1],auc[i+2]);
          printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d K(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],au[i+1],au[i+2]);
	printf("\n");}
				i=i+2;
                        }else{printf("trafoNT: Idofs ordered the wrong way...stop");
				FORTRAN(stop,());}
                   }
		   break;
  		case 2:
		 if(k==j-1){//1. untere Nebendiagonale-> speichere gleich auch Werte für Diagonale  und 1. oND
	      		e1=auc[i];
	      		e2=adc[j];
	      		e3=auc[i+1];
                        help=slavnor[3*(islavnodeentry-1)+jcol-1];
	      		au[i]=s*slavnor[3*(islavnodeentry-1)+jcol-1];
	      		t1=slavtan[6*(islavnodeentry-1)];
	      		t2=slavtan[6*(islavnodeentry-1)+1];
	      		t3=slavtan[6*(islavnodeentry-1)+2]; 
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
	      		ad[j]=t1*e1+t2*e2+t3*e3; 
	      		t1=slavtan[6*(islavnodeentry-1)+3];
	      		t2=slavtan[6*(islavnodeentry-1)+4];
	      		t3=slavtan[6*(islavnodeentry-1)+5];
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
	      		au[i+1]=t1*e1+t2*e2+t3*e3;
        if(debug==1){
	printf("trafoNT: K+d(1-3):%e,\t %e,\t %e \n ",auc[i],adc[j],auc[i+2]);
         printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d K(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],ad[j],au[i+1]);
	printf("\n");}
			i=i+1;//Achtung: i wird hochgezählt! 
                  }else{
			if(jrow==1){//->normal
				if(islavnodeentry!=jslavnodeentry){
				au[i]=0;
				}else{
				au[i]=s*slavnor[3*(islavnodeentry-1)+jcol-1];
				}
                     		e1=auc[i];
                     		e2=auc[i+1];
                     		e3=auc[i+2];
                     		t1=slavtan[6*(islavnodeentry-1)];
				t2=slavtan[6*(islavnodeentry-1)+1];
				t3=slavtan[6*(islavnodeentry-1)+2];
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
				au[i+1]=t1*e1+t2*e2+t3*e3;
                     		t1=slavtan[6*(islavnodeentry-1)+3];
				t2=slavtan[6*(islavnodeentry-1)+4];
				t3=slavtan[6*(islavnodeentry-1)+5];
				au[i+2]=t1*e1+t2*e2+t3*e3;
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){
	printf("trafoNT: K(1-3):%e,\t %e,\t %e \n ",auc[i],auc[i+1],auc[i+2]);
          printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d K(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],au[i+1],au[i+2]);
	printf("\n");}
				i=i+2;
                        }else {printf("trafoNT: Idofs ordered the wrong way...stop");
				FORTRAN(stop,());}
			
		  }
                  break;
		case 3:
		 if(k==j-2){//1. untere Nebendiagonale-> speichere gleich auch Werte für 1.uND und Diagonale
	      		e1=auc[i];
	      		e2=auc[i+1];
	      		e3=adc[j];
			help=slavnor[3*(islavnodeentry-1)+jcol-1];
	      		au[i]=s*slavnor[3*(jslavnodeentry-1)+jcol-1];
	      		t1=slavtan[6*(islavnodeentry-1)];
	      		t2=slavtan[6*(islavnodeentry-1)+1];
	      		t3=slavtan[6*(islavnodeentry-1)+2]; 
	      		au[i+1]=t1*e1+t2*e2+t3*e3; 
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);}
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);} 
	      		t1=slavtan[6*(islavnodeentry-1)+3];
	      		t2=slavtan[6*(islavnodeentry-1)+4];
	      		t3=slavtan[6*(islavnodeentry-1)+5];
	      		ad[j]=t1*e1+t2*e2+t3*e3;
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){
	printf("trafoNT: K+d(1-3):%e,\t %e,\t %e \n ",auc[i],auc[i+1],adc[j]);
        printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d K(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],au[i+1],ad[j]);
	printf("\n");}
			i=i+1;//Achtung: i wird hochgezählt! 
                  }else{
			if(jrow==1){//->normal
				if(islavnodeentry!=jslavnodeentry){
				au[i]=0;
				}else{
				au[i]=s*slavnor[3*(islavnodeentry-1)+jcol-1];
				}
                     		e1=auc[i];
                     		e2=auc[i+1];
                     		e3=auc[i+2];
                     		t1=slavtan[6*(islavnodeentry-1)];
				t2=slavtan[6*(islavnodeentry-1)+1];
				t3=slavtan[6*(islavnodeentry-1)+2];
	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
				au[i+1]=t1*e1+t2*e2+t3*e3;
                     		t1=slavtan[6*(islavnodeentry-1)+3];
				t2=slavtan[6*(islavnodeentry-1)+4];
				t3=slavtan[6*(islavnodeentry-1)+5];
				au[i+2]=t1*e1+t2*e2+t3*e3;
	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
        if(debug==1){
	printf("trafoNT: K(1-3):%e,\t %e,\t %e \n ",auc[i],auc[i+1],auc[i+2]);
         printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d hK(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],au[i+1],au[i+2]);
	printf("\n");}
				i=i+2;
                        }else {printf("trafoNT: Idofs ordered the wrong way...stop");
				FORTRAN(stop,());}
		  }
                  break;
		default: printf("problem of active set");
		  break;
            }
//	   printf("\n");
           }
          }else{// column inaktiv, row actic
            if(islavactdof[k]>0){
//        printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d \n ",j,i,k,jcol,jrow);	    
			if(jrow==1){//->normal				
				au[i]=0;
                     		e1=auc[i];
                     		e2=auc[i+1];
                     		e3=auc[i+2];
                     		t1=slavtan[6*(islavnodeentry-1)];
				t2=slavtan[6*(islavnodeentry-1)+1];
				t3=slavtan[6*(islavnodeentry-1)+2];
//	if(debug==1){  printf("trafoNT: slavtan1(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
//        if(debug==1){               printf("a= %e b= %e c= %e \n",t1*e1,t2*e2,t3*e3);}
				au[i+1]=t1*e1+t2*e2+t3*e3;
                     		t1=slavtan[6*(islavnodeentry-1)+3];
				t2=slavtan[6*(islavnodeentry-1)+4];
				t3=slavtan[6*(islavnodeentry-1)+5];
				au[i+2]=t1*e1+t2*e2+t3*e3;
//	if(debug==1){  printf("trafoNT: slavtan2(1-3):%e,\t %e,\t %e \n ",t1,t2,t3);} 
//       if(debug==1){
//	printf("trafoNT: K(1-3):%e,\t %e,\t %e \n ",auc[i],auc[i+1],auc[i+2]);
//         printf("trafoNT: j:%d, i:%d, k:%d, jcol:%d, jrow:%d hK(1-3):%e,\t %e,\t %e \n ",j,i,k,jcol,jrow, au[i],au[i+1],au[i+2]);
//	printf("\n");}
				i=i+2;
                        }else {printf("trafoNT: Idofs ordered the wrong way...stop");
				FORTRAN(stop,());}               
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
	  b[k]=gap[islavnodeentry-1];
//	  	    printf("jrow=1 %d %e\n",islavnodeentry,b[k]);
	}
	else if (jrow==2){
	  t1=slavtan[6*(islavnodeentry-1)];
	  t2=slavtan[6*(islavnodeentry-1)+1];
	  t3=slavtan[6*(islavnodeentry-1)+2];
	  e1=bhat[k-1];
	  e2=bhat[k];
	  e3=bhat[k+1];
	  b[k]=t1*e1+t2*e2+t3*e3;
//	  	    printf("jrow=2 %d %e\n",islavnodeentry,b[k]);
	  //printf("b,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	}
	else{
	  t1=slavtan[6*(islavnodeentry-1)+3];
	  t2=slavtan[6*(islavnodeentry-1)+4];
	  t3=slavtan[6*(islavnodeentry-1)+5];
	  e1=bhat[k-2];
	  e2=bhat[k-1];
	  e3=bhat[k];
	  //	    e3=b[k];
	  b[k]=t1*e1+t2*e2+t3*e3;
//	  	    printf("jrow=3 %d %e\n",islavnodeentry,b[k]);
	  //printf("b,%d,%d,%f,%f,%f,%f,%f,%f,%f\n",k+1,jrow,t1,t2,t3,e1,e2,e3,au[i]);
	}
      }
    }
/*    for (i=0;i<neq[1];i++){
  	printf("TNT:bhat(%d)= %e \t b(%d)0 %e \n",i,bhat[i],i,b[i]);
   }*/
}