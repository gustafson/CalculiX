/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h> 
#include "CalculiX.h"

void multi_rect(double *au_1,int * irow_1,int * jq_1,int n_1, int m_1,
	       double *au_2,int * irow_2,int * jq_2,int n_2, int m_2,
               double **au_rp,int **irow_rp,int * jq_r,int *nzs){ 

       /*Result fields*/
       	int *irow=NULL,ifree=1,numb,icol,i,j,k,l,m,carre=0,kflag=2,istart,icounter;
        int flag=0;
	double *au=NULL,value;
        clock_t debut;
    	clock_t fin;
        /*Perform a_1T*a_2*/

	debut=clock();
 	if (n_1!=n_2) {
		printf("Error in mutli_rec : Matrix sizes are not compatible\n");
		return;
	} 

//        nzs=n_1*m_2;
        irow=*irow_rp;
        au=*au_rp; 

        if (n_1==m_2) carre=1;
	
        jq_r[0]=1;
        for(j=0;j<m_2;j++){
                m=j+1;
		for (i=0;i<m_1;i++){
                l=i+1;
		flag=0;
                value=0.0;
		multi_scal(au_1,irow_1,jq_1,au_2,irow_2,jq_2,i,j,&value,&flag);
		if (flag!=0) insertas_ws(&irow,&l,&m,&ifree,nzs,&value,&au);
		}
           jq_r[m]=ifree;
	}

	/* Sort the column and compute jq*/ 
    
	*nzs=ifree;
	RENEW(au,double,*nzs);
	RENEW(irow,int,*nzs);
        *irow_rp=irow;*au_rp=au;   
	   fin= clock();
	printf("multi_rect : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
 return;
}
