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
#include <string.h> 
#include "CalculiX.h"

void multi_rectv(double *au_1,int * irow_1,int * jq_1,int n_1, int m_1,
		 double * b, double ** v_rp){ 

       /*Result fields*/
       	int irow,i,j;
	double *v_r=NULL,value;
        
        
        v_r=NNEW(double,m_1);   

	for(j=0;j<m_1;j++){
		for (i=jq_1[j]-1;i<jq_1[j+1]-1;i++){
//			v_r[j]-=au_1[i]*b[j];
			v_r[j]+=au_1[i]*b[irow_1[i]-1];
		}
	}

        *v_rp=v_r;

 return;
}
