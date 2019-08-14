/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */
/*     This subroutine                                                   */
/*              Copyright (C) 2013-2018 Peter A. Gustafson               */
/*                                                                       */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);             */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef EXODUSII
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include "exodusII.h"
#include "exo.h"

void exovector(double *v,ITG *iset,ITG *ntrans,char * filabl,ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3, int exoid, ITG time_step,
	       int countvar, ITG nout, ITG *node_map_inv){

  ITG nksegment;
  ITG i,j,k,l,m,n,ii,jj,kk;

  double a[9];

  /* When initializing parameter values:
     "g" (or "G")
     For global variables.
     "n" (or "N")
     For nodal variables.
     "e" (or "E")
     For element variables.
     "m" (or "M")
     For nodeset variables.
     "s" (or "S")
     For sideset variables.
  */

  int errr;

  int num_nod_vars=3;

  float *nodal_var_vals;
  nodal_var_vals = (float *) calloc (nout, sizeof(float));

  for (j=1; j<=num_nod_vars; j++){ // For each direction
    if(*iset==0){
      if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)){
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0) continue;
	  nodal_var_vals[node_map_inv[i]-1]=v[(mi[1]+1)*i+j];
	}
      }else{
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0) continue;
	  if(inotr[2*i]==0){
	    nodal_var_vals[node_map_inv[i]-1]=v[(mi[1]+1)*i+j];
	  }else{
	    ii=(mi[1]+1)*i+1;
	    jj=(mi[1]+1)*i+2;
	    kk=(mi[1]+1)*i+3;
	    FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	    nodal_var_vals[node_map_inv[i]-1]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
	  }
	}
      }
    }else{
      nksegment=(*nkcoords)/(*ngraph);
      for(k=istartset[*iset-1]-1;k<iendset[*iset-1];k++){
	if(ialset[k]>0){
	  for(l=0;l<*ngraph;l++){
	    i=ialset[k]+l*nksegment-1;
	    if(inum[i]<=0) continue;
	    if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)||(inotr[2*i]==0)){
	      nodal_var_vals[node_map_inv[i]-1]=v[(mi[1]+1)*i+j];
	    }else{
	      FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	      ii=(mi[1]+1)*i+1;
	      jj=(mi[1]+1)*i+2;
	      kk=(mi[1]+1)*i+3;
	      nodal_var_vals[node_map_inv[i]-1]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
	    }
	  }
	}else{
	  l=ialset[k-2];
	  do{
	    l-=ialset[k];
	    if(l>=ialset[k-1]) break;
	    for(m=0;m<*ngraph;m++){
	      i=l+m*nksegment-1;
	      if(inum[i]<=0) continue;
	      if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)||(inotr[2*i]==0)){
		nodal_var_vals[node_map_inv[i]-1]=v[(mi[1]+1)*i+j];
	      }else{
		FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
		ii=(mi[1]+1)*i+1;
		jj=(mi[1]+1)*i+2;
		kk=(mi[1]+1)*i+3;
		nodal_var_vals[node_map_inv[i]-1]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
	      }
	    }
	  }while(1);
	}
      }
    }

    errr = ex_put_var (exoid, time_step, EX_NODAL, j+countvar, 1, nout, nodal_var_vals);
    if (errr) printf ("ERROR storing vector data into exo file.\n");
  }

  free(nodal_var_vals);
  return;

}
#endif
