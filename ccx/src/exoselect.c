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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void exoselect(double *field1,double *field2,ITG *iset,ITG *nkcoords,ITG *inum,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ngraph,ITG *ncomp,
	       ITG *ifield,ITG *icomp,ITG *nfield,ITG *iselect,ITG exoid,
	       ITG time_step, int countvar, ITG nout, ITG *node_map_inv){

  /* storing scalars, components of vectors and tensors without additional
     transformations */

  /* number of components in field1: nfield[0]
     number of components in field2: nfield[1]

     number of entities to store: ncomp
     for each entity i, 0<=i<ncomp:
     - ifield[i]: 1=field1,2=field2
     - icomp[i]: component: 0...,(nfield[0]-1 or nfield[1]-1) */

  ITG i,j,k,l,m,n,nksegment;

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
  // ITG num_nod_vars = *ncomp;
  float *nodal_var_vals;
  nodal_var_vals = (float *) calloc (nout, sizeof(float));


  for(j=0;j<*ncomp;j++){
    if(*iset==0){
      for(i=0;i<*nkcoords;i++){
	/* check whether output is requested for solid nodes or
	   network nodes */

	if(*iselect==1){
	  if(inum[i]<=0) continue;
	}else if(*iselect==-1){
	  if(inum[i]>=0) continue;
	}else{
	  if(inum[i]==0) continue;
	}

	/* storing the entities */
	if(ifield[j]==1){
	  nodal_var_vals[node_map_inv[i]-1]=field1[i*nfield[0]+icomp[j]];
	}else{
	  nodal_var_vals[node_map_inv[i]-1]=field2[i*nfield[1]+icomp[j]];
	}
      }
    }else{
      nksegment=(*nkcoords)/(*ngraph);
      for(k=istartset[*iset-1]-1;k<iendset[*iset-1];k++){
	if(ialset[k]>0){
	
	  for(l=0;l<*ngraph;l++){
	    i=ialset[k]+l*nksegment-1;
	
	    /* check whether output is requested for solid nodes or
	       network nodes */
	
	    if(*iselect==1){
	      if(inum[i]<=0) continue;
	    }else if(*iselect==-1){
	      if(inum[i]>=0) continue;
	    }else{
	      if(inum[i]==0) continue;
	    }
	
	    /* storing the entities */
	    if(ifield[j]==1){
	      nodal_var_vals[node_map_inv[i]-1]=field1[i*nfield[0]+icomp[j]];
	    }else{
	      nodal_var_vals[node_map_inv[i]-1]=field2[i*nfield[1]+icomp[j]];
	    }
	  }
	}else{
	  l=ialset[k-2];
	  do{
	    l-=ialset[k];
	    if(l>=ialset[k-1]) break;
	    for(m=0;m<*ngraph;m++){
	      i=l+m*nksegment-1;
	
	      /* check whether output is requested for solid nodes or
		 network nodes */
	
	      if(*iselect==1){
		if(inum[i]<=0) continue;
	      }else if(*iselect==-1){
		if(inum[i]>=0) continue;
	      }else{
		if(inum[i]==0) continue;
	      }
	
	      /* storing the entities */
	      if(ifield[j]==1){
		nodal_var_vals[node_map_inv[i]-1]=field1[i*nfield[0]+icomp[j]];
	      }else{
		nodal_var_vals[node_map_inv[i]-1]=field2[i*nfield[1]+icomp[j]];
	      }
	    }
	  }while(1);
	}
      }
    }

    // Note: the results for exo tensor format are strangely ordered,
    if (*ncomp==12){
      k=j;
    }else{k=j;}

    int errr = ex_put_var (exoid, time_step, EX_NODAL, k+1+countvar, 1, nout, nodal_var_vals);
    if (errr) printf ("ERROR exoselect data for dim %i record %i.\n", j, k+1+countvar);
  }

  free(nodal_var_vals);
  return;

}

#endif
