/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                     */
/*     This subroutine                                                   */
/*              Copyright (C) 2013-2021 Peter A. Gustafson               */
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
	       int countvar){

  ITG nksegment;
  ITG i,j,k,l,m,ii,jj,kk;

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

  // Get required info out of the exodus file.  Most important number of nodes stored.
  int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets, errr;
  char title [ MAX_LINE_LENGTH +1];
  errr = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);
  
  float *nodal_var_vals = (float *) calloc(num_nodes, sizeof(float));
  ITG   *node_map       = (ITG *)   calloc(num_nodes, sizeof(ITG));
  ITG   *inode_map      = (ITG *)   calloc(num_nodes, sizeof(ITG));
  
  errr = ex_get_id_map (exoid, EX_NODE_MAP, node_map);
  if(errr)printf("*ERROR in exo: failed to get prior node map");

  for (i=0; i<num_nodes; i++){
    // Seed with NaN for non-stored output
    nodal_var_vals[i]=nanf("");
  }
  
  {j=0;
    // Pull node map and create inverse
    // WHY DO WE NEED A INVERSE NODE MAP?
    // FRD doesn't require results to be sequential be we do!
    // FRD can be out of sequence (based on set order)
    // FRD tags each result with the associated node number.
    // Hence we need inverse and FRD does not.
    // Note also that changes in output sets, element deletion, etc
    // means the initial from exo.c must be recreated each time
    for(i=0;i<*nkcoords;i++){
      if(inum[i]<=0) continue;
      while(node_map[j%*nkcoords]!=i+1){j++;}
      inode_map[i] = (j%*nkcoords);
    }
  }
  
  int num_nod_vars=3;
  
  for (j=1; j<=num_nod_vars; j++){ // For each direction
    if(*iset==0){
      if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)){
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0) continue;
	  nodal_var_vals[inode_map[i]]=v[(mi[1]+1)*(i)+j];
	}
      }else{
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0) continue;

	  if(inotr[2*i]==0){
	    nodal_var_vals[inode_map[i]]=v[(mi[1]+1)*i+j];
	  }else{
	    ii=(mi[1]+1)*i+1;
	    jj=(mi[1]+1)*i+2;
	    kk=(mi[1]+1)*i+3;
	    FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	    nodal_var_vals[inode_map[i]]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
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
	      nodal_var_vals[inode_map[i]]=v[(mi[1]+1)*i+j];
	    }else{
	      // Transform into global directions
	      FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	      ii=(mi[1]+1)*i+1;
	      jj=(mi[1]+1)*i+2;
	      kk=(mi[1]+1)*i+3;
	      nodal_var_vals[inode_map[i]]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
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
		nodal_var_vals[inode_map[i]]=v[(mi[1]+1)*i+j];
	      }else{
		FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
		ii=(mi[1]+1)*i+1;
		jj=(mi[1]+1)*i+2;
		kk=(mi[1]+1)*i+3;
		nodal_var_vals[inode_map[i]]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
	      }
	    }
	  }while(1);
	}
      }
    }

    errr = ex_put_var (exoid, time_step, EX_NODAL, j+countvar, 1, num_nodes, nodal_var_vals);
    if (errr) printf ("ERROR storing vector data into exo file.\n");
  }
    
  free(node_map);
  free(inode_map);
  free(nodal_var_vals);
  return;

}
#endif
