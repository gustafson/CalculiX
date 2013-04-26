/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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
#include "exodusII.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define rc(r,c) (r+*nkcoords*c)

void exoselect(double *field1,double *field2,int *iset,int *nkcoords,int *inum,
	       int *istartset,int *iendset,int *ialset,int *ngraph,int *ncomp,
	       int *ifield,int *icomp,int *nfield,int *iselect,int *exoid,
	       int *time_step, char *vname, int countvar){
    
  /* storing scalars, components of vectors and tensors without additional
     transformations */

  /* number of components in field1: nfield[0]
     number of components in field2: nfield[1]

     number of entities to store: ncomp
     for each entity i, 0<=i<ncomp:
     - ifield[i]: 1=field1,2=field2
     - icomp[i]: component: 0...,(nfield[0]-1 or nfield[1]-1) */
 
  int i,j,k,l,m,n,nksegment;
    
  int iw;

  float ifl;


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
  int num_nod_vars = 6;
  float *nodal_var_vals;
  nodal_var_vals = (float *) calloc (*nkcoords * num_nod_vars, sizeof(float));
  // float nodal_var_vals[*nkcoords][num_nod_vars];

  // typedef struct {
  //   int a;
  //   int b;
  //   float* data;
  // } float2d;
  
  // float2d arr2d;
  // arr2d={2,3};

  // float2d.data = malloc(arr2d.a * arr2d.b * sizeof *arr2d.data);
  
  // typedef struct {
  //   int a;
  //   int b;
  //   int* data;
  // } Int2d;

  // Int2d arr2d = { 2, 3 };
  // arr2d.data = malloc(arr2d.a * arr2d.b * sizeof *arr2d.data);

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
      
      for(n=1;n<=(int)((*ncomp+5)/6);n++){
	if(n==1){
	  // printf("%3s%10d","\nm1",i+1);
	  for(j=0;j<min(6,*ncomp);j++){
	    if(ifield[j]==1){
	      // printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
	      nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+icomp[j]];
	    }else{
	      // printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
	      nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+icomp[j]];
	    }
	  }
	}else{
	  for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
	    if(ifield[j]==1){
	      // printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
	      nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+icomp[j]];
	    }else{
	      // printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
	      nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+icomp[j]];
	    }
	  }
	}
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
	  for(n=1;n<=(int)((*ncomp+5)/6);n++){
	    if(n==1){
	      // printf("%3s%10d","\nm1",i+1);
	      for(j=0;j<min(6,*ncomp);j++){
		if(ifield[j]==1){
		  // printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
		  nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+icomp[j]];
		}else{
		  // printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
		  nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+icomp[j]];
		}
	      }
	    }else{
	      for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		if(ifield[j]==1){
		  // printf("%12.5E",field1[i*nfield[0]+j]);
		  nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+j];
		}else{
		  // printf("%12.5E",field2[i*nfield[1]+j]);
		  nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+j];
		}
	      }
	    }
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
	    
	    for(n=1;n<=(int)((*ncomp+5)/6);n++){
	      if(n==1){
		// printf("%3s%10d","\nm1",i+1);
		for(j=0;j<min(6,*ncomp);j++){
		  if(ifield[j]==1){
		    // printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
		    nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+icomp[j]];
		  }else{
		    // printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
		    nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+icomp[j]];
		  }
		}
	      }else{
		for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		  if(ifield[j]==1){
		    // printf("%12.5E",field1[i*nfield[0]+j]);
		    nodal_var_vals[rc(i,j)]=field1[i*nfield[0]+j];
		  }else{
		    // printf("%12.5E",field2[i*nfield[1]+j]);
		    nodal_var_vals[rc(i,j)]=field2[i*nfield[1]+j];
		  }
		}
	      }
	    }
	  }
	}while(1);
      }
    }
  }
    
  
  float *nodal_var_vals_out;
  nodal_var_vals_out = (float *) calloc (*nkcoords, sizeof(float));
  for (j=0; j<num_nod_vars; j++){
    for (i=0; i<*nkcoords; i++){
      nodal_var_vals_out[i]=nodal_var_vals[rc(i,j)];
    }
    int errr = ex_put_nodal_var (exoid, time_step, j+1+countvar, *nkcoords, nodal_var_vals_out);
    if (errr) printf ("ERROR storing data into exo file for dim %i record %i.\n", j, countvar+j);
  }

  free(nodal_var_vals);
  free(nodal_var_vals_out);
  return;

}
