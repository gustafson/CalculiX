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

void exoselect(double *field1,double *field2,int *iset,int *nkcoords,int *inum,
	       int *istartset,int *iendset,int *ialset,int *ngraph,int *ncomp,
	       int *ifield,int *icomp,int *nfield,int *iselect,int *exoid,
	       int *istore, char *vname){

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

  float *nodal_var_vals;

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
  nodal_var_vals = (float *) calloc (*nkcoords, sizeof(float));
  
  for (j=1; j<=num_nod_vars; j++){ // For each direction
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
	    printf("%3s%10d","\nm1",i+1);
	    for(j=0;j<min(6,*ncomp);j++){
	      if(ifield[j]==1){
		printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
	      }else{
		printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
	      }
	    }
	  }else{
	    for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
	      if(ifield[j]==1){
		printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
	      }else{
		printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
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
		printf("%3s%10d","\nm1",i+1);
		for(j=0;j<min(6,*ncomp);j++){
		  if(ifield[j]==1){
		    printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
		  }else{
		    printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
		  }
		}
	      }else{
		for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		  if(ifield[j]==1){
		    printf("%12.5E",field1[i*nfield[0]+j]);
		  }else{
		    printf("%12.5E",field2[i*nfield[1]+j]);
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
		  printf("%3s%10d","\nm1",i+1);
		  for(j=0;j<min(6,*ncomp);j++){
		    if(ifield[j]==1){
		      printf("%12.5E",field1[i*nfield[0]+icomp[j]]);
		    }else{
		      printf("%12.5E",field2[i*nfield[1]+icomp[j]]);
		    }
		  }
		}else{
		  for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		    if(ifield[j]==1){
		      printf("%12.5E",field1[i*nfield[0]+j]);
		    }else{
		      printf("%12.5E",field2[i*nfield[1]+j]);
		    }
		  }
		}
	      }
	    }
	  }while(1);
	}
      }
    }
    
    printf ("\n");
    int errr = ex_put_nodal_var (exoid, istore, j, *nkcoords, nodal_var_vals);
    if (errr) printf ("ERROR storing vector data into exo file.\n");
  }  

  free(nodal_var_vals);
  return;

}
