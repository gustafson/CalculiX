#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include "exodusII.h"

void exovector(double *v,int *iset,int *ntrans,char * filabl,int *nkcoords,
               int *inum,int *inotr,double *trab,double *co,
               int *istartset,int *iendset,int *ialset,int *mi,int *ngraph,
               int *exoid, int *istore, char *vname){
  
  int nksegment;
  int i,j,k,l,m,ii,jj,kk;
  
  float *nodal_var_vals;
  
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
  int num_nod_vars = 3;
  char *var_names[num_nod_vars];
  var_names[0] = "x";
  var_names[1] = "y";
  var_names[2] = "z";
  
  char tmpstrx[32],tmpstry[32],tmpstrz[32];
  char *x="x";
  char *y="y";
  char *z="z";
 
  x="x";
  strcpy (tmpstrx, vname);
  strcat (tmpstrx, x);
  var_names[0]=tmpstrx;
  
  y="y";
  strcpy (tmpstry, vname);
  strcat (tmpstry, y);
  var_names[1]=tmpstry;

  z="z";
  strcpy (tmpstrz, vname);
  strcat (tmpstrz, z);
  var_names[2]=tmpstrz;
  
  int errr;
  errr = ex_put_var_param (exoid, "n", num_nod_vars);
  errr = ex_put_var_names (exoid, "n", num_nod_vars, var_names);
  if (errr) printf ("ERROR in specifing the number of vars.\n"); 
  
  nodal_var_vals = (float *) calloc (*nkcoords, sizeof(float));
  
  for (j=1; j<=num_nod_vars; j++){ // For each direction
    if(*iset==0){
      if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)){
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0){
	    nodal_var_vals[i]=0;
	    continue;
	  }else{
	    nodal_var_vals[i]=v[(mi[1]+1)*i+j];
	  }
	}
      }else{
	for(i=0;i<*nkcoords;i++){
	  if(inum[i]<=0){
	    nodal_var_vals[i]=0;
	    continue;
	  }
	  if(inotr[2*i]==0){
	    nodal_var_vals[i]=v[(mi[1]+1)*i+j];
	  }else{
	    ii=(mi[1]+1)*i+1;
	    jj=(mi[1]+1)*i+2;
	    kk=(mi[1]+1)*i+3;
	    FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	    nodal_var_vals[i]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
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
	      nodal_var_vals[i]=v[(mi[1]+1)*i+j];
	    }else{
	      FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
	      ii=(mi[1]+1)*i+1;
	      jj=(mi[1]+1)*i+2;
	      kk=(mi[1]+1)*i+3;
	      nodal_var_vals[i]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
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
		// nodal_var_vals[i]=v[(mi[1]+1)*i+j];
	      }else{
		FORTRAN(transformatrix,(&trab[7*(inotr[2*i]-1)],&co[3*i],a));
		// ii=(mi[1]+1)*i+1;
		// jj=(mi[1]+1)*i+2;
		// kk=(mi[1]+1)*i+3;
		// nodal_var_vals[i]=v[ii]*a[0+(j-1)*3]+v[jj]*a[1+(j-1)*3]+v[kk]*a[2+(j-1)*3];
	      }
	    }
	  }while(1);
	}
      }
    }
    
    errr = ex_put_nodal_var (exoid, istore, j, *nkcoords, nodal_var_vals);
    if (errr) printf ("ERROR\n");
  }  

  free(nodal_var_vals);
  return;

}
