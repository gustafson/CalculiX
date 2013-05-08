#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include "exodusII.h"

void exosetfind(char *set, int *nset, int *ialset, int *istartset, int *iendset, 
		int *num_ns, int *num_ss, int *num_es, int *num_fs, int *exoid, 
		int countonly){
  
  int i,j,n;
  char tmpstr[81];
  char *space = " ";
  char *pos;
  int *node_map;

  *num_ns = 0; 
  *num_ss = 0;
  *num_es = 0;
  *num_fs = 0; 
  
  for (i=0; i<*nset; i++){
    if (countonly) {
      // Find and store the set numbers
      if (iendset[i]<0){ // GENERATE WAS USED
	n=(ialset[istartset[i]]-ialset[iendset[i]-1])/(-ialset[iendset[i]])+1;
	printf ("Size n %i\n",n);
	node_map = (int *) calloc(n, sizeof(int));
	n=0;
	for (j=ialset[istartset[i]]; j<=ialset[iendset[i]-1]; j=j-ialset[iendset[i]]){
	  node_map[n++]=j;
	}
      }else{ // GENERATE WAS NOT USED
	n=istartset[i]-iendset[i]+1;
	printf ("Size n %i\n",n);
	node_map = (int *) calloc(n, sizeof(int));
	n=0;
	for (j=ialset[istartset[i]]; j<=ialset[iendset[i]]; j++){
	  node_map[n++]=j;
	}
      }
    }
    
    strncpy(tmpstr,set+i*81,81);
    pos = strpbrk(tmpstr, space)-1;
    if(strcmp1(pos,"N")==0){
      (*num_ns)++; // printf ("Node set identified\n"); 
      // ex_put_node_set_param (exoid, i, n, 0); // CURRENTLY NO DISTRIBUTIONS ADDED
      // ex_put_node_set       (exoid, i, node_map);
      // ex_put_node_set_dist_fact (exoid, i, node_map);  // 
    }    
    if(strcmp1(pos,"E")==0) {(*num_es)++;} // printf ("Element set identified\n");}
    if(strcmp1(pos,"S")==0) {(*num_ns)++;} // printf ("Node set surface identified\n");}
    if(strcmp1(pos,"T")==0) {(*num_fs)++;} // printf ("Face set surface identified\n");} 
    if (countonly) free (node_map);
  }  

  return;
}
