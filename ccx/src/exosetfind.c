#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include "exodusII.h"

void exosetfind(char *set, int *nset, int *ialset, int *istartset, int *iendset, 
		int *num_ns, int *num_ss, int *num_es, int *num_fs, int *node_map_inv,
		int exoid, 
		int store){
  
  int i,j,k,l,n;
  char tmpstr[81];
  char *space = " ";
  char *pos;
  int *node_map;
  int errr;
  char *names[*nset];

  *num_ns = 0; 
  *num_ss = 0;
  *num_es = 0;
  *num_fs = 0; 
  

  for (i=0; i<*nset; i++){
    if (store) {
      // Find and store the set numbers
      if (ialset[iendset[i]]<0){ // GENERATE WAS USED
	n=(ialset[istartset[i]]-ialset[iendset[i]-1])/(-ialset[iendset[i]])+1;
	node_map = (int *) calloc(n, sizeof(int));
	n=0;
	for (j=ialset[istartset[i]]; j<=ialset[iendset[i]-1]; j=j-ialset[iendset[i]]){
	  node_map[n++]=node_map_inv[j-1];
	}
      }else{ // GENERATE WAS NOT USED
	n=istartset[i]-iendset[i]+1;
	node_map = (int *) calloc(n, sizeof(int));
	n=0;
	for (j=ialset[istartset[i]]; j<=ialset[iendset[i]]; j++){
	  node_map[n++]=node_map_inv[j-1];
	}
      }
    }
    
    strncpy(tmpstr,set+i*81,81);
    pos = strpbrk(tmpstr, space)-1;
    
    if(strcmp1(pos,"N")==0){
      (*num_ns)++; // printf ("Node set identified\n"); 
      if (store){ 
	errr = ex_put_node_set_param (exoid, i, n, 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	if (errr) printf ("Error writing set parameters\n");
	errr = ex_put_node_set       (exoid, i, node_map);
	if (errr) printf ("Error writing set numbers\n");
	// ex_put_node_set_dist_fact (exoid, i, node_map);  // 
      }
    }    
    if(strcmp1(pos,"E")==0) {
      (*num_es)++; // printf ("Element set identified\n");}
      /* No element set storage mechanism?
      if (store){ 
      	errr = ex_put_elem_set_param (exoid, i, n, 0);
      	if (errr) printf ("Error writing set parameters\n");
      	errr = ex_put_elem_set       (exoid, i, node_map);
      	if (errr) printf ("Error writing set numbers\n");
      	// ex_put_elem_set_dist_fact (exoid, i, node_map);  // 
	} */
    } 
    if(strcmp1(pos,"S")==0) {
      (*num_ss)++; // printf ("Node side set surface identified\n");}
      /* Side sets (node surfaces) not yet implemented
	 if (store){ 
      	errr = ex_put_side_set_param (exoid, i, n, 0);
      	if (errr) printf ("Error writing set parameters\n");
      	errr = ex_put_side_set       (exoid, i, node_map);
      	if (errr) printf ("Error writing set numbers\n");
      	// ex_put_side_set_dist_fact (exoid, i, node_map);  // 
	} */
    } 
    if(strcmp1(pos,"T")==0) {
      (*num_fs)++; // printf ("Face set surface identified\n");} 
      /* Face sets not yet implemented
	 if (store){ 
	 errr = ex_put_face_set_param (exoid, i, n, 0);
	 if (errr) printf ("Error writing set parameters\n");
	 errr = ex_put_face_set       (exoid, i, node_map);
	 if (errr) printf ("Error writing set numbers\n");
	 // ex_put_face_set_dist_fact (exoid, i, node_map);  // 
	 } */
    } 
    if (store) {free (node_map);}
    
  }  
  
  if (store){
    char *namesnset[*num_ns]; j=0;
    char *namessset[*num_ss]; k=0;
    char *nameseset[*num_es]; l=0;
    char *namesfset[*num_fs]; n=0;
    for (i=0; i<*nset; i++){
      strncpy(tmpstr,set+i*81,81);
      pos = strpbrk(tmpstr, space)-1;
      // This doesn't really work yet... don't know why.
      if(strcmp1(pos,"N")==0) {namesnset[j++]=&names[i];}
      if(strcmp1(pos,"E")==0) {nameseset[l++]=&names[i];}
      if(strcmp1(pos,"S")==0) {namessset[k++]=&names[i];}
      if(strcmp1(pos,"T")==0) {namesfset[n++]=&names[i];}
      // if(strcmp1(pos,"N")==0) {strcpy(namesnset[j++],names[i]);}
      // if(strcmp1(pos,"E")==0) {strcpy(nameseset[l++],names[i]);}
      // if(strcmp1(pos,"S")==0) {strcpy(namessset[k++],names[i]);}
      // if(strcmp1(pos,"T")==0) {strcpy(namesfset[n++],names[i]);}
    }
    
    if (*num_ns>0){
      errr = ex_put_names (exoid, EX_NODE_SET, namesnset);
      if (errr) printf ("Error writing node set names\n");
    }
    if (*num_ss>0){
      errr = ex_put_names (exoid, EX_SIDE_SET, namessset);
      if (errr) printf ("Error writing side set names\n");
    }
  }
  
  return;
}
