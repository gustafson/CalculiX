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
  
  int i,j,k,l,n,s,e,gen;
  char tmpstr[81];
  char *space = " ";
  char *pos;
  int *node_map;
  int errr;

  *num_ns = 0; 
  *num_ss = 0;
  *num_es = 0;
  *num_fs = 0; 
  
  for (i=0; i<*nset; i++){

    if (store) {
      // Find and store the set numbers
      // The pointer integers are 1 based (fortran)
      s=istartset[i]-1;
      e=iendset[i]-1;
      // Determine if generate was used
      gen=0; l=0;
      for (j=s; j<=e; j++){
	if (ialset[j]<0) {
	  gen+=(ialset[j-1]-ialset[j-2])/(-ialset[j])+1;
	  l-=3;
	}
      }
      
      // Now set the length of the set allocation
      l=e-s+1+gen+l;
      node_map = (int *) calloc(l, sizeof(int));
      
      /* Only add the generate code if there are at least
	 three points in the vector */
      n=0; j=s;
      if (l>2){
	while (j<=e-2){
	  // Account for generated ids
	  if (ialset[j+2]<0) {
	    for (k=ialset[j]; k<=ialset[j+1]; k-=ialset[j+2]){
	      node_map[n++]=node_map_inv[k];  
	    }
	    j+=3;
	  }else{
	    // Account for directly added id
	    node_map[n++]=node_map_inv[ialset[j++]];
	  }
	}
	// Must finish the last two
	if (ialset[e]>0){
	  node_map[n++]=node_map_inv[ialset[e]];
	  if (ialset[e-1]>0){
	    printf ("Tacking on the end-1\n");
	    node_map[n++]=node_map_inv[ialset[e-1]];
	  }
	}
      }else if(l>1){
	node_map[n++]=node_map_inv[ialset[s]];
	node_map[n++]=node_map_inv[ialset[e]];
      }else{
	node_map[n++]=node_map_inv[ialset[e]];
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
    // char *namessset[*num_ss]; k=0;
    // char *nameseset[*num_es]; l=0;
    // char *namesfset[*num_fs]; n=0;

    for (i=0; i<*nset; i++){
      strncpy(tmpstr,set+i*81,81);
      pos = strpbrk(tmpstr, space)-1;
      // This crashed in valgrind
      if(strcmp1(pos,"N")==0) {strncpy(namesnset[j++], tmpstr, MAX_STR_LENGTH);}
      // if(strcmp1(pos,"E")==0) {strcpy(nameseset[l++],tmpstr);}
      // if(strcmp1(pos,"S")==0) {strcpy(namessset[k++],tmpstr);}
      // if(strcmp1(pos,"T")==0) {strcpy(namesfset[n++],tmpstr);}
    }
    
    if (*num_ns>0){
      errr = ex_put_names (exoid, EX_NODE_SET, namesnset);
      if (errr) printf ("Error writing node set names\n");
    }
    /* side sets not implemented yet
       if (*num_ss>0){
       errr = ex_put_names (exoid, EX_SIDE_SET, namessset);
       if (errr) printf ("Error writing side set names\n");
       }
    */
  }
  
  return;
}
