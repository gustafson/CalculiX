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

void exosetfind(char *set, ITG *nset, ITG *ialset, ITG *istartset, ITG *iendset,
		ITG *num_ns, ITG *num_ss, ITG *num_es, ITG *num_fs, ITG *node_map_inv,
		int exoid, int store, ITG *nk){

  ITG i,j,k,l,n,s,e,gen;
  char tmpstr[81];
  char *space = " ";
  char *pos;
  ITG *set_nums;
  int errr;
  int dropped=0, unidentified=0;

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
      set_nums = (ITG *) calloc(l, sizeof(ITG));

      /* Only add the generate code if there are at least
	 three points in the vector */
      n=0; j=s;
      if (l>2){
	while (j<=e-2){
	  // Account for generated ids
	  if (ialset[j+2]<0) {
	    for (k=ialset[j]; k<=ialset[j+1]; k-=ialset[j+2]){
	      set_nums[n++]=exoset_check(k-1, node_map_inv, nk, &dropped, &unidentified);
	    }
	    j+=3;
	  }else{
	    // Account for directly added id
	    gen=ialset[j++]-1;
	    set_nums[n++]=exoset_check(gen, node_map_inv, nk, &dropped, &unidentified);
	  }
	}
	// Must finish the last two of directly added set
	if (ialset[e]>0){ // only if the last set is not a generated set
	  // 1+n++ and -1+n++ to preserve order
	  set_nums[1+n++]=exoset_check(ialset[e]-2, node_map_inv, nk, &dropped, &unidentified);
	  if (ialset[e-1]>0){
	    set_nums[-1+n++]=exoset_check(ialset[e]-3, node_map_inv, nk, &dropped, &unidentified);
	  }
	}
      }else if(l>1){
	set_nums[n++]=exoset_check(ialset[s]-1, node_map_inv, nk, &dropped, &unidentified);
	set_nums[n++]=exoset_check(ialset[e]-1, node_map_inv, nk, &dropped, &unidentified);
      }else{
	set_nums[n++]=exoset_check(ialset[e]-1, node_map_inv, nk, &dropped, &unidentified);
      }
    } // if (store)

    strncpy(tmpstr,set+i*81,81);
    pos = strpbrk(tmpstr, space)-1;

    if(strcmp1(pos,"N")==0){
      (*num_ns)++; // printf ("Node set identified\n");
      if (store){
	errr = ex_put_set_param (exoid, EX_NODE_SET, i, n, 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	if (errr) printf ("Error writing set parameters\n");
	errr = ex_put_set       (exoid, EX_NODE_SET, i, set_nums, NULL);
	if (errr) printf ("Error writing set numbers\n");
	// ex_put_set_dist_fact (exoid, EX_NODE_SET, i, set_nums);  //
      }
    }
    if(strcmp1(pos,"E")==0) {
      (*num_es)++; // printf ("Element set identified\n");}
      /* No element set storage mechanism?
      if (store){
      	errr = ex_put_set_param (exoid, EX_ELEM_SET, i, n, 0);
      	if (errr) printf ("Error writing set parameters\n");
      	errr = ex_put_set       (exoid, EX_ELEM_SET, i, set_nums, NULL);
      	if (errr) printf ("Error writing set numbers\n");
      	// ex_put_set_dist_fact (exoid, EX_ELEM_SET, i, set_nums);  //
	} */
    }
    if(strcmp1(pos,"S")==0) {
      (*num_ss)++; // printf ("Node side set surface identified\n");}
      /* Side sets (node surfaces) not yet implemented
	 if (store){
      	errr = ex_put_set_param (exoid, EX_SIDE_SET, i, n, 0);
      	if (errr) printf ("Error writing set parameters\n");
      	errr = ex_put_set       (exoid, EX_SIDE_SET, i, set_nums, side_set_side_list);
      	if (errr) printf ("Error writing set numbers\n");
      	// ex_put_set_dist_fact (exoid, EX_SIDE_SET, i, set_nums);  //
	} */
    }
    if(strcmp1(pos,"T")==0) {
      (*num_fs)++; // printf ("Face set surface identified\n");}
      /* Face sets not yet implemented
	 if (store){
	 errr = ex_put_set_param (exoid, EX_FACE_SET, i, n, 0);
	 if (errr) printf ("Error writing set parameters\n");
	 errr = ex_put_set       (exoid, EX_FACE_SET, i, set_nums);
	 if (errr) printf ("Error writing set numbers\n");
	 // ex_put_set_dist_fact (exoid, EX_FACE_SET, i, set_nums);  //
	 } */
    }
    if (store) {free (set_nums);}
  }

  if (store){
    // char *namesnset[*num_ns]; j=0;
    // char *namessset[*num_ss]; k=0;
    // char *nameseset[*num_es]; l=0;
    // char *namesfset[*num_fs]; n=0;

    for (i=0; i<*nset; i++){
      strncpy(tmpstr,set+i*81,81);
      pos = strpbrk(tmpstr, space)-1;
      // This crashed in valgrind
      // if(strcmp1(pos,"N")==0) {strncpy(namesnset[j++], tmpstr, MAX_STR_LENGTH);}
      // if(strcmp1(pos,"E")==0) {strcpy(nameseset[l++],tmpstr);}
      // if(strcmp1(pos,"S")==0) {strcpy(namessset[k++],tmpstr);}
      // if(strcmp1(pos,"T")==0) {strcpy(namesfset[n++],tmpstr);}
    }

    // if (*num_ns>0){
    //   errr = ex_put_names (exoid, EX_NODE_SET, namesnset);
    //   if (errr) printf ("Error writing node set names\n");
    // }
    /* side sets not implemented yet
       if (*num_ss>0){
       errr = ex_put_names (exoid, EX_SIDE_SET, namessset);
       if (errr) printf ("Error writing side set names\n");
       }
    */
  }

  if (dropped){
    printf ("\nExodus Output WARNING: At least one node or element is dropped from a set.\n");
    printf ("  This may be due rigid bodies or 3D expansion (beams, shells, OUTPUT=3D).\n\n");
  }

  if (unidentified){
    printf ("\nExodus Output WARNING: At least one unidentified node or element is dropped from a set.\n\n");
  }

  return;
}


ITG exoset_check(ITG n, ITG *node_map_inv, ITG *nk, int *dropped, int *unidentified){
  ITG val=0;
  // printf ("%" ITGFORMAT ", %" ITGFORMAT "\n", n, *nk);
  if (n<=*nk){
    val = node_map_inv[n]-1;
    // printf ("val node number %" ITGFORMAT "\n", val+1);
    if (val==-1) {
      *dropped = 1;
    }
  }else{
    // printf("UNIDENTIFIED %" ITGFORMAT "\n", n);
    *unidentified = 1;
  }
  return val;
}

#endif
