/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */
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

#define NAMELEN 81
#define type_ns 0
#define type_es 1
#define type_ss 2
#define type_fs 3

void exosetfind(char *set, ITG *nset, ITG *ialset, ITG *istartset, ITG *iendset,
		ITG *num_ns, ITG *num_ss, ITG *num_es, ITG *num_fs, ITG *node_map_inv,
		int exoid, int store, ITG *nk){

  // Note first call to this function counts the sets.  Second call to
  // the function (store=1) actually stores the sets.

  ITG i,j,k,l,n,s,e,gen,z;
  
  int errr;
  int dropped=0, unidentified=0;

  int settype[*nset];
  int n_in_set[*nset];

  char *names[*nset];
  // Individual set names are set after the initial count, so this
  // should work even when sized to zero in the first call
  char *names_nset[*num_ns];
  char *names_eset[*num_es];
  char *names_sset[*num_ss];
  char *names_fset[*num_fs];

  char *space = " ";
  char *pos0;
  char *pos1;

  int use_ns=0;
  int use_es=0;
  int use_ss=0;
  int use_fs=0;
    
  for (int i=0; i<*nset; i++){
    // set names are stored in set, and appear to be 80 characters in
    // length, but is deliminated by a space
    pos0 = set+i*NAMELEN;
    pos1 = strpbrk(pos0, space)-1;
    int strl = (int) (pos1-pos0);

    char* tmpstr = strndup(pos0, strl);
    names[i] = tmpstr;
    pos1 = strpbrk(pos1, space)-1;
    if(strcmp1(pos1,"N")==0){ // printf("Found node set\n");
      settype[i] = type_ns;
      char* tmpstr = strndup(pos0, strl);
      if (store){names_nset[use_ns++] = tmpstr;} else {(*num_ns)++;}
    }else if(strcmp1(pos1,"E")==0){ // printf("Found element set\n");
      settype[i] = type_es;
      // char* tmpstr = strndup(pos0, strl);
      // if (store){names_eset[use_es++] = tmpstr;} else {(*num_es)++;}
    }else if(strcmp1(pos1,"S")==0){ // printf("Found side set\n");
      settype[i] = type_ss;
      // char* tmpstr = strndup(pos0, strl);
      // if (store){names_sset[use_ss++] = tmpstr;} else {(*num_ss)++;}
    }else if(strcmp1(pos1,"T")==0){ // printf("Found face set\n");
      settype[i] = type_fs;
      // char* tmpstr = strndup(pos0, strl);
      // if (store){names_fset[use_fs++] = tmpstr;} else {(*num_fs)++;}
    }
  } // end loop i over all sets


  if (store==0){
    printf("FIRST PASS %i %i %i %i\n", *num_ns, *num_es, *num_ss, *num_fs);
    return;
  }

  // Get actual numbers for number of nodes and total number etc.
  for (i=0; i<*nset; i++){
    
    // printf("Assessing set %s with type %i\n", names[i], settype[i]);
    // ONLY WORKS FOR NSETS FOR NOW... We also need the element number inverse map
    if (settype[i] != type_ns){n_in_set[i]=0; continue;}
    
    // Find and store the set numbers
    // The pointer integers are 1 based (fortran)
    s=istartset[i]-1;
    e=iendset[i]-1;
    // Determine if generate was used
    gen=0; l=0; n=1;
    
    for (j=s; j<=e; j++){
      if (ialset[j]<0){
	// This is a generated set
	k=ialset[j-1]-ialset[j-2];
	if (k<0){
	  if (n){
	    printf("Warning: Exodus deduced a generated set with decreasing numbers.\n");
	    printf("         These numbers will be reverse. Check the input deck.\n");
	    n=0;
	  }
	  k=-k;
	}
	gen+=(k)/(-ialset[j])+1;
	l-=3;
      }
    }
    
    
    // Now set the length of the set allocation
    l=e-s+1+gen+l;
    ITG *set_nums;
    set_nums = (ITG *) calloc(l, sizeof(ITG));
    
    /* Only add the generate code if there are at least
       three points in the vector */
    n=0; j=s;
    if (l>2){
      while (j<=e-2){
	// Account for generated ids
	if (ialset[j+2]<0) {
	  if (ialset[j+1]-ialset[j]<0){
	    // Deal with reversed generated sets
	    printf("REversal\n");
	    for (k=ialset[j+1]; k<=ialset[j]; k-=ialset[j+2]){
	      z=exoset_check(k-1, node_map_inv, nk, &dropped, &unidentified);
	      if (z>=0){set_nums[n++]=z;}
	    }
	  } else {
	    for (k=ialset[j]; k<=ialset[j+1]; k-=ialset[j+2]){
	      // printf("generated set ialset[k]=%i\n", ialset[k]);
	      z=exoset_check(k-1, node_map_inv, nk, &dropped, &unidentified);
	      if (z>=0){set_nums[n++]=z;}
	    }
	  }
	  j+=3;
	} else {
	  // Account for directly added id
	  gen=ialset[j++]-1;
	  z=exoset_check(gen, node_map_inv, nk, &dropped, &unidentified);
	  if (z>=0){set_nums[n++]=z;}
	}
      }
      // Must finish the last two of directly added set
      if (ialset[e]>0){ // only if the last set is not a generated set
	// 1+n++ and -1+n++ to preserve order
	z=exoset_check(ialset[e]-2, node_map_inv, nk, &dropped, &unidentified);
	if (z>=0){set_nums[1+n++]=z;}
	if (ialset[e-1]>0){
	  z=exoset_check(ialset[e]-3, node_map_inv, nk, &dropped, &unidentified);
	  if (z>=0){set_nums[-1+n++]=z;}
	}
      }
    }else if(l>1){
      // When a generated set is only of length 2.
      z=exoset_check(ialset[s]-1, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}
      z=exoset_check(ialset[e]-1, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}
    } else {
      z=exoset_check(ialset[e]-1, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}
    }
    n_in_set[i]=n;

    // // DEBUG
    // printf("Set %s has %i members ",names[i], n);
    // for (j=0; j<n; j++){printf("%i, ",set_nums[j]);}
    // printf("\n");
    
    // Write the number of sets
    if (n_in_set[i]>0){
      if (settype[i]==type_ns){
	errr = ex_put_set_param (exoid, EX_NODE_SET, use_ns,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	if (errr) printf ("ERROR in exo: failed node set parameters\n");
	errr = ex_put_set       (exoid, EX_NODE_SET, use_ns++, set_nums, NULL);
	if (errr) printf ("ERROR in exo: failed node set\n");
      }else if (settype[i]==type_es){
	1;
	// printf("Exodus Warning: Element sets not implemented. Affected set is %s\n", names[i]);
	// I haven't figured out how to implement element sets which I think must be based on blocks.
	
	// errr = ex_put_set_param (exoid, EX_ELEM_SET, use_es,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	// if (errr) printf ("ERROR in exo: failed elem set parameters\n");
	// errr = ex_put_set       (exoid, EX_ELEM_SET, use_es++, set_nums, NULL);
	// if (errr) printf ("ERROR in exo: failed elem set\n");
      }else if (settype[i]==type_fs){
	1;
	// printf("Exodus Warning: Face sets not implemented. Affected set is %s\n", names[i]);
	
	// errr = ex_put_set_param (exoid, EX_FACE_SET, use_fs,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	// if (errr) printf ("ERROR in exo: failed face set parameters\n");
	// errr = ex_put_set       (exoid, EX_FACE_SET, use_fs++, set_nums, NULL);
	// if (errr) printf ("ERROR in exo: failed face set\n");
      }else if (settype[i]==type_ss){
	1;
	// printf("Exodus Warning: Face sets not implemented. Affected set is %s\n", names[i]);
	
	// errr = ex_put_set_param (exoid, EX_SIDE_SET, use_ss,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	// if (errr) printf ("ERROR in exo: failed side set parameters\n");
	// errr = ex_put_set       (exoid, EX_SIDE_SET, use_ss++, set_nums, NULL);
	// if (errr) printf ("ERROR in exo: failed side set\n");
      }
    }else{
      printf("Exodus Warning: Empty set skipped: %s\n", names[i]);
    }
    free(set_nums);
  } //end i loop


  /// Last thing is to issue warnings
  printf("Exodus Warning: element, face, and side sets not implemented. Affected sets are:\n");
  int warnempty=0;
  for(i=0; i<*nset; i++){
    if (settype[i]!=type_ns){
      printf(" %s", names[i]);
      if (n_in_set[i]==0){
	warnempty=1;
      }
    }
  }
  printf("\n");
  // if (warnempty){
  //   printf("Exodus Warning: empty sets not saved. Affected sets are:\n");
  //   for(i=0; i<*nset; i++){
  //     if (settype[i]!=type_ns){
  // 	if (n_in_set[i]==0){
  // 	  printf(" %s", names[i]);
  // 	}
  //     }
  //   }
  //   printf("\n");
  // }
  
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
  } else {
    *unidentified = 1;
  }
  return val;
}

#endif
