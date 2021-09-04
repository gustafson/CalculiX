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
  int dropped_set[*nk];
  
  char *names[*nset];
  // Individual set names are set after the initial count, so this
  // should work even when sized to zero in the first call
  char *names_nset[*num_ns];
  // char *names_eset[*num_es];
  // char *names_sset[*num_ss];
  // char *names_fset[*num_fs];

  char *space = " ";
  char *pos0;
  char *pos1;

  int use_ns=0;
  // int use_es=0;
  // int use_ss=0;
  // int use_fs=0;

  for (i=0; i<*nset; i++){
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
    // printf("FIRST PASS %i %i %i %i\n", *num_ns, *num_es, *num_ss, *num_fs);
    return;
  }

  int warnempty = 1;
  int warnreverse = 0;
  // Get actual numbers for number of nodes and total number etc.
  for (i=0; i<*nset; i++){
    
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
	// printf("Generating set for set %i for set %s.\n", i, names[i]);
	k=ialset[j-1]-ialset[j-2];
	if (k<0){
	  if (n){
	    warnreverse=1;
	    n=0;
	  }
	  k=-k;
	}
	// printf("k=%i\n",k);
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
	// generated ids
	if (ialset[j+2]<0) {
	  if (ialset[j+1]-ialset[j]<0){
	    // Deal with reversed generated sets
	    for (k=ialset[j+1]; k<=ialset[j]; k-=ialset[j+2]){
	      z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
	      if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
	    }
	  } else {
	    for (k=ialset[j]; k<=ialset[j+1]; k-=ialset[j+2]){ // Index arrays are fortran based (1)
	      z=exoset_check(k-1, node_map_inv, nk, &dropped, &unidentified);
	      // printf("Generated is index k=%i with resulting node index %i\n", k-1, z);
	      if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
	    }
	  }
	  j+=3;
	} else {
	  // Account for directly added id
	  k=ialset[j++];
	  z=exoset_check(gen, node_map_inv, nk, &dropped, &unidentified);
	  // printf("Direct add %i with resulting node index %i\n", gen, z);
	  if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
	}
      }
      // Must finish the last two of directly added set
      if (ialset[e]>0){ // only if the last set is not a generated set
	// 1+n++ and -1+n++ to preserve order
	k=ialset[e]-2;
	z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
	if (z>=0){set_nums[1+n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
	if (ialset[e-1]>0){
	  k=ialset[e]-3;
	  z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
	  if (z>=0){set_nums[-1+n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
	}
      }
    }else if(l>1){
      // When a generated set is only of length 2.
      k=ialset[s]-1;
      z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
      k=ialset[e]-1;
      z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
    } else {
      k=ialset[e]-1;
      z=exoset_check(k, node_map_inv, nk, &dropped, &unidentified);
      if (z>=0){set_nums[n++]=z;}else{exoset_dups(dropped_set, &dropped, k);}
    }
    n_in_set[i]=n;

    // // DEBUG
    // printf("Set %s has %i members ",names[i], n);
    // for (j=0; j<n; j++){printf("%i, ",set_nums[j]);}
    // printf("\n");
    
    // Write the number of sets
    if (n_in_set[i]>0){
      switch (settype[i])
	{
	case type_ns:
	  errr = ex_put_set_param (exoid, EX_NODE_SET, use_ns,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	  if (errr) printf ("ERROR in exo: failed node set parameters\n");
	  errr = ex_put_set       (exoid, EX_NODE_SET, use_ns++, set_nums, NULL);
	  if (errr) printf ("ERROR in exo: failed node set\n");
	  break;
	case type_es:
	  // printf("Exodus Warning: Element sets not implemented. Affected set is %s\n", names[i]);
	  // I haven't figured out how to implement element sets which I think must be based on blocks.
	  
	  // errr = ex_put_set_param (exoid, EX_ELEM_SET, use_es,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	  // if (errr) printf ("ERROR in exo: failed elem set parameters\n");
	  // errr = ex_put_set       (exoid, EX_ELEM_SET, use_es++, set_nums, NULL);
	  // if (errr) printf ("ERROR in exo: failed elem set\n");
	  // break;
	case type_fs:
	  // printf("Exodus Warning: Face sets not implemented. Affected set is %s\n", names[i]);
	  
	  // errr = ex_put_set_param (exoid, EX_FACE_SET, use_fs,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	  // if (errr) printf ("ERROR in exo: failed face set parameters\n");
	  // errr = ex_put_set       (exoid, EX_FACE_SET, use_fs++, set_nums, NULL);
	  // if (errr) printf ("ERROR in exo: failed face set\n");
	  // break
	case type_ss:
	  // printf("Exodus Warning: Face sets not implemented. Affected set is %s\n", names[i]);
	  
	  // errr = ex_put_set_param (exoid, EX_SIDE_SET, use_ss,   n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	  // if (errr) printf ("ERROR in exo: failed side set parameters\n");
	  // errr = ex_put_set       (exoid, EX_SIDE_SET, use_ss++, set_nums, NULL);
	  // if (errr) printf ("ERROR in exo: failed side set\n");
	  continue;
	}
    }else{
      if (warnempty){
	printf("\t- Empty set(s) skipped. ");
	printf("(Due to shell and beam\n\t  expansions or non-use.) ");
	printf("Affected sets are:\n");
	warnempty=0;
      }
      printf("\t\t- %s\n", names[i]);
    }
    free(set_nums);
  } //end i loop


  printf("\t- Element, and surface sets not written to .exo.\n\t  Affected sets are:\n");
  // printf("\t- Element, face, and side sets not written to .exo.\n\t  Affected sets are:\n");
  for(i=0; i<*nset; i++){
    if (settype[i]!=type_ns){
      printf("\t\t- %s\n", names[i]);
    }
  }
  
  if (warnreverse){
    printf("\t- Found a generated set with decreasing numbers.\n");
    printf("\t  These numbers will be reversed. Check the input deck.\n");
  }

  exodus_sort(dropped_set, &dropped);
  printf("\t- Inactive nodes (unused or due to shell and beam expansion):");
  for(i=0; i<dropped; i++){
    if ((i%8)==0){printf("\n\t\t");}
    printf("%i, ", dropped_set[i]+1);
  }
  printf("\n");
  ex_put_names (exoid, EX_NODE_SET, names_nset);
  return;
}

ITG exoset_check(ITG n, ITG *node_map_inv, ITG *nk, int *dropped, int *unidentified){
  // Submitted should be an index which is zero based
  // Returned should be an index which is zero based
  ITG val=0;
  // printf ("%" ITGFORMAT ", %" ITGFORMAT "\n", n, *nk);
  // printf("nk = %i\n", *nk);
  // DEBUG // int count=0;
  // DEBUG // for (ITG z=0; z<*nk; z++){
  // DEBUG //   if (node_map_inv[z]!=0){count++;
  // DEBUG //     printf("Node map input=%i, z=%i, inv=%i\n", n, z, node_map_inv[z]);
  // DEBUG //   }
  // DEBUG // };
  
  if (n<=*nk){
    val = node_map_inv[n]-1;
  } else {
    *unidentified = 1;
  }

  return val;
}

void exodus_sort(ITG *a, ITG *n){
  ITG i, j, tmp;
  
  for(i=0; i<*n; i++){
    for(j=i+1; j<*n; j++){
      if(a[j] <a[i]){
	tmp = a[i];
	a[i] = a[j];
	a[j] = tmp;
      }
    }
  }
  return;
}

void exoset_dups(ITG *set, ITG *n, ITG k){
  for (ITG i=0; i<*n; i++){
    if (set[i]==k){return;};
  }
  set[(*n)++] = k;
}
#endif
