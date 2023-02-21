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

#if defined(_WIN32) || defined(_WIN64)

#include <math.h>

// strndup() is not available on Windows
char *strndup( const char *s1, size_t n)
{
    char *copy= (char*)malloc( n+1 );
    memcpy( copy, s1, n );
    copy[n] = 0;
    return copy;
};
#endif

#define NAMELEN 81
#define type_ns 0
#define type_es 1
#define type_ss 2
#define type_fs 3

void exosetfind(char *set, ITG *nset, ITG *ialset, ITG *istartset, ITG *iendset,
		ITG *num_ns, ITG *num_ss, ITG *num_es, ITG *num_fs, ITG *node_map_inv,
		int exoid, int store, ITG *nk){

  // Paraview cannot deal with empty sets in exodus files.  Thus we
  // need to count carefully the number of sets.

  // Note first call to this function counts the sets.  Second call to
  // the function (store=1) actually stores the sets.
  ITG i,l,s,e;
  
  int errr;
  int dropped=0;

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

    // This seems redudant to below and could be optimized
    char* tmpstr = strndup(pos0, (int) (pos1-pos0));
    names[i] = tmpstr;

    // Find type
    pos1 = strpbrk(pos1, space)-1;
    if(strcmp1(pos1,"N")==0){
      settype[i] = type_ns;
    }else if(strcmp1(pos1,"E")==0){
      settype[i] = type_es;
    }else if(strcmp1(pos1,"S")==0){
      settype[i] = type_ss;
    }else if(strcmp1(pos1,"T")==0){
      settype[i] = type_fs;
    }
  } // end loop i over all sets

  int warnempty = 1;
  int warnreverse = 0;
  // Get actual numbers for number of nodes and total number etc.
  for (i=0; i<*nset; i++){
    
    // ONLY WORKS FOR NSETS FOR NOW... We also need the element number inverse map
    if (settype[i] != type_ns){n_in_set[i]=0; continue;}
    
    // Now get the length of the set allocation (stored in l)
    exoset_count_set(ialset, istartset, iendset, &warnreverse, &i, &s, &e, &l);
    if (l==0){
      exoset_warn_empty(&warnempty);
      printf("\t\t- %s\n", names[i]);
      continue;
    }
    
    // The rest is done only if there is a non-empty set!
    // Allocate an array to store the set numbers
    ITG *set_nums;
    set_nums = (ITG *) calloc(l, sizeof(ITG));

    // Actually get the set.  The actual set size (n_in_set) may be
    // different than what was allocated (l) due to node drops and
    // shell expansions etc
    exoset_get_set(ialset, istartset, iendset, &i, &s, &l, &e, 
		   set_nums, dropped_set, nk, n_in_set,
		   node_map_inv, &dropped);

    // Store the results (but only after initializing the exodus file) when store=1
    if (n_in_set[i]>0){
      switch (settype[i])
	{
	case type_ns:
	  if (store){

	    // Add to the list of names.  This seems redundant to above and could be optimized
	    pos0 = set+i*NAMELEN;
	    pos1 = strpbrk(pos0, space)-1;
	    char* tmpstr = strndup(pos0, (int) (pos1-pos0));
	    names_nset[use_ns++] = tmpstr;

	    errr = ex_put_set_param (exoid, EX_NODE_SET, *num_ns, n_in_set[i], 0); // CURRENTLY NO DISTRIBUTIONS ADDED
	    if (errr) printf ("ERROR in exo: failed node set parameters\n");
	    errr = ex_put_set       (exoid, EX_NODE_SET, *num_ns, set_nums, NULL);
	    if (errr) printf ("ERROR in exo: failed node set\n");
	  }
	  (*num_ns)++;
	  continue;
	case type_es:
	case type_fs:
	case type_ss:
	  printf("\t- Element, face, and sides sets not implemented to exodus file.\n");
	  continue;
	}
    }else{ // T
    }
    free(set_nums);
  } //end i loop

  if (!store){
    // printf("FIRST PASS %i %i %i %i\n", *num_ns, *num_es, *num_ss, *num_fs);
    return;
  }

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

  printf("\t- Inactive nodes (unused or due to shell and beam expansion):");
  for(i=0; i<dropped; i++){
    if ((i%8)==0){printf("\n\t\t");}
    printf("%i, ", dropped_set[i]+1);
  }
  printf("\n");
  ex_put_names (exoid, EX_NODE_SET, names_nset);
  return;
}

ITG exoset_check_in_set(ITG n, ITG *node_map_inv, ITG *nk, int *dropped, int *unidentified){
  // Check if the numbered item is in the set.
  // Submitted should be an index which is zero based
  // Returned should be an index which is zero based
  ITG val=0;
  if (n<=*nk){
    val = node_map_inv[n]-1;
  } else {
    *unidentified = 1;
  }
  return val;
}

void exoset_count_set(ITG *ialset, ITG *istartset, ITG *iendset, int *warnreverse, ITG *i, ITG *s, ITG *e, ITG *l){
  // Count the length of the set.  Helps to determine how much memory
  // to allocate and whether to write a set into exodus (paraview
  // cannot handle empty sets)

  ITG j,k,n,gen;
  *s=istartset[*i]-1;
  *e=iendset[*i]-1;
  
  // Determine if generate was used
  gen=0; *l=0; n=1;
    
  for (j=*s; j<=*e; j++){
    if (ialset[j]<0){
      k=ialset[j-1]-ialset[j-2];
      if (k<0){
	if (n){
	  *warnreverse=1;
	  n=0;
	}
	k=-k;
      }
      // printf("k=%i\n",k);
      gen+=(k)/(-ialset[j])+1;
      *l-=3;
    }
  }
    
  // Now set the length of the set allocation
  *l=*e-*s+1+gen+*l;
}

void exoset_get_set(ITG *ialset, ITG *istartset, ITG *iendset,
		    ITG *i, ITG *s, ITG *l, ITG *e,
		    ITG *set_nums, int *dropped_set, ITG *nk,
		    int *n_in_set, ITG *node_map_inv, int *dropped){

  // Find and store the set numbers
  // The pointer integers are 1 based (fortran)

  ITG j=*s;
  ITG n,k,gen,z;
  int unidentified=0;
  
  n=0;
  
  /* Only add the generate code if there are at least
     three points in the vector */  
  if (*l>2){
    while (j<=*e-2){
      // generated ids
      if (ialset[j+2]<0) {
	if (ialset[j+1]-ialset[j]<0){
	  // Deal with reversed generated sets
	  for (k=ialset[j+1]; k<=ialset[j]; k-=ialset[j+2]){
	    z=exoset_check_in_set(k, node_map_inv, nk, dropped, &unidentified);
	    if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=k;}
	  }
	} else {
	  for (k=ialset[j]; k<=ialset[j+1]; k-=ialset[j+2]){ // Index arrays are fortran based (1)
	    z=exoset_check_in_set(k-1, node_map_inv, nk, dropped, &unidentified);
	    // printf("Generated is index k=%i with resulting node index %i\n", k-1, z);
	    if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=k;}
	  }
	}
	j+=3;
      } else {
	// Account for directly added id
	gen=ialset[j++];
	z=exoset_check_in_set(gen, node_map_inv, nk, dropped, &unidentified);
	// printf("Direct add %i with resulting node index %i\n", gen, z);
	if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=gen;}
      }
    }
    // Must finish the last two of directly added set
    if (ialset[*e]>0){ // only if the last set is not a generated set
      // 1+n++ and -1+n++ to preserve order
      z=exoset_check_in_set(ialset[*e]-2, node_map_inv, nk, dropped, &unidentified);
      if (z>=0){set_nums[1+n++]=z;}else{dropped_set[(*dropped)++]=ialset[*e]-2;}
      if (ialset[*e-1]>0){
	z=exoset_check_in_set(ialset[*e]-3, node_map_inv, nk, dropped, &unidentified);
	if (z>=0){set_nums[-1+n++]=z;}else{dropped_set[(*dropped)++]=ialset[*e]-3;}
      }
    }
  }else if(*l>1){
    // When a generated set is only of length 2.
    z=exoset_check_in_set(ialset[*s]-1, node_map_inv, nk, dropped, &unidentified);
    if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=ialset[*s]-1;}
    z=exoset_check_in_set(ialset[*e]-1, node_map_inv, nk, dropped, &unidentified);
    if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=ialset[*e]-1;}
  } else {
    z=exoset_check_in_set(ialset[*e]-1, node_map_inv, nk, dropped, &unidentified);
    if (z>=0){set_nums[n++]=z;}else{dropped_set[(*dropped)++]=ialset[*e]-1;}
  }
  n_in_set[*i]=n;
}

void exoset_warn_empty(int *warnempty){
  if (*warnempty){
    printf("\t- Empty set(s) skipped. ");
    printf("(Could be skipped due to shell and beam\n\t  expansions or non-use.) ");
    printf("Affected sets are:\n");
    *warnempty=0;
  }
}

#endif
