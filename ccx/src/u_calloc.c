
/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2013 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <stdlib.h>
/*
 Diehl program
*/

void *u_calloc(size_t num,size_t size){

  void *a;
  if(num==0){
    a=NULL;
    return(a);
  }
      
  a=calloc(num,size);
  if(a==NULL){
    printf("*ERROR in u_calloc: error allocating memory\n");
    printf("num=%ld,size=%ld\n",num,size);
    exit(16);
  }
  else {
    return(a);
  }
}
