/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                          */

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

void remcontmpc(int *nmpc, char *labmpc, int *mpcfree, int *nodempc,
            int *ikmpc, int *ilmpc, double *coefmpc, int *ipompc){

  /*   removes the contact MPC's */

  int i;

  for(i=*nmpc;i>0;i--){
      if(strcmp1(&labmpc[20*(i-1)],"CONTACT")==0){
	  FORTRAN(mpcrem,(&i,mpcfree,nodempc,nmpc,ikmpc,ilmpc,
			  labmpc,coefmpc,ipompc));
      }
  }

  return;

}
	  
