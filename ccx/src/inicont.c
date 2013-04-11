/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

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
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

void inicont(int *ncont, int *ntie, char *tieset, int *nset, char *set,
               int *istartset, int *iendset, int *ialset, int **itietrip,
               char *lakon, int *ipkon, int *kon, int **koncontp,
	     int *ncone, double *tietol, int *ismallsliding){

  int *itietri=NULL,*koncont=NULL;

  itietri=*itietrip;koncont=*koncontp;

  FORTRAN(allocont,(ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,lakon,ncone,tietol,ismallsliding));

  if(ncont==0) return;

  itietri=NNEW(int,2**ntie);
  koncont=NNEW(int,4**ncont);

  FORTRAN(triangucont,(ncont,ntie,tieset,nset,set,istartset,iendset,
          ialset,itietri,lakon,ipkon,kon,koncont));

  *itietrip=itietri;*koncontp=koncont;
  
  return;
}
