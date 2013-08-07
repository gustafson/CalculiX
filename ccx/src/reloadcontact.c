/*     CalculiX - A 3-dimensional finite element program                 */
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
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

void reloadcontact(char *lakon, int *ipkon, int *kon,
	       int **nelemloadp, char **sideloadp, int **iamloadp, 
	       double **xloadp, int *nload, int *ne, double *t1,
	       int *iamt1, int *nam, int *ithermal, double *vold,
	       int *mi, double **xloadoldp){

  char *sideload=NULL;

  int *nelemload=NULL,*iamload=NULL;

  double *xload=NULL,*xloadold=NULL;

  nelemload=*nelemloadp;sideload=*sideloadp;iamload=*iamloadp;xload=*xloadp;
  xloadold=*xloadoldp;
  
  /* reallocating the fields for distributed loading */
  
  /* each contact face is divided into 4 new faces */

  RENEW(nelemload,int,2**nload*5);
  RENEW(sideload,char,20**nload*5);
  if(*nam>0) RENEW(iamload,int,2**nload*5);
  RENEW(xload,double,2**nload*5);
  RENEW(xloadold,double,2**nload*5);

  FORTRAN(remeshload,(ipkon,kon,lakon,nelemload,sideload,iamload,xload,
		      nload,ne,t1,iamt1,nam,ithermal,vold,mi,xloadold));

  RENEW(xloadold,double,2**nload);
/*  RENEW(nelemload,int,2**nload);
  RENEW(sideload,char,20**nload);
  RENEW(iamload,int,2**nload);
  RENEW(xload,double,2**nload);*/

  *nelemloadp=nelemload;*sideloadp=sideload;*iamloadp=iamload;*xloadp=xload;
  *xloadoldp=xloadold;

  return;
}
