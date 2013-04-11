/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <scsl_sparse.h>

void sgi_main(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, int *icol, int *irow, 
         int *neq, int *nzs, int token);

void sgi_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,int *icol, int *irow, 
                int *neq, int *nzs, int token);

void sgi_solve(double *b,int token);

void sgi_cleanup(int token);



