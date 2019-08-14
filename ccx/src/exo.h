/*     Calculix - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */
/*     This subroutine                                                   */
/*              Copyright (C) 2013-2018 Peter A. Gustafson               */
/*                                                                       */
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

void exo(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne0,
	 double *v,double *stn,ITG *inum,ITG *nmethod,ITG *kode,
	 char *filab,double *een,double *t1,double *fn,double *time,
	 double *epn,ITG *ielmat,char *matname,double *enern,
	 double *xstaten,ITG *nstate_,ITG *istep,ITG *iinc,
	 ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,
	 double *trab,ITG *inotr,ITG *ntrans,double *orab,
	 ITG *ielorien,ITG *norien,char *description,ITG *ipneigh,
	 ITG *neigh,ITG *mi,double *stx,double *vr,double *vi,
	 double *stnr,double *stni,double *vmax,double *stnmax,
	 ITG *ngraph,double *veold,double *ener,ITG *ne,double *cs,
	 char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
	 double *eenmax,double *fnr,double *fni,double *emn,
	 double *thicke,char *jobnamec,char *output,double *qfx,
         double *cdn,ITG *mortar,double *cdnr,double *cdni,ITG *nmat,
	 ITG *ielprop,double *prop);

void exosetfind(char *set, ITG *nset, ITG *ialset, ITG *istartset, ITG *iendset,
		ITG *num_ns, ITG *num_ss, ITG *num_es, ITG *num_fs, ITG *node_map_inv,
		int exoid, int store, ITG *nk);

ITG exoset_check(ITG n, ITG *node_map_inv, ITG *nk, int *dropped, int *unidentified);

void exovector(double *v,ITG *iset,ITG *ntrans,char * filabl,ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3, int exoid, ITG time_step,
	       int countvar, ITG nout, ITG *node_map_inv);

void exoselect(double *field1,double *field2,ITG *iset,ITG *nkcoords,ITG *inum,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ngraph,ITG *ncomp,
	       ITG *ifield,ITG *icomp,ITG *nfield,ITG *iselect,ITG exoid,
	       ITG time_step, int countvar, ITG nout, ITG *node_map_inv);
