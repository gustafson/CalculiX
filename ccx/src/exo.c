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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#ifdef EXODUSII
#include "exodusII.h"
#include "exo.h"
#endif

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

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
	 ITG *ielprop,double *prop){

#ifdef EXODUSII
  /* stores the results in exo format

     iselect selects which nodes are to be stored:
     iselect=-1 means only those nodes for which inum negative
     ist, i.e. network nodes
     iselect=+1 means only those nodes for which inum positive
     ist, i.e. structural nodes
     iselect=0  means both of the above */

  /* Note for frd strcmp1(output,"asc")==0 defines ascii file, otherwise binary */

  char fneig[132]="",
    material[59]="                                                          ",
    text[2]=" ";

  /* For call compatibility with frd */
  char m1[4]="";
  char m3[4]="";
  FILE *f1=NULL;

  static ITG nkcoords,nout,noutmin,noutplus,iaxial;
  ITG nterms;

  ITG i,j,k,l,m,n,o,indexe,nemax,nlayer,noutloc,iset,iselect,ncomp,nope,
    nodes,ifield[7],nfield[2],icomp[7],ifieldstate[*nstate_],icompstate[*nstate_],nelout;

  ITG mt=mi[1]+1;

  ITG ncompscalar=1,ifieldscalar[1]={1},icompscalar[1]={0},nfieldscalar[2]={1,0};
  ITG ncompvector=3,ifieldvector[3]={1,1,1},icompvector[3]={0,1,2},nfieldvector1[2]={3,0},nfieldvector0[2]={mi[1]+1,0},icompvectorlast[3]={3,4,5};
  ITG ncomptensor=6,ifieldtensor[6]={1,1,1,1,1,1},icomptensor[6]={0,1,2,3,5,4},nfieldtensor[2]={6,0};
  ITG ncompscalph=2,ifieldscalph[2]={1,2},icompscalph[2]={0,0},nfieldscalph[2]={0,0};
  ITG ncompvectph=6,ifieldvectph[6]={1,1,1,2,2,2},icompvectph[6]={1,2,3,1,2,3},nfieldvectph[2]={mi[1]+1,mi[1]+1};
  ITG ncomptensph=12,ifieldtensph[12]={1,1,1,1,1,1,2,2,2,2,2,2},icomptensph[12]={0,1,2,3,5,4,0,1,2,3,5,4},nfieldtensph[2]={6,6};

  double *vold=NULL;
 
  int errr=0, exoid=0;
  ITG num_dim, num_elem;
  ITG num_elem_blk; /* Node element blocks.  One eltype per block*/
  ITG num_ns, num_ss, num_es, num_fs; /* Node sets, side sets, element sets, face sets */
  int CPU_word_size = sizeof(float);
  int IO_word_size = sizeof(float);

  /* Filename */
  strcpy (fneig, jobnamec);
  strcat (fneig, ".exo");

  /* nkcoords is the number of nodes at the time when
     the nodal coordinates are stored in the exo file */
  nkcoords = *nk;
  ITG num_nodes = nkcoords;

  /* determining nout, noutplus and noutmin
     nout: number of structural and network nodes
     noutplus: number of structural nodes
     noutmin: number of network nodes */
  if(*nmethod!=0){
    nout=0;
    noutplus=0;
    noutmin=0;
    for(i=0;i<*nk;i++){
      if(inum[i]==0) continue;
      if(inum[i]>0) noutplus++;
      if(inum[i]<0) noutmin++;
      nout++;
    }
  }else{
    nout=*nk;
  }

  // Allocate memory for node positions
  float *x, *y, *z;
  x = (float *) calloc(nout, sizeof(float));
  y = (float *) calloc(nout, sizeof(float));
  z = (float *) calloc(nout, sizeof(float));

  // Write optional node map
  j = 0; // Counter for the exo order of the nodes
  ITG *node_map,*node_map_inv;
  node_map = (ITG *) calloc(nout, sizeof(ITG));
  node_map_inv = (ITG *) calloc(nkcoords, sizeof(ITG));
  /* storing the coordinates of the nodes */
  if(*nmethod!=0){
    for(i=0;i<*nk;i++){
      if(inum[i]==0){continue;}
      // The difference between i and j is that not all values of i
      // increment j.
      node_map[j] = i+1;
      node_map_inv[i] = j+1;
      x[j]   = co[3*i];
      y[j]   = co[3*i+1];
      z[j++] = co[3*i+2];
    }
  }else{
    for(i=0;i<*nk;i++){
      node_map[j] = i+1;
      node_map_inv[i] = j+1;
      x[j]   = co[3*i];
      y[j]   = co[3*i+1];
      z[j++] = co[3*i+2];
    }
  }

  /* first time something is written in the exo-file: store
     computational metadata, the nodal coordinates and the
     topology */
  if((*kode==1)&&((*nmethod!=5)||(*mode!=0))){
    iaxial=0.;
    num_dim = 3;
    num_elem_blk = 17;

    // Find the number of sets
    exosetfind(set, nset, ialset, istartset, iendset,
	       &num_ns, &num_ss, &num_es, &num_fs, NULL, exoid, (int) 0, nk);
    printf ("Side sets to exo file not implemented.\n");
    num_ss=0;

#ifdef LONGLONG
    // This handles LONGLONG transparently
    remove(fneig);
    exoid = ex_create (fneig, /*Filename*/
    		       EX_ALL_INT64_API,	/* create mode */
    		       &CPU_word_size,  /* CPU float word size in bytes */
    		       &IO_word_size);  /* I/O float word size in bytes */
#else
    exoid = ex_create (fneig, /*Filename*/
		       EX_CLOBBER,	/* create mode */
		       &CPU_word_size,  /* CPU float word size in bytes */
		       &IO_word_size);  /* I/O float word size in bytes */
#endif

    /* determining the number of elements */
    if(*nmethod!=0){
      nelout=0;
      for(i=0;i<*ne0;i++){
	if(ipkon[i]<=-1){
	  continue;
	}else if(strcmp1(&lakon[8*i],"ESPRNGC")==0){
	  continue;
	}else if(strcmp1(&lakon[8*i],"ESPRNGF")==0){
	  continue;
	}else if(strcmp1(&lakon[8*i],"DCOUP3D")==0){
	  continue;
	}else if(strcmp2(&lakon[8*i+6],"LC",2)==0){
	  // Count the number of layers
	  nlayer=0;
	  for(k=0;k<mi[2];k++){
	    if(ielmat[i*mi[2]+k]==0) break;
	    nlayer++;
	  }
	  // Allow an element for each layer
	  for(k=0;k<nlayer;k++){
	    nelout++;
	  }
	}else{
	  nelout++;
	}
      }
    }else{
      nelout=*ne;
    }
    num_elem = nelout;

    /* initialize file with parameters */
    printf("\nData writen to the .exo file\n");
    num_nodes=nout;
    printf("Number of nodes: %" ITGFORMAT "\n", num_nodes);
    printf("Number of elements %" ITGFORMAT "\n", num_elem);
    printf("Number of element blocks %" ITGFORMAT "\n", num_elem_blk);
    printf("Number of node sets %" ITGFORMAT "\n", num_ns);
    printf("Number of side sets %" ITGFORMAT "\n", num_ss);

    errr = ex_put_init (exoid, "CalculiX EXO File",
			num_dim, num_nodes,
			num_elem, num_elem_blk,
			num_ns, num_ss);

    // Write values to database
    errr = ex_put_coord (exoid, x, y, z);
    if(errr)printf("*ERROR in exo: failed node positions");
    errr = ex_put_id_map (exoid, EX_NODE_MAP, node_map);
    if(errr)printf("*ERROR in exo: failed node map");

    // Deallocate
    free (node_map);
    free (x);
    free (y);
    free (z);

    // Write coordinate names
    char *coord_names[3];
    coord_names[0] = "x";
    coord_names[1] = "y";
    coord_names[2] = "z";
    errr = ex_put_coord_names (exoid, coord_names);
    if(errr){printf("*ERROR in exo: failed coordinate names");}

    // Initialize enough memory to store the element numbers
    ITG *elem_map;
    elem_map = (ITG *) calloc(num_elem, sizeof(ITG));
    char curblk[6]="     ";

    ITG *blkassign;
    blkassign = (ITG *) calloc(num_elem, sizeof(ITG));

    l=0;
    for(i=0;i<*ne0;i++){ // For each element.  Composite elements are
			 // one increment in this loop and all layers
			 // have the same element number.
      if(ipkon[i]<=-1){
	continue;
      }else if(strcmp1(&lakon[8*i],"F")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"ESPRNGC")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"ESPRNGF")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"DCOUP3D")==0){
	continue;
      }else{
	indexe=ipkon[i];
      }

      elem_map[l] = i+1;

      strcpy1(curblk,&lakon[8*i],5);
      // strcpy1(curblk,&lakon[8*i],10);
      // printf ("%s\n", curblk);
      strcpy1(material,&matname[80*(ielmat[i*mi[2]]-1)],5);
      // printf ("TODO store material identifier and name.\n");

      // Identify element type
      if(strcmp1(&lakon[8*i+3],"2")==0){
	/* 20-node brick element */
	if(((strcmp1(&lakon[8*i+6]," ")==0)||
	    (strcmp1(&filab[4],"E")==0)||
	    (strcmp1(&lakon[8*i+6],"I")==0))&&
	   (strcmp2(&lakon[8*i+6],"LC",2)!=0)){
	  blkassign[l++]=1;
	}else if(strcmp2(&lakon[8*i+6],"LC",2)==0){
	  /* composite material */
	  nlayer=0;
	  for(k=0;k<mi[2];k++){
	    if(ielmat[i*mi[2]+k]==0) break;
	    nlayer++;
	  }
	  for(k=0;k<nlayer;k++){
	    nemax++;
	    elem_map[l] = i+1;
	    blkassign[l++]=2;
	  }
	}else if(strcmp1(&lakon[8*i+6],"B")==0){
	  /* 3-node beam element */
	  blkassign[l++]=3;
	}else{
	  /* 8-node 2d element */
	  blkassign[l++]=4;
	}
      }else if(strcmp1(&lakon[8*i+3],"8")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  /* 8-node brick element */
	  blkassign[l++]=5;
	}else if(strcmp1(&lakon[8*i+6],"B")==0){
	  /* 2-node 1d element */
	  blkassign[l++]=6;
	  // if(strcmp1(&lakon[8*i+4],"R")==0){
	  //   blkassign[l++]=6;
	  // }else if(strcmp1(&lakon[8*i+4],"I")==0){
	  //   blkassign[l++]=6;
	  // }
	}else{
	  /* 4-node 2d element */
	  blkassign[l++]=7;
	  /* not sure exactly what this does, probably axisymmetric? */
	  // if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  // 
	  // if(strcmp1(&lakon[8*i+4],"R")==0){
	  //   blkassign[l++]=7;
	  // }else if(strcmp1(&lakon[8*i+4],"I")==0){
	  //   blkassign[l++]=7;
	  // }	    
	}
      }else if((strcmp1(&lakon[8*i+3],"10")==0)||
	       (strcmp1(&lakon[8*i+3],"14")==0)){
	/* 10-node tetrahedral element */
	blkassign[l++]=8;
      }else if(strcmp1(&lakon[8*i+3],"4")==0){
	/* 4-node tetrahedral element */
	blkassign[l++]=9;
      }else if(strcmp1(&lakon[8*i+3],"15")==0){
	if(((strcmp1(&lakon[8*i+6]," ")==0)||
	    (strcmp1(&filab[4],"E")==0))&&
           (strcmp2(&lakon[8*i+6],"LC",2)!=0)){
	  /* 15-node wedge element */
	  blkassign[l++]=10;
	}else{
	  /* 6-node 2d element */
	  blkassign[l++]=11;
	  /* not sure exactly what this does */
	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	}
      }else if(strcmp1(&lakon[8*i+3],"6")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  /* 6-node wedge element */
	  blkassign[l++]=12;
	}else{
	  /* 3-node 2d element */ /* Shells and triangles */
	  /* not sure exactly what this does */
	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  blkassign[l++]=13;
	}
	//      }else if((strcmp1(&lakon[8*i],"D")==0)&&
	//	       (strcmp1(&lakon[8*i],"DCOUP3D")!=0)){
      }else if(strcmp1(&lakon[8*i],"D")==0){
	if(kon[indexe]==0){
	  /* 2-node 1d element (network entry element) */
	  blkassign[l++]=14;
	}else if(kon[indexe+2]==0){
	  /* 2-node 1d element (network exit element) */
	  blkassign[l++]=15;
	}else{
	  /* 3-node 1d element (genuine network element) */
	  blkassign[l++]=16;
	}
      }else if((strcmp1(&lakon[8*i],"E")==0)&&
	       (strcmp1(&lakon[8*i+6],"A")==0)){
	/* Not sure exactly what iaxial does yet */
	if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	/* 2-node 1d element (spring element) */
	blkassign[l++]=17;
      }
    }


    int num_nodes_per_elem[num_elem_blk], num_edges_per_elem[num_elem_blk], num_faces_per_elem[num_elem_blk];
    char *blknames[num_elem_blk];
    j=0;
    num_nodes_per_elem[j]=1;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="PNT";
    num_nodes_per_elem[j]=20; num_edges_per_elem[j]=12; num_faces_per_elem[j]=6;  blknames[j++]="C3D20[R]";
    num_nodes_per_elem[j]=20; num_edges_per_elem[j]=12; num_faces_per_elem[j]=6;  blknames[j++]="COMPOSITE LAYER C3D20";
    num_nodes_per_elem[j]=3;  num_edges_per_elem[j]=4;  num_faces_per_elem[j]=4;  blknames[j++]="Beam B32[R]";
    num_nodes_per_elem[j]=8;  num_edges_per_elem[j]=4;  num_faces_per_elem[j]=4;  blknames[j++]="CPS8 CPE8 CAX8 S8[R]";
    num_nodes_per_elem[j]=8;  num_edges_per_elem[j]=12; num_faces_per_elem[j]=6;  blknames[j++]="C3D8[R]";
    num_nodes_per_elem[j]=2;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="TRUSS2";
    num_nodes_per_elem[j]=4;  num_edges_per_elem[j]=4;  num_faces_per_elem[j]=4;  blknames[j++]="CPS4[RI] CPE4[RI] CAX4 S4[R]";
    num_nodes_per_elem[j]=10; num_edges_per_elem[j]=6;  num_faces_per_elem[j]=4;  blknames[j++]="C3D10";
    num_nodes_per_elem[j]=4;  num_edges_per_elem[j]=6;  num_faces_per_elem[j]=4;  blknames[j++]="C3D4";
    num_nodes_per_elem[j]=15; num_edges_per_elem[j]=9;  num_faces_per_elem[j]=5;  blknames[j++]="C3D15";
    num_nodes_per_elem[j]=6;  num_edges_per_elem[j]=3;  num_faces_per_elem[j]=3;  blknames[j++]="CPS6 CPE6 S6";
    num_nodes_per_elem[j]=6;  num_edges_per_elem[j]=3;  num_faces_per_elem[j]=3;  blknames[j++]="C3D6";
    num_nodes_per_elem[j]=3;  num_edges_per_elem[j]=3;  num_faces_per_elem[j]=3;  blknames[j++]="CPS3 CPE3 S3";
    num_nodes_per_elem[j]=2;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="2-node 1d network entry elem";
    num_nodes_per_elem[j]=2;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="2-node 1d network exit elem";
    num_nodes_per_elem[j]=3;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="2-node 1d genuine network elem";
    num_nodes_per_elem[j]=2;  num_edges_per_elem[j]=0;  num_faces_per_elem[j]=0;  blknames[j++]="2-node 1d spring elem";

    errr = ex_put_names (exoid, EX_ELEM_BLOCK, blknames);
    if(errr){printf("*ERROR in exo: cannot write block names");}


    /* write element connectivity */
    ITG *connect;
    ITG num_elem_in_blk;
    ITG blksize[num_elem_blk];

    for(l=0;l<num_elem_blk;l++){
      // First determine the size of the block
      j=0;
      m=0;
      for(i=0;i<*ne0;i++){
	if(ipkon[i]<0) continue;
	if(blkassign[m]==l){
	  if(blkassign[m]==2){//composite
	    for(k=0;k<mi[2];k++){
	      if(ielmat[i*mi[2]+k]==0) break;
	      j++;
	    }
	  }else{
	    j++;
	  }
	}
	m++;
      }

      blksize[l]=j;
      num_elem_in_blk=blksize[l];

      connect = (ITG *) calloc (num_elem_in_blk*num_nodes_per_elem[l], sizeof(ITG));
      // printf ("Size of connect %" ITGFORMAT "\n", num_elem_in_blk*num_nodes_per_elem[l]*sizeof(ITG));
      k=0; o=0;

      // Now connectivity
      for(i=0;i<*ne0;i++){
	if(ipkon[i]<0) continue;
	indexe=ipkon[i];
	if (blkassign[o]==l){
	  // printf ("block assignment %" ITGFORMAT "\n", blkassign[o]);
	  if(blkassign[o]==1){ // C3D20
	    for(m=0;m<12;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	    for(m=16;m<20;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	    for(m=12;m<16;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	  }else if(blkassign[o]==10){ // C3D15
	    for(m=0;m<9;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	    for(m=12;m<15;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	    for(m=9;m<12;m++){connect[k++] = node_map_inv[kon[indexe+m]-1];}
	  }else if (blkassign[o]==2){ // Composite
	    nlayer=0;
	    for(l=0;l<mi[2];l++){
	      if(ielmat[i*mi[2]+l]==0) break;
	      nlayer++;
	    }
	    for(n=0;n<nlayer;n++){
	      for(m=0;m<12;m++){connect[k++] = node_map_inv[kon[indexe+28+20*n+m]-1];}
	      for(m=16;m<20;m++){connect[k++] = node_map_inv[kon[indexe+28+20*n+m]-1];}
	      for(m=12;m<16;m++){connect[k++] = node_map_inv[kon[indexe+28+20*n+m]-1];}
	    }
	  }else if(blkassign[o]==4){ // 8 Node 2D elements CAX8 S8 S8R etc
	    for (j = 0; j <num_nodes_per_elem[l]; j++){
	      connect[k++] = node_map_inv[kon[indexe+20+j]-1];
	    }
	  }else if(blkassign[o]== 5 || // C3D8 or C3D8R
		   blkassign[o]== 8 || // C3D10
		   blkassign[o]== 9 || // C3D4
		   blkassign[o]==12){  // C3D6
	    for (j = 0; j <num_nodes_per_elem[l]; j++){
	      connect[k++] = node_map_inv[kon[indexe+j]-1];
	    }
	  }else if(blkassign[o]==11){ // 6-node 2D element (S6)
	    for (j = 0; j <num_nodes_per_elem[l]; j++){
	      connect[k++] = node_map_inv[kon[indexe+15+j]-1];
	    }
	  }else if(blkassign[o]==13){ // 2D triangle or 3-node shell
	    for (j = 0; j <num_nodes_per_elem[l]; j++){
	      connect[k++] = node_map_inv[kon[indexe+6+j]-1];
	    }
	  }else { // 4 node shell element?
	    for (j = 0; j <num_nodes_per_elem[l]; j++){
	      // RETAIN FOR A TIME TO SEE IF ANYTHING BREAKS.
	      // Introduced new element classifications above in 2.8
	      connect[k++] = node_map_inv[kon[indexe+8+j]-1];
	    }
	  }
	}
	o++;
      }

      int num_attr=0;
      switch (l)
	{
	case 0:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "SPHERE", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 1:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 2:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 3:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 4:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "QUAD", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 5:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 6:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 7:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "SHELL", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 8:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TETRA", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 9:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TETRA", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 10:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 11:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 12:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	case 13:
	  // KEEP FOR AWHILE CHECKING BREAKAGE
	  // errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TRIANGLE", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	default:
	  // case 14:
	  // case 15:
	  // case 16:
	  // case 17:
	  errr = ex_put_block (exoid, EX_ELEM_BLOCK, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_edges_per_elem[l], num_faces_per_elem[l], num_attr);
	  break;
	};



      if (num_elem_in_blk>0){
	errr = ex_put_conn (exoid, EX_ELEM_BLOCK, l, connect, NULL, NULL);
	if (errr)
	  printf ("ERROR in ex_put_conn %i\n", errr);
      }
      free (connect);
    }

    // Write the element map into the file
    errr = ex_put_id_map (exoid, EX_ELEM_MAP, elem_map);
    if (errr)
      printf ("ERROR in ex_put_id_map %i\n", errr);

    // Write the node sets into the file
    exosetfind(set, nset, ialset, istartset, iendset,
	       &num_ns, &num_ss, &num_es, &num_fs, node_map_inv, exoid, (int) 1, nk);

    // Free up memory which is gathering dust
    free (elem_map);
    free (blkassign);

    // Close files
    ex_update (exoid);
    ex_close (exoid);

    if(*nmethod==0){return;}

    /* End of if(*kode==1) */
  }




  /* Start the data storage section */
  float version;
  exoid = ex_open (fneig, /*Filename*/
		   EX_WRITE,	/* create mode */
		   &CPU_word_size,  /* CPU float word size in bytes */
		   &IO_word_size, /* I/O float word size in bytes */
		   &version);

  // Define an output variable for general use
  float *nodal_var_vals;

  int num_time_steps;
  float fdum;
  char cdum = 0;

  float timet;
  errr = ex_inquire (exoid, EX_INQ_TIME, &num_time_steps, &fdum, &cdum);
  errr = ex_get_time (exoid, num_time_steps, &timet);
  if (num_time_steps>0){
    printf ("\t%i Time periods in exo file, most recent at time=%f\n", num_time_steps, timet);
  } else {
    printf ("\t0 Time periods in exo file\n");
  }
  timet = (float) *time;
  ++num_time_steps;
  printf ("\tWriting new time period %" ITGFORMAT " at time=%f\n", num_time_steps, *time);
  errr = ex_put_time (exoid, num_time_steps, &timet);
  if (errr) printf ("Error storing time into exo file.\n");

  // Statically allocate an array of pointers.  Can't figure out how to do this dynamically.
  // 100 should be enough to store the variable names.
  char *var_names[100];

  int countvars=0;
  int countbool=3;

  while(countbool>0){
    // First time count
    // Second time assign names
    // Last time save data

    /*  for cyclic symmetry frequency calculations only results for
	even numbers (= odd modes, numbering starts at 0) are stored */
    if((*nmethod==2)&&(((*mode/2)*2!=*mode)&&(*noddiam>=0))){ex_close(exoid);return;}

    /* storing the displacements in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(filab,"U ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="Ux";
	  var_names[countvars++]="Uy";
	  var_names[countvars++]="Uz";
	}else{
	  iselect=1;

	  frdset(filab,set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);

	  exovector(v,&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		    exoid, num_time_steps,countvars,nout,node_map_inv);
	  countvars+=3;
	}
      }
    }
    
    /*     storing the imaginary part of displacements in the nodes
	   for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if(strcmp1(filab,"U ")==0){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="U-imag-x";
	  var_names[countvars++]="U-imag-y";
	  var_names[countvars++]="U-imag-z";
	}else{
	  exovector(&v[*nk*mt],&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		    exoid,num_time_steps,countvars,nout,node_map_inv);

	  countvars+=3;
	}
      }
    }
    
    /*     storing the imaginary part of displacements in the nodes
	   for steady state calculations */
    if((*nmethod==5)&&(*mode==0)){
      if((strcmp1(filab,"U ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="U-imag-x";
	  var_names[countvars++]="U-imag-y";
	  var_names[countvars++]="U-imag-z";
	}else{
	  iselect=1;
	  exovector(v,&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		    exoid,num_time_steps,countvars,nout,node_map_inv);
	  countvars+=3;
	}
      }
    }

    /* storing the velocities in the nodes */
    if((strcmp1(&filab[1740],"V   ")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=3;
      }else if(countbool==2){
      	var_names[countvars++]="Vx";
	var_names[countvars++]="Vy";
	var_names[countvars++]="Vz";
      }else{
	iselect=1;

	frdset(&filab[1740],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exovector(veold,&iset,ntrans,&filab[1740],&nkcoords,inum,m1,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		  exoid,num_time_steps,countvars,nout,node_map_inv);

	countvars+=3;
      }
    }

    /* storing the temperatures in the nodes */
    if(strcmp1(&filab[87],"NT  ")==0){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
      	var_names[countvars++]="NT";
      }else{
	iselect=0;

	frdset(&filab[87],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	if(*ithermal<=1){
	  exoselect(t1,t1,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		    nfieldscalar,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	}else{
	  exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		    nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	}
	countvars+=1;
      }
    }

    /* storing the electrical potential in the nodes */
    if((strcmp1(&filab[3654],"POT ")==0)&&(*ithermal==2)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
      	var_names[countvars++]="POT";
      }else{
	iselect=0;
    
	frdset(&filab[3654],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export POT to exo not tested.\n");
	countvars+=1;
      }
    }
    
    /* storing the stresses in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="Sxx";
	  var_names[countvars++]="Syy";
	  var_names[countvars++]="Szz";
	  var_names[countvars++]="Sxy";
	  var_names[countvars++]="Sxz";
	  var_names[countvars++]="Syz";
	}else{
	  iselect=1;

	  frdset(&filab[174],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }
      
    /* storing the imaginary part of the stresses in the nodes
       for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="S-imagxx";
	  var_names[countvars++]="S-imagyy";
	  var_names[countvars++]="S-imagzz";
	  var_names[countvars++]="S-imagxy";
	  var_names[countvars++]="S-imagzx";
	  var_names[countvars++]="S-imagyz";
	}else{
	  exoselect(&stn[6**nk],stn,&iset,&nkcoords,inum,istartset,iendset,
	  	    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
	  	    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }

    /* storing the imaginary part of the stresses in the nodes
       for steady state calculations */
    if((*nmethod==5)&&(*mode==0)){
      if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="S-imagxx";
	  var_names[countvars++]="S-imagyy";
	  var_names[countvars++]="S-imagzz";
	  var_names[countvars++]="S-imagxy";
	  var_names[countvars++]="S-imagzx";
	  var_names[countvars++]="S-imagyz";
	}else{
	  iselect=1;

	  frdset(&filab[174],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
      
	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=6;
	}      
      }
    }

    /* storing the electromagnetic field E in the nodes */
    if((strcmp1(&filab[3741],"EMFE")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=3;
      }else if(countbool==2){
      	var_names[countvars++]="EMF-Ex";
	var_names[countvars++]="EMF-Ey";
	var_names[countvars++]="EMF-Ez";
      }else{
	iselect=1;

	frdset(&filab[3741],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvector,ifieldvector,icompvector,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export EMF-E to exo not tested.\n");
	countvars+=3;
      }
    }

    /* storing the electromagnetic field B in the nodes */
    if((strcmp1(&filab[3828],"EMFB")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=3;
      }else if(countbool==2){
      	var_names[countvars++]="EMF-Bx";
	var_names[countvars++]="EMF-By";
	var_names[countvars++]="EMF-Bz";
      }else{
	iselect=1;

	frdset(&filab[3828],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvector,ifieldvector,icompvectorlast,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export EMF-B to exo not tested.\n");
	countvars+=3;
      }
    }
    
    /* storing the total strains in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[261],"E   ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="Exx";
	  var_names[countvars++]="Eyy";
	  var_names[countvars++]="Ezz";
	  var_names[countvars++]="Exy";
	  var_names[countvars++]="Exz";
	  var_names[countvars++]="Eyz";
	}else{
	  iselect=1;


	  frdset(&filab[261],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);

	  exoselect(een,een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=6;
	}
      }
    }
    
    /* storing the imaginary part of the total strains in the nodes
       for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if(strcmp1(&filab[261],"E   ")==0){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="E-imagxx";
	  var_names[countvars++]="E-imagyy";
	  var_names[countvars++]="E-imagzz";
	  var_names[countvars++]="E-imagxy";
	  var_names[countvars++]="E-imagxz";
	  var_names[countvars++]="E-imagyz";
	}else{
	  exoselect(&een[6**nk],een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=6;
	}
      }
    }
    
    /* storing the imaginary part of the total strains in the nodes
       for steady state calculations */
    if((*nmethod==5)&&(*mode==0)){
      if((strcmp1(&filab[261],"E   ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="E-imagxx";
	  var_names[countvars++]="E-imagyy";
	  var_names[countvars++]="E-imagzz";
	  var_names[countvars++]="E-imagxy";
	  var_names[countvars++]="E-imagxz";
	  var_names[countvars++]="E-imagyz";
	}else{
	  iselect=1;
	  
	  frdset(&filab[261],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  
	  exoselect(een,een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  printf ("Warning: export E-imag to exo not tested.\n");
	  countvars+=6;
	}
      }
    }
    
    /* storing the mechanical strains in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[2697],"ME  ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="MExx";
	  var_names[countvars++]="MEyy";
	  var_names[countvars++]="MEzz";
	  var_names[countvars++]="MExy";
	  var_names[countvars++]="MExz";
	  var_names[countvars++]="MEyz";
	}else{
	  iselect=1;
	  frdset(&filab[2697],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  exoselect(emn,emn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }

    /* storing the imaginary part of the mechanical strains in the nodes
       for the odd modes of cyclic symmetry calculations */
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if((strcmp1(&filab[2697],"ME  ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="ME-imagxx";
	  var_names[countvars++]="ME-imagyy";
	  var_names[countvars++]="ME-imagzz";
	  var_names[countvars++]="ME-imagxy";
	  var_names[countvars++]="ME-imagxz";
	  var_names[countvars++]="ME-imagyz";
	}else{
	  exoselect(&emn[6**nk],een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }
  
    /* storing the forces in the nodes */
    if((*nmethod!=5)||(*mode==-1)){    
      if((strcmp1(&filab[348],"RF  ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="RFx";
	  var_names[countvars++]="RFy";
	  var_names[countvars++]="RFz";
	}else{
	  iselect=1;
	  frdset(&filab[348],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  if((iaxial==1)&&(strcmp1(&filab[352],"I")==0)){for(i=0;i<*nk;i++){fn[1+i*mt]*=180.;fn[2+i*mt]*=180.;fn[3+i*mt]*=180.;}}
	  exovector(fn,&iset,ntrans,&filab[348],&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		    exoid,num_time_steps,countvars,nout,node_map_inv);
	  if((iaxial==1)&&(strcmp1(&filab[352],"I")==0)){for(i=0;i<*nk;i++){fn[1+i*mt]/=180.;fn[2+i*mt]/=180.;fn[3+i*mt]/=180.;}}
	  countvars+=3;
	}
      }
    }
    
    /*     storing the imaginary part of the forces in the nodes
	   for the odd modes of cyclic symmetry calculations */
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if((strcmp1(&filab[348],"RF  ")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="RF-imagx";
	  var_names[countvars++]="RF-imagy";
	  var_names[countvars++]="RF-imagz";
	}else{
	  exovector(&fn[*nk*mt],&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3,
		    exoid,num_time_steps,countvars,nout,node_map_inv);
	  countvars+=3;
	}
      }
    }

    /* storing the equivalent plastic strains in the nodes */
    if((strcmp1(&filab[435],"PEEQ")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="PEEQ";
      }else{
	iselect=1;
	frdset(&filab[435],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
	exoselect(epn,epn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		  nfieldscalar,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export PEEQ to exo not tested.\n");
	countvars+=1;
      }
    }


    /* storing the energy in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[522],"ENER")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=1;
	}else if(countbool==2){
	  var_names[countvars++]="ENER";
	}else{
	  iselect=1;
	  frdset(&filab[522],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  exoselect(enern,enern,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		    nfieldscalar,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  printf ("Warning: export ENER to exo not tested.\n");
	  countvars+=1;
	}
      }
    }

    /* storing the contact displacements and stresses at the slave nodes */
    /* node-to-face penalty */
    if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar!=1)&&(*ithermal!=2)&&(*nmethod!=2)){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="COPEN";
	var_names[countvars++]="CSLIP1";
	var_names[countvars++]="CSLIP2";
	var_names[countvars++]="CPRESS";
	var_names[countvars++]="CSHEAR1";
	var_names[countvars++]="CSHEAR2";
      }else{

	nodal_var_vals = (float *) calloc (nout, sizeof(float));

	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	}
	noutloc=*ne-i-1;

	for(j=0;j<6;j++){
	  for(i=*ne-1;i>=0;i--){
	    if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	      break;
	    strcpy1(text,&lakon[8*i+7],1);
	    nope=atoi(text)+1;
	    nodes=node_map_inv[kon[ipkon[i]+nope-1]-1]-1;
	    nodal_var_vals[nodes]=stx[6*mi[0]*i+j];
	  }

	  errr = ex_put_var (exoid, num_time_steps, EX_NODAL, 1+countvars++, 1, nout, nodal_var_vals);
	  if (errr) printf ("ERROR storing contact data into exo file.\n");
	}

	free(nodal_var_vals);
	printf ("Warning: export CONT to exo not tested.\n");
      }
    }

    /* face-to-face penalty */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==1)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  var_names[countvars++]="COPEN";
	  var_names[countvars++]="CSLIP1";
	  var_names[countvars++]="CSLIP2";
	  var_names[countvars++]="CPRESS";
	  var_names[countvars++]="CSHEAR1";
	  var_names[countvars++]="CSHEAR2";
	}else{
	  iselect=1;

	  frdset(&filab[2175],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
    
	  exoselect(cdn,cdn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  printf ("Warning: export CONT to exo not tested.\n");
	  countvars+=6;
	}
      }
    }

    /* storing imaginary part of the differential contact displacements 
       and the contact stresses for the odd modes of cyclic symmetry
       calculations */
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==1)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  var_names[countvars++]="COPEN-Imag";
	  var_names[countvars++]="CSLIP1-Imag";
	  var_names[countvars++]="CSLIP2-Imag";
	  var_names[countvars++]="CPRESS-Imag";
	  var_names[countvars++]="CSHEAR1-Imag";
	  var_names[countvars++]="CSHEAR2-Imag";
	}else{
	  iselect=1;
	  exoselect(&cdn[6**nk],cdn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }
  
    /* storing the contact energy at the slave nodes */
    if((strcmp1(&filab[2262],"CELS")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="CELS";
      }else{
	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	}
	noutloc=*ne-i-1;

	nodal_var_vals = (float *) calloc (nkcoords, sizeof(float));
	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	  nope=atoi(&lakon[8*i+7])+1;
	  nodes=node_map_inv[kon[ipkon[i]+nope-1]-1]-1;
	  nodal_var_vals[nodes]=ener[i*mi[0]];
	}

	countvars+=1;
	int errr = ex_put_var (exoid, num_time_steps, EX_NODAL, countvars, 1, nout, nodal_var_vals);
	if (errr) printf ("ERROR storing CELS data into exo file.\n");
	printf ("Warning: export CELS to exo not tested.\n");
	free(nodal_var_vals);
      }
    }

    /* storing the internal state variables in the nodes */
    if(strcmp1(&filab[609],"SDV ")==0){
      if (countbool==3){
	countvars+=*nstate_;
      }else if(countbool==2){
	for(j=1;j<=*nstate_;j++){
	  sprintf(var_names[countvars++],"SDV%" ITGFORMAT,j);
	}
      }else{
	iselect=1;

	frdset(&filab[609],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	for(i=0;i<*nstate_;i++){
	  ifieldstate[i]=1;icompstate[i]=i;
	}
	nfield[0]=*nstate_;

	exoselect(xstaten,xstaten,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,nstate_,ifieldstate,icompstate,
		  nfield,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=*nstate_;
      }
    }

    /* storing the heat flux in the nodes
       the heat flux has been extrapolated from the integration points
       in subroutine extropolate.f, taking into account whether the
       results are requested in the global system or in a local system.
       Therefore, subroutine exovector cannot be used, since it assumes
       the values are stored in the global system */
    if((strcmp1(&filab[696],"HFL ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=3;
      }else if(countbool==2){
	var_names[countvars++]="HFLx";
	var_names[countvars++]="HFLy";
	var_names[countvars++]="HFLz";
      }else{
	iselect=1;

	frdset(&filab[696],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(qfn,qfn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvector,ifieldvector,icompvector,
		  nfieldvector1,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=3;
      }
    }

    /* storing the electrical current in the nodes
       (cf. heat flux HFL above)  */
    
    if((strcmp1(&filab[3567],"ECD ")==0)&&(*ithermal==2)){
      if (countbool==3){
	countvars+=3;
      }else if(countbool==2){
	var_names[countvars++]="ECDx";
	var_names[countvars++]="ECDy";
	var_names[countvars++]="ECDz";
      }else{
	iselect=1;
    
	frdset(&filab[3567],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(qfn,qfn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvector,ifieldvector,icompvector,
		  nfieldvector1,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export ECD to exo not tested.\n");
	countvars+=3;
      }
    }
    
    /* storing the heat generation in the nodes */
    if((strcmp1(&filab[783],"RFL ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="RFL";
      }else{
	iselect=1;

	frdset(&filab[783],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
	exoselect(fn,fn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=1;
      }
    }

    /* storing the Zienkiewicz-Zhu improved stresses in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[1044],"ZZS")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="ZZSxx";
	  var_names[countvars++]="ZZSyy";
	  var_names[countvars++]="ZZSzz";
	  var_names[countvars++]="ZZSxy";
	  var_names[countvars++]="ZZSxz";
	  var_names[countvars++]="ZZSyz";
	}else{
	  FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
			   stx,&mi[0]));

	  iselect=1;

	  frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);

	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=6;
	}
      }
    }
    
    /* storing the imaginary part of the Zienkiewicz-Zhu
       improved stresses in the nodes
       for the odd modes of cyclic symmetry calculations */
    
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if((strcmp1(&filab[1044],"ZZS")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  // Note reordered relative to frd file
	  // Order must be xx  yy  zz  xy  xz  yz
	  var_names[countvars++]="ZZS-imagxx";
	  var_names[countvars++]="ZZS-imagyy";
	  var_names[countvars++]="ZZS-imagzz";
	  var_names[countvars++]="ZZS-imagxy";
	  var_names[countvars++]="ZZS-imagxz";
	  var_names[countvars++]="ZZS-imagyz";
	}else{

	  FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
			   &stx[6*mi[0]**ne],&mi[0]));

	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  printf ("Warning: export of ZZSTR imaginary to exo not tested.\n");
	  countvars+=6;
	}
      }
    }

    /* storing the error estimator in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=2;
	}else if(countbool==2){
	  var_names[countvars++]="PSTDERROR";
	  var_names[countvars++]="VMSTDERROR";
	}else{

	  nterms=6;
	  FORTRAN(errorestimator,(stx,stn,ipkon,kon,lakon,nk,ne,
				  mi,ielmat,&nterms,inum,co,vold,&filab[1048],
				  ielprop,prop));

	  iselect=1;

	  frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);

	  ncomp=2;
	  ifield[0]=1;ifield[1]=1;
	  icomp[0]=0;icomp[1]=1;

	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomp,ifield,icomp,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=2;
	}
      }
    }

    /* storing the imaginary part of the error estimator in the nodes
       for the odd modes of cyclic symmetry calculations */
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2)){
	if (countbool==3){
	  countvars+=2;
	}else if(countbool==2){
	  var_names[countvars++]="PSTD-imag-ERROR";
	  var_names[countvars++]="VMSTD-imag-ERROR";
	}else{
	  nterms=6;
	  FORTRAN(errorestimator,(&stx[6*mi[0]**ne],stn,ipkon,kon,lakon,nk,ne,
				  mi,ielmat,&nterms,inum,co,vold,&filab[1048],
				  ielprop,prop));

	  ncomp=2;
	  ifield[0]=1;ifield[1]=1;
	  icomp[0]=0;icomp[1]=1;

	  exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomp,ifield,icomp,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=2;
	}
      }
    }

    /* storing the thermal error estimator in the nodes */
    if((*nmethod!=5)||(*mode==-1)){
      if((strcmp1(&filab[2784],"HER")==0)&&(*ithermal>1)){
	if (countbool==3){
	  countvars+=1;
	}else if(countbool==2){
	  var_names[countvars++]="HFLSTD";
	}else{
	  nterms=3;
	  FORTRAN(errorestimator,(qfx,qfn,ipkon,kon,lakon,nk,ne,
				  mi,ielmat,&nterms,inum,co,vold,&filab[2788],
				  ielprop,prop));

	  iselect=1;
	  frdset(&filab[2784],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);

	  ncomp=1;
	  ifield[0]=1;
	  icomp[0]=1;

	  exoselect(qfn,qfn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomp,ifield,icomp,
		    nfieldvector1,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);
	  countvars+=1;
	}
      }
    }

    /* storing the imaginary part of the thermal error estimator in the nodes
       for the odd modes of cyclic symmetry calculations */
  
    if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
      if(strcmp1(&filab[2784],"HER")==0){
	if (countbool==3){
	  countvars+=1;
	}else if(countbool==2){
	  var_names[countvars++]="HFLSTD-IMG";
	}else{
	  nterms=3;
	  FORTRAN(errorestimator,(&qfx[3*mi[0]**ne],qfn,ipkon,kon,lakon,nk,ne,
				  mi,ielmat,&nterms,inum,co,vold,&filab[2788],
				  ielprop,prop));

	  ncomp=1;
	  ifield[0]=1;
	  icomp[0]=1;
	  
	  exoselect(qfn,qfn,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomp,ifield,icomp,
		    nfieldvector1,&iselect,exoid,num_time_steps,countvars,nout,
		    node_map_inv);

	  countvars+=1;
	}
      }
    }

    /* storing the total temperatures in the network nodes */
    if((strcmp1(&filab[1131],"TT  ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="TT";
      }else{

	iselect=-1;
	frdset(&filab[1131],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of TT to exo not tested.\n");
	countvars+=1;
      }
    }

    /* storing the total temperatures in the network nodes */
    if((strcmp1(&filab[1218],"MF  ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="MF";
      }else{

	iselect=-1;
	frdset(&filab[1218],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=1;
	if((iaxial==1)&&(strcmp1(&filab[1222],"I")==0)){for(i=0;i<*nk;i++)v[1+i*mt]*=180.;}
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	if((iaxial==1)&&(strcmp1(&filab[1222],"I")==0)){for(i=0;i<*nk;i++)v[1+i*mt]/=180.;}
	countvars+=1;
      }
    }

    /* storing the total pressure in the network nodes */
    if((strcmp1(&filab[1305],"PT  ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="PT";
      }else{

	iselect=-1;
	frdset(&filab[1305],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=2;
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=1;
      }
    }

    /* storing the static pressure in the liquid network nodes */
    if((strcmp1(&filab[1827],"PS  ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="STPRESS PS";
      }else{

	iselect=-1;
	frdset(&filab[1827],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=2;
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of PS to exo not tested.\n");
	countvars+=1;
      }
    }

    /* storing the liquid depth in the channel nodes */
    if(strcmp1(&filab[2349],"PS  ")==0){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="DEPTH";
      }else{

	iselect=-1;
	frdset(&filab[2349],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=2;
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of DEPTH to exo not tested.\n");
	countvars+=1;
      }
    }

    /* storing the critical depth in the channel nodes */
    if((strcmp1(&filab[2436],"HCRI")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="HCRIT";
      }else{

	iselect=-1;
	frdset(&filab[2436],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=3;
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of HCRIT to exo not tested.\n");
	countvars+=1;
      }
    }

    /* storing the static temperature in the network nodes */
    if((strcmp1(&filab[1392],"TS  ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=1;
      }else if(countbool==2){
	var_names[countvars++]="STTEMP";
      }else{

	iselect=-1;
	frdset(&filab[1392],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	icomp[0]=3;
	exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of STTEMP to exo not tested.\n");
	countvars+=1;
      }
    }

    /*  the remaining lines only apply to frequency calculations
	with cyclic symmetry, complex frequency and steady state calculations */

    if((*nmethod!=2)&&(*nmethod!=5)&&(*nmethod!=6)&&(*nmethod!=7)){goto WRITENAMES;}
    if((*nmethod==5)&&(*mode==-1)){goto WRITENAMES;}

    /* storing the displacements in the nodes (magnitude, phase) */
    if((strcmp1(&filab[870],"PU  ")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	// Note reordered relative to frd file
	// Order must be xx  yy  zz  xy  xz  yz
	var_names[countvars++]="PDISP-MAGx";
	var_names[countvars++]="PDISP-MAGy";
	var_names[countvars++]="PDISP-MAGz";
	var_names[countvars++]="PDISP-PHAx";
	var_names[countvars++]="PDISP-PHAz";
	var_names[countvars++]="PDISP-PHAy";
      }else{
	iselect=1;

	frdset(&filab[870],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(vr,vi,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
		  nfieldvectph,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=6;
      }
    }

    /* storing the temperatures in the nodes (magnitude, phase) */
    if((strcmp1(&filab[957],"PNT ")==0)&&(*ithermal>1)){
      if (countbool==3){
	countvars+=2;
      }else if(countbool==2){
	var_names[countvars++]="PNDTEMP MAG1";
	var_names[countvars++]="PNDTEMP PHA2";
      }else{
	iselect=1;

	frdset(&filab[957],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(vr,vi,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompscalph,ifieldscalph,icompscalph,
		  nfieldscalph,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of PNDTEMP to exo not tested.\n");
	countvars+=2;
      }
    }

    /* storing the stresses in the nodes (magnitude, phase) */
	  
    if((strcmp1(&filab[1479],"PHS ")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=12;
      }else if(countbool==2){
	// Note reordered relative to frd file
	// Order must be xx  yy  zz  xy  xz  yz
	var_names[countvars++]="PSTRESS MAGXX";
	var_names[countvars++]="PSTRESS MAGYY";
	var_names[countvars++]="PSTRESS MAGZZ";
	var_names[countvars++]="PSTRESS MAGXY";
	var_names[countvars++]="PSTRESS MAGXZ";
	var_names[countvars++]="PSTRESS MAGYZ";
	var_names[countvars++]="PSTRESS PHAXX";
	var_names[countvars++]="PSTRESS PHAYY";
	var_names[countvars++]="PSTRESS PHAZZ";
	var_names[countvars++]="PSTRESS PHAXY";
	var_names[countvars++]="PSTRESS PHAXZ";
	var_names[countvars++]="PSTRESS PHAYZ";
      }else{
	iselect=1;

	frdset(&filab[1479],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(stnr,stni,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensph,ifieldtensph,icomptensph,
		  nfieldtensph,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=12;
      }
    }

    /* storing the differential contact displacements and
       the contact stresses in the nodes (magnitude, phase)
       only for face-to-face penalty contact */
    if((strcmp1(&filab[3915],"PCON")==0)&&(*ithermal!=2)&&(*mortar==1)){
      if (countbool==3){
	countvars+=12;
      }else if(countbool==2){
	// NOT SURE IF THE ORDER IS CORRECT HERE... ORDERING GETS
	// CHANGED IN EXOSELECT
	var_names[countvars++]="MAGO  ";
	var_names[countvars++]="MAGSL1";
	var_names[countvars++]="MAGSL2";
	var_names[countvars++]="MAGP  ";
	var_names[countvars++]="MAGSH1";
	var_names[countvars++]="MAGSH2";
	var_names[countvars++]="PHAO  ";
	var_names[countvars++]="PHASL1";
	var_names[countvars++]="PHASL2";
	var_names[countvars++]="PHAP  ";
	var_names[countvars++]="PHASH1";
	var_names[countvars++]="PHASH2";
      }else{
	iselect=1;

	frdset(&filab[3915],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(cdnr,cdni,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensph,ifieldtensph,icomptensph,
		  nfieldtensph,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	countvars+=12;
	printf ("Warning: export PCON to exo not tested and likely has incorrect ordering.\n");
      }
    }

    /* storing the displacements in the nodes (magnitude, phase) */
    if(strcmp1(&filab[2610],"PRF ")==0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	// Note reordered relative to frd file
	// Order must be xx  yy  zz  xy  xz  yz
	var_names[countvars++]="PFORC MAG1";
	var_names[countvars++]="PFORC MAG2";
	var_names[countvars++]="PFORC MAG3";
	var_names[countvars++]="PFORC PHA1";
	var_names[countvars++]="PFORC PHA3";
	var_names[countvars++]="PFORC PHA2";
      }else{
	iselect=1;

	frdset(&filab[2610],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	exoselect(fnr,fni,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
		  nfieldvectph,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of PFORC to exo not tested.\n");
	countvars+=6;
      }
    }

    /* the remaining parts are for frequency calculations with cyclic symmetry only */
    if(*nmethod!=2){goto WRITENAMES;}

    /* storing the maximum displacements of the nodes in the base sector
       (components, magnitude) */
	  
    if((strcmp1(&filab[1566],"MAXU")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=4;
      }else if(countbool==2){
	var_names[countvars++]="MDISP DX";
	var_names[countvars++]="MDISP DY";
	var_names[countvars++]="MDISP DZ";
	var_names[countvars++]="MDISP ANG";
      }else{
	iselect=1;

	frdset(&filab[1566],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	ncomp=4;
	ifield[0]=1;icomp[0]=1;
	ifield[1]=1;icomp[1]=2;
	ifield[2]=1;icomp[2]=3;
	ifield[3]=1;icomp[3]=0;
	nfield[0]=4;nfield[1]=4;

	exoselect(vmax,vmax,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomp,ifield,icomp,
		  nfield,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of MDISP to exo not tested.\n");
	countvars+=4;
      }
    }

    /* storing the worst principal stress at the nodes in the basis
       sector (components, magnitude)
    
       the worst principal stress is the maximum of the absolute value
       of all principal stresses, times its original sign */
      
    if((strcmp1(&filab[1653],"MAXS")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=7;
      }else if(countbool==2){
	// Note reordered relative to frd file
	// Order must be xx  yy  zz  xy  xz  yz
	var_names[countvars++]="MSTRESS SXX";
	var_names[countvars++]="MSTRESS SYY";
	var_names[countvars++]="MSTRESS SZZ";
	var_names[countvars++]="MSTRESS SXY";
	var_names[countvars++]="MSTRESS SXZ";
	var_names[countvars++]="MSTRESS SYZ";
	var_names[countvars++]="MSTRESS MAG";
      }else{
	iselect=1;

	frdset(&filab[1653],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	ncomp=7;
	ifield[0]=1;icomp[0]=1;
	ifield[1]=1;icomp[1]=2;
	ifield[2]=1;icomp[2]=3;
	ifield[3]=1;icomp[3]=4;
	ifield[4]=1;icomp[4]=6;
	ifield[5]=1;icomp[5]=5;
	ifield[6]=1;icomp[6]=0;
	nfield[0]=7;nfield[1]=7;

	exoselect(stnmax,stnmax,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomp,ifield,icomp,
		  nfield,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of MSTRESS to exo not tested.\n");
	countvars+=7;
      }
    }

    /* storing the worst principal strain at the nodes
       in the basis sector (components, magnitude)

       the worst principal strain is the maximum of the
       absolute value of all principal strains, times
       its original sign */
    if((strcmp1(&filab[2523],"MAXE")==0)&&(*ithermal!=2)){
      if (countbool==3){
	countvars+=7;
      }else if(countbool==2){
	// Note reordered relative to frd file
	// Order must be xx  yy  zz  xy  xz  yz
	var_names[countvars++]="MAXE EXX";
	var_names[countvars++]="MAXE EYY";
	var_names[countvars++]="MAXE EZZ";
	var_names[countvars++]="MAXE EXY";
	var_names[countvars++]="MAXE EXZ";
	var_names[countvars++]="MAXE EYZ";
	var_names[countvars++]="MAXE MAG";
      }else{
	iselect=1;

	frdset(&filab[2523],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);

	ncomp=7;
	ifield[0]=1;icomp[0]=1;
	ifield[1]=1;icomp[1]=2;
	ifield[2]=1;icomp[2]=3;
	ifield[3]=1;icomp[3]=4;
	ifield[4]=1;icomp[4]=6;
	ifield[5]=1;icomp[5]=5;
	ifield[6]=1;icomp[6]=0;
	nfield[0]=7;nfield[1]=7;

	exoselect(eenmax,eenmax,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomp,ifield,icomp,
		  nfield,&iselect,exoid,num_time_steps,countvars,nout,
		  node_map_inv);
	printf ("Warning: export of MAXE to exo not tested.\n");
	countvars+=7;
      }
    }

  WRITENAMES:
    if (countbool==3){
      errr = ex_put_variable_param (exoid, EX_NODAL, countvars);
      ex_update (exoid);
      // var_names = (char *) calloc (countvars, sizeof (char));
      // var_names = (char *) calloc (countvars, MAX_STR_LENGTH);
      int ii=0; // 4 bit integer
      for (ii=0; ii<countvars; ii++)
	{
	  var_names[ii] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
    }else if(countbool==2){
      errr = ex_put_variable_names (exoid, EX_NODAL, countvars, var_names);
      if (errr) {
	printf ("Unable to update variable names.  Was output data requested?\n\n");
	printf ("  NOTE CalculiX Extras doesn't currently support\n");
	printf ("  changing output variables between steps.\n");
	printf ("  Requested output variables in this step are:\n");
	int ii=0; // 4 bit integer
	for (ii=0; ii<countvars; ii++){
	  printf("    %i %s\n", ii+1, var_names[ii]);
	}
      }

      //  for (i=0; i<countvars; i++)
      // 	{
      // 	  free(var_names[i]);
      // 	}
      ex_update (exoid);
    }

    countvars=0;
    // printf ("Decrement countbool to %i\n",--countbool);
    --countbool;
  };

  // free (var_names);
  free (node_map_inv);
  ex_update (exoid);
  ex_close(exoid);

#else
  printf ("Warning: requested exodus output support not included at compile time.\n");
#endif
  return;
}
