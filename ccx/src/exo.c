/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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
#include "exodusII.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define rc(r,c) (r+nkcoords*c)

void exo(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne0,
	 double *v,double *stn,int *inum,int *nmethod,int *kode,
	 char *filab,double *een,double *t1,double *fn,double *time,
	 double *epn,int *ielmat,char *matname,double *enern,
	 double *xstaten,int *nstate_,int *istep,int *iinc,
	 int *ithermal,double *qfn,int *mode,int *noddiam,
	 double *trab,int *inotr,int *ntrans,double *orab,
	 int *ielorien,int *norien,char *description,int *ipneigh,
	 int *neigh,int *mi,double *stx,double *vr,double *vi,
	 double *stnr,double *stni,double *vmax,double *stnmax,
	 int *ngraph,double *veold,double *ener,int *ne,double *cs,
	 char *set,int *nset,int *istartset,int *iendset,int *ialset,
	 double *eenmax,double *fnr,double *fni,double *emn,
	 double *thicke,char *jobnamec,char *output){

  /* stores the results in exo format

     iselect selects which nodes are to be stored:
     iselect=-1 means only those nodes for which inum negative
     ist, i.e. network nodes
     iselect=+1 means only those nodes for which inum positive
     ist, i.e. structural nodes
     iselect=0  means both of the above */
  
  /* Note for frd strcmp1(output,"asc")==0 defines ascii file, otherwise binary */
  
  char fneig[132]="",date[8],clock[10],newdate[21],newclock[9],
    material[6]="     ",text[2]=" ";

  static int icounter=0,nkcoords,nout,noutmin,noutplus;

  int null,one,i,j,k,indexe,nemax,nlayer,noutloc,iset,iselect,ncomp,nope,
    nodes,ifield[7],nfield[2],icomp[7],ifieldstate[*nstate_],two,three,
    icompstate[*nstate_],imaterial=0,nelout;

  int ncompscalar=1,ifieldscalar[1]={1},icompscalar[1]={0},nfieldscalar[2]={1,0};
  int ncompvector=3,ifieldvector[3]={1,1,1},icompvector[3]={0,1,2},nfieldvector1[2]={3,0},nfieldvector0[2]={mi[1]+1,0};
  int ncomptensor=6,ifieldtensor[6]={1,1,1,1,1,1},icomptensor[6]={0,1,2,3,5,4},nfieldtensor[2]={6,0};
  // int ncompscalph=2,ifieldscalph[2]={1,2},icompscalph[2]={0,0},nfieldscalph[2]={0,0};
  // int ncompvectph=6,ifieldvectph[6]={1,1,1,2,2,2},icompvectph[6]={1,2,3,1,2,3},nfieldvectph[2]={mi[1]+1,mi[1]+1};
  // int ncomptensph=12,ifieldtensph[12]={1,1,1,1,1,1,2,2,2,2,2,2},icomptensph[12]={0,1,2,3,5,4,0,1,2,3,5,4},nfieldtensph[2]={6,6};
  int iw;
  float ifl;
  double pi,oner;

  int errr, exoid;
  int num_dim, num_elem;
  int num_elem_blk; /* Node element blocks.  One eltype per block*/
  int num_ns, num_ss; /* Node sets, side sets */
  int CPU_word_size = sizeof(float);
  int IO_word_size = sizeof(float);

  /* Filename */
  strcpy (fneig, jobnamec);
  strcat (fneig, ".exo");

  /* nkcoords is the number of nodes at the time when 
     the nodal coordinates are stored in the exo file */
  nkcoords = *nk;
  int num_nodes = nkcoords;

  
  pi=4.*atan(1.);
  null=0;
  one=1;
  two=2;
  three=3;
  oner=1.;

  /* first time something is written in the exo-file: store
     computational metadata, the nodal coordinates and the
     topology */
  if(*kode==1){
    
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
	nout++;
	if(inum[i]>0) noutplus++;
	if(inum[i]<0) noutmin++;
      }
    }else{
      nout=*nk;
    }

    num_dim = 3;  
    num_elem_blk = 19;
    num_ns = 0; 
    num_ss = 0;
    exoid = ex_create (fneig, /*Filename*/
		       EX_CLOBBER,	/* create mode */
		       &CPU_word_size,  /* CPU float word size in bytes */
		       &IO_word_size);  /* I/O float word size in bytes */
    
    
    float x[*nk], y[*nk], z[*nk];
    // Write optional node map
    j = 0;
    // int node_map[*nk];
    int *node_map;
    node_map = (int *) calloc(*nk, sizeof(int));
    
    /* storing the coordinates of the nodes */
    if(*nmethod!=0){
      for(i=0;i<*nk;i++){
	// if(inum[i]==0){
	//   // Put other nodes at the origin
	//   x[i] = 0.0;
	//   y[i] = 0.0;
	//   z[i] = 0.0;
	//   continue;}
	node_map[i] = i+1;
	x[i] = co[3*i];
	y[i] = co[3*i+1];
	z[i] = co[3*i+2];
      }
    }else{
      for(i=0;i<*nk;i++){
	node_map[i] = i+1;
	x[i] = co[3*i];
	y[i] = co[3*i+1];
	z[i] = co[3*i+2];
      }
    }

    
    /* determining the number of elements */
    if(*nmethod!=0){
      nelout=0;
      for(i=0;i<*ne0;i++){
	if(ipkon[i]<0) continue;
	if(strcmp1(&lakon[8*i],"ESPRNGC")==0) continue;
	nelout++;
      }
    }else{
      nelout=*ne;
    }
    num_elem = nelout;
    
    /* initialize file with parameters */
    printf("\nData writen to the .exo file\n");
    printf("Number of nodes %i\n", num_nodes);
    printf("Number of elements %i\n", num_elem);
    printf("Number of element blocks %i\n", num_elem_blk);
    printf("Number of node sets %i\n", num_ns);
    printf("Number of side sets %i\n", num_ss);
    errr = ex_put_init (exoid, "CalculiX EXO File", num_dim,
			num_nodes, num_elem, num_elem_blk, num_ns, num_ss);

    if(errr){
      printf("*ERROR in exo: cannot open exo file for writing...");
    }
    
    /* write values to database */
    errr = ex_put_coord (exoid, x, y, z);
    errr = ex_put_node_num_map (exoid, node_map);
    if(errr){
      printf("*ERROR in exo: failed to write node map");
    }
    free (node_map);

    //Write coordinate names
    char *coord_names[3];
    coord_names[0] = "x";
    coord_names[1] = "y";
    coord_names[2] = "z";
    errr = ex_put_coord_names (exoid, coord_names);
    if(errr){
      printf("*ERROR in exo: failed to write coordinate names");
    }
    
    //SAMPLE    // Property names
    //SAMPLE    char *prop_names[num_elem];
    //SAMPLE    int fakeprop[num_elem];
    //SAMPLE    int num_props = 2;
    //SAMPLE    prop_names[0] = "TOP";
    //SAMPLE    prop_names[1] = "BOT";
    //SAMPLE    for (i=0; i<num_elem; i++)
    //SAMPLE      fakeprop[i]=1;
    //SAMPLE    errr = ex_put_prop_names (exoid, EX_ELEM_BLOCK, num_props, prop_names);
    //SAMPLE    if (errr)
    //SAMPLE      printf ("ERROR in ex_put_prop_names %i\n", errr);
    //SAMPLE    errr = ex_put_prop_array (exoid, EX_ELEM_BLOCK, prop_names[0], fakeprop);
    //SAMPLE    if (errr)
    //SAMPLE      printf ("ERROR in ex_put_prop_names %i\n", errr);
    //SAMPLE    errr = ex_put_elem_block (exoid, blkid, "SHEL", num_elem_in_blk, num_nodes_per_elem, num_attr);
    //SAMPLE    if (errr){
    //SAMPLE      printf ("ERROR in ex_put_elem_block %i\n", errr);
    //SAMPLE    
    //SAMPLE    int *idelbs;
    //SAMPLE    idelbs = (int *) calloc(10, sizeof(int));
    //SAMPLE    errr = ex_get_elem_blk_ids (exoid, idelbs);
    //SAMPLE    
    
    /* storing the topology */
    nemax=*ne0;
    num_elem=nemax;


    int num_elem_in_blk;
    
    // Initialize enough memory to store the element numbers
    int *elem_map;
    elem_map = (int *) calloc(num_elem, sizeof(int));
    char curblk[6]="     ";

    int *blkassign;
    blkassign = (int *) calloc(num_elem, sizeof(int));
    
    int num_nodes_per_elem[num_elem_blk];
    for(i=0;i<*ne0;i++){
      elem_map[i] = i+1;
      if(ipkon[i]<0){
	blkassign[i]=0; 
	continue;}
      
      strcpy1(curblk,&lakon[8*i],5);
      strcpy1(material,&matname[80*(ielmat[i*mi[2]]-1)],5);
      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i+3],"2")==0){
	if(((strcmp1(&lakon[8*i+6]," ")==0)||
	    (strcmp1(&filab[4],"E")==0)||
	    (strcmp1(&lakon[8*i+6],"I")==0))&&
	   (strcmp2(&lakon[8*i+6],"LC",2)!=0)){
	  blkassign[i]=1;
	}else if(strcmp2(&lakon[8*i+6],"LC",2)==0){
	  /* composite material */
	  nlayer=0;
	  printf ("Composite materials not implemented in exo file.\n");
	  return;
	  for(k=0;k<mi[2];k++){
	    if(ielmat[i*mi[2]+k]==0) break;
	    nlayer++;
	  }
	  for(k=0;k<nlayer;k++){
	    nemax++;
	    // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,nemax,p4,p0,material,m2);
	    // for(j=0;j<10;j++)fprintf(f1,"%10d",kon[indexe+28+20*k+j]);
	    // fprintf(f1,"\n%3s",m2);
	    // for(j=10;j<12;j++)fprintf(f1,"%10d",kon[indexe+28+20*k+j]);
	    // for(j=16;j<19;j++)fprintf(f1,"%10d",kon[indexe+28+20*k+j]);
	    // for(j=19;j<20;j++)fprintf(f1,"%10d",kon[indexe+28+20*k+j]);
	    // for(j=12;j<16;j++)fprintf(f1,"%10d",kon[indexe+28+20*k+j]);
	    // fprintf(f1,"\n");
	  }
	}else if(strcmp1(&lakon[8*i+6],"B")==0){ // Beams?
	  blkassign[i]=2;
	}else{
	  blkassign[i]=3;
	}
      }else if(strcmp1(&lakon[8*i+3],"8")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  blkassign[i]=4;
	}else if(strcmp1(&lakon[8*i+6],"B")==0){
	  if(strcmp1(&lakon[8*i+4],"R")==0){
	    blkassign[i]=5;
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    blkassign[i]=6;
	  }
	}else{
	  if(strcmp1(&lakon[8*i+4],"R")==0){
	    blkassign[i]=7;
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    blkassign[i]=8;
	  }
	}
      }else if(strcmp1(&lakon[8*i+3],"10")==0){
	blkassign[i]=9;
      }else if(strcmp1(&lakon[8*i+3],"4")==0){
	blkassign[i]=10;
      }else if(strcmp1(&lakon[8*i+3],"15")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  blkassign[i]=11;
	}else{
	  blkassign[i]=12;
	}
      }else if(strcmp1(&lakon[8*i+3],"6")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  blkassign[i]=13;
	}else{
	  blkassign[i]=14;
	}
      }else if((strcmp1(&lakon[8*i],"D")==0)&&
	       (strcmp1(&lakon[8*i],"DCOUP3D")!=0)){
	if(kon[indexe]==0){
	  blkassign[i]=15;
	}else if(kon[indexe+2]==0){
	  blkassign[i]=16;
	}else{
	  blkassign[i]=17;
	}
      }else if((strcmp1(&lakon[8*i],"E")==0)&&
	       (strcmp1(&lakon[8*i+6],"A")==0)){
	blkassign[i]=18;
      }
    }


    /* write element connectivity */
    int *connect;
    int blksize[num_elem_blk];
    int l, k;

    j=0;
    num_nodes_per_elem[j++]=1;
    num_nodes_per_elem[j++]=20;
    num_nodes_per_elem[j++]=3;
    num_nodes_per_elem[j++]=8;
    num_nodes_per_elem[j++]=8;
    num_nodes_per_elem[j++]=2;
    num_nodes_per_elem[j++]=2;
    num_nodes_per_elem[j++]=4;
    num_nodes_per_elem[j++]=4;
    num_nodes_per_elem[j++]=10;
    num_nodes_per_elem[j++]=4;
    num_nodes_per_elem[j++]=15;
    num_nodes_per_elem[j++]=6;
    num_nodes_per_elem[j++]=6;
    num_nodes_per_elem[j++]=3;
    num_nodes_per_elem[j++]=2;
    num_nodes_per_elem[j++]=2;
    num_nodes_per_elem[j++]=3;
    num_nodes_per_elem[j++]=2;

    char *blknames[num_elem_blk];
    j=0;
    blknames[j++]="PNT";
    blknames[j++]="C3D20*";
    blknames[j++]="TRUSS3?";
    blknames[j++]="C3D8*?";
    blknames[j++]="C3D8*?";
    blknames[j++]="TRUSS2";
    blknames[j++]="TRUSS2";
    blknames[j++]="SHELL";
    blknames[j++]="SHELL";
    blknames[j++]="TETRA";
    blknames[j++]="TETRA";
    blknames[j++]="WEDGE";
    blknames[j++]="HEX";
    blknames[j++]="WEDGE";
    blknames[j++]="WEDGE";
    blknames[j++]="TRUSS";
    blknames[j++]="TRUSS";
    blknames[j++]="TRUSS";
    blknames[j++]="TRUSS";
    ex_put_names (exoid, EX_ELEM_BLOCK, blknames);

    for(l=0;l<num_elem_blk;l++){
      
      // First determine the size of the block
      j=0;
      for(i=0;i<*ne0;i++){
	if(blkassign[i]==l) j++;
      }
      blksize[l]=j;
      num_elem_in_blk=blksize[l];
      
      connect = (int *) calloc (num_elem_in_blk*num_nodes_per_elem[l], sizeof(int));
      k=0;
      // Now connectivity
      for(i=0;i<*ne0;i++){
	indexe=ipkon[i];
	if (blkassign[i]==l){
	  for (j = 0; j <num_nodes_per_elem[l]; j++){
	    connect[k++] = kon[indexe+j];
	  }
	}
      }

      int num_attr=0;
      switch (l)
	{
	case 0:
	  errr = ex_put_elem_block (exoid, l, "SPHERE", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 1:
	  errr = ex_put_elem_block (exoid, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 2:
	  errr = ex_put_elem_block (exoid, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 3:
	  errr = ex_put_elem_block (exoid, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 4:
	  errr = ex_put_elem_block (exoid, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 5:
	  errr = ex_put_elem_block (exoid, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 6:
	  errr = ex_put_elem_block (exoid, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 7:
	  errr = ex_put_elem_block (exoid, l, "SHELL", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 8:
	  errr = ex_put_elem_block (exoid, l, "SHELL", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 9:
	  errr = ex_put_elem_block (exoid, l, "TETRA", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 10:
	  errr = ex_put_elem_block (exoid, l, "TETRA", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 11:
	  errr = ex_put_elem_block (exoid, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 12:
	  errr = ex_put_elem_block (exoid, l, "HEX", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 13:
	  errr = ex_put_elem_block (exoid, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	case 14:
	  errr = ex_put_elem_block (exoid, l, "WEDGE", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	default:
	  // case 15:
	  // case 16:
	  // case 17:
	  // case 18:
	  errr = ex_put_elem_block (exoid, l, "TRUSS", num_elem_in_blk, num_nodes_per_elem[l], num_attr);	  
	};
	  


      if (num_elem_in_blk>0){
	errr = ex_put_elem_conn (exoid, l, connect);
	if (errr)
	  printf ("ERROR in ex_put_elem_conn %i\n", errr);
      }
      free (connect);
    }
    
    // Write the element map into the file
    errr = ex_put_elem_num_map (exoid, elem_map); 
    if (errr)
      printf ("ERROR in ex_put_elem_num_map %i\n", errr);
    free (elem_map);

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
    
  int num_time_steps;
  float fdum;
  char cdum = 0;
  
  errr = ex_inquire (exoid, EX_INQ_TIME, &num_time_steps, &fdum, &cdum);
  printf ("Time periods existing in exo file: %i\n", num_time_steps);
  ++num_time_steps;
  float timet = (float) *time;
  printf ("Writing time period %i at time %f\n", num_time_steps, *time);
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
    if(strcmp1(filab,"U ")==0){
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
	
	exovector(v,&iset,ntrans,filab,&nkcoords,inum,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,exoid,
		  num_time_steps,countvars);
	countvars+=3;
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
	  exovector(&v[*nk*(mi[1]+1)],&iset,ntrans,filab,&nkcoords,inum,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,exoid,
		    num_time_steps,countvars);
	  printf ("Warning: export of imaginary part of displacement to exo not tested.\n");
	  countvars+=3;
	}
      }
    }
    
    /* storing the velocities in the nodes */
    if(strcmp1(&filab[1740],"V   ")==0){
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
    
	exovector(veold,&iset,ntrans,&filab[1740],&nkcoords,inum,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,exoid,
		  num_time_steps,countvars);
	printf ("Warning velocity to exo not tested.\n");
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
		    nfieldscalar,&iselect,exoid,num_time_steps,countvars);
	  printf ("Warning: export temperature to exo not tested.\n");
	}else{
	  exoselect(v,v,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
		    nfieldscalar,&iselect,exoid,num_time_steps,countvars);
	  printf ("Warning: export temperature to exo not tested.\n");
	}
	countvars+=1;	
      }
    }
    
    /* storing the stresses in the nodes */
  
    if(strcmp1(&filab[174],"S   ")==0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="S_XX";
	var_names[countvars++]="S_YY";
	var_names[countvars++]="S_ZZ";
	var_names[countvars++]="S_XY";
	var_names[countvars++]="S_YZ";
	var_names[countvars++]="S_XZ";
      }else{
	iselect=1;
    
	frdset(&filab[174],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	countvars+=6;
      }
    }

    /* storing the imaginary part of the stresses in the nodes
       for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="S-imag-xx";
	var_names[countvars++]="S-imag-yy";
	var_names[countvars++]="S-imag-zz";
	var_names[countvars++]="S-imag-xy";
	var_names[countvars++]="S-imag-yz";
	var_names[countvars++]="S-imag-zx";
      }else{
	if(strcmp1(&filab[174],"S   ")==0){      
	  exoselect(&stn[6**nk],stn,&iset,&nkcoords,inum,istartset,iendset,
	  	    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
	  	    nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	  printf ("Warning: export of imaginary part of stress to exo not tested.\n");
	  countvars+=6;
	}
      }
    }

    /* storing the total strains in the nodes */
    if(strcmp1(&filab[261],"E   ")==0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="E_XX";
	var_names[countvars++]="E_YY";
	var_names[countvars++]="E_ZZ";
	var_names[countvars++]="E_XY";
	var_names[countvars++]="E_YZ";
	var_names[countvars++]="E_XZ";
      }else{
	iselect=1;
	
    
	frdset(&filab[261],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
	
	exoselect(een,een,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	countvars+=6;
      }
    }
    
    /* storing the imaginary part of the total strains in the nodes
       for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if(strcmp1(&filab[261],"E   ")==0){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  var_names[countvars++]="E-imag_XX";
	  var_names[countvars++]="E-imag_YY";
	  var_names[countvars++]="E-imag_ZZ";
	  var_names[countvars++]="E-imag_XY";
	  var_names[countvars++]="E-imag_YZ";
	  var_names[countvars++]="E-imag_XZ";
	}else{
	  exoselect(&een[6**nk],een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	  printf ("Warning: export of imaginary part of strain to exo not tested.\n");
	  countvars+=6;
	}
      }
    }
    
    /* storing the mechanical strains in the nodes */
    if(strcmp1(&filab[2697],"ME  ")==0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="ME_XX";
	var_names[countvars++]="ME_YY";
	var_names[countvars++]="ME_ZZ";
	var_names[countvars++]="ME_XY";
	var_names[countvars++]="ME_YZ";
	var_names[countvars++]="ME_XZ";
      }else{
	iselect=1;
      
	frdset(&filab[2697],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
	exoselect(emn,emn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	printf ("Warning: export of mechanical strain to exo not tested.\n");
	countvars+=6;
      }    
    }

    /* storing the imaginary part of the mechanical strains in the nodes
       for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if(strcmp1(&filab[2697],"ME  ")==0){
	if (countbool==3){
	  countvars+=6;
	}else if(countbool==2){
	  var_names[countvars++]="ME-imag_XX";
	  var_names[countvars++]="ME-imag_YY";
	  var_names[countvars++]="ME-imag_ZZ";
	  var_names[countvars++]="ME-imag_XY";
	  var_names[countvars++]="ME-imag_YZ";
	  var_names[countvars++]="ME-imag_XZ";
	}else{
	  exoselect(&emn[6**nk],een,&iset,&nkcoords,inum,istartset,iendset,
		    ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		    nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	  printf ("Warning: export of imaginary part of mechanical strain to exo not tested.\n");
	  countvars+=6;      
	}
      }
    }

    /* storing the forces in the nodes */
    
    if(strcmp1(&filab[348],"RF  ")==0){
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
      
	exovector(fn,&iset,ntrans,&filab[348],&nkcoords,inum,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,exoid,
		  num_time_steps,countvars);
	printf ("Warning force to exo not tested.\n");
	countvars+=3;
      }
    }
    
    /*     storing the imaginary part of the forces in the nodes
	   for the odd modes of cyclic symmetry calculations */
    if(*noddiam>=0){
      if(strcmp1(&filab[348],"RF  ")==0){
	if (countbool==3){
	  countvars+=3;
	}else if(countbool==2){
	  var_names[countvars++]="RF-imagx";
	  var_names[countvars++]="RF-imagy";
	  var_names[countvars++]="RF-imagz";
	}else{
	  exovector(&fn[*nk*(mi[1]+1)],&iset,ntrans,filab,&nkcoords,inum,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,exoid,
		    num_time_steps,countvars);
	  printf ("Warning imaginary force to exo not tested.\n");
	  countvars+=3;
	}
      }
    }

    /* storing the equivalent plastic strains in the nodes */
    if(strcmp1(&filab[435],"PEEQ")==0){
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
		  nfieldscalar,&iselect,exoid,num_time_steps,countvars);
	printf ("Warning: export PEEQ to exo not tested.\n");
	countvars+=1;
      }
    }
      

  /* storing the energy in the nodes */
  
    if(strcmp1(&filab[522],"ENER")==0){
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
		  nfieldscalar,&iselect,exoid,num_time_steps,countvars);
	printf ("Warning: export ENER to exo not tested.\n");
	countvars+=1;
      }
    }
  
    /* storing the contact displacements and stresses at the slave nodes */
    if(strcmp1(&filab[2175],"CONT")==0){
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
	
	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	}
	noutloc=*ne-i-1;
	
	int num_nod_vars = 6;
	float *nodal_var_vals;
	nodal_var_vals = (float *) calloc (nkcoords * num_nod_vars, sizeof(float));
	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	  strcpy1(text,&lakon[8*i+7],1);
	  nope=atoi(text);
	  nodes=kon[ipkon[i]+nope-1];
	  // fprintf(f1,"%3s%10d",m1,nodes);
	  // for(j=0;j<6;j++)fprintf(f1,"%12.5E",stx[6*mi[0]*i+j]);
	  for(j=0;j<6;j++){
	    nodal_var_vals[rc(i,j)]=stx[6*mi[0]*i+j]; 
	  }
	}

	float *nodal_var_vals_out;
	nodal_var_vals_out = (float *) calloc (nkcoords, sizeof(float));
	for (j=0; j<num_nod_vars; j++){
	  for (i=0; i<nkcoords; i++){
	    nodal_var_vals_out[i]=nodal_var_vals[rc(i,j)];
	  }
	  int errr = ex_put_nodal_var (exoid, num_time_steps, j+1+countvars, nkcoords, nodal_var_vals_out);
	  if (errr) printf ("ERROR storing data into exo file for dim %i record %i.\n", j, countvars+j);
	}

	free(nodal_var_vals);
	free(nodal_var_vals_out);
        
	printf ("Warning: export of Contact Variables to exo not tested and not yet expected to work.\n");
	countvars+=6;
      }
    }

    /* storing the contact energy at the slave nodes */
    if(strcmp1(&filab[2262],"CELS")==0){
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
	    
	float *nodal_var_vals_out;
	nodal_var_vals_out = (float *) calloc (nkcoords, sizeof(float));
	for(i=*ne-1;i>=0;i--){
	  if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	    break;
	  nope=atoi(&lakon[8*i+7]);
	  nodes=kon[ipkon[i]+nope-1];
	  nodal_var_vals_out[i]=ener[i*mi[0]];
	}
	    
	int errr = ex_put_nodal_var (exoid, num_time_steps, countvars, nkcoords, nodal_var_vals_out);
	if (errr) printf ("ERROR storing CELS data into exo file.\n");
	printf ("Warning: export of contact energy to exo not tested and not yet expected to work.\n");
	countvars+=1;
      }
    }
  
  /* storing the internal state variables in the nodes */
    if(strcmp1(&filab[609],"SDV ")==0){
      if (countbool==3){
	countvars+=*nstate_;
      }else if(countbool==2){
	for(j=1;j<=*nstate_;j++){
	  char str[6];
	  sprintf(str, "SDV%2d",j);
	  var_names[countvars++]=str;
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
		  nfield,&iselect,exoid,num_time_steps,countvars);
	printf ("Warning: export of SDV to exo not tested and not yet expected to work.\n");
	countvars+=*nstate_;
      }
    }
  
    /* storing the heat flux in the nodes
       the heat flux has been extrapolated from the integration points
       in subroutine extropolate.f, taking into account whether the 
       results are requested in the global system or in a local system.
       Therefore, subroutine frdvector cannot be used, since it assumes
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
		  nfieldvector1,&iselect,exoid,num_time_steps,countvars);
	countvars+=3;
	printf ("Warning: export of HFL to exo not tested.\n");
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
		  nfieldvector0,&iselect,exoid,num_time_steps,countvars);
	countvars+=1;
	printf ("Warning: export of RFL to exo not tested.\n");
      }
    }
  
    /* storing the Zienkiewicz-Zhu improved stresses in the nodes */
    if(strcmp1(&filab[1044],"ZZS")==0){
      if (countbool==3){
	countvars+=6;
      }else if(countbool==2){
	var_names[countvars++]="ZZSXX";
	var_names[countvars++]="ZZSYY";
	var_names[countvars++]="ZZSZZ";
	var_names[countvars++]="ZZSXY";
	var_names[countvars++]="ZZSYZ";
	var_names[countvars++]="ZZSXZ";
      }else{
	FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
			 stx,&mi[0]));

	iselect=1;
    
	frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
	       inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	       ngraph);
    
	exoselect(stn,stn,&iset,&nkcoords,inum,istartset,iendset,
		  ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		  nfieldtensor,&iselect,exoid,num_time_steps,countvars);
	printf ("Warning: export of ZZSTR to exo not tested.\n");
	countvars+=6;
      }
    }

    //  /* storing the imaginary part of the Zienkiewicz-Zhu 
    //     improved stresses in the nodes
    //     for the odd modes of cyclic symmetry calculations */
    //  
    //  if(*noddiam>=0){
    //    if(strcmp1(&filab[1044],"ZZS")==0){
    //
    //      FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
    //		      &stx[6*mi[0]**ne],&mi[0]));
    //      
    //      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //      
    //      fprintf(f1," -4  ZZSTR       6    1\n");
    //      fprintf(f1," -5  SXX         1    4    1    1\n");
    //      fprintf(f1," -5  SYY         1    4    2    2\n");
    //      fprintf(f1," -5  SZZ         1    4    3    3\n");
    //      fprintf(f1," -5  SXY         1    4    1    2\n");
    //      fprintf(f1," -5  SYZ         1    4    2    3\n");
    //      fprintf(f1," -5  SZX         1    4    3    1\n");
    //      
    //      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
    //                nfieldtensor,&iselect,m2,f1,output,m3);
    //      
    //    }
    //  }
    //  
    //  /* storing the error estimator in the nodes */
    //  
    //  if(strcmp1(&filab[1044],"ERR")==0){
    //
    //    FORTRAN(errorestimator,(stx,stn,ipkon,inum,kon,lakon,nk,ne,
    //            mi,ielmat,thicke));
    //
    //    iselect=1;
    //    
    //    frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  ERROR       2    1\n");
    //    fprintf(f1," -5  PSTD        1    1    1    0\n");
    //    fprintf(f1," -5  VMSTD       1    2    2    0\n");
    //
    //    ncomp=2;
    //    ifield[0]=1;ifield[1]=1;
    //    icomp[0]=2;icomp[1]=4;
    //
    //    frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomp,ifield,icomp,
    //                nfieldtensor,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the imaginary part of the error estimator in the nodes
    //     for the odd modes of cyclic symmetry calculations */
    //  
    //  if(*noddiam>=0){
    //    if(strcmp1(&filab[1044],"ERR")==0){
    //
    //      FORTRAN(errorestimator,(&stx[6*mi[0]**ne],stn,ipkon,inum,kon,lakon,nk,ne,
    //			      mi,ielmat,thicke));
    //      
    //      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //      fprintf(f1," -4  ERROR       2    1\n");
    //      fprintf(f1," -5  PSTD        1    1    1    0\n");
    //      fprintf(f1," -5  VMSTD       1    2    2    0\n");
    //      
    //      ncomp=2;
    //      ifield[0]=1;ifield[1]=1;
    //      icomp[0]=2;icomp[1]=4;
    //
    //      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomp,ifield,icomp,
    //                nfieldtensor,&iselect,m2,f1,output,m3);
    //      
    //    }
    //  }
    //
    //  /* storing the total temperatures in the network nodes */
    //  
    //  if(strcmp1(&filab[1131],"TT  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[1131],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  TOTEMP      1    1\n");
    //    fprintf(f1," -5  TT          1    1    0    0\n");
    //
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the total temperatures in the network nodes */
    //  
    //  if(strcmp1(&filab[1218],"MF  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[1218],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  MAFLOW      1    1\n");
    //    fprintf(f1," -5  MF          1    1    0    0\n");
    //
    //    icomp[0]=1;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the total pressure in the network nodes */
    //  
    //  if(strcmp1(&filab[1305],"PT  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[1305],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  TOPRES      1    1\n");
    //    fprintf(f1," -5  PT          1    1    0    0\n");
    //
    //    icomp[0]=2;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the static pressure in the liquid network nodes */
    //  
    //  if(strcmp1(&filab[1827],"PS  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[1827],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  STPRES      1    1\n");
    //    fprintf(f1," -5  PS          1    1    0    0\n");
    //
    //    icomp[0]=2;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the liquid depth in the channel nodes */
    //  
    //  if(strcmp1(&filab[2349],"PS  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[2349],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  DEPTH       1    1\n");
    //    fprintf(f1," -5  DEPTH       1    1    0    0\n");
    //
    //    icomp[0]=2;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the critical depth in the channel nodes */
    //  
    //  if(strcmp1(&filab[2436],"HCRI")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[2436],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  HCRIT       1    1\n");
    //    fprintf(f1," -5  HCRIT       1    1    0    0\n");
    //
    //    icomp[0]=3;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the static temperature in the network nodes */
    //  
    //  if(strcmp1(&filab[1392],"TS  ")==0){
    //
    //    iselect=-1;
    //    frdset(&filab[1392],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  STTEMP      1    1\n");
    //    fprintf(f1," -5  TS          1    1    0    0\n");
    //
    //    icomp[0]=3;
    //    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
    //                nfieldvector0,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    /*  the remaining lines only apply to frequency calculations
	with cyclic symmetry, complex frequency and steady state calculations */
  
    // if((*nmethod!=2)&&(*nmethod!=5)&&(*nmethod!=6)){ex_close(exoid);return;}
    // if((*nmethod==5)&&(*mode==-1)){ex_close(exoid);return;}
  
    //  /* storing the displacements in the nodes (magnitude, phase) */
    //	  
    //  if(strcmp1(&filab[870],"PU  ")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[870],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  PDISP       6    1\n");
    //    fprintf(f1," -5  MAG1        1   12    1    0\n");
    //    fprintf(f1," -5  MAG2        1   12    2    0\n");
    //    fprintf(f1," -5  MAG3        1   12    3    0\n");
    //    fprintf(f1," -5  PHA1        1   12    4    0\n");
    //    fprintf(f1," -5  PHA2        1   12    5    0\n");
    //    fprintf(f1," -5  PHA3        1   12    6    0\n");
    //
    //    frdselect(vr,vi,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
    //                nfieldvectph,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the temperatures in the nodes (magnitude, phase) */
    //	  
    //  if(strcmp1(&filab[957],"PNT ")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[957],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  PNDTEMP     2    1\n");
    //    fprintf(f1," -5  MAG1        1    1    1    0\n");
    //    fprintf(f1," -5  PHA1        1    1    2    0\n");
    //
    //    frdselect(vr,vi,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompscalph,ifieldscalph,icompscalph,
    //                nfieldscalph,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the stresses in the nodes (magnitude, phase) */
    //	  
    //  if(strcmp1(&filab[1479],"PHS ")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[1479],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  PSTRESS    12    1\n");
    //    fprintf(f1," -5  MAGXX       1    4    1    1\n");
    //    fprintf(f1," -5  MAGYY       1    4    2    2\n");
    //    fprintf(f1," -5  MAGZZ       1    4    3    3\n");
    //    fprintf(f1," -5  MAGXY       1    4    1    2\n");
    //    fprintf(f1," -5  MAGYZ       1    4    2    3\n");
    //    fprintf(f1," -5  MAGZX       1    4    3    1\n");
    //    fprintf(f1," -5  PHAXX       1    4    1    1\n");
    //    fprintf(f1," -5  PHAYY       1    4    2    2\n");
    //    fprintf(f1," -5  PHAZZ       1    4    3    3\n");
    //    fprintf(f1," -5  PHAXY       1    4    1    2\n");
    //    fprintf(f1," -5  PHAYZ       1    4    2    3\n");
    //    fprintf(f1," -5  PHAZX       1    4    3    1\n");
    //
    //    frdselect(stnr,stni,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomptensph,ifieldtensph,icomptensph,
    //                nfieldtensph,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the displacements in the nodes (magnitude, phase) */
    //	  
    //  if(strcmp1(&filab[2610],"PRF ")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[2610],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  PFORC       6    1\n");
    //    fprintf(f1," -5  MAG1        1   12    1    0\n");
    //    fprintf(f1," -5  MAG2        1   12    2    0\n");
    //    fprintf(f1," -5  MAG3        1   12    3    0\n");
    //    fprintf(f1," -5  PHA1        1   12    4    0\n");
    //    fprintf(f1," -5  PHA2        1   12    5    0\n");
    //    fprintf(f1," -5  PHA3        1   12    6    0\n");
    //
    //    frdselect(fnr,fni,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
    //                nfieldvectph,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    /* the remaining parts are for frequency calculations with cyclic symmetry only */
  
    // if(*nmethod!=2){ex_close(exoid);return;}
  
    //  /* storing the maximum displacements of the nodes in the base sector
    //     (components, magnitude) */
    //	  
    //  if(strcmp1(&filab[1566],"MAXU")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[1566],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  MDISP       4    1\n");
    //    fprintf(f1," -5  DX          1    4    1    0\n");
    //    fprintf(f1," -5  DY          1    4    2    0\n");
    //    fprintf(f1," -5  DZ          1    4    3    0\n");
    //    fprintf(f1," -5  ANG         1    4    4    0\n");
    //    
    //    ncomp=4;
    //    ifield[0]=1;icomp[0]=1;
    //    ifield[1]=1;icomp[1]=2;
    //    ifield[2]=1;icomp[2]=3;
    //    ifield[3]=1;icomp[3]=0;
    //    nfield[0]=4;nfield[1]=4;
    //
    //    frdselect(vmax,vmax,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomp,ifield,icomp,
    //                nfield,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the worst principal stress at the nodes
    //     in the basis sector (components, magnitude)
    //
    //     the worst principal stress is the maximum of the
    //     absolute value of all principal stresses, times
    //     its original sign */
    //	  
    //  if(strcmp1(&filab[1653],"MAXS")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[1653],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  MSTRESS     7    1\n");
    //    fprintf(f1," -5  SXX         1    4    1    1\n");
    //    fprintf(f1," -5  SYY         1    4    2    2\n");
    //    fprintf(f1," -5  SZZ         1    4    3    3\n");
    //    fprintf(f1," -5  SXY         1    4    1    2\n");
    //    fprintf(f1," -5  SYZ         1    4    2    3\n");
    //    fprintf(f1," -5  SZX         1    4    3    1\n");
    //    fprintf(f1," -5  MAG         1    4    0    0\n");
    //    
    //    ncomp=7;
    //    ifield[0]=1;icomp[0]=1;
    //    ifield[1]=1;icomp[1]=2;
    //    ifield[2]=1;icomp[2]=3;
    //    ifield[3]=1;icomp[3]=4;
    //    ifield[4]=1;icomp[4]=6;
    //    ifield[5]=1;icomp[5]=5;
    //    ifield[6]=1;icomp[6]=0;
    //    nfield[0]=7;nfield[1]=7;
    //
    //    frdselect(stnmax,stnmax,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomp,ifield,icomp,
    //                nfield,&iselect,m2,f1,output,m3);
    //
    //  }
    //
    //  /* storing the worst principal strain at the nodes
    //     in the basis sector (components, magnitude)
    //
    //     the worst principal strain is the maximum of the
    //     absolute value of all principal strains, times
    //     its original sign */
    //	  
    //  if(strcmp1(&filab[2523],"MAXE")==0){
    //    iselect=1;
    //    
    //    frdset(&filab[2523],set,&iset,istartset,iendset,ialset,
    //	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
    //	   ngraph);
    //    
    //    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
    //	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    //
    //    fprintf(f1," -4  MSTRAIN     7    1\n");
    //    fprintf(f1," -5  EXX         1    4    1    1\n");
    //    fprintf(f1," -5  EYY         1    4    2    2\n");
    //    fprintf(f1," -5  EZZ         1    4    3    3\n");
    //    fprintf(f1," -5  EXY         1    4    1    2\n");
    //    fprintf(f1," -5  EYZ         1    4    2    3\n");
    //    fprintf(f1," -5  EZX         1    4    3    1\n");
    //    fprintf(f1," -5  MAG         1    4    0    0\n");
    //    
    //    ncomp=7;
    //    ifield[0]=1;icomp[0]=1;
    //    ifield[1]=1;icomp[1]=2;
    //    ifield[2]=1;icomp[2]=3;
    //    ifield[3]=1;icomp[3]=4;
    //    ifield[4]=1;icomp[4]=6;
    //    ifield[5]=1;icomp[5]=5;
    //    ifield[6]=1;icomp[6]=0;
    //    nfield[0]=7;nfield[1]=7;
    //
    //    frdselect(eenmax,eenmax,&iset,&nkcoords,inum,m1,istartset,iendset,
    //                ialset,ngraph,&ncomp,ifield,icomp,
    //                nfield,&iselect,m2,f1,output,m3);
    //
    //  }
    
    if (countbool==3){
      errr = ex_put_var_param (exoid, "n", countvars);
      ex_update (exoid);  
      // var_names = (char *) calloc (countvars, sizeof (char));
      // var_names = (char *) calloc (countvars, MAX_STR_LENGTH);
      for (i=0; i<countvars; i++)
	{
	  var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
      // printf ("Total count %i\n", countvars);
    }else if(countbool==2){
      // printf ("Total count %i\n", countvars);
      // printf ("Writing %i var_names\n", countvars);
      errr = ex_put_var_names (exoid, "n", countvars, var_names);
      if (errr) printf ("Unable to update variable names.");
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


  ex_update (exoid);  
  ex_close(exoid);
  return;
}
