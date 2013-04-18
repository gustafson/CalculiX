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
  

  char c[2]="C",m1[4]=" -1",m2[4]=" -2",m3[4]=" -3",
    p0[6]="    0",p1[6]="    1",p2[6]="    2",p3[6]="    3",p4[6]="    4",
    p5[6]="    5",p6[6]="    6",p7[6]="    7",p8[6]="    8",p9[6]="    9",
    p10[6]="   10",p11[6]="   11",
    p12[6]="   12", fneig[132]="",date[8],clock[10],newdate[21],newclock[9],
    material[6]="     ",text[2]=" ";

  static int icounter=0,nkcoords,nout,noutmin,noutplus;

  int null,one,i,j,k,indexe,nemax,nlayer,noutloc,iset,iselect,ncomp,nope,
    nodes,ifield[7],nfield[2],icomp[7],ifieldstate[*nstate_],two,three,
    icompstate[*nstate_],ip0=0,ip1=1,ip2=2,ip3=3,ip4=4,ip5=5,ip6=6,ip7=7,
    ip8=8,ip9=9,ip10=10,ip11=11,ip12=12,imaterial=0,nelout;

  int ncompscalar=1,ifieldscalar[1]={1},icompscalar[1]={0},nfieldscalar[2]={1,0};
  int ncompvector=3,ifieldvector[3]={1,1,1},icompvector[3]={0,1,2},nfieldvector1[2]={3,0},nfieldvector0[2]={mi[1]+1,0};
  int ncomptensor=6,ifieldtensor[6]={1,1,1,1,1,1},icomptensor[6]={0,1,2,3,5,4},nfieldtensor[2]={6,0};
  int ncompscalph=2,ifieldscalph[2]={1,2},icompscalph[2]={0,0},nfieldscalph[2]={0,0};
  int ncompvectph=6,ifieldvectph[6]={1,1,1,2,2,2},icompvectph[6]={1,2,3,1,2,3},nfieldvectph[2]={mi[1]+1,mi[1]+1};
  int ncomptensph=12,ifieldtensph[12]={1,1,1,1,1,1,2,2,2,2,2,2},icomptensph[12]={0,1,2,3,5,4,0,1,2,3,5,4},nfieldtensph[2]={6,6};
  int iw;
  float ifl;
  double pi,oner;
  
  strcpy(fneig,jobnamec);
  strcat(fneig,".exo");
  
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

    int num_nodes = nout;   
    int num_dim, num_elem;
    int num_elem_blk; /* Node element blocks.  One eltype per block*/
    int num_ns, num_ss; /* Node sets, side sets */
    int errr;
    int CPU_word_size = sizeof(float);
    int IO_word_size = sizeof(float);

    num_dim = 3;  
    num_elem_blk = 10;
    num_ns = 0; 
    num_ss = 0;
    int exoid = ex_create (fneig, /*Filename*/
			   EX_CLOBBER,	/* create mode */
			   &CPU_word_size,  /* CPU float word size in bytes */
			   &IO_word_size);  /* I/O float word size in bytes */
  
    
    float x[*nk], y[*nk], z[*nk];
    // Write optional node map
    j = 0;
    // int node_map[*nk];
    int *node_map;
    node_map = (int *) calloc(*nk, sizeof(int));
    
    /* nkcoords is the number of nodes at the time when 
       the nodal coordinates are stored in the exo file */
    nkcoords = *nk;
    num_nodes = nkcoords;

    /* storing the coordinates of the nodes */
    if(*nmethod!=0){
      for(i=0;i<*nk;i++){
	// if(inum[i]==0){
	//   // Put other nodes at the origin
	//   x[i] = 0.0;
	//   y[i] = 0.0;
	//   z[i] = 0.0;
	//   continue;}
	printf ("Node map %i\n", i+1);
	node_map[i] = i+1;
	printf ("Node map %i\n", node_map[i]);
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
      exit(0);
    }
    
    /* write values to database */
    errr = ex_put_coord (exoid, x, y, z);
    errr = ex_put_node_num_map (exoid, node_map);
    if(errr){
      printf("*ERROR in exo: failed to write node map");
      exit(0);
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
      exit(0);
    }
    
 
    //SAMPLE    // Cell eltypes = args(1).cell_value();
    //SAMPLE    // 
    //SAMPLE    
    //SAMPLE    for (i = 0; i < num_elem; i++)
    //SAMPLE      {
    //SAMPLE    	elementsval(i) = elements.elem (i); 
    //SAMPLE    	// Pull the first int for the element map
    //SAMPLE    	elem_map[i] = elementsval(i).matrix_value()(0);
    //SAMPLE    	eltypesval(i) = eltypes.elem (i);
    //SAMPLE      }
    //SAMPLE    
    //SAMPLE    // write out the connectivity array, and write out the element attributes array:
    //SAMPLE    int num_elem_in_blk, num_nodes_per_elem, num_attr;
    //SAMPLE    float *attrib;
    //SAMPLE    /* write element block parameters */
    //SAMPLE    int blkid = 10; 
    //SAMPLE    num_elem_in_blk = num_elem;
    //SAMPLE    num_nodes_per_elem = elementsval(0).columns()-1;
    //SAMPLE    num_attr = 1;
    //SAMPLE    
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
    //SAMPLE    /* write element connectivity */
    //SAMPLE    int *connect;
    //SAMPLE    connect = (int *) calloc (num_elem_in_blk*num_nodes_per_elem, sizeof(int));
    //SAMPLE    j=0; k=0;
    //SAMPLE    for (i = 0; i < num_elem_in_blk; i++){
    //SAMPLE      for (j = 1; j <= num_nodes_per_elem; j++, k++){
    //SAMPLE    	connect[k] = elementsval(i).matrix_value()(j);
    //SAMPLE      }
    //SAMPLE    }
    //SAMPLE    
    //SAMPLE    int ex_update (exoid);
    //SAMPLE    
    //SAMPLE    int *idelbs;
    //SAMPLE    idelbs = (int *) calloc(10, sizeof(int));
    //SAMPLE    errr = ex_get_elem_blk_ids (exoid, idelbs);
    //SAMPLE    
    //SAMPLE    // FIXME, this call currently does not find the blkid
    //SAMPLE    errr = ex_put_elem_conn (exoid, idelbs[0], connect);
    //SAMPLE    if (errr)
    //SAMPLE      printf ("ERROR in ex_put_elem_conn %i\n", errr);
    //SAMPLE    free (connect);
    //SAMPLE    
    //SAMPLE    int ex_update (exoid);
    //SAMPLE    
    
    /* storing the topology */
    nemax=*ne0;
    num_elem=nemax;

    // Initialize enough memory to store the element numbers
    int *elem_map;
    elem_map = (int *) calloc(num_elem, sizeof(int));
    
    for(i=0;i<*ne0;i++){
      if(ipkon[i]<0) continue;
      strcpy1(material,&matname[80*(ielmat[i*mi[2]]-1)],5);
      indexe=ipkon[i];
      elem_map[i] = i+1;
      if(strcmp1(&lakon[8*i+3],"2")==0){
	if(((strcmp1(&lakon[8*i+6]," ")==0)||
	    (strcmp1(&filab[4],"E")==0)||
	    (strcmp1(&lakon[8*i+6],"I")==0))&&
	   (strcmp2(&lakon[8*i+6],"LC",2)!=0)){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p4,p0,material,m2);
	  // for(j=0;j<10;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n%3s",m2);
	  // for(j=10;j<12;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // for(j=16;j<19;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // for(j=19;j<20;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // for(j=12;j<16;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n");
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
	}else if(strcmp1(&lakon[8*i+6],"B")==0){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p12,p0,material);
	  // fprintf(f1,"%3s%10d%10d%10d\n",m2,kon[indexe+20],kon[indexe+22],kon[indexe+21]);
	}else{
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p10,p0,material,m2);
	  // for(j=0;j<8;j++)fprintf(f1,"%10d",kon[indexe+20+j]);
	  // fprintf(f1,"\n");
	}
      }else if(strcmp1(&lakon[8*i+3],"8")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p1,p0,material,m2);
	  // for(j=0;j<8;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n");
	}else if(strcmp1(&lakon[8*i+6],"B")==0){
	  if(strcmp1(&lakon[8*i+4],"R")==0){
	    printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	    // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p11,p0,material);
	    // fprintf(f1,"%3s%10d%10d\n",m2,kon[indexe+8],kon[indexe+9]);
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	    // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p11,p0,material);
	    // fprintf(f1,"%3s%10d%10d\n",m2,kon[indexe+11],kon[indexe+12]);
	  }
	}else{
	  if(strcmp1(&lakon[8*i+4],"R")==0){
	    printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	    // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p9,p0,material,m2);
	    // for(j=0;j<4;j++)fprintf(f1,"%10d",kon[indexe+8+j]);
	    // fprintf(f1,"\n");
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	    // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p9,p0,material,m2);
	    // for(j=0;j<4;j++)fprintf(f1,"%10d",kon[indexe+11+j]);
	    // fprintf(f1,"\n");
	  }
	}
      }else if(strcmp1(&lakon[8*i+3],"10")==0){
	printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	// fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p6,p0,material,m2);
	// for(j=0;j<10;j++)fprintf(f1,"%10d",kon[indexe+j]);
	// fprintf(f1,"\n");
      }else if(strcmp1(&lakon[8*i+3],"4")==0){
	printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	// fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p3,p0,material,m2);
	// for(j=0;j<4;j++)fprintf(f1,"%10d",kon[indexe+j]);
	// fprintf(f1,"\n");
      }else if(strcmp1(&lakon[8*i+3],"15")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p5,p0,material,m2);
	  // for(j=0;j<9;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // for(j=12;j<13;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n%3s",m2);
	  // for(j=13;j<15;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // for(j=9;j<12;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n");
	}else{
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p8,p0,material,m2);
	  // for(j=0;j<6;j++)fprintf(f1,"%10d",kon[indexe+15+j]);
	  // fprintf(f1,"\n");
	}
      }else if(strcmp1(&lakon[8*i+3],"6")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p2,p0,material,m2);
	  // for(j=0;j<6;j++)fprintf(f1,"%10d",kon[indexe+j]);
	  // fprintf(f1,"\n");
	}else{
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n%3s",m1,i+1,p7,p0,material,m2);
	  // for(j=0;j<3;j++)fprintf(f1,"%10d",kon[indexe+6+j]);
	  // fprintf(f1,"\n");
	}
      }else if((strcmp1(&lakon[8*i],"D")==0)&&
	       (strcmp1(&lakon[8*i],"DCOUP3D")!=0)){
	if(kon[indexe]==0){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p11,p0,material);
	  // fprintf(f1,"%3s%10d%10d\n",m2,kon[indexe+1],kon[indexe+2]);
	}else if(kon[indexe+2]==0){
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p11,p0,material);
	  // fprintf(f1,"%3s%10d%10d\n",m2,kon[indexe],kon[indexe+1]);
	}else{
	  printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	  // fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p12,p0,material);
	  // fprintf(f1,"%3s%10d%10d%10d\n",m2,kon[indexe],kon[indexe+2],kon[indexe+1]);
	}
      }else if((strcmp1(&lakon[8*i],"E")==0)&&
	       (strcmp1(&lakon[8*i+6],"A")==0)){
	printf ("Element %i type %s material %s\n", i+1, &lakon[8*i], material);
	// fprintf(f1,"%3s%10d%5s%5s%5s\n",m1,i+1,p11,p0,material);
	// fprintf(f1,"%3s%10d%10d\n",m2,kon[indexe],kon[indexe+1]);
      }
    }

    // Write the element map into the file
    errr = ex_put_elem_num_map (exoid, elem_map); 
    if (errr)
      printf ("ERROR in ex_put_elem_num_map %i\n", errr);
    free (elem_map);

    ex_update (exoid);  
    ex_close (exoid);

    // if(*nmethod==0){fclose(f1);return;}
    
    /* End of if(*kode==1) */
  }

  return;
}
