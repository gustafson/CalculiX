/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                          */

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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void mastructcs(int *nk, int *kon, int *ipkon, char *lakon, int *ne,
	      int *nodeboun, int *ndirboun, int *nboun, int *ipompc,
	      int *nodempc, int *nmpc, int *nactdof, int *icol,
	      int *jq, int **mast1p, int **irowp, int *isolver, int *neq,
	      int *ikmpc, int *ilmpc,int *ipointer, int *nzs, 
              int *nmethod,int *ics, double *cs, char *labmpc, int *mcs, 
              int *mi,int *mortar){

  /* determines the structure of the thermo-mechanical matrices with
     cyclic symmetry;
     (i.e. the location of the nonzeros */

  char lakonl[2]=" \0";

  int i,j,k,l,jj,ll,id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,
    ist1,ist2,node1,node2,isubtract,nmast,ifree,istart,istartold,
    index1,index2,m,node,nzs_,ist,kflag,indexe,nope,isize,*mast1=NULL,
    *irow=NULL,inode,icomplex,inode1,icomplex1,inode2,
    icomplex2,kdof1,kdof2,ilength,lprev,ij,mt=mi[1]+1;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irow=*irowp;

  kflag=2;
  nzs_=nzs[1];

  /* initialisation of nactmpc */

  for(i=0;i<mt**nk;++i){nactdof[i]=0;}

  /* determining the active degrees of freedom due to elements */

  for(i=0;i<*ne;++i){
    
    if(ipkon[i]<0) continue;
    indexe=ipkon[i];
/* Bernhardi start */
    if (strcmp1(&lakon[8*i+3],"8I")==0)nope=11;
    else if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
/* Bernhardi end */
    else if(strcmp1(&lakon[8*i+3],"2")==0)nope=26;
    else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
    else if ((strcmp1(&lakon[8*i+3],"4")==0)||
	     (strcmp1(&lakon[8*i+2],"4")==0)) nope=4;
    else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else if (strcmp1(&lakon[8*i],"E")==0){
	if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
	    nope=kon[ipkon[i]-1];
	}else{
	    lakonl[0]=lakon[8*i+7];
	    nope=atoi(lakonl)+1;
	}
    }else continue;

/*    else if (strcmp1(&lakon[8*i],"E")==0){
	lakonl[0]=lakon[8*i+7];
	nope=atoi(lakonl)+1;}
	else continue;*/

    for(j=0;j<nope;++j){
      node=kon[indexe+j]-1;
      for(k=1;k<4;++k){
	nactdof[mt*node+k]=1;
      }
    }
  }

  /* determining the active degrees of freedom due to mpc's */

  for(i=0;i<*nmpc;++i){
      index=ipompc[i]-1;
      do{
	  if(nodempc[3*index+1]!=0){
//	      nactdof[mt*nodempc[3*index]+nodempc[3*index+1]-4]=1;}
	      nactdof[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]=1;}
	  index=nodempc[3*index+2];
	  if(index==0) break;
	  index--;
      }while(1);
  }
	   
  /* subtracting the SPC and MPC nodes */

  for(i=0;i<*nboun;++i){
    nactdof[mt*(nodeboun[i]-1)+ndirboun[i]]=0;
  }

  for(i=0;i<*nmpc;++i){
    index=ipompc[i]-1;
    nactdof[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]=0;
  }
  
  /* numbering the active degrees of freedom */
  
  neq[0]=0;
  for(i=0;i<*nk;++i){
    for(j=1;j<4;++j){
	if(nactdof[mt*i+j]!=0){
	++neq[0];
	nactdof[mt*i+j]=neq[0];
      }
    }
  }
  
  ifree=0;
  
    /* determining the position of each nonzero matrix element

       mast1(ipointer(i)) = first nonzero row in column i
       irow(ipointer(i))  points to further nonzero elements in 
                             column i */
      
  for(i=0;i<6**nk;++i){ipointer[i]=0;}
    
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    indexe=ipkon[i];
/*  Bernhardi start  */
    if(strcmp1(&lakon[8*i],"C3D8I")==0){nope=11;}
    else if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
/*  Bernhardi end */
    else if(strcmp1(&lakon[8*i+3],"2")==0)nope=26;
    else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else if (strcmp1(&lakon[8*i],"E")==0){
	if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
	    nope=kon[ipkon[i]-1];
	}else{
	    lakonl[0]=lakon[8*i+7];
	    nope=atoi(lakonl)+1;
	}
    }else continue;

/*    else if (strcmp1(&lakon[8*i],"E")==0){
	lakonl[0]=lakon[8*i+7];
	nope=atoi(lakonl)+1;}
	else continue;*/
      
    for(jj=0;jj<3*nope;++jj){
	
      j=jj/3;
      k=jj-3*j;
	
      node1=kon[indexe+j];
      jdof1=nactdof[mt*(node1-1)+k+1];
	
      for(ll=jj;ll<3*nope;++ll){
	  
	l=ll/3;
	m=ll-3*l;
	  
	node2=kon[indexe+l];
	jdof2=nactdof[mt*(node2-1)+m+1];
	  
	/* check whether one of the DOF belongs to a SPC or MPC */
	  
	if((jdof1!=0)&&(jdof2!=0)){
	  insert(ipointer,&mast1,&irow,&jdof1,&jdof2,&ifree,&nzs_);
	  kdof1=jdof1+neq[0];kdof2=jdof2+neq[0];
	  insert(ipointer,&mast1,&irow,&kdof1,&kdof2,&ifree,&nzs_);
	}
	else if((jdof1!=0)||(jdof2!=0)){
	  
	  /* idof1: genuine DOF
	     idof2: nominal DOF of the SPC/MPC */
	  
	  if(jdof1==0){
	    idof1=jdof2;
	    idof2=8*node1+k-7;}
	  else{
	    idof1=jdof1;
	    idof2=8*node2+m-7;}
	  
	  if(*nmpc>0){
	    
	    FORTRAN(nident,(ikmpc,&idof2,nmpc,&id));
	    if((id>0)&&(ikmpc[id-1]==idof2)){
	      
	      /* regular DOF / MPC */
	      
	      id1=ilmpc[id-1];
	      ist=ipompc[id1-1];
	      index=nodempc[3*ist-1];
	      if(index==0) continue;
	      while(1){
		inode=nodempc[3*index-3];
		icomplex=0;
		if(strcmp1(&labmpc[(id1-1)*20],"CYCLIC")==0){
                  icomplex=atoi(&labmpc[20*(id1-1)+6]);
		}
		else if(strcmp1(&labmpc[(id1-1)*20],"SUBCYCLIC")==0){
                  for(ij=0;ij<*mcs;ij++){
                    ilength=cs[17*ij+3];
                    lprev=cs[17*ij+13];
                    FORTRAN(nident,(&ics[lprev],&inode,&ilength,&id));
                    if(id>0){
                      if(ics[lprev+id-1]==inode){
                        icomplex=ij+1;
                        break;
                      }
                    }
                  }
		}
//		idof2=nactdof[mt*inode+nodempc[3*index-2]-4];
		idof2=nactdof[mt*(inode-1)+nodempc[3*index-2]];
		if(idof2!=0){
		  insert(ipointer,&mast1,&irow,&idof1,&idof2,&ifree,&nzs_);
		  kdof1=idof1+neq[0];kdof2=idof2+neq[0];
		  insert(ipointer,&mast1,&irow,&kdof1,&kdof2,&ifree,&nzs_);
		  if((icomplex!=0)&&(idof1!=idof2)){
		    insert(ipointer,&mast1,&irow,&kdof1,&idof2,&ifree,&nzs_);
		    insert(ipointer,&mast1,&irow,&idof1,&kdof2,&ifree,&nzs_);
		  }
		}
		index=nodempc[3*index-1];
		if(index==0) break;
	      }
	      continue;
	    }
	  }
	}
	
	else{
	  idof1=8*node1+k-7;
	  idof2=8*node2+m-7;
	  mpc1=0;
	  mpc2=0;
	  if(*nmpc>0){
	    FORTRAN(nident,(ikmpc,&idof1,nmpc,&id1));
	    if((id1>0)&&(ikmpc[id1-1]==idof1)) mpc1=1;
	    FORTRAN(nident,(ikmpc,&idof2,nmpc,&id2));
	    if((id2>0)&&(ikmpc[id2-1]==idof2)) mpc2=1;
	  }
	  if((mpc1==1)&&(mpc2==1)){
	    id1=ilmpc[id1-1];
	    id2=ilmpc[id2-1];
	    if(id1==id2){
	      
	      /* MPC id1 / MPC id1 */
	      
	      ist=ipompc[id1-1];
	      index1=nodempc[3*ist-1];
	      if(index1==0) continue;
	      while(1){
		inode1=nodempc[3*index1-3];
		icomplex1=0;
		if(strcmp1(&labmpc[(id1-1)*20],"CYCLIC")==0){
                  icomplex1=atoi(&labmpc[20*(id1-1)+6]);
		}
		else if(strcmp1(&labmpc[(id1-1)*20],"SUBCYCLIC")==0){
                  for(ij=0;ij<*mcs;ij++){
                    ilength=cs[17*ij+3];
                    lprev=cs[17*ij+13];
                    FORTRAN(nident,(&ics[lprev],&inode1,&ilength,&id));
                    if(id>0){
                      if(ics[lprev+id-1]==inode1){
                        icomplex1=ij+1;
                        break;
                      }
                    }
                  }
		}
//		idof1=nactdof[mt*inode1+nodempc[3*index1-2]-4];
		idof1=nactdof[mt*(inode1-1)+nodempc[3*index1-2]];
		index2=index1;
		while(1){
		  inode2=nodempc[3*index2-3];
		  icomplex2=0;
		  if(strcmp1(&labmpc[(id1-1)*20],"CYCLIC")==0){
                    icomplex2=atoi(&labmpc[20*(id1-1)+6]);
		  }
		  else if(strcmp1(&labmpc[(id1-1)*20],"SUBCYCLIC")==0){
                    for(ij=0;ij<*mcs;ij++){
                      ilength=cs[17*ij+3];
                      lprev=cs[17*ij+13];
                      FORTRAN(nident,(&ics[lprev],&inode2,&ilength,&id));
                      if(id>0){
                        if(ics[lprev+id-1]==inode2){
                          icomplex2=ij+1;
                          break;
                        }
                      }
                    }
                  }
//		  idof2=nactdof[mt*inode2+nodempc[3*index2-2]-4];
		  idof2=nactdof[mt*(inode2-1)+nodempc[3*index2-2]];
		  if((idof1!=0)&&(idof2!=0)){
		    insert(ipointer,&mast1,&irow,&idof1,&idof2,&ifree,&nzs_);
		    kdof1=idof1+neq[0];kdof2=idof2+neq[0];
		    insert(ipointer,&mast1,&irow,&kdof1,&kdof2,&ifree,&nzs_);
                    if(((icomplex1!=0)||(icomplex2!=0))&&
                       (icomplex1!=icomplex2)){
                    /*   if(((icomplex1!=0)||(icomplex2!=0))&&
                         ((icomplex1==0)||(icomplex2==0))){*/
		      insert(ipointer,&mast1,&irow,&kdof1,&idof2,&ifree,&nzs_);
		      insert(ipointer,&mast1,&irow,&idof1,&kdof2,&ifree,&nzs_);
		    }
		  }
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
	    
	    else{
	      
	      /* MPC id1 /MPC id2 */
	      
	      ist1=ipompc[id1-1];
	      index1=nodempc[3*ist1-1];
	      if(index1==0) continue;
	      while(1){
		inode1=nodempc[3*index1-3];
		icomplex1=0;
		if(strcmp1(&labmpc[(id1-1)*20],"CYCLIC")==0){
                  icomplex1=atoi(&labmpc[20*(id1-1)+6]);
		}
		else if(strcmp1(&labmpc[(id1-1)*20],"SUBCYCLIC")==0){
                  for(ij=0;ij<*mcs;ij++){
                    ilength=cs[17*ij+3];
                    lprev=cs[17*ij+13];
                    FORTRAN(nident,(&ics[lprev],&inode1,&ilength,&id));
                    if(id>0){
                      if(ics[lprev+id-1]==inode1){
                        icomplex1=ij+1;
                        break;
                      }
                    }
                  }
		}
//		idof1=nactdof[mt*inode1+nodempc[3*index1-2]-4];
		idof1=nactdof[mt*(inode1-1)+nodempc[3*index1-2]];
		ist2=ipompc[id2-1];
		index2=nodempc[3*ist2-1];
		if(index2==0){
		  index1=nodempc[3*index1-1];
		  if(index1==0){break;}
		  else{continue;}
		}
		while(1){
		  inode2=nodempc[3*index2-3];
		  icomplex2=0;
		  if(strcmp1(&labmpc[(id2-1)*20],"CYCLIC")==0){
                    icomplex2=atoi(&labmpc[20*(id2-1)+6]);
		  }
		  else if(strcmp1(&labmpc[(id2-1)*20],"SUBCYCLIC")==0){
                    for(ij=0;ij<*mcs;ij++){
                      ilength=cs[17*ij+3];
                      lprev=cs[17*ij+13];
                      FORTRAN(nident,(&ics[lprev],&inode2,&ilength,&id));
                      if(id>0){
                        if(ics[lprev+id-1]==inode2){
                          icomplex2=ij+1;
                          break;
                        }
                      }
                    }
                  }
//		  idof2=nactdof[mt*inode2+nodempc[3*index2-2]-4];
		  idof2=nactdof[mt*(inode2-1)+nodempc[3*index2-2]];
		  if((idof1!=0)&&(idof2!=0)){
		    insert(ipointer,&mast1,&irow,&idof1,&idof2,&ifree,&nzs_);
		    kdof1=idof1+neq[0];kdof2=idof2+neq[0];
		    insert(ipointer,&mast1,&irow,&kdof1,&kdof2,&ifree,&nzs_);
                    if(((icomplex1!=0)||(icomplex2!=0))&&
                       (icomplex1!=icomplex2)){
                    /*   if(((icomplex1!=0)||(icomplex2!=0))&&
                         ((icomplex1==0)||(icomplex2==0))){*/
		      insert(ipointer,&mast1,&irow,&kdof1,&idof2,&ifree,&nzs_);
		      insert(ipointer,&mast1,&irow,&idof1,&kdof2,&ifree,&nzs_);
		    }
		  }
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
	  }
	}
      }
    }
  }

  neq[0]=2*neq[0];
  neq[1]=neq[0];
  
  /* ordering the nonzero nodes in the SUPERdiagonal columns
     mast1 contains the row numbers column per column,
     irow the column numbers */
  
/*  for(i=0;i<neq[0];++i){
    itot=0;
    if(ipointer[i]==0){
      printf("*ERROR in mastructcs: zero column");
      FORTRAN(stop,());
    }
    istart=ipointer[i];
    while(1){
      ++itot;
      ikcol[itot-1]=mast1[istart-1];
      istart=irow[istart-1];
      if(istart==0) break;
    }
    FORTRAN(isortii,(ikcol,icol,&itot,&kflag));
    istart=ipointer[i];
    for(j=0;j<itot-1;++j){
      mast1[istart-1]=ikcol[j];
      istartold=istart;
      istart=irow[istart-1];
      irow[istartold-1]=i+1;
    }
    mast1[istart-1]=ikcol[itot-1];
    irow[istart-1]=i+1;
    }*/

    
  for(i=0;i<neq[0];++i){
      if(ipointer[i]==0){
	  if(i>=neq[1]) continue;
	  printf("*ERROR in mastructcs: zero column\n");
	  FORTRAN(stop,());
      }
      istart=ipointer[i];
      while(1){
	  istartold=istart;
	  istart=irow[istart-1];
	  irow[istartold-1]=i+1;
	  if(istart==0) break;
      }
  }
  
  if(neq[0]==0){
    printf("\n*WARNING: no degrees of freedom in the model\n");
    FORTRAN(stop,());
  }
  
  printf(" number of equations\n");
  printf(" %d\n",neq[0]);
  printf(" number of nonzero lower triangular matrix elements\n");
  printf(" %d\n",ifree-neq[0]);
  
  /* new meaning of icol,j1,mast1,irow:
     
     - irow is going to contain the row numbers of the SUBdiagonal
     nonzero's, column per column
     - mast1 contains the column numbers
     - icol(i)=# SUBdiagonal nonzero's in column i
     - jq(i)= location in field irow of the first SUBdiagonal
     nonzero in column i
     
     */
  
  nmast=ifree;
  
  /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
  
  FORTRAN(isortii,(mast1,irow,&nmast,&kflag));
  
  /* filtering out the diagonal elements and generating icol and jq */
  
  isubtract=0;
  for(i=0;i<neq[0];++i){icol[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
    if(mast1[i]==irow[i]){++isubtract;}
    else{
      mast1[i-isubtract]=mast1[i];
      irow[i-isubtract]=irow[i];
      if(k!=mast1[i]){
	for(l=k;l<mast1[i];++l){jq[l]=i+1-isubtract;}
	k=mast1[i];
      }
      ++icol[k-1];
    }
  }
  nmast=nmast-isubtract;
  for(l=k;l<neq[0]+1;++l){jq[l]=nmast+1;}
  
  for(i=0;i<neq[0];++i){
    if(jq[i+1]-jq[i]>0){
      isize=jq[i+1]-jq[i];
      FORTRAN(isortii,(&irow[jq[i]-1],&mast1[jq[i]-1],&isize,&kflag));
    }
  }
  
  nzs[0]=jq[neq[0]-1]-1;
  nzs[1]=nzs[0];
  nzs[2]=nzs[0];
  
  *mast1p=mast1;
  *irowp=irow;
  
  return;
  
}
