/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2024 Guido Dhondt                          */

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
#include <pthread.h>
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void mastructmatrixcs(ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *nactdof,
		      ITG *jq,ITG **mast1p,ITG *neq,ITG *ipointer, ITG *nzs_, 
		      ITG *nmethod,ITG *mi,ITG **nextp,
		      ITG *node1,ITG *k,ITG *node2,ITG *m,ITG *ifree,
		      char *labmpc,ITG *mcs,double *cs,ITG *ics){

  /* determines the structure of the thermo-mechanical matrices;
     (i.e. the location of the nonzeros */

  ITG id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,ist1,ist2,
    index1,index2,ist,*mast1=NULL,mt=mi[1]+1,*next=NULL,ij,
    inode1,inode2,icomplex1,icomplex2,kdof1,kdof2,inode,icomplex,
    ilength,lprev;

  mast1=*mast1p;next=*nextp;

  /* caveat: k and m take values 1..3 for dof in x,y,z
     (FORTRAN convention) */
  
  jdof1=nactdof[mt*(*node1-1)+(*k)];
  jdof2=nactdof[mt*(*node2-1)+(*m)];

  //  printf("%d,%d,%d,%d,%d,%d\n",*node1,*k,*node2,*m,jdof1,jdof2);
  
  /* check whether one of the DOF belongs to a SPC or MPC */
  	  
  if((jdof1>0)&&(jdof2>0)){
    insert(ipointer,&mast1,&next,&jdof1,&jdof2,ifree,nzs_);
    kdof1=jdof1+neq[0];kdof2=jdof2+neq[0];
    insert(ipointer,&mast1,&next,&kdof1,&kdof2,ifree,nzs_);
  }
  else if((jdof1>0)||(jdof2>0)){
	  
    /* idof1: genuine DOF
       idof2: nominal DOF of the SPC/MPC */
	  
    if(jdof1<=0){
      idof1=jdof2;
      idof2=jdof1;}
    else{
      idof1=jdof1;
      idof2=jdof2;}
	  
    if(*nmpc>0){
	    
      if(idof2!=2*(idof2/2)){
	      
	/* regular DOF / MPC */
	      
	id1=(-idof2+1)/2;
	ist=ipompc[id1-1];
	index=nodempc[3*ist-1];
	if(index==0) return;
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
	  idof2=nactdof[mt*(inode-1)+nodempc[3*index-2]];
	  if(idof2>0){
	    insert(ipointer,&mast1,&next,&idof1,&idof2,ifree,nzs_);
	    kdof1=idof1+neq[0];kdof2=idof2+neq[0];
	    insert(ipointer,&mast1,&next,&kdof1,&kdof2,ifree,nzs_);
	    if((icomplex!=0)&&(idof1!=idof2)){
	      insert(ipointer,&mast1,&next,&kdof1,&idof2,ifree,nzs_);
	      insert(ipointer,&mast1,&next,&idof1,&kdof2,ifree,nzs_);
	    }
	  }
	  index=nodempc[3*index-1];
	  if(index==0) break;
	}
	return;
      }
    }
  }
	
  else{
    idof1=jdof1;
    idof2=jdof2;
    mpc1=0;
    mpc2=0;
    if(*nmpc>0){
      if(idof1!=2*(idof1/2)) mpc1=1;
      if(idof2!=2*(idof2/2)) mpc2=1;
    }
    if((mpc1==1)&&(mpc2==1)){
      id1=(-idof1+1)/2;
      id2=(-idof2+1)/2;
      if(id1==id2){
	      
	/* MPC id1 / MPC id1 */
	      
	ist=ipompc[id1-1];
	index1=nodempc[3*ist-1];
	if(index1==0) return;
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
	    idof2=nactdof[mt*(inode2-1)+nodempc[3*index2-2]];
	    if((idof1>0)&&(idof2>0)){
	      insert(ipointer,&mast1,&next,&idof1,&idof2,ifree,nzs_);
	      kdof1=idof1+neq[0];kdof2=idof2+neq[0];
	      insert(ipointer,&mast1,&next,&kdof1,&kdof2,ifree,nzs_);
	      if(((icomplex1!=0)||(icomplex2!=0))&&
		 (icomplex1!=icomplex2)){
		insert(ipointer,&mast1,&next,&kdof1,&idof2,ifree,nzs_);
		insert(ipointer,&mast1,&next,&idof1,&kdof2,ifree,nzs_);
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
	if(index1==0) return;
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
	  idof1=nactdof[mt*(inode1-1)+nodempc[3*index1-2]];
	  ist2=ipompc[id2-1];
	  index2=nodempc[3*ist2-1];
	  if(index2==0){
	    index1=nodempc[3*index1-1];
	    if(index1==0){break;}
	    else{return;}
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
	    idof2=nactdof[mt*(inode2-1)+nodempc[3*index2-2]];
	    if((idof1>0)&&(idof2>0)){
	      insert(ipointer,&mast1,&next,&idof1,&idof2,ifree,nzs_);
	      kdof1=idof1+neq[0];kdof2=idof2+neq[0];
	      insert(ipointer,&mast1,&next,&kdof1,&kdof2,ifree,nzs_);
	      if(((icomplex1!=0)||(icomplex2!=0))&&
		 (icomplex1!=icomplex2)){
		insert(ipointer,&mast1,&next,&kdof1,&idof2,ifree,nzs_);
		insert(ipointer,&mast1,&next,&idof1,&kdof2,ifree,nzs_);
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
  
  *mast1p=mast1;*nextp=next;
  
  return;
}
