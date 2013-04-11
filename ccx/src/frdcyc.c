/*     CalculiX - A 3-dimensional finite element program                   */
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
#include <string.h>
#include "CalculiX.h"

void frdcyc(double *co,int *nk,int *kon,int *ipkon,char *lakon,int *ne,double *v,
	    double *stn,int *inum,int *nmethod,int *kode,char *filab,
	    double *een,double *t1,double *fn,double *time,double *epn,
	    int *ielmat,char *matname, double *cs, int *mcs, int *nkon,
            double *enern, double *xstaten, int *nstate_, int *istep,
            int *iinc, int *iperturb, double *ener, int *mint_, char *output,
            int *ithermal, double *qfn, int *ialset, int *istartset,
            int *iendset, double *trab, int *inotr, int *ntrans,
	    double *orab, int *ielorien, int *norien, double *sti,double *veold){

  /* duplicates fields for static cyclic symmetric calculations */

  char *lakont=NULL,description[13]="            ";
  int nkt,icntrl,*kont=NULL,*ipkont=NULL,*inumt=NULL,*ielmatt=NULL,net,i,l,
     imag=0,mode=-1,noddiam=0,ngraph,*inocs=NULL,*ielcs=NULL,l1,l2,is,
      jj,node,i1,i2,nope,iel,indexe,j,ielset,*inotrt=NULL;
  double *vt=NULL,*fnt=NULL,*stnt=NULL,*eent=NULL,*cot=NULL,*t1t=NULL,
         *epnt=NULL,*enernt=NULL,*xstatent=NULL,theta,pi,t[3],*qfnt=NULL,
         *vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL;

  int *ipneigh=NULL,*neigh=NULL;

  pi=4.*atan(1.);

  /* determining the maximum number of sectors to be plotted */

  ngraph=1;
  for(j=0;j<*mcs;j++){
    if(cs[17*j+4]>ngraph) ngraph=cs[17*j+4];
  }

  /* assigning nodes and elements to sectors */

  inocs=NNEW(int,*nk);
  ielcs=NNEW(int,*ne);
  ielset=cs[12];
  if((*mcs!=1)||(ielset!=0)){
    for(i=0;i<*nk;i++) inocs[i]=-1;
    for(i=0;i<*ne;i++) ielcs[i]=-1;
  }

  for(i=0;i<*mcs;i++){
    is=cs[17*i+4];
    if(is==1) continue;
    ielset=cs[17*i+12];
    if(ielset==0) continue;
    for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
      if(ialset[i1]>0){
        iel=ialset[i1]-1;
        if(ipkon[iel]<0) continue;
        ielcs[iel]=i;
        indexe=ipkon[iel];
        if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
        else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
        else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
        else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
        else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
        else {nope=6;}
        for(i2=0;i2<nope;++i2){
          node=kon[indexe+i2]-1;
          inocs[node]=i;
        }
      }
      else{
        iel=ialset[i1-2]-1;
        do{
          iel=iel-ialset[i1];
          if(iel>=ialset[i1-1]-1) break;
          if(ipkon[iel]<0) continue;
          ielcs[iel]=i;
          indexe=ipkon[iel];
          if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
          else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
          else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
          else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
          else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
          else {nope=6;}
          for(i2=0;i2<nope;++i2){
            node=kon[indexe+i2]-1;
            inocs[node]=i;
          }
        }while(1);
      }
    } 
  }

  cot=NNEW(double,3**nk*ngraph);
  if(*ntrans>0)inotrt=NNEW(int,2**nk*ngraph);

  if((strcmp1(&filab[0],"U   ")==0)||
     ((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal>=2)))
    vt=NNEW(double,5**nk*ngraph);
  if((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal<2))
    t1t=NNEW(double,*nk*ngraph);
  if(strcmp1(&filab[12],"S   ")==0)
    stnt=NNEW(double,6**nk*ngraph);
  if(strcmp1(&filab[18],"E   ")==0)
    eent=NNEW(double,6**nk*ngraph);
  if((strcmp1(&filab[24],"RF  ")==0)||(strcmp1(&filab[54],"RFL ")==0))
    fnt=NNEW(double,4**nk*ngraph);
  if(strcmp1(&filab[30],"PEEQ")==0)
    epnt=NNEW(double,*nk*ngraph);
  if(strcmp1(&filab[36],"ENER")==0)
    enernt=NNEW(double,*nk*ngraph);
  if(strcmp1(&filab[42],"SDV ")==0)
    xstatent=NNEW(double,*nstate_**nk*ngraph);
  if(strcmp1(&filab[48],"HFL ")==0)
    qfnt=NNEW(double,3**nk*ngraph);

  /* the topology only needs duplication the first time it is
     stored in the frd file (*kode=1) */

  if(*kode==1){
    kont=NNEW(int,*nkon*ngraph);
    ipkont=NNEW(int,*ne*ngraph);
    lakont=NNEW(char,8**ne*ngraph);
    ielmatt=NNEW(int,*ne*ngraph);
  }
  inumt=NNEW(int,*nk*ngraph);
  
  nkt=ngraph**nk;
  net=ngraph**ne;

  /* copying the coordinates of the first sector */
  
  for(l=0;l<3**nk;l++){cot[l]=co[l];}
  if(*ntrans>0){for(l=0;l<*nk;l++){inotrt[2*l]=inotr[2*l];}}

  /* copying the topology of the first sector */
  
  if(*kode==1){
      for(l=0;l<*nkon;l++){kont[l]=kon[l];}
      for(l=0;l<*ne;l++){ipkont[l]=ipkon[l];}
      for(l=0;l<8**ne;l++){lakont[l]=lakon[l];}
      for(l=0;l<*ne;l++){ielmatt[l]=ielmat[l];}
  }  

  /* generating the coordinates for the other sectors */
  
  icntrl=1;
  
  FORTRAN(rectcyl,(cot,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag));
  
  for(jj=0;jj<*mcs;jj++){
    is=cs[17*jj+4];
    for(i=1;i<is;i++){
      
      theta=i*2.*pi/cs[17*jj];
      
      for(l=0;l<*nk;l++){
        if(inocs[l]==jj){
	  cot[3*l+i*3**nk]=cot[3*l];
	  cot[1+3*l+i*3**nk]=cot[1+3*l]+theta;
	  cot[2+3*l+i*3**nk]=cot[2+3*l];
        }
      }
      
      if(*ntrans>0){
	  for(l=0;l<*nk;l++){
	      if(inocs[l]==jj){
		  inotrt[2*l+i*2**nk]=inotrt[2*l];
	      }
	  }
      }
      
      if(*kode==1){
        
        for(l=0;l<*nkon;l++){kont[l+i**nkon]=kon[l]+i**nk;}
        for(l=0;l<*ne;l++){
          if(ielcs[l]==jj){
            if(ipkon[l]>=0){
              ipkont[l+i**ne]=ipkon[l]+i**nkon;
              ielmatt[l+i**ne]=ielmat[l];
              for(l1=0;l1<8;l1++){
                l2=8*l+l1;
                lakont[l2+i*8**ne]=lakon[l2];
              }
            }
            else ipkont[l+i**ne]=-1;
	  }
        }
      }
    }
  }

  icntrl=-1;
    
  FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
      &imag));
  
  /* mapping the results to the other sectors */
  
  for(l=0;l<*nk;l++){inumt[l]=inum[l];}
  
  icntrl=2;
  
  FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag));
  
  if((strcmp1(&filab[0],"U   ")==0)||
     ((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal>=2)))
    for(l=0;l<5**nk;l++){vt[l]=v[l];};
  if((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal<2))
    for(l=0;l<*nk;l++){t1t[l]=t1[l];};
  if(strcmp1(&filab[12],"S   ")==0)
    for(l=0;l<6**nk;l++){stnt[l]=stn[l];};
  if(strcmp1(&filab[18],"E   ")==0)
    for(l=0;l<6**nk;l++){eent[l]=een[l];};
  if((strcmp1(&filab[24],"RF  ")==0)||(strcmp1(&filab[54],"RFL ")==0))
    for(l=0;l<4**nk;l++){fnt[l]=fn[l];};
  if(strcmp1(&filab[30],"PEEQ")==0)
    for(l=0;l<*nk;l++){epnt[l]=epn[l];};
  if(strcmp1(&filab[36],"ENER")==0)
    for(l=0;l<*nk;l++){enernt[l]=enern[l];};
  if(strcmp1(&filab[42],"SDV ")==0)
    for(l=0;l<*nstate_**nk;l++){xstatent[l]=xstaten[l];};
  if(strcmp1(&filab[48],"HFL ")==0)
    for(l=0;l<3**nk;l++){qfnt[l]=qfn[l];};
  
  for(jj=0;jj<*mcs;jj++){
    is=cs[17*jj+4];
    for(i=1;i<is;i++){
    
      for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
    
      if((strcmp1(&filab[0],"U   ")==0)||
         ((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal>=2))){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<4;l2++){
              l=5*l1+l2;
              vt[l+5**nk*i]=v[l];
            }
          }
        }
      }
    
      if((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal<2)){
        for(l=0;l<*nk;l++){
          if(inocs[l]==jj) t1t[l+*nk*i]=t1[l];
        }
      }
    
      if(strcmp1(&filab[12],"S   ")==0){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<6;l2++){
              l=6*l1+l2;
              stnt[l+6**nk*i]=stn[l];
            }
          }
        }
      }
    
      if(strcmp1(&filab[18],"E   ")==0){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<6;l2++){
              l=6*l1+l2;
              eent[l+6**nk*i]=een[l];
            }
          }
        }
      }
    
      if((strcmp1(&filab[24],"RF  ")==0)||(strcmp1(&filab[54],"RFL ")==0)){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<4;l2++){
              l=4*l1+l2;
              fnt[l+4**nk*i]=fn[l];
            }
          }
        }
      }
    
      if(strcmp1(&filab[30],"PEEQ")==0){
        for(l=0;l<*nk;l++){
          if(inocs[l]==jj) epnt[l+*nk*i]=epn[l];
        }
      } 
    
      if(strcmp1(&filab[36],"ENER")==0){
        for(l=0;l<*nk;l++){
          if(inocs[l]==jj) enernt[l+*nk*i]=enern[l];
        }
      } 
    
      if(strcmp1(&filab[42],"SDV ")==0){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<*nstate_;l2++){
              l=*nstate_*l1+l2;
              xstatent[l+*nstate_**nk*i]=xstaten[l];
            }
          } 
        }
      }
    
      if(strcmp1(&filab[48],"HFL ")==0){
        for(l1=0;l1<*nk;l1++){
          if(inocs[l1]==jj){
            for(l2=0;l2<3;l2++){
              l=3*l1+l2;
              qfnt[l+3**nk*i]=qfn[l];
            }
          }
        }
      }
    }
  }
  
  icntrl=-2;
  
  FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
     &imag));
  
    ipneigh=NNEW(int,nkt);neigh=NNEW(int,40*net);
  FORTRAN(out,(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,kode,
	       filab,eent,t1t,fnt,time,epnt,ielmatt,matname,enernt,
               xstatent,nstate_,istep,iinc,iperturb,ener,mint_,output,
               ithermal,qfnt,&mode,&noddiam,trab,inotrt,ntrans,orab,ielorien,
               norien,description,ipneigh,neigh,sti,vr,vi,stnr,stni,
               vmax,stnmax,&ngraph,veold,&net,cs));
  free(ipneigh);free(neigh);
  
  if((strcmp1(&filab[0],"U   ")==0)||
     ((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal>=2))) free(vt);
  if((strcmp1(&filab[6],"NT  ")==0)&&(*ithermal<2)) free(t1t);
  if(strcmp1(&filab[12],"S   ")==0) free(stnt);
  if(strcmp1(&filab[18],"E   ")==0) free(eent);
  if((strcmp1(&filab[24],"RF  ")==0)||(strcmp1(&filab[54],"RFL ")==0))
        free(fnt);
  if(strcmp1(&filab[30],"PEEQ")==0) free(epnt);
  if(strcmp1(&filab[36],"ENER")==0) free(enernt);
  if(strcmp1(&filab[42],"SDV ")==0) free(xstatent);
  if(strcmp1(&filab[48],"HFL ")==0) free(qfnt);

  if(*kode==1){
    free(kont);free(ipkont);free(lakont);free(ielmatt);
  }
  free(inumt);free(cot);if(*ntrans>0)free(inotrt);
  return;
}

