/*     CalculiX - A 3-Dimensional finite element program                   */
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
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif

void expand(double *co, int *nk, int *kon, int *ipkon, char *lakon,
	     int *ne, int *nodeboun, int *ndirboun, double *xboun, int *nboun, 
	     int *ipompc, int *nodempc, double *coefmpc, char *labmpc,
             int *nmpc, int *nodeforc, int *ndirforc,double *xforc, 
             int *nforc, int *nelemload, char *sideload, double *xload,
             int *nload, int *nactdof, int *neq, 
	     int *nmethod, int *ikmpc, int *ilmpc, int *ikboun, int *ilboun,
	     double *elcon, int *nelcon, double *rhcon, int *nrhcon,
	     double *alcon, int *nalcon, double *alzero, int *ielmat,
	     int *ielorien, int *norien, double *orab, int *ntmat_,
	     double *t0,int *ithermal,double *prestr, int *iprestr, 
	     double *vold,int *iperturb, double *sti, int *nzs,  
	     double *adb, double *aub,char *filab, double *eme,
             double *plicon, int *nplicon, double *plkcon,int *nplkcon,
             double *xstate, int *npmat_, char *matname, int *mint_,
	     int *ics, double *cs, int *mpcend, int *ncmat_,
             int *nstate_, int *mcs, int *nkon, double *ener,
             char *jobnamec, char *output, char *set, int *nset,int *istartset,
             int *iendset, int *ialset, int *nprint, char *prlab,
             char *prset, int *nener, double *trab, 
             int *inotr, int *ntrans, double *ttime, double *fmpc,
	     int *nev, double *z, int *iamboun, double *xbounold,
             int *nsectors, int *nm,int *icol,int *irow,int *nzl, int *nam,
             int *ipompcold, int *nodempcold, double *coefmpcold,
             char *labmpcold, int *nmpcold, double *xloadold, int *iamload,
             double *t1old,double *t1,int *iamt1, double *xstiff){

  /* calls the Arnoldi Package (ARPACK) for cyclic symmetry calculations */
  
    char *filabt;
    int *inum=NULL,k,idir,lfin,j,lint,iout=0,index,inode,id,i,idof,
      ielas,icmd,kk,l,nkt,icntrl,imag=1,icomplex,kk4,kk5,kk6,iterm,
	lprev,ilength,ij,i1,i2,iel,ielset,node,indexe,nope,ml1,
	*inocs=NULL,*ielcs=NULL,jj,l1,l2,is,*ncocon=NULL,nlabel,
        new,nodeleft,noderight[3],numnodes,ileft,kflag=2,itr,locdir,
        neqh,j1,nodenew;
    double *stn=NULL,*v=NULL,*temp_array=NULL,*vini=NULL,
	*een=NULL,cam[3],*f=NULL,*fn=NULL,qa[2],*epn=NULL,*stiini=NULL,
	*xstateini=NULL,theta,pi,*coefmpcnew=NULL,t[3],ctl,stl,
	*stx=NULL,*enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*cocon=NULL,
	*qfx=NULL,*qfn=NULL,xreal,ximag,*vt=NULL,sum,*aux=NULL,coefright[3],
        coef,a[9],ratio;
    
    /* dummy arguments for the results call */
    
    double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

    pi=4.*atan(1.);
    
    v=NNEW(double,10**nk);
    vt=NNEW(double,5**nk**nsectors);
    
    fn=NNEW(double,8**nk);
    stn=NNEW(double,12**nk);
    inum=NNEW(int,*nk);
    stx=NNEW(double,6**mint_**ne);
    
    nlabel=15;
    filabt=NNEW(char,6*nlabel);
    for(i=1;i<6*nlabel;i++) filabt[i]=' ';
    filabt[0]='U';
    
    temp_array=NNEW(double,neq[1]);
    coefmpcnew=NNEW(double,*mpcend);
    
    nkt=*nsectors**nk;
    
    /* assigning nodes and elements to sectors */
    
    inocs=NNEW(int,*nk);
    ielcs=NNEW(int,*ne);
    ielset=cs[12];
    if((*mcs!=1)||(ielset!=0)){
	for(i=0;i<*nk;i++) inocs[i]=-1;
	for(i=0;i<*ne;i++) ielcs[i]=-1;
    }
    
    for(i=0;i<*mcs;i++){
	is=cs[17*i];
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
    
    /* generating the coordinates for the other sectors */
    
    icntrl=1;
    
    FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filabt,&imag));
    
    for(jj=0;jj<*mcs;jj++){
	is=cs[17*jj];
	for(i=1;i<is;i++){
	    
	    theta=i*2.*pi/cs[17*jj];
	    
	    for(l=0;l<*nk;l++){
		if(inocs[l]==jj){
		    co[3*l+i*3**nk]=co[3*l];
		    co[1+3*l+i*3**nk]=co[1+3*l]+theta;
		    co[2+3*l+i*3**nk]=co[2+3*l];
		    if(*ntrans>0) inotr[2*l+i*2**nk]=inotr[2*l];
		}
	    }
	    for(l=0;l<*nkon;l++){kon[l+i**nkon]=kon[l]+i**nk;}
	    for(l=0;l<*ne;l++){
		if(ielcs[l]==jj){
		    if(ipkon[l]>=0){
			ipkon[l+i**ne]=ipkon[l]+i**nkon;
			ielmat[l+i**ne]=ielmat[l];
			if(*norien>0) ielorien[l+i**ne]=ielorien[l];
			for(l1=0;l1<8;l1++){
			    l2=8*l+l1;
			    lakon[l2+i*8**ne]=lakon[l2];
			}
		    }
		}
	    }
	}
    }
    
    icntrl=-1;
    
    FORTRAN(rectcyl,(co,vt,fn,stn,qfn,een,cs,&nkt,&icntrl,t,filabt,&imag));
    
/* copying the boundary conditions
   (SPC's must be defined in cylindrical coordinates) */
    
    for(i=1;i<*nsectors;i++){
	for(j=0;j<*nboun;j++){
	    nodeboun[i**nboun+j]=nodeboun[j]+i**nk;
	    ndirboun[i**nboun+j]=ndirboun[j];
	    xboun[i**nboun+j]=xboun[j];
	    xbounold[i**nboun+j]=xbounold[j];
	    if(*nam>0) iamboun[i**nboun+j]=iamboun[j];
	    ikboun[i**nboun+j]=ikboun[j]+8*i**nk;
	    ilboun[i**nboun+j]=ilboun[j]+i**nboun;
	}
    }
    
/* distributed loads */
    
    for(i=0;i<*nload;i++){
	if(nelemload[2*i+1]<*nsectors){
	    nelemload[2*i]+=*ne*nelemload[2*i+1];
	}else{
	    nelemload[2*i]+=*ne*(nelemload[2*i+1]-(*nsectors));
	}
    }

  /*  sorting the elements with distributed loads */

    if(*nload>0){
	if(*nam>0){
	    FORTRAN(isortiddc2,(nelemload,iamload,xload,xloadold,sideload,nload,&kflag));
	}else{
	    FORTRAN(isortiddc1,(nelemload,xload,xloadold,sideload,nload,&kflag));
	}
    }
    
/* point loads */
    
    for(i=0;i<*nforc;i++){
	if(nodeforc[2*i+1]<*nsectors){
	    nodeforc[2*i]+=*nk*nodeforc[2*i+1];
	}else{
	    nodeforc[2*i]+=*nk*(nodeforc[2*i+1]-(*nsectors));
	}
    }

    neqh=neq[1]/2;

/* expand nactdof */

    for(i=1;i<*nsectors;i++){
	lint=i*4**nk;
	for(j=0;j<4**nk;j++){
	    if(nactdof[j]!=0){
		nactdof[lint+j]=nactdof[j]+i*neqh;
	    }else{
		nactdof[lint+j]=0;
	    }
	}
    }

    /*zcopy=NNEW(double,*nev*neq[1]);
      for(i=0;i<*nev*neq[1];i++) zcopy[i]=z[i];*/
    
    /* loop over all eigenvalues; the loop starts from the highest eigenvalue
       so that the reuse of z is not a problem
       z before: real and imaginary part for a segment for all eigenvalues
       z after: real part for all segments for all eigenvalues */

    lfin=0;
    for(j=*nev-1;j>-1;--j){
	lint=2*j*neqh;

	/* calculating the cosine and sine of the phase angle */

	for(jj=0;jj<*mcs;jj++){
	    theta=nm[j]*2.*pi/cs[17*jj];
	    cs[17*jj+14]=cos(theta);
	    cs[17*jj+15]=sin(theta);
	}
	
	/* generating the cyclic MPC's (needed for nodal diameters
	   different from 0 */
	
	eei=NNEW(double,6**mint_**ne);
	
	for(k=0;k<2*neqh;k+=neqh){
	    
	    for(i=0;i<6**mint_**ne;i++){eme[i]=0.;}
	    
	    if(k==0) {kk=0;kk4=0;kk5=0;kk6=0;}
	    else {kk=*nk;kk4=4**nk;kk5=5**nk;kk6=6**nk;}
	    for(i=0;i<*nmpc;i++){
		index=ipompc[i]-1;
		/* check whether thermal mpc */
		if(nodempc[3*index+1]==0) continue;
		coefmpcnew[index]=coefmpc[index];
		while(1){
		    index=nodempc[3*index+2];
		    if(index==0) break;
		    index--;
		    
		    icomplex=0;
		    inode=nodempc[3*index];
		    if(strcmp1(&labmpc[20*i],"CYCLIC")==0){
			icomplex=atoi(&labmpc[20*i+6]);}
		    else if(strcmp1(&labmpc[20*i],"SUBCYCLIC")==0){
			for(ij=0;ij<*mcs;ij++){
			    lprev=cs[ij*17+13];
			    ilength=cs[ij*17+3];
			    FORTRAN(nident,(&ics[lprev],&inode,&ilength,&id));
			    if(id!=0){
				if(ics[lprev+id-1]==inode){icomplex=ij+1;break;}
			    }
			}
		    }
		    
		    if(icomplex!=0){
			idir=nodempc[3*index+1];
			idof=nactdof[4*(inode-1)+idir]-1;
			if(idof==-1){xreal=1.;ximag=1.;}
			else{xreal=z[lint+idof];ximag=z[lint+idof+neqh];}
			if(k==0) {
			    if(fabs(xreal)<1.e-30)xreal=1.e-30;
			    coefmpcnew[index]=coefmpc[index]*
				(cs[17*(icomplex-1)+14]+
                                 ximag/xreal*cs[17*(icomplex-1)+15]);}
			else {
			    if(fabs(ximag)<1.e-30)ximag=1.e-30;
			    coefmpcnew[index]=coefmpc[index]*
				(cs[17*(icomplex-1)+14]-
                                 xreal/ximag*cs[17*(icomplex-1)+15]);}
		    }
		    else{coefmpcnew[index]=coefmpc[index];}
		}
	    }
	    
	    FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,&v[kk5],&stn[kk6],inum,
	      stx,elcon,
	      nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	      norien,orab,ntmat_,t0,t0,ithermal,
	      prestr,iprestr,filab,eme,&een[kk6],iperturb,
	      f,&fn[kk4],nactdof,&iout,qa,vold,&z[lint+k],
	      nodeboun,ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neqh,veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,&icmd,
	      ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],sti,
	      xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload));
	    
	}
	free(eei);
	
	sum=0.;

	/* mapping the results to the other sectors */
	
	icntrl=2;imag=1;
	
	FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filabt,&imag));
	
	/* basis sector */

	for(l=0;l<5**nk;l++){vt[l]=v[l];}

        for(i=0;i<neq[1];++i) temp_array[i]=0.;
	FORTRAN(op,(&neq[1],aux,&z[j*neq[1]],temp_array,adb,aub,icol,irow,&neq[1]));
	for(i=0;i<neq[1];++i){
	    sum+=z[j*neq[1]+i]*temp_array[i];
	}
	/*	printf("j=%d,i=0,sumnew=%e\n",j,sum);*/

	/* other sectors */

	for(jj=0;jj<*mcs;jj++){
	    ilength=cs[17*jj+3];
	    lprev=cs[17*jj+13];
	    for(i=1;i<*nsectors;i++){
		
		theta=i*nm[j]*2.*pi/cs[17*jj];
		ctl=cos(theta);
		stl=sin(theta);
		
		for(l1=0;l1<*nk;l1++){
		    if(inocs[l1]==jj){
			
			/* check whether node lies on axis */
			
			ml1=-l1-1;
			FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
			if(id!=0){
			    if(ics[lprev+id-1]==ml1){
				for(l2=0;l2<4;l2++){
				    l=5*l1+l2;
				    vt[l+5**nk*i]=v[l];
				}
				continue;
			    }
			}
			for(l2=0;l2<4;l2++){
			    l=5*l1+l2;
			    vt[l+5**nk*i]=ctl*v[l]-stl*v[l+5**nk];
			}
		    }
		}
	    }
	}
	
	icntrl=-2;imag=0;
	
	FORTRAN(rectcyl,(co,vt,fn,stn,qfn,een,cs,&nkt,&icntrl,t,filabt,
			 &imag));
	
/* storing the displacements into the expanded eigenvectors */
	
	for(i=0;i<*nk;i++){
	    for(j1=0;j1<4;j1++){
		if(nactdof[4*i+j1]!=0){
		    for(k=0;k<*nsectors;k++){
			z[j**nsectors*neqh+k*neqh+nactdof[4*i+j1]-1]=
			    vt[k*5**nk+5*i+j1];
		    }
		}
	    }	    
	}

/* normalizing the eigenvectors with the mass */

	if (nm[j]==0||(nm[j]==(int)((cs[0]/2))&&(fmod(cs[0],2.)==0.))){sum=sqrt(cs[0]);}
	else{sum=sqrt(cs[0]/2);}
	lint=j**nsectors*neqh;
	for(i=0;i<*nsectors*neqh;++i) z[lint+i]=z[lint+i]/sum;
    }

/* copying the multiple point constraints */

    *nmpc=0;
    *mpcend=0;
    for(i=0;i<*nsectors;i++){
	if(i==0){
	    ileft=*nsectors-1;
	}else{
	    ileft=i-1;
	}
	for(j=0;j<*nmpcold;j++){
	    ipompc[*nmpc]=*mpcend+1;
	    ikmpc[*nmpc]=ikmpc[j]+8*i**nk;
	    ilmpc[*nmpc]=ilmpc[j]+i**nmpcold;
	    strcpy1(&labmpc[20**nmpc],&labmpcold[20*j],20);
	    if(strcmp1(&labmpcold[20*j],"CYCLIC")==0){
		index=ipompcold[j]-1;
		nodeleft=nodempcold[3*index];
		idir=nodempcold[3*index+1];
		index=nodempcold[3*index+2]-1;
		numnodes=0;
		do{
		    node=nodempcold[3*index];
		    new=1;
		    for(k=0;k<numnodes;k++){
			if(node==noderight[k]){
			    new=0;
			    break;
			}
		    }
		    if(new==1){
// 			if(numnodes==3){
// 			    printf("*ERROR in expand: more than three independent nodes\n");
// 			    FORTRAN(stop,());
// 			}
			noderight[numnodes]=node;
			coefright[numnodes]=coefmpcold[index];
			numnodes++;
		    }
		    index=nodempcold[3*index+2]-1;
		    if(index==-1) break;
		}while(1);
		if(numnodes>0){
		    sum=0.;
		    for(k=0;k<numnodes;k++){
			sum+=coefright[k];
		    }
		    for(k=0;k<numnodes;k++){
			coefright[k]/=sum;
		    }
		}else{coefright[0]=1.;}
		nodempc[3**mpcend]=nodeleft+i**nk;
		nodempc[3**mpcend+1]=idir;
		nodempc[3**mpcend+2]=*mpcend+2;
		coefmpc[*mpcend]=1.;
		for(k=0;k<numnodes;k++){
		    (*mpcend)++;
		    nodempc[3**mpcend]=noderight[k]+ileft**nk;
		    nodempc[3**mpcend+1]=idir;
		    nodempc[3**mpcend+2]=*mpcend+2;
		    coefmpc[*mpcend]=-coefright[k];
		}
		nodempc[3**mpcend+2]=0;
		(*mpcend)++;
	    }else{
		index=ipompcold[j]-1;
		iterm=0;
		do{
		    iterm++;
		    node=nodempcold[3*index];
		    idir=nodempcold[3*index+1];
		    coef=coefmpcold[index];

		    /* check whether SUBCYCLIC MPC: if the current node
                       is an independent node of a CYCLIC MPC, the
                       node in the new MPC should be the cylic previous
                       one */

		    nodenew=node+i**nk;
		    if(strcmp1(&labmpcold[20*j],"SUBCYCLIC")==0){
			for(ij=0;ij<*mcs;ij++){
			    lprev=cs[ij*17+13];
			    ilength=cs[ij*17+3];
			    FORTRAN(nident,(&ics[lprev],&node,&ilength,&id));
			    if(id!=0){
				if(ics[lprev+id-1]==node){
				    nodenew=node+ileft**nk;
				    break;
				}
			    }
			}
		    }
		    
		    /* modification of the MPC coefficients if
                       cylindrical coordinate system is active
                       it is assumed that all terms in the MPC are
                       either in the radial, the circumferential
                       or axial direction  */
		    
		    if(*ntrans<=0){itr=0;}
		    else if(inotr[2*node-2]==0){itr=0;}
		    else{itr=inotr[2*node-2];}
		    
		    if(iterm==1) locdir=-1;
		    
		    if((itr!=0)&&(idir!=0)){
			if(trab[7*itr-1]<0){
			    FORTRAN(transformatrix,(&trab[7*itr-7],
						    &co[3*node-3],a));
			    if(iterm==1){
				for(k=0;k<3;k++){
				    if(fabs(a[3*k+idir-1]-coef)<1.e-10){
					FORTRAN(transformatrix,(&trab[7*itr-7],
								&co[3*nodenew-3],a));
					coef=a[3*k+idir-1];
					locdir=k;
					break;
				    }
				    if(fabs(a[3*k+idir-1]+coef)<1.e-10){
					FORTRAN(transformatrix,(&trab[7*itr-7],
								&co[3*nodenew-3],a));
					coef=-a[3*k+idir-1];
					locdir=k;
					break;
				    }
				}
			    }else{
				if(locdir!=-1){
				    if(fabs(a[3*locdir+idir-1])>1.e-10){
					ratio=coef/a[3*locdir+idir-1];
				    }else{ratio=0.;}
				    FORTRAN(transformatrix,(&trab[7*itr-7],
							    &co[3*nodenew-3],a));
				    coef=ratio*a[3*locdir+idir-1];
				}
			    }		
			}
		    }
		    
		    nodempc[3**mpcend]=nodenew;
		    nodempc[3**mpcend+1]=idir;
		    coefmpc[*mpcend]=coef;
		    index=nodempcold[3*index+2]-1;
		    if(index==-1) break;
		    nodempc[3**mpcend+2]=*mpcend+2;
		    (*mpcend)++;
		}while(1);
		nodempc[3**mpcend+2]=0;
		(*mpcend)++;
	    }
	    (*nmpc)++;
	}
    }
    
    /*      for(i=0;i<*nmpc;i++){
	j=i+1;
	FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
	}*/

    /* copying the temperatures */

    if(*ithermal!=0){
	for(i=1;i<*nsectors;i++){
	    lint=i**nk;
	    for(j=0;j<*nk;j++){
		t0[lint+j]=t0[j];
		t1old[lint+j]=t1old[j];
		t1[lint+j]=t1[j];
	    }
	}
	if(*nam>0){
	    for(i=1;i<*nsectors;i++){
		lint=i**nk;
		for(j=0;j<*nk;j++){
		    iamt1[lint+j]=iamt1[j];
		}
	    }
	}
    }
  
    *nk=nkt;
    (*ne)*=(*nsectors);
    (*nkon)*=(*nsectors);
    (*nboun)*=(*nsectors);
    (neq[1])*=(*nsectors/2);
    
    free(temp_array);free(coefmpcnew);
    free(v);free(vt);free(fn);free(stn);free(inum);free(stx);
    free(inocs);free(ielcs);free(filabt);
    
    return;
}
