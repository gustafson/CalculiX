/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

  /* calculating the initial acceleration at the start of the step
     for dynamic calculations */

  if((*nmethod==4)&&(*ithermal!=2)){
    bet=(1.-*alpha)*(1.-*alpha)/4.;
    gam=0.5-*alpha;

    /* calculating the stiffness and mass matrix 
       the stress must be determined to calculate the 
       stiffness matrix*/

    reltime=0.;
    time=0.;
    dtime=0.;

    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload));
    
    time=0.;
    dtime=1.;
    
    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's)
        if contact arises the number of MPC's can also change */

    cam[0]=0.;cam[1]=0.;
    if(*ithermal>1){radflowload(itg,ieg,&ntg,&ntr,&ntm,
       ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,vold,
       shcon,nshcon,ipkon,kon,co,pmid,e1,e2,e3,iptri,
       kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,fij,
       dist,idist,area,nflow,ikboun,xboun,nboun,ithermal,&iinc,&iit,
       cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
       ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
       nactdog,nacteq,nodeboun,ndirboun,&network,rhcon,
       nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,ctrl,
       xloadold,&reltime,nmethod,set);}

    if((icascade==2)||(ncont!=0)){
	/**nmpc=nmpcref;*/
	memmpc_=memmpcref_;mpcfree=mpcfreeref;
	RENEW(nodempc,int,3*memmpcref_);
	for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
	RENEW(coefmpc,double,memmpcref_);
	for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
    }

    if(ncont!=0){
	*ne=ne0;*nkon=nkon0;
	contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
             ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	     co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
	     ifcont1,ifcont2,&ne0,vini,nmethod,nmpc,&mpcfree,&memmpc_,
	     &ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
             iperturb,ikboun,nboun);
	icascade=1;
    }
    newstep=0;
    FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		       nmpc,ikboun,ilboun,nboun,xbounold,aux,iaux,
		       &maxlenmpc,ikmpc,ilmpc,&icascade,
		       kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
                       &iit,&idiscon,&ncont,trab,ntrans,ithermal));
    if((icascade==2)||(ncont!=0)){
	/*nmpcref=*nmpc;*/
	memmpcref_=memmpc_;mpcfreeref=mpcfree;
	RENEW(nodempcref,int,3*memmpc_);
	for(k=0;k<3*memmpc_;k++){nodempcref[k]=nodempc[k];}
	RENEW(coefmpcref,double,memmpc_);
	for(k=0;k<memmpc_;k++){coefmpcref[k]=coefmpc[k];}
    }

    if(icascade>0) remastruct(ipompc,&coefmpc,&nodempc,nmpc,
              &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
              labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
              kon,ipkon,lakon,ne,nnn,nactdof,icol,jq,&irow,isolver,
	      neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	      &adb,&aub,ithermal,iperturb,mass);
    
    iout=-1;
    ielas=1;

    fn=NNEW(double,4**nk);
    stx=NNEW(double,6**mint_**ne);
    if(*ithermal>1) qfx=NNEW(double,3**mint_**ne);

    FORTRAN(results,(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	  prestr,iprestr,filab,eme,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounold,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,&bet,
          &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mint_,&ielas,&icmd,
          ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,sti,xstaten,
          eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
          ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
          nelemload,nload,ikmpc,ilmpc,istep,&iinc));
    
    free(fn);free(stx);if(*ithermal>1)free(qfx);
    
    iout=0;
    ielas=0;

    reltime=0.;
    time=0.;
    dtime=0.;

    if(*iexpl<=1){intscheme=1;}

    /* in mafillsm the stiffness and mass matrix are computed;
       The primary aim is to calculate the mass matrix (not 
       lumped for an implicit dynamic calculation, lumped for an
       explicit dynamic calculation). However:
       - for an implicit calculation the mass matrix is "doped" with
       a small amount of stiffness matrix, therefore the calculation
       of the stiffness matrix is needed.
       - for an explicit calculation the stiffness matrix is not 
       needed at all. Since the calculation of the mass matrix alone
       is not possible in mafillsm, the determination of the stiffness
       matrix is taken as unavoidable "ballast". */

    ad=NNEW(double,neq[1]);
    au=NNEW(double,nzs[1]);

    FORTRAN(mafillsm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounold,nboun,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	      nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	      nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
	      nmethod,ikmpc,ilmpc,ikboun,ilboun,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	      ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	      xstiff,npmat_,&dtime,matname,mint_,
              ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
	      physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		      &coriolis,ibody,xloadold,&reltime,veold));

    if(nmethod==0){
      
      /* error occurred in mafill: storing the geometry in frd format */

      ++*kode;
      ipneigh=NNEW(int,*nk);neigh=NNEW(int,40**ne);
      FORTRAN(out,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,
           een,t1,fn,ttime,epn,ielmat,matname,enern,xstaten,nstate_,istep,
           &iinc,iperturb,ener,mint_,output,ithermal,qfn,&mode,&noddiam,
           trab,inotr,ntrans,orab,ielorien,norien,description,
	   ipneigh,neigh,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ne,cs));
      free(ipneigh);free(neigh);
      
      FORTRAN(stop,());
      
    }

    /* calculating the acceleration at the start of the step.
       This can be different from the acceleration at the end
       of the last step due to a discontinuous loading increase */

    /*   reltime=0.;
    time=0.;
    dtime=0.;

    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload));*/

    /* determining the external loading vector */

    /*   FORTRAN(rhs,(co,nk,kon,ipkon,lakon,ne,
	 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	 nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	 nbody,cgr,fext,nactdof,&neq[1],
	 nmethod,ikmpc,ilmpc,
	 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	 ielmat,ielorien,norien,orab,ntmat_,
	 t0,t1act,ithermal,iprestr,vold,iperturb,
	 iexpl,plicon,nplicon,plkcon,nplkcon,
	 npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody));*/

    /* mass x acceleration = f(external)-f(internal) 
       only for the mechanical loading*/

    for(k=0;k<neq[0];++k){
      b[k]=fext[k]-f[k];
    }

    if(*iexpl<=1){

      /* a small amount of stiffness is added to the mass matrix
       otherwise the system leads to huge accelerations in 
       case of discontinuous load changes at the start of the step */

	dtime=*tinc/10.;
/*	dtime=*tinc/5.; */
	scal1=bet*dtime*dtime*(1.+*alpha);
	for(k=0;k<neq[0];++k){
	    /*ad[k]=adb[k];*/
	    ad[k]=adb[k]+scal1*ad[k];
	}
	for(k=0;k<nzs[0];++k){
	    /*au[k]=aub[k];*/
	    au[k]=aub[k]+scal1*au[k];
	}
        if(*isolver==0){
#ifdef SPOOLES
          spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
                  &symmetryflag,&inputformat);
#else
          printf("*ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
          FORTRAN(stop,());
#endif
        }
        else if((*isolver==2)||(*isolver==3)){
          preiter(ad,&au,b,&icol,&irow,&neq[0],&nzs[0],isolver,iperturb);
        }
        else if(*isolver==4){
#ifdef SGI
          token=1;
          sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
#else
          printf("*ERROR in nonlingeo: the SGI library is not linked\n\n");
          FORTRAN(stop,());
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[0],&nzs[0]);
#else
          printf("*ERROR in nonlingeo: the TAUCS library is not linked\n\n");
          FORTRAN(stop,());
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0]);
#else
          printf("*ERROR in nonlingeo: the PARDISO library is not linked\n\n");
          FORTRAN(stop,());
#endif
        }
    }

    else{
	for(k=0;k<neq[0];++k){
	    b[k]=(fext[k]-f[k])/adb[k];
	}
    }
    
    /* for thermal loading the acceleration is set to zero */

    for(k=neq[0];k<neq[1];++k){
	b[k]=0.;
    }

    /* storing the acceleration accold */
    
    for(k=0;k<4**nk;++k){
	if(nactdof[k]!=0){accold[k]=b[nactdof[k]-1];}
    }

    free(ad);free(au);

    /* the mass matrix is kept for subsequent calculations, therefore,
       no new mass calculation is necessary for the remaining iterations
       in the present step */

    mass[0]=0;intscheme=0;

  }
