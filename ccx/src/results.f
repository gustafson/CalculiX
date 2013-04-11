!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,filab,eme,een,iperturb,f,fn,
     &  nactdof,iout,qa,vold,b,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
     &  veold,accold,bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,
     &  nplkcon,
     &  xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,sti,
     &  xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
     &  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
     &  nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea)
!
!     calculates and prints the displacements, temperatures and forces 
!     at the nodes and the stress and  strain  at the reduced integration 
!     points and at the nodes
!
!     iout=-2: v is assumed to be known and is used to
!              calculate strains, stresses..., no result output
!              corresponds to iout=-1 with in addition the
!              calculation of the internal energy density
!     iout=-1: v is assumed to be known and is used to
!              calculate strains, stresses..., no result output;
!              is used to take changes in SPC's and MPC's at the
!              start of a new increment or iteration into account
!     iout=0: v is calculated from the system solution
!             and strains, stresses.. are calculated, no result output
!     iout=1:  v is calculated from the system solution and strains,
!              stresses.. are calculated, requested results output
!     iout=2: v is assumed to be known and is used to 
!             calculate strains, stresses..., requested results output
!
      implicit none
!
      logical calcul_fn,calcul_f,calcul_cauchy,calcul_qa,cauchy,
     &  fluid,force,intpointvar
!
      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*),lakonl
      character*20 labmpc(*)
      character*80 amat,matname(*)
      character*81 set(*),prset(*)
      character*87 filab(*)
!
      integer kon(*),konl(20),inum(*),iperm(20),ikmpc(*),ilmpc(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(*),ielorien(*),
     &  ntmat_,ipkon(*),mi(2),
     &  nactdof(0:mi(2),*),nodeboun(*),nelemload(2,*),
     &  ndirboun(*),ipompc(*),nodempc(3,*),ikboun(*),ilboun(*),
     &  ncocon(2,*),inotr(2,*),iorienglob,iflag,nload,nshcon,
     &  istep,iinc,mt,nk,ne,mattyp,ithermal(2),iprestr,i,j,k,m1,m2,jj,
     &  i1,m3,m4,kk,iener,indexe,nope,norien,iperturb(*),iout,
     &  nal,icmd,ihyper,nboun,nmpc,nmethod,ist,ndir,node,index,
     &  neq,kode,imat,mint3d,nfield,ndim,iorien,ielas,
     &  istiff,ncmat_,nstate_,incrementalmpc,jmin,jmax,
     &  nset,istartset(*),iendset(*),ialset(*),nprint,ntrans,ikin,
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,20),stiini(6,mi(1),*),
     &  stx(6,mi(1),*),stn(6,*),xl(3,20),vl(0:mi(2),20),stre(6),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),vini(0:mi(2),*),qfx(3,mi(1),*),qfn(3,*),
     &  alzero(*),orab(7,*),elas(21),rho,f(*),fn(0:mi(2),*),fnl(3,9),
     &  skl(3,3),beta(6),q(0:mi(2),20),vkl(0:3,3),cam(5),
     &  t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),een(6,*),ckl(3,3),
     &  vold(0:mi(2),*),b(*),xboun(*),coefmpc(*),eloc(9),
     &  veold(0:mi(2),*),springarea(*),
     &  accold(0:mi(2),*),elconloc(21),eth(6),xkl(3,3),
     &  voldl(0:mi(2),20),epn(*),
     &  xikl(3,3),ener(mi(1),*),enern(*),sti(6,mi(1),*),emec(6),
     &  eei(6,mi(1),*),enerini(mi(1),*),cocon(0:6,ntmat_,*),emec0(6),
     &  fmpc(*),shcon,sph,c1,vel(1:3,20),veoldl(0:mi(2),20),
     &  e,un,al,um,am1,xi,et,ze,tt,exx,eyy,ezz,exy,exz,eyz,
     &  xsj,qa(3),vj,t0l,t1l,bet,gam,dtime,forcempc,scal1,scal2,bnac,
     &  fixed_disp,weight,pgauss(3),vij,coconloc(6),qflux(3),time,ttime,
     &  t1lold,plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(82),
     &  vokl(3,3),xstateini(nstate_,mi(1),*),vikl(3,3),trab(7,*),
     &  xstaten(nstate_,*)
!
      include "gauss.f"
!
      data iflag /3/
      data iperm /5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/
!
      mt=mi(2)+1
      fluid=.false.
      intpointvar=.true.
!
      if(ithermal(1).le.1) then
         jmin=1
         jmax=3
      elseif(ithermal(1).eq.2) then
         jmin=0
         jmax=min(mi(2),2)
      else
         jmin=0
         jmax=3
      endif
!
      if((iout.ne.2).and.(iout.gt.-1)) then
!
         if((nmethod.ne.4).or.(iperturb(1).le.1)) then
            if(ithermal(1).ne.2) then
               do i=1,nk
                  do j=1,3
                     if(nactdof(j,i).ne.0) then
                        bnac=b(nactdof(j,i))
                     else
                        bnac=0.d0
                     endif
                     v(j,i)=v(j,i)+bnac
                     if((iperturb(1).ne.0).and.(nmethod.eq.1)) then
                        if(dabs(bnac).gt.cam(1)) then
                           cam(1)=dabs(bnac)
c                           cam(4)=i+0.5d0
                           cam(4)=nactdof(j,i)-0.5d0
                        endif
                     endif
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
                  if(nactdof(0,i).ne.0) then
                     bnac=b(nactdof(0,i))
                  else
                     bnac=0.d0
                  endif
                  v(0,i)=v(0,i)+bnac
                  if((iperturb(1).ne.0).and.(nmethod.eq.1)) then
                     if(dabs(bnac).gt.cam(2)) then
                        cam(2)=dabs(bnac)
c                        cam(5)=i+0.5d0
                        cam(5)=nactdof(0,i)-0.5d0
                     endif
                  endif
               enddo
            endif
c!
c!     extracting the displacement information from the solution
c!
c            do i=1,nk
c               do j=jmin,jmax
c                  if(nactdof(j,i).ne.0) then
c                     v(j,i)=b(nactdof(j,i))
c                  else
c                     v(j,i)=0.d0
c                  endif
c               enddo
c            enddo
c!     
c!     for static perturbation steps v represents the incremental
c!     displacements. For the total displacement vold must be added.
c!
c            if((iperturb(1).ne.0).and.(nmethod.eq.1)) then
c               if(ithermal(1).ne.2) then
c                  do i=1,nk
c                     do j=1,3
c                        if(dabs(v(j,i)).gt.cam(1)) then
c                           cam(1)=dabs(v(j,i))
c                           cam(4)=i+0.5d0
c                        endif
c                        v(j,i)=v(j,i)+vold(j,i)
c                     enddo
c                  enddo
c               endif
c               if(ithermal(1).gt.1) then
c                  do i=1,nk
c                     if(dabs(v(0,i)).gt.cam(2)) then
c                        cam(2)=dabs(v(0,i))
c                        cam(5)=i+0.5d0
c                     endif
c                     v(0,i)=v(0,i)+vold(0,i)
c                  enddo
c               endif
c!
c!              copy pressure and mass flow values
c!
c               if(ithermal(1).eq.2) then
c                  do i=1,nk
c                     do j=1,min(2,mi(2))
c                        v(j,i)=vold(j,i)
c                     enddo
c                  enddo
c               endif
c            endif
!
         else
!
!           direct integration dynamic step
!           b contains the acceleration increment
!
            if(ithermal(1).ne.2) then
               scal1=bet*dtime*dtime
               scal2=gam*dtime
               do i=1,nk
                  do j=1,3
                     if(nactdof(j,i).ne.0) then
                        bnac=b(nactdof(j,i))
                     else
                        bnac=0.d0
                     endif
                     v(j,i)=v(j,i)+scal1*bnac
                     if(dabs(scal1*bnac).gt.cam(1)) then
                        cam(1)=dabs(scal1*bnac)
c                        cam(4)=i+0.5d0
                        cam(4)=nactdof(j,i)-0.5d0
                     endif
                     veold(j,i)=veold(j,i)+scal2*bnac
                     accold(j,i)=accold(j,i)+bnac
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
                  if(nactdof(0,i).ne.0) then
                     bnac=b(nactdof(0,i))
                  else
                     bnac=0.d0
                  endif
                  v(0,i)=v(0,i)+bnac
                  if(dabs(bnac).gt.cam(2)) then
                     cam(2)=dabs(bnac)
c                     cam(5)=i+0.5d0
                     cam(5)=nactdof(0,i)-0.5d0
                  endif
                  if(nactdof(0,i).ne.0) then
                     cam(3)=max(cam(3),dabs(v(0,i)-vini(0,i)))
                  endif
                  veold(0,i)=0.d0
               enddo
            endif
c!
c!           copy pressure and mass flow values
c!
c            if(ithermal(1).eq.2) then
c               do i=1,nk
c                  do j=1,min(mi(2),2)
c                     v(j,i)=vold(j,i)
c                  enddo
c               enddo
c            endif
         endif
!
      endif
!
!        initialization
!
      calcul_fn=.false.
      calcul_f=.false.
      calcul_qa=.false.
      calcul_cauchy=.false.
!
!     determining which quantities have to be calculated
!
      if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.0))) 
     &      then
         if((iout.lt.1).and.(iout.gt.-2)) then
            calcul_fn=.true.
            calcul_f=.true.
            calcul_qa=.true.
         elseif((iout.ne.-2).and.(iperturb(2).eq.1)) then    
            calcul_cauchy=.true.
         endif
      endif
!
      if(iout.gt.0) then
         if((filab(5)(1:4).eq.'RF  ').or.
     &      (filab(10)(1:4).eq.'RFL ')) then
            calcul_fn=.true.
         else
            do i=1,nprint
               if((prlab(i)(1:4).eq.'RF  ').or.
     &            (prlab(i)(1:4).eq.'RFL ')) then
                  calcul_fn=.true.
                  exit
               endif
            enddo
         endif
      endif
!
!     initializing fn
!
      if(calcul_fn) then
         do i=1,nk
            do j=0,mi(2)
               fn(j,i)=0.d0
            enddo
         enddo
      endif
!
!     initializing f
!
      if(calcul_f) then
         do i=1,neq
            f(i)=0.d0
         enddo
      endif
!
!     SPC's and MPC's have to be taken into account for 
!     iout=0,1 and -1
!
      if(abs(iout).lt.2) then
!
!     inserting the boundary conditions
!
      do i=1,nboun
         if(ndirboun(i).gt.3) cycle
         fixed_disp=xboun(i)
         if((nmethod.eq.4).and.(iperturb(1).gt.1)) then
            ndir=ndirboun(i)
            node=nodeboun(i)
            if(ndir.gt.0) then
               accold(ndir,node)=(xboun(i)-v(ndir,node))/
     &              (bet*dtime*dtime)
               veold(ndir,node)=veold(ndir,node)+
     &              gam*dtime*accold(ndir,node)
            else
               veold(ndir,node)=(xboun(i)-v(ndir,node))/dtime
            endif
         endif
         v(ndirboun(i),nodeboun(i))=fixed_disp
      enddo
!
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!
c
c      to be checked: should replace the lines underneath do i=1,nmpc
c
c      incrementalmpc=iperturb(2)
      do i=1,nmpc
         if((labmpc(i)(1:20).eq.'                    ').or.
     &      (labmpc(i)(1:7).eq.'CONTACT').or.
     &      (labmpc(i)(1:6).eq.'CYCLIC').or.
     &      (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
            incrementalmpc=0
         else
            if((nmethod.eq.2).or.(nmethod.eq.3).or.
     &            ((iperturb(1).eq.0).and.(nmethod.eq.1)))
     &         then
               incrementalmpc=0
            else
               incrementalmpc=1
            endif
         endif
         ist=ipompc(i)
         node=nodempc(1,ist)
         ndir=nodempc(2,ist)
         if(ndir.eq.0) then
            if(ithermal(1).lt.2) cycle
         elseif(ndir.gt.3) then
            cycle
         else
            if(ithermal(1).eq.2) cycle
         endif
         index=nodempc(3,ist)
         fixed_disp=0.d0
         if(index.ne.0) then
            do
               if(incrementalmpc.eq.0) then
                  fixed_disp=fixed_disp-coefmpc(index)*
     &                 v(nodempc(2,index),nodempc(1,index))
               else
                  fixed_disp=fixed_disp-coefmpc(index)*
     &                 (v(nodempc(2,index),nodempc(1,index))-
     &                  vold(nodempc(2,index),nodempc(1,index)))
               endif
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
         fixed_disp=fixed_disp/coefmpc(ist)
         if(incrementalmpc.eq.1) then
            fixed_disp=fixed_disp+vold(ndir,node)
         endif
         if((nmethod.eq.4).and.(iperturb(1).gt.1)) then
            if(ndir.gt.0) then
               accold(ndir,node)=(fixed_disp-v(ndir,node))/
     &              (bet*dtime*dtime)
               veold(ndir,node)=veold(ndir,node)+
     &              gam*dtime*accold(ndir,node)
            else
               veold(ndir,node)=(fixed_disp-v(ndir,node))/dtime
            endif
         endif
         v(ndir,node)=fixed_disp
      enddo
      endif
!
!     check whether there are any strain output requests
!
      iener=0
      ikin=0
      if((filab(7)(1:4).eq.'ENER').or.(filab(27)(1:4).eq.'CELS')) then
         iener=1
      endif

      do i=1,nprint
         if((prlab(i)(1:4).eq.'ENER').or.(prlab(i)(1:4).eq.'ELSE').or.
     &        (prlab(i)(1:4).eq.'CELS')) then
            iener=1
         elseif(prlab(i)(1:4).eq.'ELKE') then
            ikin=1
         endif
      enddo
!     
      qa(1)=0.d0
      nal=0
!
!     check whether integration point variables are needed in
!     modal dynamics calculations
!
      if((nmethod.eq.4).and.(iperturb(1).lt.2)) then
         intpointvar=.false.
         if((filab(3)(1:4).eq.'S   ').or.
     &      (filab(4)(1:4).eq.'E   ').or.
     &      (filab(5)(1:4).eq.'RF  ').or.
     &      (filab(6)(1:4).eq.'PEEQ').or.
     &      (filab(7)(1:4).eq.'ENER').or.
     &      (filab(8)(1:4).eq.'SDV ').or.
     &      (filab(13)(1:4).eq.'ZZS ').or.
     &      (filab(18)(1:4).eq.'PHS ').or.
     &      (filab(20)(1:4).eq.'MAXS').or.
     &      (filab(26)(1:4).eq.'CONT').or.
     &      (filab(27)(1:4).eq.'CELS')) intpointvar=.true.
         do i=1,nprint
            if((prlab(i)(1:4).eq.'S   ').or.
     &           (prlab(i)(1:4).eq.'E   ').or.
     &           (prlab(i)(1:4).eq.'PEEQ').or.
     &           (prlab(i)(1:4).eq.'ENER').or.
     &           (prlab(i)(1:4).eq.'ELKE').or.
     &           (prlab(i)(1:4).eq.'CELS').or.
     &           (prlab(i)(1:4).eq.'SDV ').or.
     &           (prlab(i)(1:4).eq.'RF  ')) then
               intpointvar=.true.
               exit
            endif
         enddo
      endif
!     
!     calculation of the stresses in the integration points
! 
      if(((ithermal(1).le.1).or.(ithermal(1).ge.3)).and.
     &     (intpointvar)) then
!
c         do i=1,nk
c            write(*,*) 'results v ',i,(v(j,i),j=1,3)
c         enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         imat=ielmat(i)
         amat=matname(imat)
         if(norien.gt.0) then
            iorien=ielorien(i)
         else
            iorien=0
         endif
!
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(1:1).eq.'E') then
            read(lakon(i)(8:8),'(i1)') nope
!
!           local contact spring number
!
            if(lakon(i)(7:7).eq.'C') konl(nope+1)=kon(indexe+nope+1)
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            mint3d=8
         elseif(lakon(i)(4:4).eq.'2') then
            mint3d=27
         elseif(lakon(i)(4:5).eq.'10') then
            mint3d=4
         elseif(lakon(i)(4:4).eq.'4') then
            mint3d=1
         elseif(lakon(i)(4:5).eq.'15') then
            mint3d=9
         elseif(lakon(i)(4:4).eq.'6') then
            mint3d=2
         elseif(lakon(i)(1:1).eq.'E') then
            mint3d=0
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
               voldl(k,j)=vold(k,konl(j))
            enddo
c            write(*,*) 'noeie',i,konl(j),(vl(k,j),k=1,3)
         enddo
!
!        check for hyperelastic material
!
         if(nelcon(1,imat).lt.0) then
            ihyper=1
         else
            ihyper=0
         endif
!
!        q contains the nodal forces per element; initialisation of q
!
         if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &      then
            do m1=1,nope
               do m2=0,mi(2)
                  q(m2,m1)=fn(m2,konl(m1))
               enddo
            enddo
         endif
!
!        calculating the forces for the contact elements
!
         if(mint3d.eq.0) then
!
            lakonl=lakon(i)
!
!           "normal" spring and dashpot elements
!
            if(lakonl(7:7).eq.'A') then
               kode=nelcon(1,imat)
               t0l=0.d0
               t1l=0.d0
               if(ithermal(1).eq.1) then
                  t0l=(t0(konl(1))+t0(konl(2)))/2.d0
                  t1l=(t1(konl(1))+t1(konl(2)))/2.d0
               elseif(ithermal(1).ge.2) then
                  t0l=(t0(konl(1))+t0(konl(2)))/2.d0
                  t1l=(vold(0,konl(1))+vold(0,konl(2)))/2.d0
               endif
            endif
!
!           spring elements (including contact springs)
!     
            if(lakonl(2:2).eq.'S') then
!
!              velocity may be needed for contact springs
!
               if(lakonl(7:7).eq.'C') then
                  do j=1,nope
                     do k=1,3
                        veoldl(k,j)=veold(k,konl(j))
                     enddo
                  enddo
               endif
               call springforc(xl,konl,vl,imat,elcon,nelcon,elas,
     &              fnl,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &              plicon,nplicon,npmat_,veoldl,ener(1,i),iener,
     &              stx(1,1,i),mi,springarea(konl(nope+1)),nmethod)
               do j=1,nope
                  do k=1,3
                     fn(k,konl(j))=fn(k,konl(j))+fnl(k,j)
                  enddo
               enddo
!
!              dashpot elements (including contact dashpots)
!
            elseif((nmethod.eq.4).or.
     &             ((nmethod.eq.1).and.(iperturb(1).ge.2))) then
               do j=1,nope
                  konl(j)=kon(indexe+j)
                  do k=1,3
                     vel(k,j)=veold(k,konl(j))
                  enddo
               enddo
               call dashforc(xl,konl,vl,imat,elcon,nelcon,
     &              elas,fn,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,
     &              elconloc,plicon,nplicon,npmat_,vel,time,nmethod,mi)
            endif
         elseif(ikin.eq.1) then
            do j=1,nope
               do k=1,3
                  veoldl(k,j)=veold(k,konl(j))
               enddo
            enddo            
         endif
!
         do jj=1,mint3d
            if(lakon(i)(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,jj)
c               if(nope.eq.20) xi=xi+1.d0
               et=gauss3d2(2,jj)
               ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            elseif(lakon(i)(4:4).eq.'2') then
c               xi=gauss3d3(1,jj)+1.d0
               xi=gauss3d3(1,jj)
               et=gauss3d3(2,jj)
               ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif(lakon(i)(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif(lakon(i)(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakon(i)(4:5).eq.'15') then
               xi=gauss3d8(1,jj)
               et=gauss3d8(2,jj)
               ze=gauss3d8(3,jj)
               weight=weight3d8(jj)
            elseif(lakon(i)(4:4).eq.'6') then
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif
!
            if(nope.eq.20) then
               if(lakon(i)(7:7).eq.'A') then
                  call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
               elseif((lakon(i)(7:7).eq.'E').or.
     &                (lakon(i)(7:7).eq.'S')) then
                  call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               endif
            elseif(nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m2=1,3
               do m3=1,3
                  vkl(m2,m3)=0.d0
               enddo
            enddo
!
            do m1=1,nope
               do m2=1,3
                  do m3=1,3
                     vkl(m2,m3)=vkl(m2,m3)+shp(m3,m1)*vl(m2,m1)
                  enddo
c                  write(*,*) 'vnoeie',i,konl(m1),(vkl(m2,k),k=1,3)
               enddo
            enddo
!
!           for frequency analysis or buckling with preload the
!           strains are calculated with respect to the deformed
!           configuration
!           for a linear iteration within a nonlinear increment:
!           the tangent matrix is calculated at strain at the end
!           of the previous increment
!
            if((iperturb(1).eq.1).or.(iperturb(1).eq.-1))then
               do m2=1,3
                  do m3=1,3
                     vokl(m2,m3)=0.d0
                  enddo
               enddo
!
               do m1=1,nope
                  do m2=1,3
                     do m3=1,3
                        vokl(m2,m3)=vokl(m2,m3)+
     &                       shp(m3,m1)*voldl(m2,m1)
                     enddo
                  enddo
               enddo
            endif
!
            kode=nelcon(1,imat)
!
!           calculating the strain
!
!           attention! exy,exz and eyz are engineering strains!
!
            exx=vkl(1,1)
            eyy=vkl(2,2)
            ezz=vkl(3,3)
            exy=vkl(1,2)+vkl(2,1)
            exz=vkl(1,3)+vkl(3,1)
            eyz=vkl(2,3)+vkl(3,2)
!
!           for frequency analysis or buckling with preload the
!           strains are calculated with respect to the deformed
!           configuration
!
            if(iperturb(1).eq.1) then
               exx=exx+vokl(1,1)*vkl(1,1)+vokl(2,1)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,1)
               eyy=eyy+vokl(1,2)*vkl(1,2)+vokl(2,2)*vkl(2,2)+
     &              vokl(3,2)*vkl(3,2)
               ezz=ezz+vokl(1,3)*vkl(1,3)+vokl(2,3)*vkl(2,3)+
     &              vokl(3,3)*vkl(3,3)
               exy=exy+vokl(1,1)*vkl(1,2)+vokl(1,2)*vkl(1,1)+
     &              vokl(2,1)*vkl(2,2)+vokl(2,2)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,2)+vokl(3,2)*vkl(3,1)
               exz=exz+vokl(1,1)*vkl(1,3)+vokl(1,3)*vkl(1,1)+
     &              vokl(2,1)*vkl(2,3)+vokl(2,3)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,3)+vokl(3,3)*vkl(3,1)
               eyz=eyz+vokl(1,2)*vkl(1,3)+vokl(1,3)*vkl(1,2)+
     &              vokl(2,2)*vkl(2,3)+vokl(2,3)*vkl(2,2)+
     &              vokl(3,2)*vkl(3,3)+vokl(3,3)*vkl(3,2)
            endif
!
c            if(iperturb(1).ge.2) then
            if(iperturb(2).eq.1) then
!     
!                 Lagrangian strain
!     
               exx=exx+(vkl(1,1)**2+vkl(2,1)**2+vkl(3,1)**2)/2.d0
               eyy=eyy+(vkl(1,2)**2+vkl(2,2)**2+vkl(3,2)**2)/2.d0
               ezz=ezz+(vkl(1,3)**2+vkl(2,3)**2+vkl(3,3)**2)/2.d0
               exy=exy+vkl(1,1)*vkl(1,2)+vkl(2,1)*vkl(2,2)+
     &              vkl(3,1)*vkl(3,2)
               exz=exz+vkl(1,1)*vkl(1,3)+vkl(2,1)*vkl(2,3)+
     &              vkl(3,1)*vkl(3,3)
               eyz=eyz+vkl(1,2)*vkl(1,3)+vkl(2,2)*vkl(2,3)+
     &              vkl(3,2)*vkl(3,3)
!     
            endif
!
!              storing the local strains
!
            if(iperturb(1).ne.-1) then
               eloc(1)=exx
               eloc(2)=eyy
               eloc(3)=ezz
               eloc(4)=exy/2.d0
               eloc(5)=exz/2.d0
               eloc(6)=eyz/2.d0
            else
!
!              linear iteration within a nonlinear increment:
!
               eloc(1)=vokl(1,1)+
     &           (vokl(1,1)**2+vokl(2,1)**2+vokl(3,1)**2)/2.d0
               eloc(2)=vokl(2,2)+
     &           (vokl(1,2)**2+vokl(2,2)**2+vokl(3,2)**2)/2.d0
               eloc(3)=vokl(3,3)+
     &           (vokl(1,3)**2+vokl(2,3)**2+vokl(3,3)**2)/2.d0
               eloc(4)=(vokl(1,2)+vokl(2,1)+vokl(1,1)*vokl(1,2)+
     &              vokl(2,1)*vokl(2,2)+vokl(3,1)*vokl(3,2))/2.d0
               eloc(5)=(vokl(1,3)+vokl(3,1)+vokl(1,1)*vokl(1,3)+
     &              vokl(2,1)*vokl(2,3)+vokl(3,1)*vokl(3,3))/2.d0
               eloc(6)=(vokl(2,3)+vokl(3,2)+vokl(1,2)*vokl(1,3)+
     &              vokl(2,2)*vokl(2,3)+vokl(3,2)*vokl(3,3))/2.d0
            endif
!
!                 calculating the deformation gradient (needed to
!                 convert the element stiffness matrix from spatial
!                 coordinates to material coordinates
!                 deformation plasticity)
!
            if((kode.eq.-50).or.(kode.le.-100)) then
!
!                    calculating the deformation gradient
!
               xkl(1,1)=vkl(1,1)+1
               xkl(2,2)=vkl(2,2)+1.
               xkl(3,3)=vkl(3,3)+1.
               xkl(1,2)=vkl(1,2)
               xkl(1,3)=vkl(1,3)
               xkl(2,3)=vkl(2,3)
               xkl(2,1)=vkl(2,1)
               xkl(3,1)=vkl(3,1)
               xkl(3,2)=vkl(3,2)
!
!                    calculating the Jacobian
!
               vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &              -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &              +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
!
!              inversion of the deformation gradient (only for
!              deformation plasticity)
!
               if(kode.eq.-50) then
!
                  ckl(1,1)=(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))/vj
                  ckl(2,2)=(xkl(1,1)*xkl(3,3)-xkl(1,3)*xkl(3,1))/vj
                  ckl(3,3)=(xkl(1,1)*xkl(2,2)-xkl(1,2)*xkl(2,1))/vj
                  ckl(1,2)=(xkl(1,3)*xkl(3,2)-xkl(1,2)*xkl(3,3))/vj
                  ckl(1,3)=(xkl(1,2)*xkl(2,3)-xkl(2,2)*xkl(1,3))/vj
                  ckl(2,3)=(xkl(2,1)*xkl(1,3)-xkl(1,1)*xkl(2,3))/vj
                  ckl(2,1)=(xkl(3,1)*xkl(2,3)-xkl(2,1)*xkl(3,3))/vj
                  ckl(3,1)=(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))/vj
                  ckl(3,2)=(xkl(3,1)*xkl(1,2)-xkl(1,1)*xkl(3,2))/vj
!
!                 converting the Lagrangian strain into Eulerian
!                 strain (only for deformation plasticity)
!
                  cauchy=.false.
                  call str2mat(eloc,ckl,vj,cauchy)
               endif
!                        
            endif
!
!                 calculating fields for incremental plasticity
!
            if(kode.le.-100) then
!
!              calculating the deformation gradient at the
!              start of the increment
!
!              calculating the displacement gradient at the
!              start of the increment
!
               do m2=1,3
                  do m3=1,3
                     vikl(m2,m3)=0.d0
                  enddo
               enddo
!
               do m1=1,nope
                  do m2=1,3
                     do m3=1,3
                        vikl(m2,m3)=vikl(m2,m3)
     &                       +shp(m3,m1)*vini(m2,konl(m1))
                     enddo
                  enddo
               enddo
!
!              calculating the deformation gradient of the old
!              fields
!
               xikl(1,1)=vikl(1,1)+1
               xikl(2,2)=vikl(2,2)+1.
               xikl(3,3)=vikl(3,3)+1.
               xikl(1,2)=vikl(1,2)
               xikl(1,3)=vikl(1,3)
               xikl(2,3)=vikl(2,3)
               xikl(2,1)=vikl(2,1)
               xikl(3,1)=vikl(3,1)
               xikl(3,2)=vikl(3,2)
!
!              calculating the Jacobian
!
               vij=xikl(1,1)*(xikl(2,2)*xikl(3,3)
     &              -xikl(2,3)*xikl(3,2))
     &              -xikl(1,2)*(xikl(2,1)*xikl(3,3)
     &              -xikl(2,3)*xikl(3,1))
     &              +xikl(1,3)*(xikl(2,1)*xikl(3,2)
     &              -xikl(2,2)*xikl(3,1))
!
!              stresses at the start of the increment
!
               do m1=1,6
                  stre(m1)=stiini(m1,jj,i)
               enddo
!
            endif
!
!                 prestress values
!
            if(iprestr.ne.1) then
               do kk=1,6
                  beta(kk)=0.d0
               enddo
            else
               do kk=1,6
                  beta(kk)=-prestr(kk,jj,i)
               enddo
            endif
!
            if(ithermal(1).ge.1) then
!
!              calculating the temperature difference in
!              the integration point
!
               t0l=0.d0
               t1l=0.d0
               if(ithermal(1).eq.1) then
                  if(lakon(i)(4:5).eq.'8 ') then
                     do i1=1,nope
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+t1(konl(i1))/8.d0
                     enddo
                  elseif(lakon(i)(4:6).eq.'20 ') then
                     call lintemp(t0,t1,konl,nope,jj,t0l,t1l)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*t1(konl(i1))
                     enddo
                  endif
               elseif(ithermal(1).ge.2) then
                  if(lakon(i)(4:5).eq.'8 ') then
                     do i1=1,nope
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+vold(0,konl(i1))/8.d0
                     enddo
                  elseif(lakon(i)(4:6).eq.'20 ') then
                     call lintemp_th(t0,vold,konl,nope,jj,t0l,t1l,mi)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                     enddo
                  endif
               endif
               tt=t1l-t0l
            endif
!
!                 calculating the coordinates of the integration point
!                 for material orientation purposes (for cylindrical
!                 coordinate systems)
!
            if((iorien.gt.0).or.(kode.le.-100)) then
               do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                     pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data; for linear elastic materials
!                 this includes the calculation of the stiffness
!                 matrix
!
            istiff=0
!     
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &           nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &           elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &           istiff,elconloc,eth,kode,plicon,nplicon,
     &           plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,i,jj,
     &           xstiff,ncmat_)
!
!           determining the mechanical strain
!
            if(ithermal(1).ne.0) then
               do m1=1,6
                  emec(m1)=eloc(m1)-eth(m1)
                  emec0(m1)=eme(m1,jj,i)
               enddo
            else
               do m1=1,6
                  emec(m1)=eloc(m1)
                  emec0(m1)=eme(m1,jj,i)
               enddo
            endif
!
!           subtracting the plastic initial strains
!
            if(iprestr.eq.2) then
               do m1=1,6
                  emec(m1)=emec(m1)-prestr(m1,jj,i)
               enddo
            endif
!
!           calculating the local stiffness and stress
!
            call mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &           icmd,beta,stre,xkl,ckl,vj,xikl,vij,
     &           plconloc,xstate,xstateini,ielas,
     &           amat,t1l,dtime,time,ttime,i,jj,nstate_,mi(1),
     &           iorien,pgauss,orab,eloc,mattyp,qa(3),istep,iinc)
!
            do m1=1,21
               xstiff(m1,jj,i)=elas(m1)
            enddo
!
            if(iperturb(1).eq.-1) then
!
!                    if the forced displacements were changed at
!                    the start of a nonlinear step, the nodal
!                    forces due do this displacements are 
!                    calculated in a purely linear way, and
!                    the first iteration is purely linear in order
!                    to allow the displacements to redistribute
!                    in a quasi-static way (only applies to
!                    quasi-static analyses (*STATIC))
!
               eloc(1)=exx-vokl(1,1)
               eloc(2)=eyy-vokl(2,2)
               eloc(3)=ezz-vokl(3,3)
               eloc(4)=exy-(vokl(1,2)+vokl(2,1))
               eloc(5)=exz-(vokl(1,3)+vokl(3,1))
               eloc(6)=eyz-(vokl(2,3)+vokl(3,2))
!
               if(mattyp.eq.1) then
                  e=elas(1)
                  un=elas(2)
                  um=e/(1.d0+un)
                  al=un*um/(1.d0-2.d0*un)
                  um=um/2.d0
                  am1=al*(eloc(1)+eloc(2)+eloc(3))
                  stre(1)=am1+2.d0*um*eloc(1)
                  stre(2)=am1+2.d0*um*eloc(2)
                  stre(3)=am1+2.d0*um*eloc(3)
                  stre(4)=um*eloc(4)
                  stre(5)=um*eloc(5)
                  stre(6)=um*eloc(6)
               elseif(mattyp.eq.2) then
                  stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)
     &                 +eloc(3)*elas(4)
                  stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)
     &                 +eloc(3)*elas(5)
                  stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)
     &                 +eloc(3)*elas(6)
                  stre(4)=eloc(4)*elas(7)
                  stre(5)=eloc(5)*elas(8)
                  stre(6)=eloc(6)*elas(9)
               elseif(mattyp.eq.3) then
                  stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)+
     &                 eloc(3)*elas(4)+eloc(4)*elas(7)+
     &                 eloc(5)*elas(11)+eloc(6)*elas(16)
                  stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)+
     &                 eloc(3)*elas(5)+eloc(4)*elas(8)+
     &                 eloc(5)*elas(12)+eloc(6)*elas(17)
                  stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)+
     &                 eloc(3)*elas(6)+eloc(4)*elas(9)+
     &                 eloc(5)*elas(13)+eloc(6)*elas(18)
                  stre(4)=eloc(1)*elas(7)+eloc(2)*elas(8)+
     &                 eloc(3)*elas(9)+eloc(4)*elas(10)+
     &                 eloc(5)*elas(14)+eloc(6)*elas(19)
                  stre(5)=eloc(1)*elas(11)+eloc(2)*elas(12)+
     &                 eloc(3)*elas(13)+eloc(4)*elas(14)+
     &                 eloc(5)*elas(15)+eloc(6)*elas(20)
                  stre(6)=eloc(1)*elas(16)+eloc(2)*elas(17)+
     &                 eloc(3)*elas(18)+eloc(4)*elas(19)+
     &                 eloc(5)*elas(20)+eloc(6)*elas(21)
               endif
            endif
! 
!           updating the internal energy
!
            if((iout.gt.0).or.(iout.eq.-2)) then
               if(ithermal(1).eq.0) then
                  do m1=1,6
                     eth(m1)=0.d0
                  enddo
               endif
               if(iener.eq.1) then
                  ener(jj,i)=enerini(jj,i)+
     &                 ((eloc(1)-eth(1)-eme(1,jj,i))*
     &                  (stre(1)+stiini(1,jj,i))+
     &                  (eloc(2)-eth(2)-eme(2,jj,i))*
     &                  (stre(2)+stiini(2,jj,i))+
     &                  (eloc(3)-eth(3)-eme(3,jj,i))*
     &                  (stre(3)+stiini(3,jj,i)))/2.d0+
     &            (eloc(4)-eth(4)-eme(4,jj,i))*(stre(4)+stiini(4,jj,i))+
     &            (eloc(5)-eth(5)-eme(5,jj,i))*(stre(5)+stiini(5,jj,i))+
     &            (eloc(6)-eth(6)-eme(6,jj,i))*(stre(6)+stiini(6,jj,i))

               endif
!
               eme(1,jj,i)=eloc(1)-eth(1)
               eme(2,jj,i)=eloc(2)-eth(2)
               eme(3,jj,i)=eloc(3)-eth(3)
               eme(4,jj,i)=eloc(4)-eth(4)
               eme(5,jj,i)=eloc(5)-eth(5)
               eme(6,jj,i)=eloc(6)-eth(6)
!
               eei(1,jj,i)=eloc(1)
               eei(2,jj,i)=eloc(2)
               eei(3,jj,i)=eloc(3)
               eei(4,jj,i)=eloc(4)
               eei(5,jj,i)=eloc(5)
               eei(6,jj,i)=eloc(6)
            endif
!
!     updating the kinetic energy
!
            if(ikin.eq.1) then
               
               call materialdata_rho(rhcon,nrhcon,imat,rho,t1l,
     &              ntmat_)
               do m1=1,3
                  vel(m1,1)=0.d0
                  do i1= 1,nope
                     vel(m1,1)=vel(m1,1)+shp(4,i1)*veoldl(m1,i1)
                  enddo
               enddo
               ener(jj,i+ne)=rho*(vel(1,1)*vel(1,1)+
     &              vel(2,1)*vel(2,1)+ vel(3,1)*vel(3,1))/2.d0
            endif
!
            skl(1,1)=stre(1)
            skl(2,2)=stre(2)
            skl(3,3)=stre(3)
            skl(2,1)=stre(4)
            skl(3,1)=stre(5)
            skl(3,2)=stre(6)
!
            stx(1,jj,i)=skl(1,1)
            stx(2,jj,i)=skl(2,2)
            stx(3,jj,i)=skl(3,3)
            stx(4,jj,i)=skl(2,1)
            stx(5,jj,i)=skl(3,1)
            stx(6,jj,i)=skl(3,2)
!
            skl(1,2)=skl(2,1)
            skl(1,3)=skl(3,1)
            skl(2,3)=skl(3,2)
!
!                 calculation of the nodal forces
!
c            if(iperturb(2).eq.0) then
c               do m1=1,3
c                  do m2=1,3
c                     vkl(m1,m2)=0.d0
c                  enddo
c               enddo
c            endif
!
            if(calcul_fn)then
!
!                    calculating fn using skl
!
               do m1=1,nope
                  do m2=1,3
!
!                          linear elastic part
!                           
                     do m3=1,3
                        fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                       xsj*skl(m2,m3)*shp(m3,m1)*weight
                     enddo
!
!                          nonlinear geometric part
!
c                     if(iperturb(1).ge.2) then
                     if(iperturb(2).eq.1) then
                        do m3=1,3
                           do m4=1,3
                              fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                             xsj*skl(m4,m3)*weight*
     &                             (vkl(m2,m4)*shp(m3,m1)+
     &                             vkl(m2,m3)*shp(m4,m1))/2.d0
                           enddo
                        enddo
                     endif
!
                  enddo
               enddo
            endif
!
!           calculation of the Cauchy stresses
!
            if(calcul_cauchy) then
!
!              changing the displacement gradients into
!              deformation gradients
!
c               if(kode.ne.-50) then
               if((kode.ne.-50).and.(kode.gt.-100)) then
                  xkl(1,1)=vkl(1,1)+1
                  xkl(2,2)=vkl(2,2)+1.
                  xkl(3,3)=vkl(3,3)+1.
                  xkl(1,2)=vkl(1,2)
                  xkl(1,3)=vkl(1,3)
                  xkl(2,3)=vkl(2,3)
                  xkl(2,1)=vkl(2,1)
                  xkl(3,1)=vkl(3,1)
                  xkl(3,2)=vkl(3,2)
!
                  vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &                 -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &                 +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
               endif
!
               do m1=1,3
                  do m2=1,m1
                     ckl(m1,m2)=0.d0
                     do m3=1,3
                        do m4=1,3
                           ckl(m1,m2)=ckl(m1,m2)+
     &                          skl(m3,m4)*xkl(m1,m3)*xkl(m2,m4)
                        enddo
                     enddo
                     ckl(m1,m2)=ckl(m1,m2)/vj
                  enddo
               enddo
!
               stx(1,jj,i)=ckl(1,1)
               stx(2,jj,i)=ckl(2,2)
               stx(3,jj,i)=ckl(3,3)
               stx(4,jj,i)=ckl(2,1)
               stx(5,jj,i)=ckl(3,1)
               stx(6,jj,i)=ckl(3,2)
            endif
!
         enddo
!
!        q contains the contributions to the nodal force in the nodes
!        belonging to the element at stake from other elements (elements
!        already treated). These contributions have to be
!        subtracted to get the contributions attributable to the element
!        at stake only
!
         if(calcul_qa) then
            do m1=1,nope
               do m2=1,3
                  qa(1)=qa(1)+dabs(fn(m2,konl(m1))-q(m2,m1))
               enddo
            enddo
            nal=nal+3*nope
         endif
      enddo
!
      if(calcul_qa) then
         if(nal.gt.0) then
            qa(1)=qa(1)/nal
         endif
      endif
!
      endif
!
!     calculation of temperatures and thermal flux
!
      qa(2)=0.d0
      nal=0
!
!     check whether integration point variables are needed in
!     modal dynamics calculations
!
      if((nmethod.eq.4).and.(iperturb(1).lt.2)) then
         intpointvar=.false.
         if((filab(9)(1:4).eq.'HFL ').or.
     &      (filab(10)(1:4).eq.'RFL ')) intpointvar=.true.
         do i=1,nprint
            if((prlab(i)(1:4).eq.'HFL ').or.
     &         (prlab(i)(1:4).eq.'RFL ')) intpointvar=.true.
         enddo
      endif
!
      if((ithermal(1).ge.2).and.(intpointvar)) then
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         imat=ielmat(i)
         amat=matname(imat)
         if(norien.gt.0) then
            iorien=ielorien(i)
         else
            iorien=0
         endif
!
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            if(lakon(i)(6:7).eq.'RA') then
               mint3d=4
            else
               mint3d=8
            endif
         elseif(lakon(i)(4:4).eq.'2') then
            mint3d=27
         elseif(lakon(i)(4:5).eq.'10') then
            mint3d=4
         elseif(lakon(i)(4:4).eq.'4') then
            mint3d=1
         elseif(lakon(i)(4:5).eq.'15') then
            mint3d=9
         elseif(lakon(i)(4:4).eq.'6') then
            mint3d=2
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
            enddo
            vl(0,j)=v(0,konl(j))
            voldl(0,j)=vold(0,konl(j))
         enddo
!
!        q contains the nodal forces per element; initialisation of q
!
         if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &          then
            do m1=1,nope
               q(0,m1)=fn(0,konl(m1))
            enddo
         endif
!
         do jj=1,mint3d
            if(lakon(i)(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,jj)
               et=gauss3d2(2,jj)
               ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            elseif(lakon(i)(4:4).eq.'2') then
               xi=gauss3d3(1,jj)
               et=gauss3d3(2,jj)
               ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif(lakon(i)(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif(lakon(i)(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakon(i)(4:5).eq.'15') then
               xi=gauss3d8(1,jj)
               et=gauss3d8(2,jj)
               ze=gauss3d8(3,jj)
               weight=weight3d8(jj)
            elseif(lakon(i)(4:4).eq.'6') then
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif
!
            if(nope.eq.20) then
               if(lakon(i)(7:7).eq.'A') then
                  call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
               elseif((lakon(i)(7:7).eq.'E').or.
     &                (lakon(i)(7:7).eq.'S')) then
                  call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               endif
            elseif(nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
            c1=xsj*weight
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m3=1,3
               vkl(0,m3)=0.d0
            enddo
!
            do m1=1,nope
               do m3=1,3
                  vkl(0,m3)=vkl(0,m3)+shp(m3,m1)*vl(0,m1)
               enddo
            enddo
!
            kode=ncocon(1,imat)
!
!              calculating the temperature difference in
!              the integration point
!
            t1lold=0.d0
            t1l=0.d0
            if(lakon(i)(4:5).eq.'8 ') then
               do i1=1,nope
                  t1lold=t1lold+vold(0,konl(i1))/8.d0
                  t1l=t1l+v(0,konl(i1))/8.d0
               enddo
            elseif(lakon(i)(4:6).eq.'20 ') then
               call lintemp_th(t0,vold,konl,nope,jj,t0l,t1lold,mi)
               call lintemp_th(t0,v,konl,nope,jj,t0l,t1l,mi)
            else
               do i1=1,nope
                  t1lold=t1lold+shp(4,i1)*vold(0,konl(i1))
                  t1l=t1l+shp(4,i1)*v(0,konl(i1))
               enddo
            endif
!
!           calculating the coordinates of the integration point
!           for material orientation purposes (for cylindrical
!           coordinate systems)
!
            if((iorien.gt.0).or.(kode.le.-100)) then
               do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                     pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data; for linear elastic materials
!                 this includes the calculation of the stiffness
!                 matrix
!
            istiff=0
!
            call materialdata_th(cocon,ncocon,imat,iorien,pgauss,orab,
     &           ntmat_,coconloc,mattyp,t1l,rhcon,nrhcon,rho,shcon,
     &           nshcon,sph,xstiff,jj,i,istiff,mi(1))
!
            call thermmodel(amat,i,jj,kode,coconloc,vkl,dtime,
     &           time,ttime,mi(1),nstate_,xstateini,xstate,qflux,xstiff,
     &           iorien,pgauss,orab,t1l,t1lold,vold,co,lakon(i),konl,
     &           ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc)
! 
            qfx(1,jj,i)=qflux(1)
            qfx(2,jj,i)=qflux(2)
            qfx(3,jj,i)=qflux(3)
            if(lakon(i)(6:7).eq.'RA') then
               qfx(1,jj+4,i)=qflux(1)
               qfx(2,jj+4,i)=qflux(2)
               qfx(3,jj+4,i)=qflux(3)
            endif
!
!           calculation of the nodal flux
!
            if(calcul_fn)then
!
!                    calculating fn using skl
!
               if(lakon(i)(6:7).eq.'RA') then
                  do m1=1,nope
                     fn(0,konl(m1))=fn(0,konl(m1))
     &                  -c1*(qflux(1)*(shp(1,m1)+shp(1,iperm(m1)))
     &                      +qflux(2)*(shp(2,m1)+shp(2,iperm(m1)))
     &                      +qflux(3)*(shp(3,m1)+shp(3,iperm(m1))))
                  enddo
               else
                  do m1=1,nope
                     do m3=1,3
                        fn(0,konl(m1))=fn(0,konl(m1))-
     &                       c1*qflux(m3)*shp(m3,m1)
                     enddo
                  enddo
               endif
            endif
         enddo
!
!        q contains the contributions to the nodal force in the nodes
!        belonging to the element at stake from other elements (elements
!        already treated). These contributions have to be
!        subtracted to get the contributions attributable to the element
!        at stake only
!
         if(calcul_qa) then
            do m1=1,nope
               qa(2)=qa(2)+dabs(fn(0,konl(m1))-q(0,m1))
            enddo
            nal=nal+nope
         endif
      enddo
!
      endif
!
      if(calcul_qa) then
         if(nal.gt.0) then
            qa(2)=qa(2)/nal
         endif
      endif
!
!     subtracting the mpc force (for each linear mpc there is one
!     force; the actual force in a node belonging to the mpc is
!     obtained by multiplying this force with the nodal coefficient.
!     The force has to be subtracted from f, since it does not
!     appear on the rhs of the equations system
!
      if(calcul_fn)then
        do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.gt.3) cycle
            forcempc=fn(ndir,node)/coefmpc(ist)
            fmpc(i)=forcempc
            fn(ndir,node)=0.d0
            index=nodempc(3,ist)
            if(index.eq.0) cycle
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               fn(ndir,node)=fn(ndir,node)-coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         enddo
      endif
!
!     calculating the system force vector
!
      if(calcul_f) then
         do i=1,nk
            do j=0,mi(2)
               if(nactdof(j,i).ne.0) then
                  f(nactdof(j,i))=fn(j,i)
               endif
            enddo
         enddo
      endif
!
!     adding the mpc force again to fn
!
      if(calcul_fn)then
         do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.gt.3) cycle
            forcempc=fmpc(i)
            fn(ndir,node)=forcempc*coefmpc(ist)
            index=nodempc(3,ist)
!
!           nodes not belonging to the structure have to be
!           taken out
!
            if(labmpc(i)(1:7).eq.'MEANROT') then
               if(nodempc(3,nodempc(3,index)).eq.0) cycle
            elseif(labmpc(i)(1:10).eq.'PRETENSION') then
               if(nodempc(3,index).eq.0) cycle
            elseif(labmpc(i)(1:5).eq.'RIGID') then
               if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,inde
     &x))))).eq.0) cycle
            else
               if(index.eq.0) cycle
            endif
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               fn(ndir,node)=fn(ndir,node)+coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(labmpc(i)(1:7).eq.'MEANROT') then
                  if(nodempc(3,nodempc(3,index)).eq.0) exit
               elseif(labmpc(i)(1:10).eq.'PRETENSION') then
                  if(nodempc(3,index).eq.0) exit
               elseif(labmpc(i)(1:5).eq.'RIGID') then
                  if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,i
     &ndex))))).eq.0) exit
               else
                  if(index.eq.0) exit
               endif
            enddo
         enddo
      endif
!
!     no print requests
!
      if(iout.le.0) return
!
!     output in dat file (with *NODE PRINT or *EL PRINT)
!
      call printout(set,nset,istartset,iendset,ialset,nprint,
     &  prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &  mi(1),nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &  orab,ielorien,norien,nk,ne,inum,filab,vold,ikin)
!
!     interpolation in the original nodes of 1d and 2d elements
!     this operation has to be performed in any case since
!     the interpolated values may be needed as boundary conditions
!     in the next step (e.g. the temperature in a heat transfer
!     calculation as boundary condition in a subsequent static
!     step)
!
      if(filab(1)(5:5).ne.' ') then
         nfield=mt
         cflag=filab(1)(5:5)
         force=.false.
         call map3dto1d2d(v,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,cflag,co,vold,force,mi)
      endif
!
!     user defined output
!
      call uout(v,mi)
!
      if((filab(2)(1:4).eq.'NT  ').and.(ithermal(1).le.1)) then
         if(filab(2)(5:5).eq.'I') then
            nfield=1
            cflag=filab(2)(5:5)
            force=.false.
            call map3dto1d2d(t1,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,force,mi)
         endif
      endif
!
!     determining the stresses in the nodes for output in frd format
!
      if((filab(3)(1:4).eq.'S   ').or.(filab(18)(1:4).eq.'PHS ').or.
     &   (filab(20)(1:4).eq.'MAXS')) then
         nfield=6
         ndim=6
         if((norien.gt.0).and.(filab(3)(6:6).eq.'L')) then
            iorienglob=1
         else
            iorienglob=0
         endif
         cflag=filab(3)(5:5)
!
         call extrapolate(stx,stn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
!
      endif
!
!     determining the strains in the nodes for output in frd format
!
      if(filab(4)(1:4).eq.'E   ') then
         nfield=6
         ndim=6
         if((norien.gt.0).and.(filab(4)(6:6).eq.'L')) then
            iorienglob=1
         else
            iorienglob=0
         endif
         cflag=filab(4)(5:5)
         call extrapolate(eei,een,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
!     determining the plastic equivalent strain in the nodes 
!     for output in frd format
!
      if(filab(6)(1:4).eq.'PEEQ') then
         nfield=1
         ndim=nstate_
         iorienglob=0
         cflag=filab(6)(5:5)
         call extrapolate(xstate,epn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
!     determining the total energy in the nodes 
!     for output in frd format
!
      if(filab(7)(1:4).eq.'ENER') then
         nfield=1
         ndim=1
         iorienglob=0
         cflag=filab(7)(5:5)
         call extrapolate(ener,enern,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
!     determining the internal state variables in the nodes 
!     for output in frd format
!
      if(filab(8)(1:4).eq.'SDV ') then
         nfield=nstate_
         ndim=nstate_
         if((norien.gt.0).and.(filab(9)(6:6).eq.'L')) then
            write(*,*) '*WARNING in results: SDV variables cannot'
            write(*,*) '         be stored in a local frame;'
            write(*,*) '         the global frame will be used'
         endif
         iorienglob=0
         cflag=filab(8)(5:5)
         call extrapolate(xstate,xstaten,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
!     determining the heat flux in the nodes for output in frd format
!
      if((filab(9)(1:4).eq.'HFL ').and.(ithermal(1).gt.1)) then
         nfield=3
         ndim=3
         if((norien.gt.0).and.(filab(9)(6:6).eq.'L')) then
            iorienglob=1
         else
            iorienglob=0
         endif
         cflag=filab(9)(5:5)
         call extrapolate(qfx,qfn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
!     if no element quantities requested in the nodes: calculate
!     inum if nodal quantities are requested: used in subroutine frd
!     to determine which nodes are active in the model 
!
c      if((filab(3)(1:4).ne.'S   ').and.(filab(4)(1:4).ne.'E   ').and.
c     &    (filab(6)(1:4).ne.'PEEQ').and.(filab(7)(1:4).ne.'ENER').and.
c     &    (filab(8)(1:4).ne.'SDV ').and.(filab(9)(1:4).ne.'HFL ').and.
c     &    (iinc.le.1)) then
      if((filab(3)(1:4).ne.'S   ').and.(filab(4)(1:4).ne.'E   ').and.
     &    (filab(6)(1:4).ne.'PEEQ').and.(filab(7)(1:4).ne.'ENER').and.
     &    (filab(8)(1:4).ne.'SDV ').and.(filab(9)(1:4).ne.'HFL ')) then
!
         nfield=0
         ndim=0
         iorienglob=0
         cflag=filab(1)(5:5)
         call extrapolate(stx,stn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &        nelemload,nload,nodeboun,nboun,fluid,ndirboun,vold,
     &        ithermal,force)
      endif
!
      if(fluid) then
         call fluidextrapolate(v,ipkon,inum,kon,lakon,ne,mi)
      endif
!
      return
      end
