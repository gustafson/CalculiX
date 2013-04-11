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
      subroutine mafillsmas(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,
     &  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &  ad,au,bb,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &  ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &  nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,prestr,
     &  iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &  physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &  coriolis,ibody,xloadold,reltime,veold,springarea)
!
!     filling the stiffness matrix in spare matrix format (sm)
!     asymmetric contributions
!
      implicit none
!
      logical mass,stiffness,buckling,rhsi,coriolis
!
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),ikmpc(*),
     &  ilmpc(*),ikboun(*),ilboun(*),mi(2),
     &  nactdof(0:mi(2),*),konl(20),irow(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(*),ielorien(*),
     &  ipkon(*),intscheme,ncocon(2,*),nshcon(*),ipobody(2,*),nbody,
     &  ibody(3,*)
!
      integer nk,ne,nboun,nmpc,nforc,nload,neq,nzl,nmethod,
     &  ithermal,iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &  ll,jdof1,jdof2,node1,node2,
     &  ntmat_,indexe,nope,norien,iexpl,ncmat_,istep,iinc
!
      integer nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &  p2(3),ad(*),au(*),bodyf(3),bb(*),xloadold(2,*),
     &  t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),ff(60),
     &  sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),reltime,
     &  alcon(0:6,ntmat_,*),physcon(*),cocon(0:6,ntmat_,*),
     &  shcon(0:3,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &  springarea(*)
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),veold(0:mi(2),*)
!
      real*8 om,dtime,ttime,time
!
!     storing the symmetric matrix in asymmetric format
!
      do i=1,nzs(3)
         au(nzs(3)+i)=au(i)
      enddo
!
!     mechanical analysis: asymmetric contributions
!
      do i=1,ne
!
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'E') cycle
        indexe=ipkon(i)
        nope=4
!
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
!
        call e_c3d(co,nk,konl,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,i,
     &          nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &          alzero,ielmat,ielorien,norien,orab,ntmat_,
     &          t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,
     &          nload,idist,sti,stx,iexpl,plicon,
     &          nplicon,plkcon,nplkcon,xstiff,npmat_,
     &          dtime,matname,mi(1),ncmat_,mass,stiffness,buckling,rhsi,
     &          intscheme,ttime,time,istep,iinc,coriolis,xloadold,
     &          reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &          springarea)
!
        do jj=1,3*nope
!
          j=(jj-1)/3+1
          k=jj-3*(j-1)
!
          node1=kon(indexe+j)
          jdof1=nactdof(k,node1)
!
          do ll=1,3*nope
!
            l=(ll-1)/3+1
            m=ll-3*(l-1)
!
            node2=kon(indexe+l)
            jdof2=nactdof(m,node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.ne.0).and.(jdof2.ne.0)) then
               call add_sm_st_as(au,ad,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),jj,ll,nzs)
            endif
          enddo
        enddo
      enddo
!
      return
      end
