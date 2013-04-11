!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine mafilltrhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,
     &  b,nactdoh,neqt,nmethod,ikmpc,ilmpc,ikboun,
     &  ilboun,rhcon,nrhcon,ielmat,ntmat_,t0,ithermal,vold,voldcon,nzst,
     &  dtl,matname,mi,ncmat_,physcon,shcon,nshcon,ttime,time,
     &  istep,iinc,ibody,xloadold,reltimef,cocon,ncocon,nelemface,
     &  sideface,nface,compressible,v,voldtu,yy,turbulent,nea,neb,
     &  dtimef,ipvar,var,ipvarf,varf)
!
!     filling the rhs b of the velocity equations (step 1)
!
      implicit none
!
      character*1 sideface(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),nelemface(*),nface,
     &  ikmpc(*),ilmpc(*),ikboun(*),ilboun(*),nactdoh(0:4,*),konl(20),
     &  nrhcon(*),mi(*),ielmat(mi(3),*),ipkon(*),nshcon(*),ipobody(2,*),
     &  nbody,ibody(3,*),ncocon(2,*),compressible,nea,neb,ipvar(*),
     &  ipvarf(*)
!
      integer nk,ne,nboun,nmpc,nforc,nload,neqt,nmethod,
     &  ithermal,nzst,i,j,idist,jj,id,ist,index,jdof1,idof1,
     &  node1,kflag,ntmat_,indexe,nope,i0,ncmat_,istep,iinc,
     &  turbulent
!
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &  p2(3),bodyf(3),b(*),xloadold(2,*),reltimef,cocon(0:6,ntmat_,*),
     &  t0(*),vold(0:mi(2),*),voldcon(0:4,*),ff(60),v(0:mi(2),*),yy(*),
     &  rhcon(0:1,ntmat_,*),physcon(*),voldtu(2,*),
     &  shcon(0:3,ntmat_,*),xbody(7,*),var(*),varf(*)
!
      real*8 om,dtimef,ttime,time,dtl(*)
!
      kflag=2
      i0=0
!
      do i=1,neqt
         b(i)=0.d0
      enddo
!     
!     distributed forces (body forces or thermal loads or
!     residual stresses or distributed face loads)
!     
      if((nbody.ne.0).or.(ithermal.ne.0).or.
     &     (nload.ne.0)) then
         idist=1
      else
         idist=0
      endif
!     
      do i=nea,neb
!     
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'F') cycle
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
c         do j=1,nope
c            konl(j)=kon(indexe+j) 
c         enddo
!     
         om=0.d0
!     
         if(nbody.gt.0) then
!     
!     assigning centrifugal forces
!     
            bodyf(1)=0.
            bodyf(2)=0.
            bodyf(3)=0.
!     
            index=i
            do
               j=ipobody(1,index)
               if(j.eq.0) exit
               if(ibody(1,j).eq.1) then
                  om=xbody(1,j)
                  p1(1)=xbody(2,j)
                  p1(2)=xbody(3,j)
                  p1(3)=xbody(4,j)
                  p2(1)=xbody(5,j)
                  p2(2)=xbody(6,j)
                  p2(3)=xbody(7,j)
!     
!     assigning gravity forces
!     
               elseif(ibody(1,j).eq.2) then
                  bodyf(1)=bodyf(1)+xbody(1,j)*xbody(2,j)
                  bodyf(2)=bodyf(2)+xbody(1,j)*xbody(3,j)
                  bodyf(3)=bodyf(3)+xbody(1,j)*xbody(4,j)
               endif
               index=ipobody(2,index)
               if(index.eq.0) exit
            enddo
         endif
!     
         call e_c3d_trhs(co,nk,kon(indexe+1),lakon(i),p1,p2,om,bodyf,
     &        nbody,ff,i,nmethod,rhcon,nrhcon,
     &        ielmat,ntmat_,vold,voldcon,nelemload,
     &        sideload,xload,nload,idist,dtimef,matname,mi(1),
     &        ttime,time,istep,iinc,xloadold,reltimef,shcon,nshcon,
     &        cocon,ncocon,physcon,nelemface,sideface,nface,
     &        ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,compressible,v,
     &        voldtu,yy,turbulent,ipvar,var,ipvarf,varf)
!     
         do jj=1,nope
!     
            node1=kon(indexe+jj)
c            ff(jj)=ff(jj)*dtl(node1)/dtimef
            jdof1=nactdoh(0,node1)
!     
!     distributed forces
!     
            if(jdof1.eq.0) then
               if(nmpc.ne.0) then
                  idof1=(node1-1)*8
                  call nident(ikmpc,idof1,nmpc,id)
                  if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                     id=ilmpc(id)
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        jdof1=nactdoh(nodempc(2,index),
     &                       nodempc(1,index))
                        if(jdof1.ne.0) then
                           b(jdof1)=b(jdof1)
     &                          -coefmpc(index)*ff(jj)
     &                          /coefmpc(ist)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                  endif
               endif
               cycle
            endif
            b(jdof1)=b(jdof1)+ff(jj)
!     
         enddo
      enddo
!
!     nonlocal time stepping for compressible steady state calculations
!     
c      if((compressible.eq.1).and.(nmethod.eq.1)) then
c         do i=1,nk
c            if(nactdoh(0,i).gt.0) then
c               b(nactdoh(0,i))=b(nactdoh(0,i))*dtl(i)/dtimef
c            endif
c         enddo
c      endif
!     
      return
      end
      
