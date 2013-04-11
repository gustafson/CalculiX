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
      subroutine createddentry(itie,ipkon,kon,nodes,
     &  lakon,islavsurf,itiefac,contribution,co,vold,
     &  iponoels,inoels,mi,pslavdual)
!
!     compute the Dd[p,q] matrix entry for contact problems
!     Author: Li, Yang
!
      implicit none
!
      character*8 lakon(*)
!
      logical contained
!
      integer itie,ipkon(*),kon(*),konl(20),iflag,l,j,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),
     &  nope,islavsurf(2,*),itiefac(2,*),ifaces,nelems,jfaces,
     &  mint2d,m,nopes,iq,ii,nodes,indexe,iponoels(*),
     &  inoels(3,*),index1,islavsurfentry,locs,mi(*),ns
!
      real*8 ets,xis,weights,xl2s(3,8),xsj2s(3),co(3,*),
     &  vold(0:mi(2),*),shp2s(7,8),xs2s(3,2),contribution,
     &  shp2d(4,8),xsj2d(3),xs2d(3,2),
     &  pslavdual(16,*)
!
      include "gauss.f"
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             4,6,3,1/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      contribution = 0.d0
      itie = itie + 1
      index1=iponoels(nodes)
      do
         if(index1.eq.0) exit
         islavsurfentry=inoels(1,index1)
         locs=inoels(2,index1)
c      WRITE(*,*) "islavsurfentry",islavsurfentry,"index1",index1
         if((islavsurfentry.lt.itiefac(1,itie)).or.
     &        (islavsurfentry.gt.itiefac(2,itie))) exit
         ifaces=islavsurf(1,islavsurfentry)
         nelems = int(ifaces/10)
         jfaces = ifaces - nelems*10
         indexe = ipkon(nelems)
!     
!     Decide the max integration points number, junt consider 2D situation 
!     
         if(lakon(nelems)(4:5).eq.'8R') then
            mint2d=1
         elseif((lakon(nelems)(4:4).eq.'8').or.
     &           (lakon(nelems)(4:6).eq.'20R')) then
            mint2d=4
         elseif(lakon(nelems)(4:4).eq.'2') then
            mint2d=9
         elseif(lakon(nelems)(4:5).eq.'10') then
            mint2d=3
         elseif(lakon(nelems)(4:4).eq.'4') then
            mint2d=1
         endif
!     
!     treatment of wedge faces
!     
         if(lakon(nelems)(4:4).eq.'6') then
            mint2d=1
            if(jfaces.le.2) then
               nopes=3
            else
               nopes=4
            endif
         endif
         if(lakon(nelems)(4:5).eq.'15') then
            if(jfaces.le.2) then
               mint2d=3
               nopes=6
            else
               mint2d=4
               nopes=8
            endif
         endif
!     
         if(lakon(nelems)(4:4).eq.'2') then
            nope=20
            nopes=8
         elseif(lakon(nelems)(4:4).eq.'8') then
            nope=8
            nopes=4
         elseif(lakon(nelems)(4:5).eq.'10') then
            nope=10
            nopes=6
         elseif(lakon(nelems)(4:4).eq.'4') then
            nope=4
            nopes=3
         elseif(lakon(nelems)(4:5).eq.'15') then
            nope=15
         elseif(lakon(nelems)(4:4).eq.'6') then
            nope=6
         endif
!     
!     determining the nodes belonging to the slave face
!     and their coordinates    
!     
         do j=1,nope
            konl(j)=kon(ipkon(nelems)+j)
         enddo
!     
         if((nope.eq.20).or.(nope.eq.8)) then
            do m=1,nopes
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifaceq(m,jfaces)))+
     &                 vold(j,konl(ifaceq(m,jfaces)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do m=1,nopes
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &                 vold(j,konl(ifacet(m,jfaces)))
               enddo
            enddo
         else
            do m=1,nopes
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &                 vold(j,konl(ifacew1(m,jfaces)))
               enddo
            enddo
         endif
!     
         do m = 1,mint2d
            if((lakon(nelems)(4:5).eq.'8R').or.
     &           ((lakon(nelems)(4:4).eq.'6').and.
     &           (nopes.eq.4))) then
               xis=gauss2d1(1,m)
               ets=gauss2d1(2,m)
               weights=weight2d1(m)
            elseif((lakon(nelems)(4:4).eq.'8').or.
     &              (lakon(nelems)(4:6).eq.'20R').or.
     &              ((lakon(nelems)(4:5).eq.'15').and.
     &              (nopes.eq.8))) then
               xis=gauss2d2(1,m)
               ets=gauss2d2(2,m)
               weights=weight2d2(m)
            elseif(lakon(nelems)(4:4).eq.'2') then
               xis=gauss2d3(1,m)
               ets=gauss2d3(2,m)
               weights=weight2d3(m)
            elseif((lakon(nelems)(4:5).eq.'10').or.
     &              ((lakon(nelems)(4:5).eq.'15').and.
     &              (nopes.eq.6))) then
               xis=gauss2d5(1,m)
               ets=gauss2d5(2,m)
               weights=weight2d5(m)
            elseif((lakon(nelems)(4:4).eq.'4').or.
     &              ((lakon(nelems)(4:4).eq.'6').and.
     &              (nopes.eq.3))) then
               xis=gauss2d4(1,m)
               ets=gauss2d4(2,m)
               weights=weight2d4(m)
            endif
!     
            iflag = 2
            ns=islavsurfentry
            if(nopes.eq.8) then
               call shape8q(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
               call dualshape8q(xis,ets,xl2s,xsj2d,xs2d,shp2d,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
               call dualshape4q(xis,ets,xl2s,xsj2d,xs2d,shp2d,ns,
     &                pslavdual,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
               call dualshape6tri(xis,ets,xl2s,xsj2d,xs2d,shp2d,iflag)
            else
               call shape3tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
               call dualshape3tri(xis,ets,xl2s,xsj2d,xs2d,shp2d,iflag)
            endif
!     
            contribution=contribution+shp2s(4,locs)*
     &           weights*dsqrt(xsj2s(1)**2+xsj2s(2)**2+xsj2s(3)**2)
         enddo
         index1=inoels(3,index1)
      enddo
      itie = itie - 1
!     
      return 
      end
