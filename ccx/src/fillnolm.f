!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
!
c> subroutine to fill in values of LM for slavenodes without LM
c> 
c> 

      subroutine fillnolm(islav,node,itie,ipkon,kon,lakon,
     &  islavsurf,itiefac,iponoels,inoels,mi,
     &  pslavdual,nslavnode,islavnode,cdisp)
!         
!
!
!     Author: Saskia Sitzmann
!
      implicit none
!
      character*8 lakon(*)
!
      logical debug
!
      integer ipkon(*),kon(*),konl(20),iflag,m,l,j,jj,node,islav,
     &  indexe,islavsurf(2,*),iponoels(*),inoels(3,*),
     &  itiefac(2,*),ifaces,nelemens,jfaces,ifacem,
     &  mint2d,indexf,nopes1,nopes2,nodem,nodesf,nopes,
     &  locs,locm,mi(*),ns,mint2dloc1,mint2dloc2,
     &  ifs,ifm,nope1,nope2,getiface,
     &  jfacem,nelemenm,icounter,idummy,ifac,
     &  icounter2,nslavnode(*),islavnode(*),ict,id,itie,index2,
     &  index1,node2,nelems,islavsurfentry,nfaces
!
      real*8 cdisp(6,*),
     &  weight,dx,help,
     &  ets,xis,etm,xim,xl2s(3,8),xsj2s(3),xsj2s2(3),
     &  shp2s(4,8),xs2s(3,2),xl2m(3,8),xsj2m(3),shp2m(7,8),xs2m(3,2),
     &  contribution(3),pslavdual(16,*),
     &  shp2s2(7,8),xs2s2(3,2),lcoordq(2,8),lcoordt(2,6)
!

      data lcoordq / -1.0, -1.0,
     & 1.0,-1.0,
     & 1.0, 1.0,
     & -1.0, 1.0,
     & 0.0, -1.0,
     & 1.0, 0.0,
     & 0.0, 1.0,
     & -1.0, 0.0 /

      data lcoordt / 0.0, 0.0,
     & 1.0, 0.0,
     & 0.0, 1.0,
     & 0.5, 0.0,
     & 0.5, 0.5,
     & 0.0, 0.5 /

      debug=.false.
c      if(node.eq.4769)debug=.true.
      nfaces=0
      iflag=2
      do j=1,3
       contribution(j) = 0.d0
      enddo
      itie = itie + 1
      index1=iponoels(node)
       if(debug)then
       write(*,*) "*********node",node,"***************"
       endif
      do
         if(index1.eq.0) exit
         nfaces=nfaces+1
         islavsurfentry=inoels(1,index1)
         locs=inoels(2,index1)
         if(debug)then
      WRITE(*,*) "node",node,"islavsurfentry",islavsurfentry,islav+1
         endif
         if((islavsurfentry.lt.itiefac(1,itie)).or.
     &        (islavsurfentry.gt.itiefac(2,itie))) exit
         ifaces=islavsurf(1,islavsurfentry)
         nelems = int(ifaces/10)
         jfaces = ifaces - nelems*10
         indexe = ipkon(nelems)
         call getnumberofnodes(nelems,jfaces,lakon,nope1,
     &      nopes1,idummy)
         do j=1,nope1
            konl(j)=kon(ipkon(nelems)+j)
         enddo
c         do m=1,nopes1
c              ifac=getiface(m,jfaces,nope1)
c              do j=1,3
c                     xl2s(j,m)=co(j,konl(ifac))+
c     &                vold(j,konl(ifac))       
c              enddo
c         enddo
        if(nopes1.eq.4 .or. nopes1.eq.8)then
         xis=lcoordq(1,locs)
         ets=lcoordq(2,locs)
        else
         xis=lcoordq(1,locs)
         ets=lcoordq(2,locs)
        endif
c        write(*,*)"locs",locs, "xi,eta",xis,ets
        ns=islavsurfentry
        if(nopes1.eq.8) then
               call dualshape8q(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
        elseif(nopes1.eq.4) then
               call dualshape4q(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &              pslavdual,iflag)
        elseif(nopes1.eq.6) then
               call dualshape6tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &                    pslavdual,iflag)
        else
               call dualshape3tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &                    pslavdual,iflag)
        endif
         do m=1,nopes1
            ifac=getiface(m,jfaces,nope1)
            node2=konl(ifac)
            call nident(islavnode(nslavnode(itie)+1), node2, 
     &              nslavnode(itie+1)-nslavnode(itie), id)
               index2=nslavnode(itie)+id
            if(debug)then
            write(*,*) "     m",m,"node2",node2,index2
            write(*,*) contribution(1),cdisp(4,index2),shp2s(4,m)
            endif           
            do j=1,3
             contribution(j)=contribution(j)
     &           +cdisp(j+3,index2)*shp2s(4,m)
            enddo
         enddo
         
         index1=inoels(3,index1)
      enddo
      
      do j=1,3
       cdisp(j+3,islav+1)=contribution(j)/nfaces
      enddo
        if(debug)then
        WRITE(*,*) "node",node,"contr_total",cdisp(4,islav+1)
        endif
      itie = itie - 1
      return
      end