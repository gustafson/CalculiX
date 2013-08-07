!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
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
      subroutine checktriaedge(node1,node2,ipe,ime,iactiveline,
     &  nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &  ifreeintersec,xn,co,nopes,xl2,itri,idin,vold,mi)
!
!     check whether triangular master edge cuts the slave surface
!     edges
!
      implicit none
!
      logical invert,active
!
      integer node1,node2,ipe(*),ime(4,*),indexl,iactiveline(3,*),
     &  nactiveline,id,indexi,intersec(2,*),nvertex,lvertex(13),i,
     &  ifreeintersec,nopes,nintersec,j,itri,index1,index2,
     &  k,node,ithree,idin,mi(*)
!
      real*8 xntersec(3,*),pvertex(3,*),pr(3),xm(3),xn(3),co(3,*),
     &  dd,xl2(3,*),rc(3),dc(3),al,ratio(8),dist,xil,
     &  etl,al2,inter(3),err,vold(0:mi(2),*)
!
      data ithree /3/
!
!     check whether the first node of the edge has a lower number
!     than the second node. If not, the line is stored in reverse
!     order in field ime and the invert flag is set to true
!
      err=1d-6
      invert=.false.
      if(node2.lt.node1) then
         node=node1
         node1=node2
         node2=node
         invert=.true.
      endif
!
!     retrieving the number of the line in field ime: indexl
!
      indexl=ipe(node1)
      do
         if(ime(1,indexl).eq.node2) exit
         indexl=ime(4,indexl)
         if(indexl.eq.0) then
            write(*,*) '*ERROR in checktriaedge: line was not'
            write(*,*) itri,"node1",node1, "node2",node2
            write(*,*) '       properly catalogued'
            stop
         endif
      enddo
!
!     check whether line is active (i.e. lies on the progressing
!     front)
!
      active=.false.
      call nidentk(iactiveline,indexl,nactiveline,id,ithree)
      if(id.gt.0) then
         if(iactiveline(1,id).eq.indexl) then
            active=.true.
         endif
      endif
!
      if(active) then
!
!        retrieving the intersection points and storing them in 
!        pvertex...
!         
         indexi=iactiveline(3,id)
!
!        check whether there is at least one intersection
!
         if(indexi.gt.0) then
            nvertex=nvertex+1
            lvertex(nvertex)=intersec(1,indexi)
            do i=1,3
               pvertex(i,nvertex)=xntersec(i,indexi)
            enddo
            indexi=intersec(2,indexi)
!
!           check whether there is a second intersection
!
            if(indexi.ne.0) then
               nvertex=nvertex+1
!
!              for two intersections the orientation of the line
!              is important
!
               if(invert) then
                  lvertex(nvertex)=lvertex(nvertex-1)
                  lvertex(nvertex-1)=intersec(1,indexi)
                  do i=1,3
                     pvertex(i,nvertex)=pvertex(i,nvertex-1)
                     pvertex(i,nvertex-1)=xntersec(i,indexi)
                  enddo
               else
                  lvertex(nvertex)=intersec(1,indexi)
                  do i=1,3
                     pvertex(i,nvertex)=xntersec(i,indexi)
                  enddo
               endif
            endif
         endif
!
!        remove the line from the active stack
!
c!        restore intersec/ifreeintersec
c!
c         indexi=iactiveline(3,id)
c         do
c!
c! Inversion
c!
c            idummy=indexi
c            indexi=intersec(2,indexi)
c            intersec(2,idummy)=0
c            if(indexi.eq.0) exit
c         enddo
!
!        restore iactiveline/nactiveline
!
         nactiveline=nactiveline-1
         do i=id,nactiveline
            do k=1,3
               iactiveline(k,i)=iactiveline(k,i+1)
            enddo
         enddo
      else
!
!        line was not active: check for intersections
!
         do i=1,3
             pr(i)=co(i,node2)+vold(i,node2)-
     &         co(i,node1)-vold(i,node1)
         enddo
!
!        normal on a plane through the line and vector xn
!
         xm(1)=xn(2)*pr(3)-xn(3)*pr(2)
         xm(2)=xn(3)*pr(1)-xn(1)*pr(3)
         xm(3)=xn(1)*pr(2)-xn(2)*pr(1)
         dd=dsqrt(xm(1)**2+xm(2)**2+xm(3)**2)
         do i=1,3
            xm(i)=xm(i)/dd
         enddo
!
!        check for intersections with the slave edges
!
         nintersec=0
         do j=1,nopes
            if(j.ne.nopes) then
               do i=1,3
                  rc(i)=co(i,node1)+vold(i,node1)-xl2(i,j)
                  dc(i)=xl2(i,j+1)-xl2(i,j)
               enddo
            else
               do i=1,3
                  rc(i)=co(i,node1)+vold(i,node1)-xl2(i,j)
                  dc(i)=xl2(i,1)-xl2(i,j)
               enddo
            endif
            al=(xm(1)*rc(1)+xm(2)*rc(2)+xm(3)*rc(3))/
     &         (xm(1)*dc(1)+xm(2)*dc(2)+xm(3)*dc(3))
!
!           the intersection point must lie in between the
!           triangular vertices
!
            if((al.ge.1.d0).or.(al.le.0.d0)) cycle
!           intersection found: catalogueing the line as active
!           and storing the intersection
!
            do i=1,3
               inter(i)=xl2(i,j)+al*dc(i)
            enddo
!
!
!
            al2=(((inter(1)-co(1,node1)-vold(1,node1))*pr(1)+
     &        (inter(2)-co(2,node1)-vold(2,node1))*pr(2)+
     &        (inter(3)-co(3,node1)-vold(3,node1))*pr(3))-
     &        (pr(1)*xn(1)+pr(2)*xn(2)+pr(3)*xn(3))*
     &        ((inter(1)-co(1,node1)-vold(1,node1))*xn(1)+
     &        (inter(2)-co(2,node1)-vold(2,node1))*xn(2)+
     &        (inter(3)-co(3,node1)-vold(3,node1))*xn(3)))/
     &        ((pr(1)**2+pr(2)**2+pr(3)**2)-
     &    (pr(1)*xn(1)+pr(2)*xn(2)+pr(3)*xn(3))**2)
!
            if((al2.ge.1.0d0).or.(al2.le.0.0d0)) cycle 
!
            if(nintersec.eq.0) then
               nactiveline=nactiveline+1
               ifreeintersec=ifreeintersec+1
               do k=nactiveline,id+2,-1
                  do i=1,3
                     iactiveline(i,k)=iactiveline(i,k-1)
                  enddo
               enddo
               iactiveline(1,id+1)=indexl
               iactiveline(2,id+1)=itri
               iactiveline(3,id+1)=ifreeintersec
               nintersec=nintersec+1
            elseif(nintersec.eq.1) then
               ifreeintersec=ifreeintersec+1
               intersec(2,iactiveline(3,id+1))=ifreeintersec
               nintersec=nintersec+1
            else
               write(*,*) '*ERROR in checktriaedge: no more'
               write(*,*) '       than two intersections allowed'
               stop
            endif
!
!           update intersec and xntersec
!
            intersec(1,ifreeintersec)=j
            do i=1,3
               xntersec(i,ifreeintersec)=inter(i)
            enddo
            call attachline(xl2,xntersec(1,ifreeintersec),nopes,
     &           ratio,dist,xil,etl,xn)
c            ifreeintersec=intersec(2,ifreeintersec)
            intersec(2,ifreeintersec)=0
         enddo
!
!        if there are two intersections, their order has to be
!        checked
!                  
         if(nintersec.eq.2) then
!
!           check order of crossings
!            
            index1=iactiveline(3,id+1)
            index2=intersec(2,index1)
!
!           measuring the distance from node1
!
            if(((xntersec(1,index1)-co(1,node1)-vold(1,node1))**2+
     &           (xntersec(2,index1)-co(2,node1)-vold(2,node1))**2+
     &           (xntersec(3,index1)-co(3,node1)-vold(3,node1))**2).gt.
     &           ((xntersec(1,index2)-co(1,node1)-vold(1,node1))**2+
     &           (xntersec(2,index2)-co(2,node1)-vold(2,node1))**2+
     &           (xntersec(3,index2)-co(3,node1)-vold(3,node1))**2) )
     &               then
!     
               	   iactiveline(3,id+1)=index2
                   intersec(2,index2)=index1
                   intersec(2,index1)=0
            endif
         endif
!     
c         indexi=iactiveline(3,id+1)
c         if((indexi.gt.0).and.(nintersec.gt.0)) then
         if(nintersec.gt.0) then
            indexi=iactiveline(3,id+1)
            nvertex=nvertex+1
            lvertex(nvertex)=intersec(1,indexi)
            do i=1,3
               pvertex(i,nvertex)=xntersec(i,indexi)
            enddo
            indexi=intersec(2,indexi)
!     
!     check whether there is a second intersection
!     
            if((indexi.ne.0)) then
               nvertex=nvertex+1
!     
!     for two intersections the orientation of the line
!     is important
!     
               if(invert) then
                  lvertex(nvertex)=lvertex(nvertex-1)
                  lvertex(nvertex-1)=intersec(1,indexi)
                  do i=1,3
                     pvertex(i,nvertex)=pvertex(i,nvertex-1)
                     pvertex(i,nvertex-1)=xntersec(i,indexi)
                  enddo
               else
                  lvertex(nvertex)=intersec(1,indexi)
                  do i=1,3
                     pvertex(i,nvertex)=xntersec(i,indexi)
                  enddo
               endif
            endif
         endif
!     
!        if there are no intersections the line has to be set
!        active if node1 lies inside
!
         if((idin.gt.0).and.(nintersec.eq.0)) then
            nactiveline=nactiveline+1
            do k=nactiveline,id+2,-1
               do i=1,3
                  iactiveline(i,k)=iactiveline(i,k-1)
               enddo
            enddo
            iactiveline(1,id+1)=indexl
            iactiveline(2,id+1)=itri
            iactiveline(3,id+1)=0
         endif
      endif
! 
      if(invert) then
         node=node1
         node1=node2
         node2=node
      endif
      return
      end
