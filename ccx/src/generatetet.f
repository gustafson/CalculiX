!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine generatetet(kontet,ifatet,ielement,inodfa,
     &  ifreefa,planfa,ipofa,nodes,cotet)
!
      implicit none
!
      integer nodes(4),nodef(3),kontet(4,*),ifatet(4,*),inodfa(4,*),
     &  ipofa(*),ifreetet,ifreefa,ifree,ig(3,4),j
!
      integer n1,n2,n3,n4,nxa,nxb,nxc,nxd,nya,nyb,nyc,nyd,nza,nzb,nzc,
     &  nzd,nx1,nx2,ny1,ny2,nz1,nz2,index,j1,j2,j3,i,n,kflag,idum,
     &  node,indexold,ielement,ndx,ndy,ndz
!
      real*8 planfa(4,*),cotet(3,*),dd
!
c      real*8 p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z,p4x,p4y,p4z,a11,
c     &  a12,a13,a21,a22,a23,a31,a32,a33,d1,d2,d3,det,x,y,z,r,xmin,
c     &  ymin,zmin,charlen,dd
!
      data ig /2,3,4,3,4,1,4,1,2,1,2,3/
c      data ig /2,3,4,1,3,4,1,2,4,1,2,3/
!
c      n=4
c      kflag=1
c      call isortii(nodes,idum,n,kflag)
!
!     updating the node per element relationship
!
      do i=1,4
         kontet(i,ielement)=nodes(i)
      enddo
!
!     creating faces
!
      do i=1,4
         nodef(1)=nodes(ig(1,i))
         nodef(2)=nodes(ig(2,i))
         nodef(3)=nodes(ig(3,i))
!
         n=3
         kflag=1
         call isortii(nodef,idum,n,kflag)
!
!        check whether face already exists
!
         node=nodef(1)
         index=ipofa(node)
!
         do
            if(index.eq.0) exit
            if((inodfa(2,index).eq.nodef(2)).and.
     &         (inodfa(3,index).eq.nodef(3))) exit
            indexold=index
            index=inodfa(4,index)
         enddo
!
         if(index.eq.0) then
            index=ifreefa
            ifreefa=inodfa(4,ifreefa)
            if(ifreefa.eq.0) then
               write(*,*) '*ERROR in generatet: increase the dimension'
               write(*,*) '       of inodfa'
            endif
            inodfa(1,index)=nodef(1)
            inodfa(2,index)=nodef(2)
            inodfa(3,index)=nodef(3)
            inodfa(4,index)=0
            if(ipofa(node).eq.0) then
               ipofa(node)=index
            else
               inodfa(4,indexold)=index
            endif
!
            call planeeq(cotet,nodef,planfa(1,index))
c            write(*,*) index,(planfa(j,index),j=1,4)
!
         endif
!
!        the face number in ifatet is negative, if the equation
!        of the face plane is such, that its value in the 
!        remaining node of the tetrahedron is negative
!
         dd=planfa(1,index)*cotet(1,nodes(i))+
     &      planfa(2,index)*cotet(2,nodes(i))+
     &      planfa(3,index)*cotet(3,nodes(i))+
     &      planfa(4,index)
         if(dabs(dd).lt.1.d-10) then
            write(*,*) '*WARNING in generatetet: element ',ielement
            write(*,*) '         is extremely flat'
            write(*,*) '         the element is deleted'
            ielement=ielement-1
            return
c            write(*,*) 'dd= ',dd
         endif
         if(dd.ge.0.d0) then
            ifatet(i,ielement)=index
         else
            ifatet(i,ielement)=-index
         endif
      enddo
!
!     finding the center and radius of the circumscribed sphere
!
c      p1x=cotet(1,nodes(1))
c      p1y=cotet(2,nodes(1))
c      p1z=cotet(3,nodes(1))
c!
c      p2x=cotet(1,nodes(2))
c      p2y=cotet(2,nodes(2))
c      p2z=cotet(3,nodes(2))
c!
c      p3x=cotet(1,nodes(3))
c      p3y=cotet(2,nodes(3))
c      p3z=cotet(3,nodes(3))
c!
c      p4x=cotet(1,nodes(4))
c      p4y=cotet(2,nodes(4))
c      p4z=cotet(3,nodes(4))
c!
c      a11=p1x-p2x
c      a12=p1y-p2y
c      a13=p1z-p2z
c      d1=((p1x+p2x)*a11+(p1y+p2y)*a12+(p1z+p2z)*a13)/2.d0
c!
c      a21=p1x-p3x
c      a22=p1y-p3y
c      a23=p1z-p3z
c      d2=((p1x+p3x)*a21+(p1y+p3y)*a22+(p1z+p3z)*a23)/2.d0
c!
c      a31=p1x-p4x
c      a32=p1y-p4y
c      a33=p1z-p4z
c      d3=((p1x+p4x)*a31+(p1y+p4y)*a32+(p1z+p4z)*a33)/2.d0
c!
c      det=a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+
c     &     a13*(a21*a32-a31*a22)
c      x=(d1*(a22*a33-a23*a32)-d2*(a12*a33-a32*a13)+
c     &     d3*(a12*a23-a22*a13))/det
c      y=(-d1*(a21*a33-a31*a23)+d2*(a11*a33-a31*a13)-
c     &     d3*(a11*a23-a21*a13))/det
c      z=(d1*(a21*a32-a31*a22)-d2*(a11*a32-a31*a12)+
c     &       d3*(a11*a22-a21*a12))/det
!
      return
      end
