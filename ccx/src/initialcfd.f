!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine initialcfd(ne,ipkon,kon,lakon,co,coel,cofa,nface,
     &  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa)
!
!     calculating geometric variables of the cells and their faces
!
      implicit none
!
      character*8 lakon(*)
!
      integer ne,ipkon(*),kon(*),nface,ielfa(4,*),ipnei(*),neiel(*),
     &  ifaceq(8,6),i,j,k,indexe,kflag,index1,index2,iloc1,iloc2,
     &  nodes(4),iel1,iel2,iel3
!
      real*8 co(3,*),coel(3,*),cofa(3,*),area(*),xxn(3,*),xxi(3,*),
     &  xle(*),xlen(*),xlet(*),xrlfa(3,*),cosa(*),xsj2(3),xi,et,
     &  shp2(7,4),xs2(3,7),xl2(3,8),x,y,z,xl13
!
!     nodes belonging to the cell faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     coordinates of the center of the cells
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'F') cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,3
               do k=1,8
                  coel(j,i)=coel(j,i)+co(j,kon(indexe+k))
               enddo
               coel(j,i)=coel(j,i)/8.d0
            enddo
         endif
      enddo
!
      kflag=2
!
!     loop over all faces
!
      do i=1,nface
         iel1=ielfa(1,i)
         indexe=ipkon(iel1)
         iloc1=ielfa(4,i)
         if(lakon(i)(4:4).eq.'8') then
!
!           coordinates of the face centers
!
            do j=1,4
               nodes(j)=kon(indexe+ifaceq(j,iloc1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/4.d0
            enddo
!
            xi=0.d0
            et=0.d0
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
!
!           area of the face
!
            area(i)=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
!
!           normal and xi-vector on face viewed from cell 1
!
            index1=ipnei(iel1)+iloc1
            do k=1,3
               xxn(k,index1)=xsj2(k)/area(i)
               xxi(k,index1)=cofa(k,i)-coel(k,iel1)
            enddo
!
!           distance from face center to the center of cell 1
!
            xle(index1)=dsqrt(xxi(1,index1)**2+xxi(2,index1)**2+
     &                         xxi(3,index1)**2)
            do k=1,3
               xxi(k,index1)=xxi(k,index1)/xle(index1)
            enddo
!
            iel2=ielfa(2,i)
!
!           check whether there is an adjacent cell
!
            if(iel2.ne.0) then
               index2=ipnei(iel2)
               do iloc2=1,6
                  index2=index2+1
                  if(neiel(index2).eq.iel1) exit
               enddo
!
!              normal and xi-vector on face viewed from cell 2
!
               do k=1,3
                  xxi(k,index2)=cofa(k,i)-coel(k,iel2)
                  xxn(k,index2)=-xxn(k,index1)
               enddo
               xle(index1)=dsqrt(xxi(1,index2)**2+xxi(2,index2)**2+
     &              xxi(3,index2)**2)
               do k=1,3
                  xxi(k,index2)=xxi(k,index2)/xle(index2)
               enddo
!
!              distance from the face center to the center of the
!              adjacent cell
!
               xlen(index1)=xle(index2)
               xlen(index2)=xle(index1)
!
               x=coel(1,iel1)-coel(1,iel2)
               y=coel(2,iel1)-coel(2,iel2)
               z=coel(3,iel1)-coel(3,iel2)
!
!              distance between the cell center and the center of the
!              adjacent cell
!
               xlet(index1)=dsqrt(x*x+y*y+z*z)
               xlet(index2)=xlet(index1)
               xrlfa(1,i)=xle(index1)/(xle(index1)+xle(index2))
               xrlfa(2,i)=xle(index2)/(xle(index1)+xle(index2))
!
!              angle between the normal and the xi-vector
!
               cosa(index1)=xxn(1,index1)*xxi(1,index1)+
     &                      xxn(2,index1)*xxi(2,index1)+
     &                      xxn(3,index1)*xxi(3,index1)
               cosa(index2)=xxn(1,index2)*xxi(1,index2)+
     &                      xxn(2,index2)*xxi(2,index2)+
     &                      xxn(3,index2)*xxi(3,index2)
            else
!
!              external face: determining the cell next to the
!              adjacent cell
!
               iel3=ielfa(3,i)
               xl13=dsqrt((coel(1,iel1)-coel(1,iel3))**2+
     &                    (coel(2,iel1)-coel(2,iel3))**2+
     &                    (coel(3,iel1)-coel(3,iel3)))
               xrlfa(1,i)=(xl13+xle(index1))/xl13
               xrlfa(3,i)=1.d0-xrlfa(1,i)
            endif
         endif
      enddo
!     
      return
      end
