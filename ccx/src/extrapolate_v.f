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
      subroutine extrapolate_v(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xboun,sumfix,sumfree,xxn,area)
!
!     inter/extrapolation of v at the center of the elements
!     to the center of the faces
!
      implicit none
!
      integer nface,ielfa(3,*),ifabou(*),iel1,iel2,iel3,i,j,ipointer
!
      real*8 xrlfa(3,*),vel(0:5,*),vfa(0:5,*),xboun(*),xl1,xl2,
     &  sumfix,sumfree,xxn(3,*),area(*)
!  
      sumfix=0.d0
      sumfree=0.d0
!
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            xl2=xrlfa(2,i)
            do j=1,3
               vfa(j,i)=xl1*vel(j,iel1)+xl2*vel(j,iel2)
            enddo
         else
            iel3=ielfa(3,i)
            if(iel2.lt.0) then
               ipointer=-iel2
!
!              global x-direction
!
               if(ifabou(ipointer+1).ne.0) then
!     
!     v_1 given
!     
                  vfa(1,i)=xboun(ifabou(ipointer+1))
                  sumfix=sumfix+vfa(1,i)*area(i)*vfa(5,i)*xxn(1,iel1)
               elseif(ifabou(ipointer+4).ne.0) then
!     
!     p given; linear interpolation
!     
                  vfa(1,i)=xl1*vel(1,iel1)+xrlfa(3,i)*vel(1,iel3)
                  sumfree=sumfree+vfa(1,i)*area(i)*vfa(5,i)*xxn(1,iel1)
               else
!     
!     constant extrapolation
!     
                  vfa(1,i)=vel(1,iel1)
                  sumfree=sumfree+vfa(1,i)*area(i)*vfa(5,i)*xxn(1,iel1)
               endif
!     
!     global y-direction
!     
               if(ifabou(ipointer+2).ne.0) then
!     
!     v_2 given
!     
                  vfa(2,i)=xboun(ifabou(ipointer+2))
                  sumfix=sumfix+vfa(2,i)*area(i)*vfa(5,i)*xxn(2,iel1)
               elseif(ifabou(ipointer+4).ne.0) then
!     
!     p given; linear interpolation
!     
                  vfa(2,i)=xl1*vel(2,iel1)+xrlfa(3,i)*vel(2,iel3)
                  sumfree=sumfree+vfa(2,i)*area(i)*vfa(5,i)*xxn(2,iel1)
               else
!     
!     constant extrapolation
!     
                  vfa(2,i)=vel(2,iel1)
                  sumfree=sumfree+vfa(2,i)*area(i)*vfa(5,i)*xxn(2,iel1)
               endif
!     
!     global z-direction
!     
               if(ifabou(ipointer+3).ne.0) then
!     
!     v_3 given
!     
                  vfa(3,i)=xboun(ifabou(ipointer+3))
                  sumfix=sumfix+vfa(3,i)*area(i)*vfa(5,i)*xxn(3,iel1)
               elseif(ifabou(ipointer+4).ne.0) then
!     
!     p given; linear interpolation
!     
                  vfa(3,i)=xl1*vel(3,iel1)+xrlfa(3,i)*vel(3,iel3)
                  sumfree=sumfree+vfa(3,i)*area(i)*vfa(5,i)*xxn(3,iel1)
               else
!     
!     constant extrapolation
!     
                  vfa(3,i)=vel(3,iel1)
                  sumfree=sumfree+vfa(3,i)*area(i)*vfa(5,i)*xxn(3,iel1)
               endif
!     
            else
!     
!     constant extrapolation
!     
               do j=1,3
                  vfa(j,i)=vel(j,iel1)
                  sumfree=sumfree+vfa(j,i)*area(i)*vfa(5,i)*xxn(j,iel1)
               enddo
            endif
         endif
      enddo
!            
      return
      end
