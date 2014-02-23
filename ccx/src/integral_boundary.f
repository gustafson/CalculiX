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
      subroutine integral_boundary(sumfix,sumfree,ifaext,nfaext,ielfa,
     &  ifabou,vfa)
!
!     modifying the velocity values along the boundary which are not fixed
!     by boundary conditions such that the mass conservation holds
!
      implicit none
!
      integer i,j,k,nfaext,ifaext(*),ielfa(4,*),ipointer,ifabou(*)
!
      real*8 xi,sumfix,sumfree,vfa(0:5,*)
!
      xi=-sumfix/sumfree
!
!     loop over all external faces
!
      do j=1,nfaext
!
!        i is the row corresponding to the face in field ielfa
!
         i=ifaext(j)
         ipointer=-ielfa(2,i)
         if(ipointer.ne.0) then
            do k=1,3
               if(ifabou(ipointer+k).eq.0) then
                  vfa(k,i)=vfa(k,i)*xi
               endif
            enddo
         else
            do k=1,3
               vfa(k,i)=vfa(k,i)*xi
            enddo
         endif
      enddo
!            
      return
      end
