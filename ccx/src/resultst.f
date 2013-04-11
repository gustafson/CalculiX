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
      subroutine resultst(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc)
!
!     calculates the energy correction (STEP 4) in the nodes
!
      implicit none
!
      integer ipompc(*),nodempc(3,*),nmpc,nk,nactdoh(0:4,*),i,j,ist,
     &  node,ndir,index
!
      real*8 coefmpc(*),sol(*),v(0:4,*),fixed_disp
!
!     extracting the energy correction from the solution
!
      do i=1,nk
         if(nactdoh(0,i).ne.0) then
            v(0,i)=sol(nactdoh(0,i))
         else
            v(0,i)=0.d0
         endif
      enddo
c      write(*,*) 'sol307',v(0,307)
!     
!     inserting the mpc information
!     
      do i=1,nmpc
         ist=ipompc(i)
         node=nodempc(1,ist)
         ndir=nodempc(2,ist)
         if(ndir.ne.0) cycle
         index=nodempc(3,ist)
         fixed_disp=0.d0
         if(index.ne.0) then
            do
               fixed_disp=fixed_disp-coefmpc(index)*
     &              v(nodempc(2,index),nodempc(1,index))
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
         fixed_disp=fixed_disp/coefmpc(ist)
         v(ndir,node)=fixed_disp
      enddo
!
      return
      end
