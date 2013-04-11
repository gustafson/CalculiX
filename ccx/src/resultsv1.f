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
      subroutine resultsv1(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc)
!
!     calculates the velocity correction (STEP 1) in the nodes
!
      implicit none
!
      integer ipompc(*),nodempc(3,*),nmpc,nk,nactdoh(0:4,*),i,j,ist,
     &  node,ndir,index
!
      real*8 coefmpc(*),sol(*),v(0:4,*),fixed_disp
!
!     extracting the 1st velocity correction from the solution (STEP 1)
!
      do i=1,nk
         do j=1,3
            if(nactdoh(j,i).ne.0) then
               v(j,i)=sol(nactdoh(j,i))
            else
               v(j,i)=0.d0
            endif
         enddo
c         write(*,*) 'sollll ',i,(v(j,i),j=1,3)
      enddo
c      write(*,*) 'sol307',v(1,307),v(2,307),v(3,307)
!     
!     inserting the mpc information
!     
      do i=1,nmpc
         ist=ipompc(i)
         node=nodempc(1,ist)
         ndir=nodempc(2,ist)
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
