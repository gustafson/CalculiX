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
      subroutine fillknotmpc(co,ipompc,nodempc,coefmpc,labmpc,
     &  nmpc,nmpcold)
!
!     updates the coefficients in nonlinear MPC's
!
      implicit none
!
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),irefnode,irotnode,idir,
     &  nmpc,index,ii,inode,nmpcold,iexpnode
!
      real*8 co(3,*),coefmpc(*),e(3,3,3)
!
      data e /0.,0.,0.,0.,0.,-1.,0.,1.,0.,
     &        0.,0.,1.,0.,0.,0.,-1.,0.,0.,
     &        0.,-1.,0.,1.,0.,0.,0.,0.,0./
!
      do ii=nmpcold+1,nmpc
         if(labmpc(ii)(1:4).eq.'KNOT') then
!
!           dependent node
!
            index=ipompc(ii)
            inode=nodempc(1,index)
            idir=nodempc(2,index)
!
!           translation node
!
            index=nodempc(3,index)
            irefnode=nodempc(1,index)
!
!           expansion node
!
            index=nodempc(3,index)
            iexpnode=nodempc(1,index)
            coefmpc(index)=co(idir,irefnode)-co(idir,inode)
!
!           rotation node
!
            index=nodempc(3,index)
            irotnode=nodempc(1,index)
!
!           determining the coefficients of the rotational degrees
!           of freedom
!
            coefmpc(index)=e(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &           e(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &           e(idir,3,1)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=e(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &           e(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &           e(idir,3,2)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=e(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &           e(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &           e(idir,3,3)*(co(3,irefnode)-co(3,inode))
!
         endif
      enddo
!
      return
      end
