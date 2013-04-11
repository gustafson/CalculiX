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
      subroutine loadaddt(nelement,label,valfilm,valtemp,nelemload,
     &  sideload,xload,nload,nload_,iamload,iamptemp,
     &  iampfilm,nam,node)
!
!     adds a thermal dload condition to the data base
!
      implicit none
!
      character*20 label,sideload(*)
!
      integer nelemload(2,*),iamload(2,*),
     &  nelement,nload,nload_,j,iamptemp,nam,iampfilm,node
!
      real*8 xload(2,*),valfilm,valtemp
!
      do j=1,nload
         if((nelemload(1,j).eq.nelement).and.
     &      (sideload(j).eq.label)) then
            xload(1,j)=valfilm
            xload(2,j)=valtemp
            nelemload(2,j)=node
            if(nam.gt.0) then
               iamload(1,j)=iampfilm
               iamload(2,j)=iamptemp
            endif
            return
         endif
      enddo
      nload=nload+1
      if(nload.gt.nload_) then
         write(*,*) '*ERROR in loadadd: increase nload_'
         stop
      endif
      nelemload(1,nload)=nelement
      sideload(nload)=label
      xload(1,nload)=valfilm
      xload(2,nload)=valtemp
      nelemload(2,nload)=node
      if(nam.gt.0) then
         iamload(1,nload)=iampfilm
         iamload(2,nload)=iamptemp
      endif
!
      return
      end

