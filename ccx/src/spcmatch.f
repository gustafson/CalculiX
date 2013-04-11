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
      subroutine spcmatch(xboun,nodeboun,ndirboun,nboun,xbounold,
     &   nodebounold,ndirbounold,nbounold,ikboun,ilboun,vold,reorder,
     &   nreorder)
!
!     matches SPC boundary conditions of one step with those of
!     the previous step
!
      implicit none
!
      integer nodeboun(*),ndirboun(*),nboun,nodebounold(*),ilboun(*),
     &  ndirbounold(*),nbounold,i,kflag,idof,id,nreorder(*),ikboun(*)
!
      real*8 xboun(*),xbounold(*),vold(0:3,*),reorder(*)
!
      kflag=2
!
      do i=1,nboun
         nreorder(i)=0
      enddo
!
c      do i=1,nboun
c         ikboun(i)=7*(nodeboun(i)-1)+ndirboun(i)
c         ilboun(i)=i
c      enddo
c!
c      if(nboun.gt.0) call isortii(ikboun,ilboun,nboun,kflag)
!
      do i=1,nbounold
         idof=7*(nodebounold(i)-1)+ndirbounold(i)
         if(nboun.gt.0) then
            call nident(ikboun,idof,nboun,id)
         else
            id=0
         endif
         if((id.gt.0).and.(ikboun(id).eq.idof)) then
            reorder(ilboun(id))=xbounold(i)
            nreorder(ilboun(id))=1
         endif
      enddo
!
      do i=1,nboun
         if(nreorder(i).eq.0) then
            reorder(i)=vold(ndirboun(i),nodeboun(i))
         endif
      enddo
!
      do i=1,nboun
         xbounold(i)=reorder(i)
      enddo
!
      return
      end

