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
      subroutine writeinput(inpc,ipoinp,inp,nline,ninp,ipoinpc)
!
      implicit none
!
      character*1 inpc(*)
      character*20 nameref(12)
!
      integer nline,i,j,ninp,ipoinp(2,12),inp(3,ninp),ipoinpc(0:*)
!
      data nameref /'RESTART,READ','NODE','ELEMENT','NSET',
     &              'ELSET','TRANSFORM','MATERIAL','ORIENTATION',
     &              'SECTION','SURFACE','TIE','REST'/
!
      open(16,file='input.inpc',status='unknown')
      do i=1,nline
         write(16,'(1x,i6,1x,1320a1)') i,
     &       (inpc(j),j=ipoinpc(i-1)+1,ipoinpc(i))
      enddo
      close(16)
!
      open(16,file='input.ipoinp',status='unknown')
      do i=1,12
         write(16,'(1x,a20,1x,i6,1x,i6)') nameref(i),(ipoinp(j,i),j=1,2)
      enddo
      close(16)
!
      open(16,file='input.inp',status='unknown')
      do i=1,ninp
         write(16,'(1x,i3,1x,i6,1x,i6,1x,i6)') i,(inp(j,i),j=1,3)
      enddo
      close(16)
!
      return
      end
