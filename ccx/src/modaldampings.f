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
      subroutine modaldampings(inpc,textpart,nmethod,xmodal,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *MODAL DAMPING
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*)
!
      real*8 xmodal(9)
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in modaldampings: *MODAL DAMPING can only'
         write(*,*) '  be used within a STEP'
         stop
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR in modaldampings: definition not complete'
         write(*,*) '       '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmodal(1)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) xmodal(2)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

