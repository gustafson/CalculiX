!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &  noutloc,description,kode,nmethod,fmat)
!
!     stores the results header in frd format
!
      implicit none
!
      character*8 fmat
      character*12 description
      character*132 text
!
      integer  icounter,noddiam,null,mode,noutloc,kode,nmethod
!
      real*8 oner,time,pi,cs(17,*)
!
      text='    1PSTEP'
      icounter=icounter+1
      write(text(25:36),'(i12)') icounter
      write(7,'(a132)') text
      if(nmethod.eq.2) then
         text='    1PGM'
         write(text(25:36),'(e12.6)') oner
         write(7,'(a132)') text
         text='    1PGK'
         write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
         write(7,'(a132)') text
         text='    1PHID'
         write(text(25:36),'(i12)') noddiam
         write(7,'(a132)') text
         if(noddiam.ge.0) then
            text='    1PAX'
            write(text(25:36),'(1p,e12.5)') cs(6,1)
            write(text(37:48),'(1p,e12.5)') cs(7,1)
            write(text(49:60),'(1p,e12.5)') cs(8,1)
            write(text(61:72),'(1p,e12.5)') cs(9,1)
            write(text(73:84),'(1p,e12.5)') cs(10,1)
            write(text(85:96),'(1p,e12.5)') cs(11,1)
            write(7,'(a132)') text
         endif
         text='    1PSUBC'
         write(text(25:36),'(i12)') null
         write(7,'(a132)') text
         text='    1PMODE'
         write(text(25:36),'(i12)') mode+1
         write(7,'(a132)') text
      endif
!     
      if(abs(nmethod).eq.1) then
         text=
     & '  100CL       .00000E+00                                 0    1'
      elseif(nmethod.eq.2) then
         text=
     & '  100CL       .00000E+00                                 2    1'
      elseif(nmethod.eq.3) then
         text=
     & '  100CL       .00000E+00                                 4    1'
      elseif((nmethod.eq.4).or.(nmethod.eq.5)) then
         text=
     & '  100CL       .00000E+00                                 1    1'
      else
         text=
     & '  100CL       .00000E+00                                 3    1'
      endif
      write(text(25:36),'(i12)') noutloc
      text(37:48)=description
      if(nmethod.eq.2) text(64:68)='MODAL'
      text(75:75)='1'
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(7,'(a132)') text
!     
      return
      end
      
      
      
