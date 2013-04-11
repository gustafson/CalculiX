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
      subroutine openfile(jobname,output)
!
      implicit none
!
      logical exi
      character*3 output
      character*132 jobname,fnin,fndat,fnfrd,fnsta,fnonf
      integer i
!
!     opening the input  and output file
!
      do i=1,132
         if(jobname(i:i).eq.' ') exit
      enddo
      i=i-1
      if(i.gt.128) then
         write(*,*) '*ERROR in openfile: input file name is too long:'
         write(*,'(a132)') jobname(1:132)
         write(*,*) '       exceeds 128 characters'
         stop
      endif
!
      fnin=jobname(1:i)//'.inp'
      inquire(file=fnin,exist=exi)
      if(exi) then
         open(1,file=fnin,status='old',err=1)
      else
         write(*,*) '*ERROR in openfile: input file ',fnin
         write(*,*) 'does not exist'
         stop
      endif
!
      fndat=jobname(1:i)//'.dat'
      open(5,file=fndat,status='unknown',err=51)
      close(5,status='delete',err=52)
      open(5,file=fndat,status='unknown',err=51)
c      rewind(5)
!
      if(output.ne.'onf') then
         fnfrd=jobname(1:i)//'.frd'
         open(7,file=fnfrd,status='unknown',err=71)
         close(7,status='delete',err=72)
         open(7,file=fnfrd,status='unknown',err=71)
c         rewind(7)
      endif
!
      fnsta=jobname(1:i)//'.sta'
      open(8,file=fnsta,status='unknown',err=81)
      close(8,status='delete',err=82)
      open(8,file=fnsta,status='unknown',err=81)
c      rewind(8)
      write(8,100)
      write(8,101)
 100  format('SUMMARY OF JOB INFORMATION')
 101  format('  STEP      INC     ATT  ITRS     TOT TIME     STEP TIME      
     &    INC TIME')
!
      if(output.eq.'onf') then
         fnonf=jobname(1:i)//'.onf'
         open(11,file=fnonf,status='unknown',err=111)
         close(11,status='delete',err=112)
         open(11,file=fnonf,status='new',err=111)
      endif
!
      return
!
 1    write(*,*) '*ERROR in openfile: could not open file ',fnin
      stop
 51   write(*,*) '*ERROR in openfile: could not open file ',fndat
      stop
 52   write(*,*) '*ERROR in openfile: could not delete file ',fndat
      stop
 71   write(*,*) '*ERROR in openfile: could not open file ',fnfrd
      stop
 72   write(*,*) '*ERROR in openfile: could not delete file ',fnfrd
      stop
 81   write(*,*) '*ERROR in openfile: could not open file ',fnsta
      stop
 82   write(*,*) '*ERROR in openfile: could not delete file ',fnsta
      stop
 111  write(*,*) '*ERROR in openfile: could not open file ',fnonf
      stop
 112  write(*,*) '*ERROR in openfile: could not delete file ',fnonf
      stop
      end
