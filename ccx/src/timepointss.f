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
      subroutine timepointss(inpc,textpart,amname,amta,namta,nam,
     &  nam_,namtot_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc)
!
!     reading the input deck: *AMPLITUDE
!
      implicit none
!
      character*1 inpc(*)
      character*80 amname(*)
      character*132 textpart(16)
!
      integer namta(3,*),nam,nam_,istep,istat,n,key,i,namtot,
     &  namtot_,irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),ipos,
     &  ipoinpc(0:*)
!
      real*8 amta(2,*),x
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in amplitudes: *AMPLITUDE should be'
         write(*,*) '  placed before all step definitions'
         stop
      endif
!
      nam=nam+1
      if(nam.gt.nam_) then
         write(*,*) '*ERROR in amplitudes: increase nam_'
         stop
      endif
      namta(3,nam)=nam
      amname(nam)='
     &                           '
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            amname(nam)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*)'*ERROR in amplitudes: amplitude name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       amplitude name:',textpart(i)(1:132)
               stop
            endif
         elseif(textpart(i)(1:14).eq.'TIME=TOTALTIME') then
            namta(3,nam)=-nam
         endif
      enddo
!
      if(amname(nam).eq.'                                               
     &                                 ') then
         write(*,*) '*ERROR in amplitudes: Amplitude has no name'
         call inputerror(inpc,ipoinpc,iline)
      endif
!
      if(nam.eq.1) then
         namtot=0
      else
         namtot=namta(2,nam-1)
      endif
      namta(1,nam)=namtot+1
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         do i=1,8
            if(textpart(i)(1:1).ne.' ') then  
               namtot=namtot+1
               if(namtot.gt.namtot_) then
                  write(*,*) '*ERROR in amplitudes: increase namtot_'
                  stop
               endif
               read(textpart(i),'(f20.0)',iostat=istat) x
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               amta(1,namtot)=x
               namta(2,nam)=namtot
            else
               exit
            endif
         enddo
      enddo
!
      if(namta(1,nam).gt.namta(2,nam)) then
         ipos=index(amname(nam),' ')
         write(*,*) '*WARNING in timepointss: *TIME POINTS definition ',
     &        amname(nam)(1:ipos-1)
         write(*,*) '         has no data points'
         nam=nam-1
      endif
!
      return
      end

