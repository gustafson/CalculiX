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
      subroutine elastics(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc)
!
!     reading the input deck: *ELASTIC
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ityp,ncmat_,irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in elastics: *ELASTIC should be placed'
         write(*,*) '  before all step definitions'
         stop
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR in elastics: *ELASTIC should be preceded'
         write(*,*) '  by a *MATERIAL card'
         stop
      endif
!
      ityp=2
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
               ityp=2
            elseif(textpart(i)(6:10).eq.'ORTHO') then
               ityp=9
            elseif(textpart(i)(6:10).eq.'ANISO') then
               ityp=21
            endif
            exit
         endif
      enddo
!
      nelcon(1,nmat)=ityp
!
      if(ityp.eq.2) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR in elastics: increase ntmat_'
               stop
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 elcon(i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
            if(textpart(3)(1:1).ne.' ') then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(0,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      elseif(ityp.eq.9) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR in elastics: increase ntmat_'
               stop
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 elcon(i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) '*ERROR in elastics: orthotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            do i=1,1
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                 elcon(8+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
            if(textpart(2)(1:1).ne.' ') then
               read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &                       elcon(0,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      elseif(ityp.eq.21) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR in elastics: increase ntmat_'
               stop
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) '*ERROR in elastics: anisotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &            elcon(8+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) '*ERROR in elastics: anisotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            do i=1,5
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(16+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
            if(textpart(6)(1:1).ne.' ') then
               read(textpart(6)(1:20),'(f20.0)',iostat=istat)
     &                  elcon(0,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      endif
!
      return
      end

