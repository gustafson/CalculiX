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
      subroutine dashpots(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &        plicon,nplicon,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &        nset,ielmat,ielorien,ipoinpc)
!
!     reading the input deck: *DASHPOT
!
      implicit none
!
      logical linear
!
      character*1 inpc(*)
      character*80 matname(*)
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &  n,key,i,nplicon(0:ntmat_,*),ncmat_,istat,istartset(*),
     &  iendset(*),irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),nmat_,
     &  ialset(*),ipos,nset,j,k,ielmat(*),ielorien(*),ipoinpc(0:*)  
!
      real*8 plicon(0:2*npmat_,ntmat_,*),
     &  elcon(0:ncmat_,ntmat_,*)
!
      linear=.true.
!
      ntmat=0
      npmat=0
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in springs: *SPRING should be placed'
         write(*,*) '  before all step definitions'
         stop
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR in materials: increase nmat_'
         stop
      endif
      matname(nmat)(1:7)='DASHPOT'
      do i=8,80
         matname(nmat)(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:9).eq.'NONLINEAR') then
c            linear=.false.
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         endif
      enddo
!
      if(linear) then
         nelcon(1,nmat)=2
!
!        linear spring
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR in springs: increase ntmat_'
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
      else
c         nelcon(1,nmat)=-51
c!
c!        kinematic hardening coefficients
c!
c         do
c            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
c     &           ipoinp,inp,ipoinpc)
c            if((istat.lt.0).or.(key.eq.1)) exit
c            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
c            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
c!
c!           first temperature
c!
c            if(ntmat.eq.0) then
c               npmat=0
c               ntmat=ntmat+1
c               if(ntmat.gt.ntmat_) then
c                  write(*,*) '*ERROR in springs: increase ntmat_'
c                  stop
c               endif
c               nplicon(0,nmat)=ntmat
c               plicon(0,ntmat,nmat)=temperature
c!
c!           new temperature
c!
c            elseif(plicon(0,ntmat,nmat).ne.temperature) then
c               npmat=0
c               ntmat=ntmat+1
c               if(ntmat.gt.ntmat_) then
c                  write(*,*) '*ERROR in springs: increase ntmat_'
c                  stop
c               endif
c               nplicon(0,nmat)=ntmat
c               plicon(0,ntmat,nmat)=temperature
c            endif
c            do i=1,2
c               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
c     &              plicon(2*npmat+i,ntmat,nmat)
c               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
c            enddo
c            npmat=npmat+1
c            if(npmat.gt.npmat_) then
c               write(*,*) '*ERROR in springs: increase npmat_'
c               stop
c            endif
c            nplicon(ntmat,nmat)=npmat
c         enddo
      endif
!
      if(ntmat.eq.0) then
         write(*,*) '*ERROR in springs: *SPRING card without data'
         stop
      endif
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR in springs: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
!
!     assigning the elements of the set the appropriate material
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            ielmat(ialset(j))=nmat
            ielorien(ialset(j))=0
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               ielmat(k)=nmat
               ielorien(k)=0
            enddo
         endif
      enddo
!
      return
      end

