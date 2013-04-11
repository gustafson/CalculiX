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
      subroutine steps(inpc,textpart,iperturb,iprestr,nbody,
     &  nforc,nload,ithermal,t0,t1,nk,irstrt,istep,istat,n,jmax,ctrl,
     &  iline,ipol,inl,ipoinp,inp,newstep,ipoinpc,physcon)
!
!     reading the input deck: *STEP
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer iperturb(*),nforc,nload,ithermal,nk,istep,istat,n,key,
     &  i,j,iprestr,jmax(2),irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  newstep,nbody,ipoinpc(0:*)
!
      real*8 t0(*),t1(*),ctrl(*),physcon(*)
!
      if(newstep.eq.1) then
         write(*,*) '*ERROR in steps: *STEP statement detected'
         write(*,*) '       within step ',istep
         stop
      else
         newstep=1
      endif
!
      if(iperturb(1).lt.2) iperturb(1)=0
      if(irstrt.lt.0) irstrt=0
      istep=istep+1
      jmax(1)=100
      jmax(2)=10000
      physcon(9)=0.5d0
!
      do i=2,n
         if(textpart(i)(1:12).eq.'PERTURBATION') then
            iperturb(1)=1
            iperturb(2)=1
!
!           removing the present loading (check!!)
!
            nforc=0
            iprestr=0
            if((ithermal.eq.1).or.(ithermal.eq.3)) then
               do j=1,nk
                  t1(j)=t0(j)
               enddo
            endif
!
         elseif((textpart(i)(1:6).eq.'NLGEOM').and.
     &          (textpart(i)(7:9).ne.'=NO')) then
!
!           geometrically nonlinear calculations
!
            iperturb(2)=1
            if(iperturb(1).eq.0) then
               iperturb(1)=2
            elseif(iperturb(1).eq.1) then
               write(*,*) '*ERROR in steps: PERTURBATION and NLGEOM'
               write(*,*) '       are mutually exclusive; '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
c!
c!          to ensure linear calculations for 1d and 2d elements and 
c!          for nonlinear MPCs, the
c!          convergence criteria were set extremely high. If nonlinear
c!          calculations are requested, these criteria must be reset
c!
c            if(ctrl(19).eq.1.d+30) then
c               ctrl(19)=0.005
c               ctrl(20)=0.01
c            endif
!
         elseif(textpart(i)(1:4).eq.'INC=') then
!
!           maximum number of increments
!
            read(textpart(i)(5:14),'(i10)',iostat=istat) jmax(1)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
         elseif(textpart(i)(1:5).eq.'INCF=') then
!
!           maximum number of fluid increments
!
            read(textpart(i)(6:15),'(i10)',iostat=istat) jmax(2)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:16).eq.'TURBULENCEMODEL=') then
!
!           turbulence model
!
            if(textpart(i)(17:25).eq.'NONE') then
               physcon(9)=0.5d0
            elseif(textpart(i)(17:25).eq.'K-EPSILON') then
               physcon(9)=1.5d0
            elseif(textpart(i)(17:23).eq.'K-OMEGA') then
               physcon(9)=2.5d0
            elseif(textpart(i)(17:19).eq.'SST') then
               physcon(9)=3.5d0
            endif
         endif
      enddo
c!
c!     to ensure linear calculations for 1d and 2d elements and 
c!     for nonlinear MPCs, the
c!     convergence criteria were set extremely high. If nonlinear
c!     calculations are requested, these criteria must be reset
c!
c      if((iperturb(1).eq.3).and.(ctrl(19).eq.1.d+30)) then
c         ctrl(19)=0.005
c         ctrl(20)=0.01
c      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end




