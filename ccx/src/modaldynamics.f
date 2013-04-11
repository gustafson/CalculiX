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
      subroutine modaldynamics(inpc,textpart,nmethod,tinc,tper,iexpl,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,iperturb,isolver,
     &  cs,mcs,ipoinpc)
!
!     reading the input deck: *MODAL DYNAMIC
!
      implicit none
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,iexpl,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),iperturb,isolver,i,mcs,ipoinpc(0:*)
!
      real*8 tinc,tper,cs(17,*)
!
      iexpl=0
      iperturb=0
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in modaldynamics: *MODAL DYNAMIC can only'
         write(*,*) '  be used within a STEP'
         stop
      endif
!
!     default solver
!
      if(isolver.eq.0) then
         solver(1:7)='SPOOLES                          '
      elseif(isolver.eq.1) then
         solver(1:7)='PROFILE                          '
      elseif(isolver.eq.2) then
         solver(1:16)='ITERATIVESCALING                '
      elseif(isolver.eq.3) then
         solver(1:17)='ITERATIVECHOLESKY               '
      elseif(isolver.eq.4) then
         solver(1:3)='SGI                              '
      elseif(isolver.eq.5) then
         solver(1:5)='TAUCS                            '
      endif
!
      do i=2,n
         if(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:7).eq.'PROFILE') then
         write(*,*) '*WARNING in modaldynamics: the profile solver is'
         write(*,*) '         not allowed in modal dynamics'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) '*WARNING in modaldynamics: the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) '*WARNING in modaldynamics: the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:13).eq.'MATRIXSTORAGE') then
         isolver=6
      else
         write(*,*) '*WARNING in modaldynamics: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      if((isolver.eq.2).or.(isolver.eq.3)) then
         write(*,*) '*ERROR in modaldynamics: the default solver ',
     & solver
         write(*,*) '       cannot be used for modal dynamic'
         write(*,*) '       calculations '
         stop
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR in modaldynamics: definition not complete'
         write(*,*) '       '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
      read(textpart(1)(1:20),'(f20.0)',iostat=istat)tinc
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      read(textpart(2)(1:20),'(f20.0)',iostat=istat)tper
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
      nmethod=4
!
!     correction for cyclic symmetric structures:
!       if the present step was not preceded by a frequency step
!       no nodal diameter has been selected. To make sure that
!       mastructcs is called instead of mastruct a fictitious
!       minimum nodal diameter is stored
!
      if((mcs.ne.0).and.(cs(2,1)<0.d0)) cs(2,1)=0.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

