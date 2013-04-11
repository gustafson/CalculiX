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
      subroutine statics(inpc,textpart,nmethod,iperturb,isolver,istep,
     &  istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,inp,
     &  ithermal,cs,ics,tieset,istartset,
     &  iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &  ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc)
!
!     reading the input deck: *STATIC
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!
      implicit none
!
      character*1 inpc(*)
      character*20 labmpc(*),solver
      character*81 set(*),tieset(3,*)
      character*132 textpart(16)
!
      integer nmethod,iperturb,isolver,istep,istat,n,key,i,idrct,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal,ics(*),
     &  istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,ikmpc(*),ilmpc(*),mpcfree,nset,mcs,ipoinpc(0:*)
!
      real*8 tinc,tper,tmin,tmax,cs(17,*),coefmpc(*)
!
      idrct=0
      tmin=0.d0
      tmax=0.d0
!
      if((iperturb.eq.1).and.(istep.gt.1)) then
         write(*,*) '*ERROR in statics: perturbation analysis is'
         write(*,*) '       not provided in a *STATIC step. Perform'
         write(*,*) '       a genuine nonlinear geometric calculation'
         write(*,*) '       instead (parameter NLGEOM)'
         stop
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in statics: *STATIC can only be used'
         write(*,*) '  within a STEP'
         stop
      endif
!
!     no heat transfer analysis
!
      if(ithermal.gt.1) then
         ithermal=1
      endif
!
!     default solver
!
      if(isolver.eq.0) then
         solver(1:7)='SPOOLES'
      elseif(isolver.eq.2) then
         solver(1:16)='ITERATIVESCALING'
      elseif(isolver.eq.3) then
         solver(1:17)='ITERATIVECHOLESKY'
      elseif(isolver.eq.4) then
         solver(1:3)='SGI'
      elseif(isolver.eq.5) then
         solver(1:5)='TAUCS'
      endif
!
      do i=2,n
         if(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         elseif(textpart(i)(1:6).eq.'DIRECT') then
            idrct=1
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         isolver=2
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         isolver=3
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      else
         write(*,*) '*WARNING in statics: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      nmethod=1
!
!     check for nodes on a cyclic symmetry axis
!
      if(mcs.eq.0) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      else
         n=3
         textpart(2)='NMIN=0
     &
     &                '
         textpart(3)='NMAX=0
     &
     &                '
         nmethod=2
         call selcycsymmods(inpc,textpart,cs,ics,tieset,istartset,
     &        iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &        ilmpc,mpcfree,mcs,set,nset,labmpc,istep,istat,n,iline,
     &        ipol,inl,ipoinp,inp,nmethod,key,ipoinpc)
         nmethod=1
         do i=1,mcs
            cs(2,i)=-0.5
            cs(3,i)=-0.5
         enddo
      endif
!
      if((istat.lt.0).or.(key.eq.1)) then
         if(iperturb.ge.2) then
            write(*,*) '*WARNING in statics: a nonlinear geometric analy
     &sis is requested'
            write(*,*) '         but no time increment nor step is speci
     &fied'
            write(*,*) '         the defaults (1,1) are used'
            tinc=1.d0
            tper=1.d0
            tmin=1.d-5
            tmax=1.d+30
         endif
         return
      endif
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
      if(tper.lt.0.d0) then
         write(*,*) '*ERROR in statics: step size is negative'
         stop
      elseif(tper.le.0.d0) then
         tper=1.d0
      endif
      if(tinc.lt.0.d0) then
         write(*,*) '*ERROR in statics: initial increment size is negati
     &ve'
         stop
      elseif(tinc.le.0.d0) then
         tinc=tper
      endif
      if(tinc.gt.tper) then
         write(*,*) '*ERROR in statics: initial increment size exceeds s
     &tep size'
         stop
      endif
!      
      if(idrct.ne.1) then
         if(dabs(tmin).lt.1.d-10) then
            tmin=min(tinc,1.d-5*tper)
         endif
         if(dabs(tmax).lt.1.d-10) then
            tmax=1.d+30
         endif
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

