!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
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
      subroutine surfacebehaviors(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,npmat_,plicon,nplicon)
!
!     reading the input deck: *SURFACE BEHAVIOR
!
      implicit none
!
      character*1 inpc(*),pressureoverclosure
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ncmat_,irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ntmat,npmat,npmat_,nplicon(0:ntmat_,*)
!
      real*8 elcon(0:ncmat_,ntmat_,*),plicon(0:2*npmat_,ntmat_,*),
     &  temperature
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR reading *SURFACE BEHAVIOR:'
         write(*,*) '       *SURFACE BEHAVIOR should be placed'
         write(*,*) '       before all step definitions'
         stop
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *SURFACE BEHAVIOR:'
         write(*,*) '       *SURFACE BEHAVIOR should be preceded'
         write(*,*) '       by a *SURFACE INTERACTION card'
         stop
      endif
      pressureoverclosure=' '
!
      do i=2,n
         if(textpart(i)(1:27).eq.'PRESSURE-OVERCLOSURE=LINEAR') then
            pressureoverclosure='L'
         elseif(textpart(i)(1:32).eq.'PRESSURE-OVERCLOSURE=EXPONENTIAL') 
     &      then
            pressureoverclosure='E'
         elseif(textpart(i)(1:38).eq.'PRESSURE-OVERCLOSURE=TABULAR') 
     &      then
            pressureoverclosure='T'
         else
            write(*,*) 
     &   '*WARNING reading *SURFACE BEHAVIOR: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
      if(pressureoverclosure.eq.' ') then
         write(*,*) '*ERROR reading *SURFACE BEHAVIOR:'
         write(*,*) '       no PRESSURE-OVERCLOSURE defined on the'
         write(*,*) '       *SURFACE BEHAVIOR card'
         stop
      endif
!
      nelcon(1,nmat)=2
      nelcon(2,nmat)=1
!
      if(pressureoverclosure.eq.'E') then
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR: data'
            write(*,*) '       line is lacking for exponential'
            write(*,*) '       behavior'
            call inputerror(inpc,ipoinpc,iline)
         endif
!     
         elcon(3,1,nmat)=1.5d0
!     
!     exponential overclosure
!     
         do i=1,2
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(i,1,nmat)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         enddo
!     
!     checking the values
!     
         if(elcon(1,1,nmat).le.0.d0) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR: c_0 must'
            write(*,*) '       exceed zero'
            stop
         endif
         if(elcon(2,1,nmat).lt.0.d0) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR: p_0 must'
            write(*,*) '       not be smaller than zero'
            stop
         endif
!     
!     transforming the parameters c_0 into
!     beta such that the pressure p satisfies:
!     p=p_0*dexp(-beta*distance)
!     where n is the normal to the master surface, and
!     distance is the distance between slave node and
!     master surface (negative for penetration)
!     
         elcon(1,1,nmat)=dlog(100.d0)/elcon(1,1,nmat)
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!     
      elseif(pressureoverclosure.eq.'L') then
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR: data'
            write(*,*) '       line is lacking for linear'
            write(*,*) '       behavior'
            call inputerror(inpc,ipoinpc,iline)
         endif
!     
         elcon(3,1,nmat)=2.5d0
!     
!     linear overclosure
!     
!     linear spring stiffness
!     
         read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &        elcon(2,1,nmat)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!     
         if(elcon(2,1,nmat).le.0.d0) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR: K must'
            write(*,*) '       be strictly positive'
            stop
         endif
!     
!     tension at large clearances
!     
         read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &        elcon(1,1,nmat)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!     
!     if nonpositive: error
!     elcon(1,1,nmat)<0 indicates linear spring stiffness
!     
         if(elcon(1,1,nmat).lt.1.d-30) then
            write(*,*) '*ERROR reading *SURFACE BEHAVIOR:'
            write(*,*) '       tension at large clearances'
            write(*,*) '       < 1.e-30'
            call inputerror(inpc,ipoinpc,iline)
            stop
         else
            elcon(1,1,nmat)=-elcon(1,1,nmat)
         endif
!     
!     value of c0coef. If the clearance is inferior to 
!     c0coef*sqrt(slave_area) a contact spring element
!     is generated
!     
         read(textpart(3)(1:20),'(f20.0)',iostat=istat)
     &        elcon(4,1,nmat)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!     
!     default value
!     
         if(elcon(4,1,nmat).le.0.d0) then
            elcon(4,1,nmat)=1.d-3
         endif
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!     
      elseif(pressureoverclosure.eq.'T') then
         elcon(3,1,nmat)=3.5d0
         elcon(4,1,nmat)=1.d-3
         nelcon(1,nmat)=-51
!     
!     tabular
!     
         ntmat=0
         npmat=0
!     
         do
!     
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat)temperature
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!     
!     first temperature
!     
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *SURFACE BEHAVIOR:'
                  write(*,*) '       increase ntmat_'
                  stop
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
!     
!     new temperature
!     
            elseif(plicon(0,ntmat,nmat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *SURFACE BEHAVIOR:' 
                  write(*,*) '       increase ntmat_'
                  stop
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plicon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
               write(*,*) 
     &              '*ERROR reading *SURFACE BEHAVIOR: increase npmat_'
               stop
            endif
            nplicon(ntmat,nmat)=npmat
         enddo
!     
         if(ntmat.eq.0) then
            write(*,*) 
     &           '*ERROR reading *SURFACE BEHAVIOR: *SURFACE BEHAVIOR'
            write(*,*) '       card without data'
            stop
         endif
!
!        check whether the difference between the overclosure data
!        points is at least equal to the smallest difference in 
!        double precision numbers (no vertical slope in the pressure
!        versus overclosure curve allowed)
!
         do i=1,npmat-1
            if(plicon(2*i+2,1,nmat)-plicon(2*i,1,nmat).lt.1.d-10) then
               plicon(2*i+2,1,nmat)=plicon(2*i,1,nmat)+1.d-10
            endif
         enddo
!     
      endif
!     
      elcon(0,1,nmat)=0.d0
!     
      return
      end

