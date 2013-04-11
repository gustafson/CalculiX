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
      subroutine surfacebehaviors(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc)
!
!     reading the input deck: *SURFACE BEHAVIOR
!
      implicit none
!
      character*1 inpc(*),pressureoverclosure
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ncmat_,irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in surfacebehaviors:'
         write(*,*) '       *SURFACE BEHAVIOR should be placed'
         write(*,*) '       before all step definitions'
         stop
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR in surfacebehaviors:'
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
         else
            write(*,*) 
     &        '*WARNING in surfacebehaviors: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
      if(pressureoverclosure.eq.' ') then
         write(*,*) '*ERROR in surfacebehaviors:'
         write(*,*) '       no PRESSURE-OVERCLOSURE defined on the'
         write(*,*) '       *SURFACE BEHAVIOR card'
         stop
      endif
!
      nelcon(1,nmat)=2
      nelcon(2,nmat)=1
!
!     no temperature dependence allowed; last line is decisive
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         if(pressureoverclosure.eq.'E') then
!
!           exponential overclosure
!
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &              elcon(i,1,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            enddo
!     
!     checking the values
!     
            if(elcon(1,1,nmat).le.0.d0) then
               write(*,*) '*ERROR in surfacebehaviors: c_0 must'
               write(*,*) '       exceed zero'
               stop
            endif
            if(elcon(2,1,nmat).lt.0.d0) then
               write(*,*) '*ERROR in surfacebehaviors: p_0 must'
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
         else
!
!           linear overclosure
!
!
!           linear spring stiffness
!
            read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &           elcon(2,1,nmat)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
!           pressure at large clearances
!
            read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &           elcon(1,1,nmat)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
!           if nonpositive: default value
!           elcon(1,1,nmat)<0 indicates linear spring stiffness
!
            if(elcon(1,1,nmat).lt.1.d-30) then
               elcon(1,1,nmat)=-10.d0/(4.d0*datan(1.d0))
            else
               elcon(1,1,nmat)=-elcon(1,1,nmat)
            endif
!     
!     checking the values
!     
            if(elcon(2,1,nmat).lt.0.d0) then
               write(*,*) '*ERROR in surfacebehaviors: K must'
               write(*,*) '       not be smaller than zero'
               stop
            endif
         endif
!     
         elcon(0,1,nmat)=0.d0
      enddo
!     
      return
      end

