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
      subroutine modelchanges(inpc,textpart,tieset,istat,n,iline,
     &           ipol,inl,ipoinp,inp,ntie,ipoinpc,istep)
!
!     reading the input deck: *MODEL CHANGE
!
      implicit none
!
      logical contactpair,add,remove 
!
      character*1 inpc(*)
      character*81 tieset(3,*)
      character*132 textpart(16)
!
      integer istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ipoinpc(0:*),iposslave,iposmaster,itie,istep
!
      if(istep.eq.0) then
         write(*,*) '*ERROR reading *MODEL CHANGE: *MODEL CHANGE'
         write(*,*) '       cannot be used before the first step'
         stop
      endif
!
      contactpair=.false.
      add=.false.
      remove=.false.
!
      do i=2,n
         if(textpart(i)(1:16).eq.'TYPE=CONTACTPAIR') then
            contactpair=.true.
         elseif(textpart(i)(1:3).eq.'ADD') then
            add=.true.
         elseif(textpart(i)(1:6).eq.'REMOVE') then
            remove=.true.
         else
            write(*,*) 
     &       '*WARNING reading *MODEL CHANGE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
!     checking the validity of the input
!
      if(.not.contactpair) then
         write(*,*) '*ERROR reading *MODEL CHANGE: model change can'
         write(*,*) '       only be used for contact pairs'
         stop
      endif
!
      if((.not.add).and.(.not.remove)) then
         write(*,*) '*ERROR reading *MODEL CHANGE: at least ADD or'
         write(*,*) '        REMOVE has to be selected'
         stop
      endif
!
      if(add.and.remove) then
         write(*,*) '*ERROR reading *MODEL CHANGE: ADD and REMOVE'
         write(*,*) '       cannot both be selected'
         stop
      endif
!
!     reading the slave and the master surface
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR reading *MODEL CHANGE: definition of the '
         write(*,*) '      contact pair is not complete.'
         stop
      endif
!
!     selecting the appropriate action
!
      iposslave=index(textpart(1)(1:80),' ')
      iposmaster=index(textpart(2)(1:80),' ')
      do i=1,ntie
         if((tieset(1,i)(81:81).ne.'C').and.
     &      (tieset(1,i)(81:81).ne.'-')) cycle
         ipos=index(tieset(2,i),' ')-1
         if(ipos.ne.iposslave) cycle
         if(tieset(2,i)(1:ipos-1).ne.textpart(1)(1:ipos-1)) cycle
         ipos=index(tieset(3,i),' ')-1
         if(ipos.ne.iposmaster) cycle
         if(tieset(3,i)(1:ipos-1).ne.textpart(2)(1:ipos-1)) cycle
         itie=i
         exit
      enddo
!
      if(add) then
         tieset(1,itie)(81:81)='C'
      else
         tieset(1,itie)(81:81)='-'
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



