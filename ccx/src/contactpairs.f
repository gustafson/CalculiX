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
      subroutine contactpairs(inpc,textpart,tieset,cs,istep,
     &                istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,
     &                iperturb,matname,nmat,ipoinpc,tietol)
!
!     reading the input deck: *CONTACT PAIR
!
      implicit none
!
      character*1 inpc(*)
      character*81 tieset(3,*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,iperturb,nmat,ipoinpc(0:*)
!
      real*8 cs(17,*),tietol(*)
      character*80 matname(*),material
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in contactpairs: *CONTACT PAIR should'
         write(*,*) '  be placed before all step definitions'
         stop
      endif
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR in contactpairs: increase ntie_'
         stop
      endif
      tietol(ntie)=0.d0
!
      do i=2,n
         if(textpart(i)(1:12).eq.'INTERACTION=') then
            material=textpart(i)(13:92)
         elseif(textpart(i)(1:12).eq.'SMALLSLIDING') then
            tietol(ntie)=-1.d0
         elseif(textpart(i)(1:14).eq.'IN-FACESLIDING') then
            tietol(ntie)=-2.d0
         endif
      enddo
!
!     check for the existence of the surface interaction
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR in contactpairs: nonexistent surface'
         write(*,*) '       interaction; '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
      cs(1,ntie)=i+0.5d0
!
      tieset(1,ntie)(81:81)='C'
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR in contactpairs: definition of the '
         write(*,*) '      contact pair is not complete.'
         stop
      endif
!
      tieset(2,ntie)(1:80)=textpart(1)(1:80)
      tieset(2,ntie)(81:81)=' '
      ipos=index(tieset(2,ntie),' ')
      tieset(2,ntie)(ipos:ipos)='S'
!
      tieset(3,ntie)(1:80)=textpart(2)(1:80)
      tieset(3,ntie)(81:81)=' '
      ipos=index(tieset(3,ntie),' ')
      tieset(3,ntie)(ipos:ipos)='T'
!
!     the definition of a contact pair triggers a geometrically 
!     nonlinear calculation
!
      if(iperturb.eq.0) then
         iperturb=2
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



