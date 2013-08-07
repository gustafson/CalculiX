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
      subroutine contactpairs(inpc,textpart,tieset,cs,istep,
     &                istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,
     &                iperturb,matname,nmat,ipoinpc,tietol,set,nset,
     &                mortar)
!
!     reading the input deck: *CONTACT PAIR
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*),material
      character*81 tieset(3,*),noset,set(*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,iperturb(2),nmat,ipoinpc(0:*),nset,j,
     &  mortar
!
      real*8 cs(17,*),tietol(2,*),adjust
!
!     tietol contains information on:
!            - small (tietol<0) or large (tietol>0) sliding
!            - the adjust value (only if dabs(tietol)>=2,
!                 adjust=dabs(tietol)-2
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in contactpairs: *CONTACT PAIR should'
         write(*,*) '  be placed before all step definitions'
         stop
      endif
!
      mortar=0
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR in contactpairs: increase ntie_'
         stop
      endif
      tietol(1,ntie)=1.d0
      do j=1,80
         tieset(1,ntie)(j:j)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:12).eq.'INTERACTION=') then
            material=textpart(i)(13:92)
         elseif(textpart(i)(1:12).eq.'SMALLSLIDING') then
            tietol(1,ntie)=-tietol(1,ntie)
         elseif(textpart(i)(1:7).eq.'ADJUST=') then
            read(textpart(i)(8:25),'(f20.0)',iostat=istat) adjust
            if(istat.gt.0) then
               noset(1:80)=textpart(i)(8:87)
               noset(81:81)=' '
               ipos=index(noset,' ')
               noset(ipos:ipos)='N'
               do j=1,nset
                  if(set(j).eq.noset) exit
               enddo
               if(j.gt.nset) then
                  noset(ipos:ipos)=' '
                  write(*,*) '*ERROR in contactpairs: adjust node set',
     &                   noset
                  write(*,*) '       has not been defined'
                  call inputerror(inpc,ipoinpc,iline)
                  stop
               endif
               do j=1,ipos-1
                  tieset(1,ntie)(j:j)=noset(j:j)
               enddo
               do j=ipos,80
                  tieset(1,ntie)(j:j)=' '
               enddo
            else
               tietol(1,ntie)=dsign(1.d0,tietol(1,ntie))*(2.d0+adjust)
            endif
         elseif(textpart(i)(1:21).eq.'TYPE=SURFACETOSURFACE') then
            mortar=1
         else
            write(*,*) 
     &        '*WARNING in contactpairs: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
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
      tietol(2,ntie)=i+0.5d0
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
!
!     storing the slave surface
!
      if(mortar.eq.1) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='T'
      else
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='S'
      endif
!
      tieset(3,ntie)(1:80)=textpart(2)(1:80)
      tieset(3,ntie)(81:81)=' '
      ipos=index(tieset(3,ntie),' ')
      tieset(3,ntie)(ipos:ipos)='T'
!
!     the definition of a contact pair triggers a call to
!     nonlingeo (for static calculations) but not automatically
!     to the nonlinear calculation of strains (i.e.
!     iperturb(2) should be zero unless NLGEOM is activated)
!
      if(iperturb(1).eq.0) then
         iperturb(1)=2
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



