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
      subroutine materialdata_fl(imat,ntmat_,t1l,shcon,nshcon,sph,r,
     &  dvi,rhcon,nrhcon,rho,cocon,ncocon,cond)
!
      implicit none
!
!     determines the following gas properties:
!     specific heat, the dynamic viscosity and the specific gas constant
!
      integer imat,ntmat_,id,nshcon(*),two,four,nrhcon(*),ncocon(2,*),
     &  ncoconst,seven
!
      real*8 t1l,shcon(0:3,ntmat_,*),sph,r,dvi,rhcon(0:1,ntmat_,*),
     &  rho,cocon(0:6,ntmat_,*),cond
!
      two=2
      four=4
      seven=7
!
!     calculating the density (needed for liquids)
!
      call ident2(rhcon(0,1,imat),t1l,nrhcon(imat),two,id)
      if(nrhcon(imat).eq.0) then
         continue
      elseif(nrhcon(imat).eq.1) then
         rho=rhcon(1,1,imat)
      elseif(id.eq.0) then
         rho=rhcon(1,1,imat)
      elseif(id.eq.nrhcon(imat)) then
         rho=rhcon(1,id,imat)
      else
         rho=rhcon(1,id,imat)+
     &        (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &        (t1l-rhcon(0,id,imat))/
     &        (rhcon(0,id+1,imat)-rhcon(0,id,imat))
      endif
!     
!     calculating the specific heat and the dynamic viscosity
!
      call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
      if(nshcon(imat).eq.0) then
         continue
      elseif(nshcon(imat).eq.1) then
         sph=shcon(1,1,imat)
         dvi=shcon(2,1,imat)
      elseif(id.eq.0) then
         sph=shcon(1,1,imat)
         dvi=shcon(2,1,imat)
      elseif(id.eq.nshcon(imat)) then
         sph=shcon(1,id,imat)
         dvi=shcon(2,id,imat)
      else
         sph=shcon(1,id,imat)+
     &        (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
         dvi=shcon(2,id,imat)+
     &        (shcon(2,id+1,imat)-shcon(2,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
      endif
!
!     specific gas constant
!
      r=shcon(3,1,imat)
!
!     calculating the conductivity coefficients
!
      ncoconst=ncocon(1,imat)
      if(ncoconst.ne.1) then
         write(*,*) '*ERROR in materialdata_fl'
         write(*,*) 
     &        '       conductivity for fluids must be isotropic'
         stop
      endif
!     
      call ident2(cocon(0,1,imat),t1l,ncocon(2,imat),seven,id)
      if(ncocon(2,imat).eq.0) then
         cond=0.d0
         continue
      elseif(ncocon(2,imat).eq.1) then
         cond=cocon(1,1,imat)
      elseif(id.eq.0) then
         cond=cocon(1,1,imat)
      elseif(id.eq.ncocon(2,imat)) then
         cond=cocon(1,id,imat)
      else
         cond=(cocon(1,id,imat)+
     &        (cocon(1,id+1,imat)-cocon(1,id,imat))*
     &        (t1l-cocon(0,id,imat))/
     &        (cocon(0,id+1,imat)-cocon(0,id,imat)))
     &        
      endif
!     
      return
      end







