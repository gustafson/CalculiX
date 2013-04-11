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
      subroutine materialdata_tg_sec(imat,ntmat_,t1l,shcon,nshcon,sph,r,
     &  dvi,rhcon,nrhcon,rho,physcon)
!
      implicit none
!
!     determines the following gas properties:
!     specific heat, the dynamic viscosity and the specific gas constant
!
!     the difference with materialdata_tg is that the specific gas
!     constant cp is here the secant value and not the differential
!     value. For the differential value we have:
!            dh=cp*dT
!     and consequently
!            h=int_from_0_to_T cp*dT cp*dT
!     For the secant value one has:
!            h=cp_secant*T
!
      integer imat,ntmat_,id,nshcon(*),two,four,nrhcon(*),i
!
      real*8 t1l,shcon(0:3,ntmat_,*),sph,r,dvi,rhcon(0:1,ntmat_,*),
     &  rho,physcon(3)
!
      two=2
      four=4
!
!     calculating the density (needed for liquids)
!
      call ident2(rhcon(0,1,imat),t1l,nrhcon(imat),two,id)
      if(nrhcon(imat).eq.0) then
         rho=0.d0
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
         sph=(shcon(0,1,imat)-physcon(1))*shcon(1,1,imat)
         do i=2,nshcon(imat)
            sph=sph+(shcon(0,i,imat)-shcon(0,i-1,imat))*
     &              (shcon(1,i,imat)+shcon(1,i-1,imat))/2.d0
         enddo
         sph=sph+(t1l-shcon(0,nshcon(imat),imat))*
     &           (shcon(1,nshcon(imat),imat))/(t1l-physcon(1))
         dvi=shcon(2,id,imat)
      else
         sph=shcon(1,id,imat)+
     &        (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
         sph=(t1l-shcon(0,id,imat))*(sph+shcon(1,id,imat))/2.d0
         do i=2,id
            sph=sph+(shcon(0,i,imat)-shcon(0,i-1,imat))*
     &              (shcon(1,i,imat)+shcon(1,i-1,imat))/2.d0
         enddo
         sph=sph+(shcon(0,1,imat)-physcon(1))*shcon(1,1,imat)
         sph=sph/(t1l-physcon(1))
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
      return
      end







