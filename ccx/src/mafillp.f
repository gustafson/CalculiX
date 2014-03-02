!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine mafillp(ne,lakon,nactdoh,ipnei,neifa,neiel,vfa,area,
     &  adfa,xlet,cosa,volume,au,ad,jq,irow,nzs,ap,ielfa,ifabou,xle,
     &  b,xxn,compressible)
!
!     filling the lhs and rhs to calculate the first correction to the
!     pressure p'
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,ne,jdof1,nactdoh(*),indexf,ipnei(*),j,neifa(*),
     &  neiel(*),iel,ifa,jdof2,irow(*),nzs,ielfa(4,*),compressible,
     &  ifabou(*)
!
      real*8 coef,vfa(0:5,*),volume(*),area(*),adfa(*),xlet(*),
     &  cosa(*),ad(*),au(*),jq(*),xle(*),xxn(3,*),ap(*),b(*)
!
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         jdof1=nactdoh(i)
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,6
!
!              diffusion
!
               indexf=indexf+1
               ifa=neifa(indexf)
               iel=neiel(indexf)
               if(iel.ne.0) then
                  jdof2=nactdoh(iel)
                  coef=vfa(5,ifa)*(volume(i)+volume(iel))*area(ifa)*
     &                  adfa(ifa)/
     &                  (2.d0*xlet(indexf)*cosa(indexf))
                  call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,
     &                   coef,nzs)
                  call add_sm_fl(au,ad,jq,irow,jdof1,jdof2,
     &                   -coef,nzs)
               else
                  if(ielfa(2,ifa).lt.0) then
                     if(ifabou(abs(ielfa(2,ifa))+4).ne.0) then
!
!                       pressure given
!                        
                        coef=vfa(5,ifa)*volume(i)*area(ifa)*
     &                     adfa(ifa)/(2.d0*xle(indexf)*cosa(indexf))
                        call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,
     &                       coef,nzs)
                     endif
                  endif
               endif
!
!              save the coefficient for the rhs of p''
!
               ap(ifa)=coef
!
               b(jdof1)=b(jdof1)-vfa(5,ifa)*area(ifa)+
     &                  (vfa(1,ifa)*xxn(1,indexf)+
     &                   vfa(2,ifa)*xxn(2,indexf)+
     &                   vfa(3,ifa)*xxn(3,indexf))
!
!              convection
!
               if(compressible.eq.1) then
c                  flux=(vfa(1,ifa)*xxn(1,indexf)+
c     &                  vfa(2,ifa)*xxn(2,indexf)+
c     &                  vfa(3,ifa)*xxn(3,indexf))*
c     &                 area(ifa)/(rr*vfa(0,ifa))
c                  if(flux.gt.0.d0) then
c                     call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,
c     &                   flux,nzs)
c                     b(jdof1)=b(jdof1)+(ppfa(ifa)-pp(i))*flux
c                  else
c                     if(iel.ne.0) then
c                     call add_sm_fl(au,ad,jq,irow,jdof1,jdof2,
c     &                   flux,nzs)
c                     b(jdof1)=b(jdof1)+(ppfa(ifa)-pp(iel))*flux
c                  else
c!
c!                    inlet: no correction needed
c!
c                  endif
               endif
            enddo
         endif
      enddo
!     
      return
      end
