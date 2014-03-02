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
      subroutine mafillv(ne,nactdoh,ipnei,neifa,neiel,vfa,xxn,area,
     &  au,ad,jq,irow,nzs,b,vel,cosa,umfa,xlet,xle,gradvfa,xxi,
     &  body,volume,compressible,ielfa,lakon,ifabou,nbody)
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,ne,jdof1,nactdoh(*),indexf,ipnei(*),j,ifa,iel,neifa(*),
     &  neiel(*),jdof2,jq(*),irow(*),nzs,iwall,compressible,ielfa(4,*),
     &  ipointer,ifabou(*),nbody
!
      real*8 flux,vfa(0:5,*),xxn(3,*),area(*),au(*),ad(*),b(ne,3),
     &  vel(0:5,*),cosa(*),umfa(*),xlet(*),xle(*),coef,gradvfa(3,3,*),
     &  xxi(3,*),body(3,*),volume(*),alphav
!
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         jdof1=nactdoh(i)
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,6
!
!              convection
!
               indexf=indexf+1
               ifa=neifa(indexf)
               iel=neiel(indexf)
               jdof2=nactdoh(iel)
               flux=(vfa(1,ifa)*xxn(1,indexf)+
     &               vfa(2,ifa)*xxn(2,indexf)+
     &               vfa(3,ifa)*xxn(3,indexf))
     &               *vfa(5,ifa)*area(ifa)
!
               if(flux.ge.0.d0) then
                  call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,flux,nzs)
                  b(jdof1,1)=b(jdof1,1)-(vfa(1,ifa)-vel(1,i))*flux
                  b(jdof1,2)=b(jdof1,2)-(vfa(2,ifa)-vel(2,i))*flux
                  b(jdof1,3)=b(jdof1,3)-(vfa(3,ifa)-vel(3,i))*flux
               else
                  if(iel.ne.0) then
                     call add_sm_fl(au,ad,jq,irow,jdof1,jdof2,flux,nzs)
                     b(jdof1,1)=b(jdof1,1)-(vfa(1,ifa)-vel(1,iel))*flux
                     b(jdof1,2)=b(jdof1,2)-(vfa(2,ifa)-vel(2,iel))*flux
                     b(jdof1,3)=b(jdof1,3)-(vfa(3,ifa)-vel(3,iel))*flux
                  else
                     b(jdof1,1)=-vfa(1,ifa)*flux
                     b(jdof1,2)=-vfa(2,ifa)*flux
                     b(jdof1,3)=-vfa(3,ifa)*flux
                  endif
               endif
!
!              diffusion
!
               if(iel.ne.0) then
                  coef=umfa(ifa)*area(ifa)/xlet(indexf)
                  call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,coef,nzs)
                  call add_sm_fl(au,ad,jq,irow,jdof1,jdof2,-coef,nzs)
                  b(jdof1,1)=b(jdof1,1)+umfa(ifa)*area(ifa)*
     &             (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &              gradvfa(1,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &              gradvfa(1,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  b(jdof1,2)=b(jdof1,2)+umfa(ifa)*area(ifa)*
     &             (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &              gradvfa(2,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &              gradvfa(2,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  b(jdof1,3)=b(jdof1,3)+umfa(ifa)*area(ifa)*
     &             (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &              gradvfa(3,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &              gradvfa(3,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
               else
                  iwall=0
                  ipointer=abs(ielfa(2,ifa))
                  if(ipointer.gt.0) then
                     iwall=ifabou(ipointer+5)
                  endif
                  if(iwall.eq.0) then
!
!                    boundary, but no wall
!
                     coef=umfa(ifa)*area(ifa)/xle(indexf)
                     call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,coef,nzs)
                     b(jdof1,1)=b(jdof1,1)+coef*vfa(1,ifa)
                     b(jdof1,2)=b(jdof1,2)+coef*vfa(2,ifa)
                     b(jdof1,3)=b(jdof1,3)+coef*vfa(3,ifa)
                     b(jdof1,1)=b(jdof1,1)+umfa(ifa)*area(ifa)*
     &                  (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                  gradvfa(1,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                  gradvfa(1,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                     b(jdof1,2)=b(jdof1,2)+umfa(ifa)*area(ifa)*
     &                  (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                  gradvfa(2,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                  gradvfa(2,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                     b(jdof1,3)=b(jdof1,3)+umfa(ifa)*area(ifa)*
     &                  (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                  gradvfa(3,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                  gradvfa(3,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  else
!     
!                    wall
!     
                     coef=umfa(ifa)*area(ifa)/(xle(indexf)*cosa(indexf))
                     call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,coef,nzs)
                     coef=(vel(1,i)*xxn(1,indexf)+
     &                     vel(2,i)*xxn(2,indexf)+
     &                     vel(3,i)*xxn(3,indexf))*umfa(ifa)*area(ifa)/
     &                     (xle(indexf)*cosa(indexf))
                     b(jdof1,1)=b(jdof1,1)+coef*xxn(1,indexf)
                     b(jdof1,2)=b(jdof1,2)+coef*xxn(2,indexf)
                     b(jdof1,3)=b(jdof1,3)+coef*xxn(3,indexf)
                  endif
               endif
!     
!     compressible
!     
               if(compressible.eq.1) then
                  b(jdof1,1)=b(jdof1,1)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(1,1,ifa)*xxn(1,indexf)+
     &                 gradvfa(2,1,ifa)*xxn(2,indexf)+
     &                 gradvfa(3,1,ifa)*xxn(3,indexf)-
     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
     &                 gradvfa(3,3,ifa))*xxn(1,indexf)/3.d0)
                  b(jdof1,2)=b(jdof1,2)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(1,2,ifa)*xxn(1,indexf)+
     &                 gradvfa(2,2,ifa)*xxn(2,indexf)+
     &                 gradvfa(3,2,ifa)*xxn(3,indexf)-
     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
     &                 gradvfa(3,3,ifa))*xxn(2,indexf)/3.d0)
                  b(jdof1,3)=b(jdof1,3)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(1,3,ifa)*xxn(1,indexf)+
     &                 gradvfa(2,3,ifa)*xxn(2,indexf)+
     &                 gradvfa(3,3,ifa)*xxn(3,indexf)-
     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
     &                 gradvfa(3,3,ifa))*xxn(3,indexf)/3.d0)
               endif
!     
!     pressure
!     
               b(jdof1,1)=b(jdof1,1)-vfa(4,ifa)*xxn(1,indexf)*area(ifa)
               b(jdof1,2)=b(jdof1,2)-vfa(4,ifa)*xxn(2,indexf)*area(ifa)
               b(jdof1,3)=b(jdof1,3)-vfa(4,ifa)*xxn(3,indexf)*area(ifa)
            enddo
!     
!     body force
!     
            if(nbody.gt.0) then
               b(jdof1,1)=b(jdof1,1)+vel(5,i)*body(1,i)*volume(i)
               b(jdof1,2)=b(jdof1,2)+vel(5,i)*body(2,i)*volume(i)
               b(jdof1,3)=b(jdof1,3)+vel(5,i)*body(3,i)*volume(i)
            endif
!     
!     underrelaxation      
!     
            ad(jdof1)=ad(jdof1)/alphav
            b(jdof1,1)=b(jdof1,1)+
     &                 (1.d0-alphav)*ad(jdof1)*vel(1,i)/alphav
            b(jdof1,2)=b(jdof1,2)+
     &                 (1.d0-alphav)*ad(jdof1)*vel(2,i)/alphav
            b(jdof1,3)=b(jdof1,3)+
     &                 (1.d0-alphav)*ad(jdof1)*vel(3,i)/alphav
         endif
      enddo
!     
      return
      end
