!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine extrapol_dpel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
     &  xle,xxi,icyclic,xxn,ipnei,ifatie,
     &  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
     &  iflag,xxj,xlet)
!
!     iflag=0: calculate the pressure gradient at the element center
!              only
!     iflag=1: calculated the pressure gradient at the element center
!              and the face center
!
!
!     calculation of the pressure correction gradient at the
!     element centers from the pressure correction values 
!     at the element centers
!
      implicit none
!
      character*20 labmpc(*)
!
      integer nface,ielfa(4,*),ifabou(*),i,iel1,iel2,nef,ibou,
     &  neifa(*),icyclic,ifa,indexf,l,m,ipnei(*),ifatie(*),
     &  is,ie,nmpc,ipompc(*),nodempc(3,*),ifaext(*),nfaext,nactdoh(*),
     &  iel,index,mpc,ipointer,k,ielorig,iface,iflag
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xbounact(*),xl1,xl2,
     &   vfap(0:7,nface),gradpcel(3,*),gradpcfa(3,*),rf(3),area(*),
     &   volume(*),xle(*),xxi(3,*),c(3,3),gradnor,xxn(3,*),coefmpc(*),
     &   coefnorm,sum,xxj(3,*),xlet(*),dd
!
      intent(in) nface,ielfa,xrlfa,vel,
     &  ifabou,xbounact,nef,neifa,rf,area,volume,
     &  xle,xxi,icyclic,xxn,ipnei,ifatie,
     &  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh
!
      intent(inout) vfa,gradpcel,gradpcfa
!     
!     extrapolating the pressure correction from the element
!     centers to the faces
!
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements: interpolation
!
            vfap(4,i)=xl1*vel(iel1,4)+xrlfa(2,i)*vel(iel2,4)
         elseif(ielfa(3,i).ne.0) then
!
!           boundary face; more than one layer
!            
            ibou=0
            if(iel2.lt.0) then
               if(ifabou(-iel2+4).gt.0) then
                  ibou=ifabou(-iel2+4)
               endif
            endif
!
            if(ibou.gt.0) then
!
!              pressure boundary condition
!
               vfap(4,i)=0.d0
            else
!
!              extrapolation
!
               vfap(4,i)=xl1*vel(iel1,4)
     &                 +xrlfa(3,i)*vel(abs(ielfa(3,i)),4)
           endif
         else
!
!           boundary face; one layer
!
            vfap(4,i)=vel(iel1,4)
         endif
      enddo
!
!     Multiple point constraints
!
      if(nmpc.gt.0) then
         do k=1,nfaext
            i=ifaext(k)
            if(ielfa(2,i).ge.0) cycle
            ipointer=-ielfa(2,i)
            if(ifabou(ipointer+4).ge.0) cycle
            mpc=-ifabou(ipointer+4)
!     
            index=ipompc(mpc)
            sum=0.d0
            coefnorm=0.d0
            do
               if(index.eq.0) exit
               if(nodempc(1,index).lt.0) then
!
!                 a negative number refers to a boundary
!                 condition (fields nodeboun, ndirboun..)
!                 resulting from a SPC in local coordinates
!                  
                  sum=sum+0.d0
               else
!
!                 face term
!
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  sum=sum+coefmpc(index)
     &                 *vfa(nodempc(2,index),neifa(ipnei(iel)+iface))
                  coefnorm=coefnorm+coefmpc(index)**2
               endif
               index=nodempc(3,index)
            enddo
!
!           distribute the sum across all terms which are not
!           fixed by boundary conditions
!
            index=ipompc(mpc)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).gt.0) then
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  vfa(nodempc(2,index),neifa(ipnei(iel)+iface))=
     &                 vfa(nodempc(2,index),neifa(ipnei(iel)+iface))-
     &                 sum*coefmpc(index)/coefnorm
               endif
               index=nodempc(3,index)
            enddo
         enddo
      endif
!
!     calculate the gradient of the pressure correction at the center of
!     the elements
!
      do i=1,nef
!
!        initialization
!     
         do l=1,3
            gradpcel(l,i)=0.d0
         enddo
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do l=1,3
               gradpcel(l,i)=gradpcel(l,i)+
     &              vfap(4,ifa)*area(ifa)*xxn(l,indexf)
            enddo
         enddo
!     
!        dividing by the volume of the element
!     
         do l=1,3
            gradpcel(l,i)=gradpcel(l,i)/volume(i)
         enddo
      enddo
! 
!     interpolate/extrapolate the pressure correction gradient from the
!     center of the elements to the center of the faces
!   
      if(iflag.eq.1) then
         do i=1,nface
            iel1=ielfa(1,i)
            xl1=xrlfa(1,i)
            iel2=ielfa(2,i)
            if(iel2.gt.0) then
!
!              face in between two elements
!
               xl2=xrlfa(2,i)
               if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
                  do l=1,3
                     gradpcfa(l,i)=xl1*gradpcel(l,iel1)+
     &                    xl2*gradpcel(l,iel2)
                  enddo
               elseif(ifatie(i).gt.0) then
                  do l=1,3
                     gradpcfa(l,i)=xl1*gradpcel(l,iel1)+xl2*
     &                    (gradpcel(1,iel2)*c(l,1)+
     &                    gradpcel(2,iel2)*c(l,2)+
     &                    gradpcel(3,iel2)*c(l,3))
                  enddo
               else
                  do l=1,3
                     gradpcfa(l,i)=xl1*gradpcel(l,iel1)+xl2*
     &                    (gradpcel(1,iel2)*c(1,l)+
     &                    gradpcel(2,iel2)*c(2,l)+
     &                    gradpcel(3,iel2)*c(3,l))
                  enddo
               endif
            elseif(ielfa(3,i).ne.0) then
!     
!              boundary face; more than one layer; extrapolation
!     
               do l=1,3
                  gradpcfa(l,i)=xl1*gradpcel(l,iel1)+
     &                 xrlfa(3,i)*gradpcel(l,abs(ielfa(3,i)))
               enddo
            else
!     
!              boundary face; one layer
!     
               indexf=ipnei(iel1)+ielfa(4,i)
               gradnor=gradpcel(1,iel1)*xxi(1,indexf)+
     &              gradpcel(2,iel1)*xxi(2,indexf)+
     &              gradpcel(3,iel1)*xxi(3,indexf)
               do l=1,3
                  gradpcfa(l,i)=gradpcel(l,iel1)
     &                 -gradnor*xxi(l,indexf)
               enddo
            endif
         enddo
!
!     correct the facial pressure gradients:
!     Moukalled et al. p 289
!
         do i=1,nface
            iel2=ielfa(2,i)
            if(iel2.gt.0) then
               iel1=ielfa(1,i)
               indexf=ipnei(iel1)+ielfa(4,i)
               dd=(vel(iel2,4)-vel(iel1,4))/xlet(indexf)
     &              -gradpcfa(1,i)*xxj(1,indexf)
     &              -gradpcfa(2,i)*xxj(2,indexf)
     &              -gradpcfa(3,i)*xxj(3,indexf)
               do k=1,3
                  gradpcfa(k,i)=gradpcfa(k,i)+dd*xxj(k,indexf)
               enddo
            endif
         enddo
!
      endif
!            
      return
      end
