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
      subroutine hrv_mod_smart(nface,ielfa,vel,gradvel,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef,flux,vfa)
!
!     use the modified smart scheme to determine the facial
!     velocity
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:7),gradvel(3,3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*),vfa(0:7,*)
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,ipnei,vel,vfa,flux,gradvel,xxj,xlet)
c$omp& private(i,iel2,iel1,j,indexf,vcd,vud,phic)
c$omp do
      do i=1,nface
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
            do j=1,3
               vcd=vel(iel1,j)-vel(iel2,j)
               if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,j))) vcd=0.d0
!     
               vud=2.d0*xlet(indexf)*
     &              (gradvel(j,1,iel1)*xxj(1,indexf)+
     &              gradvel(j,2,iel1)*xxj(2,indexf)+
     &              gradvel(j,3,iel1)*xxj(3,indexf))
!     
               if(dabs(vud).lt.1.d-20) then
!     
!     upwind difference
!     
                  vfa(j,i)=vel(iel1,j)
                  cycle
               endif
!     
               phic=1.d0+vcd/vud
!     
               if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!     
!     upwind difference
!     
                  vfa(j,i)=vel(iel1,j)
               elseif(phic.le.1.d0/6.d0) then
                  vfa(j,i)=3.d0*vel(iel1,j)-2.d0*vel(iel2,j)+2.d0*vud
               elseif(phic.le.0.7d0) then
                  vfa(j,i)=3.d0*vel(iel1,j)/4.d0+vel(iel2,j)/4.d0
     &                    +vud/8.d0
               else
                  vfa(j,i)=vel(iel1,j)/3.d0+2.d0*vel(iel2,j)/3.d0
               endif
            enddo
         else
!
            do j=1,3
               vcd=vel(iel2,j)-vel(iel1,j)
               if(dabs(vcd).lt.1.d-3*dabs(vel(iel2,j))) vcd=0.d0
!     
               vud=-2.d0*xlet(indexf)*
     &              (gradvel(j,1,iel2)*xxj(1,indexf)+
     &              gradvel(j,2,iel2)*xxj(2,indexf)+
     &              gradvel(j,3,iel2)*xxj(3,indexf))
!     
               if(dabs(vud).lt.1.d-20) then
!
!     upwind difference
!     
                  vfa(j,i)=vel(iel2,j)
                  cycle
               endif
!     
               phic=1.d0+vcd/vud
!     
               if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!     
!     upwind difference
!     
                  vfa(j,i)=vel(iel2,j)
               elseif(phic.le.1.d0/6.d0) then
                  vfa(j,i)=3.d0*vel(iel2,j)-2.d0*vel(iel1,j)+2.d0*vud
               elseif(phic.le.0.7d0) then
                  vfa(j,i)=3.d0*vel(iel2,j)/4.d0+vel(iel1,j)/4.d0
     &                    +vud/8.d0
               else
                  vfa(j,i)=vel(iel2,j)/3.d0+2.d0*vel(iel1,j)/3.d0
               endif
            enddo
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
