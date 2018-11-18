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
      subroutine correctflux(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
     &  vel,xlet,ielfa,xle,ifabou,xxnj,gradpcfa)
!
!     correction of v due to the balance of mass
!     the correction is in normal direction to the face
!
      implicit none
!
      integer i,nef,indexf,ipnei(*),neifa(*),neiel(*),iel2,ielfa(4,*),
     &  iel,ifa,ifabou(*),indexb
!
      real*8 flux(*),vfa(0:7,*),advfa(*),area(*),vel(nef,0:7),xlet(*),
     &  xle(*),totflux,maxflux,xxnj(3,*),gradpcfa(3,*)
!
c      maxflux=0.d0
c$omp parallel default(none)
c$omp& shared(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,vel,xlet,ielfa,
c$omp&        xle,ifabou,xxnj,gradpcfa)
c$omp& private(i,indexf,ifa,iel,iel2,totflux,maxflux,indexb)
c$omp do
      do i=1,nef
c         totflux=0.d0
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.gt.0) then
!
!              internal face
!
               flux(indexf)=flux(indexf)+vfa(5,ifa)*advfa(ifa)*area(ifa)
     &                      *(vel(i,4)-vel(iel,4))/xlet(indexf)
            else
               indexb=-ielfa(2,ifa)
               if(indexb.gt.0) then
                  if(((ifabou(indexb+1).eq.0).or.
     &                 (ifabou(indexb+2).eq.0).or.
     &                 (ifabou(indexb+3).eq.0)).and.
     &                 (ifabou(indexb+4).ne.0)) then
c               iel2=ielfa(2,ifa)
c               if(iel2.lt.0) then
c                  if(ifabou(-iel2+4).ne.0) then
!
!                    external face with pressure boundary conditions
!
                     flux(indexf)=flux(indexf)
     &                      +vfa(5,ifa)*advfa(ifa)*area(ifa)
     &                      *vel(i,4)/xle(indexf)
                  endif
               endif
            endif
c            totflux=totflux+flux(indexf)
         enddo
c         if(dabs(totflux).gt.maxflux)maxflux=dabs(totflux)
c         write(*,*) 'correctvfa mass check ',i,totflux,maxflux
      enddo
c$omp end do
c$omp end parallel
! 
      return
      end
