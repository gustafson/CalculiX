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
      subroutine calcdimknot(idepnodes,ndepnodes,co,irotnode,iexpnode,
     & idimknot)
!
!     determines the dimensionality of the cloud
!     of nodes belonging to a knot (by calculating the second
!     order moments).
!
      implicit none
!
      integer idepnodes(*),ndepnodes,irotnode,iexpnode,i,node,n,
     &  matz,ier,idimknot
!
      real*8 co(3,*),tx,ty,tz,txx,tyy,tzz,txy,txz,tyz,fv1(3),fv2(3),
     &  c(3,3),w(3),z(3,3)
!
      tx=0.d0
      ty=0.d0
      tz=0.d0
      txx=0.d0
      tyy=0.d0
      tzz=0.d0
      txy=0.d0
      txz=0.d0
      tyz=0.d0
!
!     building the first and second order moments
!
      do i=1,ndepnodes
         node=idepnodes(i)
         tx=tx+co(1,node)
         ty=ty+co(2,node)
         tz=tz+co(3,node)
         txx=txx+co(1,node)**2
         tyy=tyy+co(2,node)**2
         tzz=tzz+co(3,node)**2
         txy=txy+co(1,node)*co(2,node)
         txz=txz+co(1,node)*co(3,node)
         tyz=tyz+co(2,node)*co(3,node)
      enddo
!
      c(1,1)=txx-tx*tx/ndepnodes
      c(2,2)=tyy-ty*ty/ndepnodes
      c(3,3)=tzz-tz*tz/ndepnodes
      c(1,2)=txy-tx*ty/ndepnodes
      c(1,3)=txz-tx*tz/ndepnodes
      c(2,3)=tyz-ty*tz/ndepnodes
      c(2,1)=c(1,2)
      c(3,1)=c(1,3)
      c(3,2)=c(2,3)
!
!     calculating the eigenvalues and eigenvectors
!
      n=3
      matz=1
      call rs(n,n,c,w,matz,z,fv1,fv2,ier)
      if(ier.ne.0) then
         write(*,*) '*ERROR in createknot: calculation of the'
         write(*,*) '       eigenvalues was not successfull'
         stop
      endif
!
!     determining the dimensionality
!
      idimknot=3
      if(dabs(w(2)).lt.1.d-10) then
         idimknot=1
      elseif(dabs(w(1)).lt.1.d-10) then
         idimknot=2
      endif
!
      if(idimknot.eq.2) then
         do i=1,3
            co(i,irotnode)=z(i,2)
            co(i,iexpnode)=z(i,3)
         enddo
      endif
!
      return
      end
