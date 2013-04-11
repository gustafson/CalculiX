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
      subroutine dasol(al,au,ad,b,jp,neq,energy)
      implicit none
c....solution of symmetric equations stored in profile form
c....coefficient matrix must be decomposed into its triangular
c....factors using datri before using dasol.
      integer jp(*),ior,iow,is,neq,j,jr,jh
      real*8 al(*),au(*),ad(*),b(*),zero,energy,bd,dot
      common /iofile/ ior,iow
      data zero/0.0d0/
c....find the first nonzero entry in the right hand side
      do 100 is=1,neq
        if(b(is).ne.zero) go to 200
  100 continue
      write(iow,2000) 
      if(ior.lt.0) write(*,2000)
      return
  200 if(is.lt.neq) then
c....reduce the right hand side
        do 300 j = is+1,neq
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
            b(j) = b(j) - dot(al(jr+1),b(j-jh),jh)
          endif
  300   continue
      endif
c....multiply by inverse of diagonal elements
      energy = zero
      do 400 j = is,neq
        bd = b(j)
        b(j) = b(j)*ad(j)
        energy=energy+bd*b(j)
  400 continue
c....backsubstitution
      if(neq.gt.1)then
        do 500 j = neq,2,-1
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
            call saxpb(au(jr+1),b(j-jh),-b(j),jh,b(j-jh))
          endif
  500   continue
      endif
      return
 2000 format(' ***DASOL WARNING 1*** Zero right-hand-side vector')
      end
