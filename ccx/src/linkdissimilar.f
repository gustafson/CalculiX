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
      subroutine linkdissimilar(co,nk,ics,csab,xn,yn,zn,ncsnodes,
     &  rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,straight,
     &  lcs,nodef,ratio,nterms,rp,zp,netri,
     &  nodesonaxis,nodel,noder,ifacetet,inodface)
!
!     links dissimilar meshes for cyclic symmetry
!
      implicit none
!
      logical nodesonaxis
!
      integer nneigh,jcs(*),j,nk,ics(*),nodef(8),noder,nodel,
     &  nrcg(*),ncsnodes,nzcg(*),nterms,
     &  ineigh(20),itrimax,lcs(*),i,
     &  netri,kontri(3,*),itri,ifirst,ilast,
     &  ifacetet(*),inodface(*)
!
      real*8 co(3,*),csab(7),xn,yn,zn,rp,zp,xi,et,
     &  straight(9,*),zcscg(*),rcscg(*),zcs0cg(*),rcs0cg(*),distmax,
     &  dist,ratio(8),pneigh(0:3,8),pnode(3)
!
!     finding for each node on the right side a corresponding
!     triangle
!
      nneigh=10
!
      call near2d(rcs0cg,zcs0cg,rcscg,zcscg,nrcg,nzcg,rp,zp,
     &     netri,ineigh,nneigh)
!     
      distmax=1.d30
!     
      do j=1,nneigh
         itri=ineigh(j)
         dist=max(rp*straight(1,itri)+zp*straight(2,itri)+
     &        straight(3,itri),0.d0)+
     &        max(rp*straight(4,itri)+zp*straight(5,itri)+
     &        straight(6,itri),0.d0)+
     &        max(rp*straight(7,itri)+zp*straight(8,itri)+
     &        straight(9,itri),0.d0)
         if(dist.le.0.d0) then
            distmax=0.d0
            itrimax=itri
            exit
         endif
         if(dist.lt.distmax) then
            itrimax=itri
            distmax=dist
         endif
      enddo
!
      itri=itrimax
!
      ilast=ifacetet(itri)
      do
         itri=itri-1
         if(itri.eq.0) then
            ifirst=0
            exit
         endif
         if(ifacetet(itri).ne.ilast) then
            ifirst=ifacetet(itri)
            exit
         endif
      enddo
      nterms=ilast-ifirst
      do i=1,nterms
         nodef(i)=inodface(ifirst+i)
         do j=1,3
            pneigh(j,i)=co(j,nodef(i))
         enddo
      enddo
      do j=1,3
         pnode(j)=co(j,noder)
      enddo
!
      call attach(pneigh,pnode,nterms,ratio,dist,xi,et)
!
      do j=1,3
         co(j,noder)=pnode(j)
      enddo
!
      write(*,*) 'distance moved: ',dsqrt(dist)
!     
      return
      end
      
