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
      subroutine linkdissimilar(co,csab,
     &  rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,straight,
     &  nodef,ratio,nterms,rp,zp,netri,
     &  nodei,ifacetet,inodface,noded,xn,yn,zn,ier)
!
!     links dissimilar meshes for cyclic symmetry
!
      implicit none
!
      integer nneigh,j,nodef(8),nodei,nrcg(*),
     &  nzcg(*),nterms,ineigh(20),itrimax,i,netri,
     &  itri,ifirst,ilast,ifacetet(*),inodface(*),noded,ier
!
      real*8 co(3,*),csab(7),rp,zp,xi,et,xn,yn,zn,rp1,zp1,rp2,zp2,
     &  straight(9,*),zcscg(*),rcscg(*),zcs0cg(*),rcs0cg(*),distmax,
     &  dist,ratio(8),pneigh(3,8),pnode(3),xap,yap,zap,
     &  x12,y12,z12,x13,y13,z13,area,typdist
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
         pnode(j)=co(j,nodei)
      enddo
!
      xap=pnode(1)-csab(1)
      yap=pnode(2)-csab(2)
      zap=pnode(3)-csab(3)
!     
      zp1=xap*xn+yap*yn+zap*zn
      rp1=dsqrt((xap-zp1*xn)**2+(yap-zp1*yn)**2+(zap-zp1*zn)**2)
!
      call attach(pneigh,pnode,nterms,ratio,dist,xi,et)
!
      xap=pnode(1)-csab(1)
      yap=pnode(2)-csab(2)
      zap=pnode(3)-csab(3)
!     
      zp2=xap*xn+yap*yn+zap*zn
      rp2=dsqrt((xap-zp2*xn)**2+(yap-zp2*yn)**2+(zap-zp2*zn)**2)
!
      do j=1,3
         co(j,nodei)=pnode(j)
      enddo
!
!     check whether this distance is inferior to the tolerance
!     only for projections on the exterior border of a face
!
      if((((nterms.eq.4).or.(nterms.eq.8)).and.
     &    ((xi.le.-1.d0).or.(xi.ge.1.d0).or.
     &     (et.le.-1.d0).or.(et.ge.1.d0))).or.
     &   (((nterms.eq.3).or.(nterms.eq.6)).and.
     &    ((xi.le.0.d0).or.(et.le.0.d0).or.(xi+et.ge.1.d0)))) then
!
!        calculating a typical distance of the face
!
         x12=pneigh(1,2)-pneigh(1,1)
         y12=pneigh(2,2)-pneigh(2,1)
         z12=pneigh(3,2)-pneigh(3,1)
         x13=pneigh(1,3)-pneigh(1,1)
         y13=pneigh(2,3)-pneigh(2,1)
         z13=pneigh(3,3)-pneigh(3,1)
         area=dsqrt((y12*z13-y13*z12)**2+
     &              (x12*z13-x13*z12)**2+
     &              (x12*y13-x13*y12)**2)
         typdist=dsqrt(area)
        
         dist=dsqrt((rp2-rp1)**2+(zp2-zp1)**2)
         if(dist.ge.typdist/10.d0) then
            write(*,*) '*WARNING in linkdissimilar: no suitable partner'
            write(*,*) '         face found for node', noded,'.'
            write(*,*) 
     &      '         Nodes belonging to the best partner face:'
            write(*,*) (nodef(i),i=1,nterms)
            write(*,*) '         3-D Euclidean distance: ',dist
            write(*,*) 
            ier=-1
         endif
      endif
!     
      return
      end
      
