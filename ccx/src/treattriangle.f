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
      subroutine treattriangle(inodesin,nnodesin,inodesout,
     &  nnodesout,nopes,slavstraight,xn,co,xl2s,ipe,ime,iactiveline,
     &  nactiveline,intersec,xntersec,ifreeintersec,itri,koncont,
     &  itriacornerl,nintpoint,pslavsurf,ncont,imastsurf,pmastsurf,
     &  xl2m,nnodelem,vold,mi,pnodesin,straight,gapmints,ssurf)
!
!     cuts a triangle of the master surface with a slave surface
!
      integer inodesin(*),nnodesin,nvertex,lvertex(13),inodesout(*),
     &  nnodesout,nopes,ipe(*),ime(4,*),iactiveline(3,*),nactiveline,
     &  intersec(2,*),ifreeintersec,itri,koncont(4,*),itriacornerl(4),
     &  i,j,k,nintpoint,ncont,idin,imastsurf(*),nnodelem,mi(2),ijk,
     &  ninsert,nodel,ssurf
!
      real*8 pvertex(3,13),pnodesin(3,*),slavstraight(20),xn(3),
     &  co(3,*),xilm,etlm,del,straight(16,*),spin,
     &  xl2s(3,*),xntersec(3,*),p(3,7),p1(3),p2(3),pslavsurf(3,*),
     &  ratio(8),dist,xil,etl,area,areax,areay,areaz,pmastsurf(2,*),
     &  xl2m(3,8),vold(0:mi(2),*),al,gapmints(*)
!
      include "gauss.f"
!
      data ijk /0/
      save ijk
!
      nvertex=0
!
      node1=koncont(1,itri)
      node2=koncont(2,itri)
      node3=koncont(3,itri)
!
!     check whether node 1 lies inside S
!
      call checktriavertex(inodesin,nnodesin,node1,nvertex,pvertex,
     &  lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &  xn,co,xl2s,vold,mi)
!
!     intersections of line node1-node2 with the edges of S
!
!     test pour idin
!
      if (node1.lt.node2) then
         call nident(inodesin,node1,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node1) then
               idin=0
            endif
         endif
      else 
         call nident(inodesin,node2,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node2) then
               idin=0
            endif
         endif
      endif
!     
!     
      call checktriaedge(node1,node2,ipe,ime,iactiveline,
     &     nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &     ifreeintersec,xn,co,nopes,xl2s,itri,idin,vold,mi)
!     
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!     
!     
!     check whether node 2 lies inside S
!     
      call checktriavertex(inodesin,nnodesin,node2,nvertex,pvertex,
     &     lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &     xn,co,xl2s,vold,mi)
!     
!     intersections of line node2-node3 with the edges of S
!     
!     test pour idin
!     
      if (node2.lt.node3) then
         call nident(inodesin,node2,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node2) then
               idin=0
            endif
         endif
      else 
         call nident(inodesin,node3,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node3) then
               idin=0
            endif
         endif
      endif
!     
      call checktriaedge(node2,node3,ipe,ime,iactiveline,
     &     nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &     ifreeintersec,xn,co,nopes,xl2s,itri,idin,vold,mi)
!     
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!     
!     check whether node 3 lies inside S
!     
      call checktriavertex(inodesin,nnodesin,node3,nvertex,pvertex,
     &     lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &     xn,co,xl2s,vold,mi)
!     
!     intersections of line node3-node1 with the edges of S
!     
!     test pour idin
!     
      if (node3.lt.node1) then
         call nident(inodesin,node3,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node3) then
               idin=0
            endif
         endif  
      else 
         call nident(inodesin,node1,nnodesin,idin)
         if (idin.gt.0) then 
            if (inodesin(idin).ne.node1) then
               idin=0
            endif
         endif
      endif
!     
      call checktriaedge(node3,node1,ipe,ime,iactiveline,
     &     nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &     ifreeintersec,xn,co,nopes,xl2s,itri,idin,vold,mi)
!     
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!     
!     check if there is always a slave node with value 1
!     
      do nodel=4,1,-1
         if(itriacornerl(nodel).eq.1) then
            nvertex=nvertex+1
          if (nvertex.lt.4) then
              do i=1,3
                  pvertex(i,nvertex)=xl2s(i,nodel)
               enddo
               itriacornerl(nodel)=2
            else
!
!           The convexity of the polygone has to be checked
!
            ninsert=nvertex-1
         spin=((pvertex(2,1)-xl2s(2,nodel))*(pvertex(3,ninsert)-
     & xl2s(3,nodel))-
     &   (pvertex(3,1)-xl2s(3,nodel))*(pvertex(2,ninsert)-
     & xl2s(2,nodel)))*xn(1)-
     &        ((pvertex(1,1)-xl2s(1,nodel))*(pvertex(3,ninsert)-
     & xl2s(3,nodel))-
     &   (pvertex(3,1)-xl2s(3,nodel))*(pvertex(1,ninsert)-
     & xl2s(1,nodel)))*xn(2)+
     &        ((pvertex(1,1)-xl2s(1,nodel))*(pvertex(2,ninsert)-
     &  xl2s(2,nodel))-
     &   (pvertex(2,1)-xl2s(2,nodel))*(pvertex(1,ninsert)-
     &  xl2s(1,nodel)))*xn(3)
10            if ((spin.gt.0.d0).and.(ninsert.gt.0)) then 
               pvertex(1,ninsert+1)=pvertex(1,ninsert)
               pvertex(2,ninsert+1)=pvertex(2,ninsert)
               pvertex(3,ninsert+1)=pvertex(3,ninsert)
         spin=((pvertex(2,ninsert)-xl2s(2,nodel))*(pvertex(3,ninsert-1)-
     & xl2s(3,nodel))-
     &   (pvertex(3,ninsert)-xl2s(3,nodel))*(pvertex(2,ninsert-1)-
     & xl2s(2,nodel)))*xn(1)-
     &        ((pvertex(1,ninsert)-xl2s(1,nodel))*(pvertex(3,ninsert-1)-
     & xl2s(3,nodel))-
     &   (pvertex(3,ninsert)-xl2s(3,nodel))*(pvertex(1,ninsert-1)-
     & xl2s(1,nodel)))*xn(2)+
     &        ((pvertex(1,ninsert)-xl2s(1,nodel))*(pvertex(2,ninsert-1)-
     & xl2s(2,nodel))-
     &   (pvertex(2,ninsert)-xl2s(2,nodel))*(pvertex(1,ninsert-1)-
     & xl2s(1,nodel)))*xn(3)
               ninsert=ninsert-1
            goto 10
          endif 
               pvertex(1,ninsert+1)=xl2s(1,nodel)
               pvertex(2,ninsert+1)=xl2s(2,nodel)
               pvertex(3,ninsert+1)=xl2s(3,nodel)
               itriacornerl(nodel)=2
          endif
         endif
      enddo
!     
!     generating integration points on the slave surface S
!     
      do k=1,nvertex-2
         p1(1)=pvertex(1,1+k)-pvertex(1,1)
         p1(2)=pvertex(2,1+k)-pvertex(2,1)
         p1(3)=pvertex(3,1+k)-pvertex(3,1)
         p2(1)=pvertex(1,2+k)-pvertex(1,1)
         p2(2)=pvertex(2,2+k)-pvertex(2,1)
         p2(3)=pvertex(3,2+k)-pvertex(3,1)
         areax=((p1(2)*p2(3))-(p2(2)*p1(3)))**2
         areay=(-(p1(1)*p2(3))+(p2(1)*p1(3)))**2
         areaz=((p1(1)*p2(2))-(p2(1)*p1(2)))**2
         area=dsqrt(areax+areay+areaz)/2.
!
!        storing the triangulation of the slave surfaces
!
         ijk=ijk+1
         write(40,100) ijk,(pvertex(i,1),i=1,3)
         ijk=ijk+1
         write(40,100) ,ijk,(pvertex(i,k+1),i=1,3)
         ijk=ijk+1
         write(40,100) ijk,(pvertex(i,k+2),i=1,3)
         write(40,101) ijk-2,ijk-2,ijk-1
         write(40,101) ijk-1,ijk-1,ijk
         write(40,101) ijk,ijk,ijk-2
 100     format('PNT ',i10,'P',3(1x,e15.8))
 101     format('LINE ',i10,'L',i10,'P ',i10,'P')
!
!     7 points scheme
!     
         do i=1,7
            do j=1,3
               p(j,i)=pvertex(j,1)*gauss2d6(1,i)+
     &              pvertex(j,1+k)*gauss2d6(2,i)+
     &              pvertex(j,2+k)*(1.d0-gauss2d6(1,i)-gauss2d6(2,i))
!     
            enddo
!     
            nintpoint=nintpoint+1
!
!           projection of the integration point in the mean
!           slave plane onto the slave surface
!
            call attachline(xl2s,p(1,i),nopes,ratio,dist,xil,etl,xn)
!
!           Calculation of the gap function at the integration point
!
            al=-(straight(16,itri)+straight(13,itri)*p(1,i)
     &      +straight(14,itri)*p(2,i)+straight(15,itri)*p(3,i))/
     &      (straight(13,itri)*xn(1)+straight(14,itri)*xn(2)
     &      +straight(15,itri)*xn(3))            
            gapmints(nintpoint)=al
!
!           calculation of the intersection with the master triangle
!
            p(1,i)=p(1,i)+al*xn(1)
            p(2,i)=p(2,i)+al*xn(2)
            p(3,i)=p(3,i)+al*xn(3)
!
!           projection of the master integration point onto the
!           master surface in order to get the local coordinates
!
            call attach(xl2m,p(1,i),nnodelem,ratio,dist,xilm,etlm)
!
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
            pslavsurf(3,nintpoint)=area*weight2d6(i)
            pmastsurf(1,nintpoint)=xilm
            pmastsurf(2,nintpoint)=etlm
            imastsurf(nintpoint)=koncont(4,itri)
         enddo
      enddo
!
      return
      end
      
