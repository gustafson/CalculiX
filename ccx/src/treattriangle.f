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
     &  nnodesout,nopes,slavstraight,xn,co,xl2,ipe,ime,iactiveline,
     &  nactiveline,intersec,xntersec,ifreeintersec,itri,koncont,
     &  itriacornerl,nintpoint,pslavsurf,ncont,imastsurf,pmastsurf,
     &  pneigh,nnodelem,vold,mi,pnodesin)
!
!     cuts a triangle of the master surface with a slave surface
!
!
      integer inodesin(*),nnodesin,nvertex,lvertex(13),inodesout(*),
     &  nnodesout,nopes,ipe(*),ime(4,*),iactiveline(3,*),nactiveline,
     &  intersec(2,*),ifreeintersec,itri,koncont(4,*),itriacornerl(4),
     &  i,j,k,nintpoint,ncont,idin,imastsurf(*),nnodelem,mi(2)
!
      real*8 pvertex(3,13),pnodesin(3,*),slavstraight(20),xn(3),
     &  co(3,*),xilm,etlm,
     &  xl2(3,*),xntersec(3,*),p(3,3),p1(3),p2(3),pslavsurf(3,*),
     &  ratio(8),dist,xil,etl,area,areax,areay,areaz,pmastsurf(2,*),
     &  pneigh(3,8),vold(0:mi(2),*)
!
      include "gauss.f"
!
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
     &  xn,co,xl2,vold,mi)
!
!     intersections of line node1-node2 with the edges of S
!
!	test pour idin
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
     &  nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &  ifreeintersec,xn,co,nopes,xl2,itri,idin,vold,mi)
!
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!
      call checkslavevertex(lvertex,nvertex,pvertex,
     &  itriacornerl,xl2)
!
!     check whether node 2 lies inside S
!

      call checktriavertex(inodesin,nnodesin,node2,nvertex,pvertex,
     &  lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &  xn,co,xl2,vold,mi)
!
!     intersections of line node2-node3 with the edges of S
!
!	test pour idin
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
     &  nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &  ifreeintersec,xn,co,nopes,xl2,itri,idin,vold,mi)
!
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!
      call checkslavevertex(lvertex,nvertex,pvertex,
     &  itriacornerl,xl2)
!
!
!     check whether node 3 lies inside S
!
      call checktriavertex(inodesin,nnodesin,node3,nvertex,pvertex,
     &  lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &  xn,co,xl2,vold,mi)
!
!     intersections of line node3-node1 with the edges of S
!
!	test pour idin
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
     &  nactiveline,intersec,xntersec,nvertex,pvertex,lvertex,
     &  ifreeintersec,xn,co,nopes,xl2,itri,idin,vold,mi)
!
!     if there are intersections, check whether the S-vertex at
!     the end of an intersected S-edge must be included
!
      call checkslavevertex(lvertex,nvertex,pvertex,
     &  itriacornerl,xl2)
!
!	check if all the intern salve vertexes have been treated
!
	ilast=itriacornerl(1)
	do nodel=4,1,-1 
	   if(itriacornerl(nodel).eq.1) then
	      if(ilast.eq.2) then
		nvertex=nvertex+1
		  do i=1,3
		    pvertex(i,nvertex)=xl2(i,nodel)
		  enddo
		  lvertex(nvertex)=0
		  itriacornerl(nodel)=2
	      else
		ilast=itriacornerl(nodel)
	      endif
	   else
		ilast=itriacornerl(nodel)		
	   endif
	enddo
!
!	check if there is always a slave node with value 1
!
	do nodel=4,1,-1
	 if(itriacornerl(nodel).eq.1) then
		nvertex=nvertex+1
		  do i=1,3
		    pvertex(i,nvertex)=xl2(i,nodel)
		  enddo
		  lvertex(nvertex)=0
		  itriacornerl(nodel)=2	 
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
       if ((areax.gt.0d0).or.(areay.gt.0d0).or.(areaz.gt.0d0)) then
          area=dsqrt(areax+areay+areaz)/2.
c       WRITE(*,*) "A ",area, "itri",itri,"boucle k",k
          if (area.lt.1d-4) cycle
       else
          cycle
       endif
!      7 points scheme
!
!      center
! 
         do i=1,1
            do j=1,3
               p(j,i)=pvertex(j,1)*1/3+
     &              pvertex(j,1+k)*1/3+
     &              pvertex(j,2+k)*1/3
!
            enddo
!
            call attach(xl2,p(1,i),nopes,ratio,dist,xil,etl)
	    call attach(pneigh,p(1,i),nnodelem,ratio,dist,xilm,etlm)
            nintpoint=nintpoint+1
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
	    pslavsurf(3,nintpoint)=area*weight2d7(1)
	    pmastsurf(1,nintpoint)=xilm
	    pmastsurf(2,nintpoint)=etlm
 	    imastsurf(nintpoint)=koncont(4,itri)
!
         enddo
!
!      first 3 points
!
         do i=1,3
            do j=1,3
               p(j,i)=pvertex(j,1)*gauss2d71(1,i)+
     &              pvertex(j,1+k)*gauss2d71(2,i)+
     &              pvertex(j,2+k)*(1.d0-gauss2d71(1,i)-gauss2d71(2,i))
!
            enddo
!
            call attach(xl2,p(1,i),nopes,ratio,dist,xil,etl)
	    call attach(pneigh,p(1,i),nnodelem,ratio,dist,xilm,etlm)
            nintpoint=nintpoint+1
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
	    pslavsurf(3,nintpoint)=area*weight2d7(2)
	    pmastsurf(1,nintpoint)=xilm
	    pmastsurf(2,nintpoint)=etlm
 	    imastsurf(nintpoint)=koncont(4,itri)
         enddo
!       last three points
         do i=1,3
            do j=1,3
               p(j,i)=pvertex(j,1)*gauss2d72(1,i)+
     &              pvertex(j,1+k)*gauss2d72(2,i)+
     &              pvertex(j,2+k)*(1.d0-gauss2d72(1,i)-gauss2d72(2,i))
!
            enddo
!
            call attach(xl2,p(1,i),nopes,ratio,dist,xil,etl)
	    call attach(pneigh,p(1,i),nnodelem,ratio,dist,xilm,etlm)
            nintpoint=nintpoint+1
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
	    pslavsurf(3,nintpoint)=area*weight2d7(3)
	    pmastsurf(1,nintpoint)=xilm
	    pmastsurf(2,nintpoint)=etlm
 	    imastsurf(nintpoint)=koncont(4,itri)
         enddo
       enddo
       return
       end
