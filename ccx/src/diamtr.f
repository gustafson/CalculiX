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
      subroutine diamtr(n,e2,adj,xadj,mask,ls,xls,hlevel,snode,nc)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer nc,j,snode,degree,mindeg,istrt,istop,hsize,node,jstrt,
     & jstop,ewidth,i,width,depth,enode,n,sdepth,e2,xadj(n+1),adj(e2),
     & xls(n+1),ls(n),mask(n),hlevel(n)
!
      mindeg=n
      do 10 i=1,n
         if(mask(i).eq.0) then
            degree=xadj(i+1)-xadj(i)
            if(degree.lt.mindeg) then
               snode=i
               mindeg=degree
            endif
         endif
 10   continue
!
      call rootls(n,snode,n+1,e2,adj,xadj,mask,ls,xls,sdepth,width)
!
      nc=xls(sdepth+1)-1
!
 15   continue
!
      hsize=0
      istrt=xls(sdepth)
      istop=xls(sdepth+1)-1
      do 20 i=istrt,istop
         node=ls(i)
         hsize=hsize+1
         hlevel(hsize)=node
         xls(node)=xadj(node+1)-xadj(node)
 20   continue
!
      if(hsize.gt.1) call isorti(hsize,hlevel,n,xls)
!
      istop=hsize
      hsize=1
      degree=xls(hlevel(1))
      do 25 i=2,istop
         node=hlevel(i)
         if(xls(node).ne.degree) then
            degree=xls(node)
            hsize=hsize+1
            hlevel(hsize)=node
         endif
 25   continue
!
      ewidth=nc+1
      do 30 i=1,hsize
        node=hlevel(i)
!
        call rootls(n,node,ewidth,e2,adj,xadj,mask,ls,xls,depth,width)
        if(width.lt.ewidth) then
!
           if(depth.gt.sdepth) then
!
              snode=node
              sdepth=depth
              go to 15
           endif
!
           enode=node
           ewidth=width
        endif
 30   continue
!
      if(node.ne.enode) then
        call rootls(n,enode,nc+1,e2,adj,xadj,mask,ls,xls,depth,width)
      endif
!
      do 50 i=1,depth
         jstrt=xls(i)
         jstop=xls(i+1)-1
         do 40 j=jstrt,jstop
            mask(ls(j))=i-1
 40      continue
 50   continue
      end
         
