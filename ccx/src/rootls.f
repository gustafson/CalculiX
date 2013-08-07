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
      subroutine rootls(n,root,maxwid,e2,adj,xadj,mask,ls,xls,depth,
     &                  width)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer root,depth,nbr,maxwid,lstrt,lstop,lwdth,node,nc,width,n,
     & jstrt,jstop,i,j,e2,xadj(n+1),adj(e2),mask(n),xls(n+1),ls(n)
!
      mask(root)=1
      ls(1)=root
      nc=1
      width=1
      depth=0
      lstop=0
      lwdth=1
 10   if(lwdth.gt.0) then
!
         lstrt=lstop+1
         lstop=nc
         depth=depth+1
         xls(depth)=lstrt
!
         do 30 i=lstrt,lstop
            node=ls(i)
            jstrt=xadj(node)
            jstop=xadj(node+1)-1
            do 20 j=jstrt,jstop
               nbr=adj(j)
               if(mask(nbr).eq.0) then
                  nc=nc+1
                  ls(nc)=nbr
                  mask(nbr)=1
               endif
 20         continue
 30      continue
!
         lwdth=nc-lstop
         width=max(lwdth,width)
!
         if(width.ge.maxwid) go to 35
         go to 10
      endif
      xls(depth+1)=lstop+1
!
 35   continue
      do 40 i=1,nc
         mask(ls(i))=0
 40   continue
      end
