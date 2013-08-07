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
      subroutine graph(n,ne,inpn,npn,xnpn,iadj,adj,xadj)
!
!     Sloan routine (Int.J.Num.Meth.Engng 28, 2651-2679(1989))
!
      implicit none
!
      integer n,ne,nodej,nodek,mstrt,iadj,i,j,k,jstrt,jstop,lstrt,
     & lstop,l,nen1,mstop,m,inpn,xnpn(ne+1),npn(inpn),adj(iadj),
     & xadj(n+1)
!
      do 5 i=1,iadj
         adj(i)=0
 5    continue
      do 10 i=1,n
         xadj(i)=0
 10   continue
!
      do 30 i=1,ne
         jstrt=xnpn(i)
         jstop=xnpn(i+1)-1
         nen1=jstop-jstrt
         do 20 j=jstrt,jstop
            nodej=npn(j)
            xadj(nodej)=xadj(nodej)+nen1
 20      continue
 30   continue
!
      l=1
      do 40 i=1,n
         l=l+xadj(i)
         xadj(i)=l-xadj(i)
 40   continue
      xadj(n+1)=l
!
      do 90 i=1,ne
         jstrt=xnpn(i)
         jstop=xnpn(i+1)-1
         do 80 j=jstrt,jstop-1
            nodej=npn(j)
            lstrt=xadj(nodej)
            lstop=xadj(nodej+1)-1
            do 70 k=j+1,jstop
               nodek=npn(k)
               do 50 l=lstrt,lstop
                  if(adj(l).eq.nodek) go to 70
                  if(adj(l).eq.0) go to 55
 50            continue
               write(6,1000)
               stop
 55            continue
               adj(l)=nodek
               mstrt=xadj(nodek)
               mstop=xadj(nodek+1)-1
               do 60 m=mstrt,mstop
                  if(adj(m).eq.0) go to 65
 60            continue
               write(6,1000)
               stop
 65            continue
               adj(m)=nodej
 70         continue
 80      continue
 90   continue
!
      k=0
      jstrt=1
      do 110 i=1,n
         jstop=xadj(i+1)-1
         do 100 j=jstrt,jstop
            if(adj(j).eq.0) go to 105
            k=k+1
            adj(k)=adj(j)
 100     continue
 105     continue
         xadj(i+1)=k+1
         jstrt=jstop+1
 110  continue
!
 1000 format(//,1x,'***error in graph***',
     &       //,1x,'cannot assemble node adjacency list',
     &       //,1x,'check npn and xnpn arrays')
      end
