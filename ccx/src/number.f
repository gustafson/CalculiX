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
      subroutine number(n,nc,snode,lstnum,e2,adj,xadj,s,q,p)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer nc,lstnum,jstrt,jstop,istop,nbr,nabor,i,j,next,addres,
     & nn,node,snode,istrt,maxprt,prty,n,w1,w2,e2,q(nc),xadj(n+1),
     & adj(e2),p(n),s(n)
!
      parameter(w1=1,w2=2)
!
      do 10 i=1,nc
         node=q(i)
         p(node)=w1*s(node)-w2*(xadj(node+1)-xadj(node)+1)
         s(node)=-2
 10   continue
!
      nn=1
      q(nn)=snode
      s(snode)=-1
!
 30   if(nn.gt.0) then
!
         addres=1
         maxprt=p(q(1))
         do 35 i=2,nn
            prty=p(q(i))
            if(prty.gt.maxprt) then
               addres=i
               maxprt=prty
            endif
 35      continue
!
         next=q(addres)
!
         q(addres)=q(nn)
         nn=nn-1
         istrt=xadj(next)
         istop=xadj(next+1)-1
         if(s(next).eq.-1) then
!
            do 50 i=istrt,istop
!
               nbr=adj(i)
               p(nbr)=p(nbr)+w2
!
               if(s(nbr).eq.-2) then
                  nn=nn+1
                  q(nn)=nbr
                  s(nbr)=-1
               endif
 50         continue
         endif
!
         lstnum=lstnum+1
         s(next)=lstnum
!
         do 80 i=istrt,istop
            nbr=adj(i)
            if(s(nbr).eq.-1) then
!
               p(nbr)=p(nbr)+w2
               s(nbr)=0
!
               jstrt=xadj(nbr)
               jstop=xadj(nbr+1)-1
               do 60 j=jstrt,jstop
                  nabor=adj(j)
!
                  p(nabor)=p(nabor)+w2
                  if(s(nabor).eq.-2) then
!
                     nn=nn+1
                     q(nn)=nabor
                     s(nabor)=-1
                  endif
 60            continue
            endif
 80      continue
         go to 30
      endif
      end
