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
      subroutine profil(n,nnn,e2,adj,xadj,oldpro,newpro,
     &  oldpro_exp,newpro_exp)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer newpro,i,j,n,jstrt,jstop,oldpro,newmin,oldmin,e2,nnn(n),
     & xadj(n+1),adj(e2),inc_oldpro,inc_newpro,oldpro_exp,newpro_exp
!
      oldpro=0
      newpro=0
      oldpro_exp=0
      newpro_exp=0
      do 20 i=1,n
         jstrt=xadj(i)
         jstop=xadj(i+1)-1
         if(jstrt.gt.jstop) cycle
         oldmin=adj(jstrt)
         newmin=nnn(adj(jstrt))
!
         do 10 j=jstrt+1,jstop
            oldmin=min(oldmin,adj(j))
            newmin=min(newmin,nnn(adj(j)))
 10      continue
!
         inc_oldpro=dim(i,oldmin)
         if(2147483647-oldpro.lt.inc_oldpro) then
            oldpro_exp=oldpro_exp+1
            inc_oldpro=inc_oldpro-2147483647
         endif
         oldpro=oldpro+inc_oldpro
!
         inc_newpro=dim(nnn(i),newmin)
         if(2147483647-newpro.lt.inc_newpro) then
            newpro_exp=newpro_exp+1
            inc_newpro=inc_newpro-2147483647
         endif
         newpro=newpro+inc_newpro
 20   continue
!
      inc_oldpro=n
      if(2147483647-oldpro.lt.inc_oldpro) then
         oldpro_exp=oldpro_exp+1
         inc_oldpro=inc_oldpro-2147483647
      endif
      oldpro=oldpro+inc_oldpro
!
      inc_newpro=n
      if(2147483647-newpro.lt.inc_newpro) then
         newpro_exp=newpro_exp+1
         inc_newpro=inc_newpro-2147483647
      endif
      newpro=newpro+inc_newpro
!
      return
      end
