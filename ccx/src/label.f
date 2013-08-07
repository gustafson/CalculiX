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
      subroutine label(n,e2,adj,xadj,nnn,iw,oldpro,newpro,
     &  oldpro_exp,newpro_exp)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer n,i1,i2,i3,i,snode,lstnum,nc,oldpro,newpro,e2,xadj(n+1),
     & adj(e2),nnn(n),iw(3*n+1),oldpro_exp,newpro_exp
!
      do 10 i=1,n
         nnn(i)=0
 10   continue
!
      i1=1
      i2=i1+n
      i3=i2+n+1
!
      lstnum=0
 20   if(lstnum.lt.n) then
!
        call diamtr(n,e2,adj,xadj,nnn,iw(i1),iw(i2),iw(i3),snode,nc)
!
        call number(n,nc,snode,lstnum,e2,adj,xadj,nnn,iw(i1),iw(i2))
        go to 20
      endif
!
      call profil(n,nnn,e2,adj,xadj,oldpro,newpro,oldpro_exp,
     &  newpro_exp)
!
      if((oldpro_exp.lt.newpro_exp).or.
     &   ((oldpro_exp.eq.newpro_exp).and.(oldpro.lt.newpro))) then
         do 30 i=1,n
            nnn(i)=i
 30      continue
         newpro=oldpro
         newpro_exp=oldpro_exp
      endif
      end
