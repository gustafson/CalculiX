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
      subroutine dredu(al,au,ad,jh,flg  ,dj)
      implicit real*8 (a-h,o-z)
c....reduce diagonal element in triangular decomposition
      logical flg
      real*8 al(jh),au(jh),ad(jh)
      do 100 j = 1,jh
        ud = au(j)*ad(j)
        dj = dj - al(j)*ud
        au(j) = ud
  100 continue
c....finish computation of column of al for unsymmetric matrices
      if(flg) then
        do 200 j = 1,jh
          al(j) = al(j)*ad(j)
  200   continue
      endif
      return 
      end
