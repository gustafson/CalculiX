!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2013 Guido Dhondt
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
!     this subroutine enables to calculate the discharge coefficient of a stepped labyrinth seal
!     as a function of the reynolds number, the ratios s/l ,r/b and p1/p2
!
!     the related data can be found in
!     "Some aerodynamic Aspects of Engine Secondary air systems"
!     H. Zimmermann
!     ASME 89-GT-209 
!     Table p 7
!
      subroutine cd_lab_reynolds(reynolds,cd_reynolds)
!
      implicit none
!
      integer id,n6
!
      real*8 reynolds , cd_reynolds
!
      real*8 tab_reynolds(6)
      data tab_reynolds
     &     /220.d0,630.d0,1260d0,2300d0,3200d0,4300d0/
!
      real*8 tab_cd(6)
      data tab_cd
     &    / 0.32d0,0.39d0,0.44d0,0.49d0,0.25d0,0.54d0/
!
      data n6 /6/
!
      call ident(tab_reynolds,reynolds,n6,id)
      
      if(id.eq.1) then
         cd_reynolds=tab_cd(1)
      elseif(id.eq.18) then
         cd_reynolds=tab_cd(6)
      else
         cd_reynolds=tab_cd(id)+(tab_cd(id+1)-tab_cd(id))
     &        *(reynolds-tab_reynolds(id))
     &        /(tab_reynolds(id+1)-tab_reynolds(id))
      endif
!     
      return
!     
      end
      
