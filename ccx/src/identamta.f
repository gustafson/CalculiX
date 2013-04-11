!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
!                                                                              
!     identifies the position id of reftime in an ordered array
!     amta(1,istart...iend) of real numbers; amta is defined as amta(2,*)
!  
!     id is such that amta(1,id).le.reftime and amta(1,id+1).gt.reftime
!                                                                             
      SUBROUTINE identamta(amta,reftime,istart,iend,ID)                                               
      IMPLICIT none
!
      integer id,istart,iend,n2,m
      real*8 amta(2,*),reftime
      id=istart-1
      if(iend.lt.istart) return
      N2=iend+1
      DO                                                         
         M=(N2+ID)/2
         IF(reftime.GE.amta(1,M)) then
            ID=M     
         else
            N2=M  
         endif
         IF((N2-ID).EQ.1) return
      enddo
      END
