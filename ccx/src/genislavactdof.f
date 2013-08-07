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
      subroutine genislavactdof(ntie,tieset,neq,nactdof,nslavnode,
     &     islavact,islavactdof,islavnode,mi)
!
!     Author : Samoela Rakotonanahary, Saskia Sitzmann
!     genislavactdof get the field islavactdof in order to 
!     help calculating the tangential matrices.
!
!     islavactdof is the inverse of nactdof for active slave nodes:
!     it links an active slave degree of freedom to the 
!     corresponding slave node position in field islavnode and the
!     global (x-y-z) degree of freedom
!     
      integer i,j,k,ntie,neq(*),node,nslavnode(*),
     &     mi(*),nactdof(0:mi(2),*),
     &     islavact(*),islavactdof(*),islavnode(*)
      character*81 tieset(3,*)
!     
!     close the contact.fbd file
!
      close(20)
      close(30)
      close(40)
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do j = nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(j)
            if(islavact(j).gt.-1) then
               do k=1,3
                  if (nactdof(k,node).eq.0) cycle
                  islavactdof(nactdof(k,node))=10*j+k
               enddo
            endif         
         enddo
      enddo
!
      return
      end
