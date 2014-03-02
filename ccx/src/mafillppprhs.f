!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine mafillppprhs(ne,lakon,nactdoh,ipnei,neifa,neiel,neij,
     &    gradpp,xxi,cosa,xxn,ap,xle,xlen,b,ielfa,ifabou)
!
!     calculating the right hand side of the equations to determine p''
!
      implicit none
!
      character*8 lakon(*)
!
      integer ne,nactdoh(*),ipnei(*),neifa(*),neiel(*),neij(*),
     &  i,j,jdof1,indexf,ifa,iel,j2,indexf2,ielfa(4,*),ifabou(*)
!
      real*8 gradpp(3,*),xxi(3,*),cosa(*),xxn(3,*),ap(*),xle(*),xlen(*),
     &  rhs,b(*)
!
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         jdof1=nactdoh(i)
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,6
!
!              diffusion
!
               indexf=indexf+1
               ifa=neifa(indexf)
               iel=neiel(indexf)
               if(iel.ne.0) then
                  j2=neij(indexf)
                  indexf2=ipnei(iel)+j2
                  rhs=((gradpp(1,iel)*(xxi(1,indexf2)
     &                            -cosa(indexf2)*xxn(1,indexf2))+
     &                     gradpp(2,iel)*(xxi(2,indexf2)
     &                            -cosa(indexf2)*xxn(2,indexf2))+
     &                     gradpp(3,iel)*(xxi(3,indexf2)
     &                            -cosa(indexf2)*xxn(3,indexf2)))
     &                            *xlen(indexf2)
     &                   -(gradpp(1,i)*(xxi(1,indexf)
     &                            -cosa(indexf)*xxn(1,indexf))+
     &                     gradpp(2,i)*(xxi(2,indexf)
     &                            -cosa(indexf)*xxn(2,indexf))+
     &                     gradpp(3,i)*(xxi(3,indexf)
     &                            -cosa(indexf)*xxn(3,indexf)))
     &                            *xle(indexf))*ap(ifa)
                  b(jdof1)=b(jdof1)+rhs
               else
                  if(ielfa(2,ifa).lt.0) then
                     if(ifabou(-ielfa(2,ifa))+4.ne.0) then
!
!                        pressure given
!
                        rhs=(-(gradpp(1,i)*(xxi(1,indexf)
     &                            -cosa(indexf)*xxn(1,indexf))+
     &                     gradpp(2,i)*(xxi(2,indexf)
     &                            -cosa(indexf)*xxn(2,indexf))+
     &                     gradpp(3,i)*(xxi(3,indexf)
     &                            -cosa(indexf)*xxn(3,indexf)))
     &                            *xle(indexf))*ap(ifa)
                        b(jdof1)=b(jdof1)+rhs
                     endif
                  endif
               endif
            enddo
         endif
      enddo
!     
      return
      end
