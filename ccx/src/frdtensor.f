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
      subroutine frdtensor(stn,iset,nkcoords,inum,m1,istartset,iendset,
     &             ialset,ngraph)
!
!     stores a tensor result (2nd order) in frd format
!
      implicit none
!
      character*3 m1
!
      integer iset,nkcoords,inum(*),ngraph,nksegment,
     &  istartset(*),iendset(*),ialset(*),i,j,k,l,m,kal(2,6)
!
      real*8 stn(6,*)
!
      data kal /1,1,2,2,3,3,1,2,1,3,2,3/
!
      if(iset.eq.0) then
         do i=1,nkcoords
            if(inum(i).le.0) cycle
            write(7,101) m1,i,(stn(j,i),j=1,4),
     &           stn(6,i),stn(5,i)
         enddo
      else
         nksegment=nkcoords/ngraph
         do k=istartset(iset),iendset(iset)
            if(ialset(k).gt.0) then
               do l=0,ngraph-1
                  i=ialset(k)+l*nksegment
                  if(inum(i).le.0) cycle
                  write(7,101) m1,i,(stn(j,i),j=1,4),
     &                 stn(6,i),stn(5,i)
               enddo
            else
               l=ialset(k-2)
               do
                  l=l-ialset(k)
                  if(l.ge.ialset(k-1)) exit
                  do m=0,ngraph-1
                     i=l+m*nksegment
                     if(inum(i).le.0) cycle
                     write(7,101) m1,i,(stn(j,i),j=1,4),
     &                    stn(6,i),stn(5,i)
                  enddo
               enddo
            endif
         enddo
      endif
!     
 101  format(a3,i10,1p,6e12.5)
!     
      return
      end
      

