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
      subroutine frdscalar(epn,iset,nkcoords,inum,m1,
     &  istartset,iendset,ialset,ngraph,iselect)
!
!     stores a scalar result in frd format
!
      implicit none
!
      character*3 m1
!
      integer iset,nkcoords,inum(*),nksegment,ngraph,
     &  istartset(*),iendset(*),ialset(*),i,j,k,l,m,iselect
!
      real*8 epn(*)
!
      if(iset.eq.0) then
         do i=1,nkcoords
            if(iselect.eq.1) then
               if(inum(i).le.0) cycle
            elseif(iselect.eq.0) then
               if(inum(i).eq.0) cycle
            endif
            write(7,101) m1,i,epn(i)
         enddo
      else
         nksegment=nkcoords/ngraph
         do k=istartset(iset),iendset(iset)
            if(ialset(k).gt.0) then
               do l=0,ngraph-1
                  i=ialset(k)+l*nksegment
                  if(iselect.eq.1) then
                     if(inum(i).le.0) cycle
                  elseif(iselect.eq.0) then
                     if(inum(i).eq.0) cycle
                  endif
                  write(7,101) m1,i,epn(i)
               enddo
            else
               l=ialset(k-2)
               do
                  l=l-ialset(k)
                  if(l.ge.ialset(k-1)) exit
                  do m=0,ngraph-1
                     i=l+m*nksegment
                     if(iselect.eq.1) then
                        if(inum(i).le.0) cycle
                     elseif(iselect.eq.0) then
                        if(inum(i).eq.0) cycle
                     endif
                     write(7,101) m1,i,epn(i)
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


