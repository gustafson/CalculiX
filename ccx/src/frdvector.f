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
      subroutine frdvector(v,iset,ntrans,filabl,nkcoords,inum,m1,inotr,
     &  trab,co,istartset,iendset,ialset,mi,ngraph)
!
!     stores a vector result in frd format
!
      implicit none
!
      character*3 m1
      character*87 filabl
!
      integer mi(2),iset,ntrans,nkcoords,inum(*),inotr(2,*),
     &  istartset(*),iendset(*),ialset(*),i,j,k,l,m,ngraph,
     &  nksegment
!
      real*8 v(0:mi(2),*),trab(7,*),co(3,*),a(3,3)
!
      if(iset.eq.0) then
         if((ntrans.eq.0).or.(filabl(6:6).eq.'G')) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               write(7,101) m1,i,(v(j,i),j=1,3)
            enddo
         else
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               if(inotr(1,i).eq.0) then
                  write(7,101) m1,i,(v(j,i),j=1,3)
               else
                  call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                  write(7,101) m1,i,
     &                 v(1,i)*a(1,1)+v(2,i)*a(2,1)+v(3,i)*a(3,1),
     &                 v(1,i)*a(1,2)+v(2,i)*a(2,2)+v(3,i)*a(3,2),
     &                 v(1,i)*a(1,3)+v(2,i)*a(2,3)+v(3,i)*a(3,3)
               endif
            enddo
         endif
      else
         nksegment=nkcoords/ngraph
         do k=istartset(iset),iendset(iset)
            if(ialset(k).gt.0) then
               do l=0,ngraph-1
                  i=ialset(k)+l*nksegment
                  if(inum(i).le.0) cycle
                  if((ntrans.eq.0).or.(filabl(6:6).eq.'G').or.
     &                 (inotr(1,i).eq.0)) then
                     write(7,101) m1,i,(v(j,i),j=1,3)
                  else
                     call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                     write(7,101) m1,i,
     &                    v(1,i)*a(1,1)+v(2,i)*a(2,1)+v(3,i)*a(3,1),
     &                    v(1,i)*a(1,2)+v(2,i)*a(2,2)+v(3,i)*a(3,2),
     &                    v(1,i)*a(1,3)+v(2,i)*a(2,3)+v(3,i)*a(3,3)
                  endif
               enddo
            else
               l=ialset(k-2)
               do
                  l=l-ialset(k)
                  if(l.ge.ialset(k-1)) exit
                  do m=0,ngraph-1
                     i=l+m*nksegment
                     if(inum(i).le.0) cycle
                     if((ntrans.eq.0).or.(filabl(6:6).eq.'G').or.
     &                    (inotr(1,i).eq.0)) then
                        write(7,101) m1,i,(v(j,i),j=1,3)
                     else
                        call transformatrix(trab(1,inotr(1,i)),
     &                       co(1,i),a)
                        write(7,101) m1,i,
     &                       v(1,i)*a(1,1)+v(2,i)*a(2,1)+v(3,i)*a(3,1),
     &                       v(1,i)*a(1,2)+v(2,i)*a(2,2)+v(3,i)*a(3,2),
     &                       v(1,i)*a(1,3)+v(2,i)*a(2,3)+v(3,i)*a(3,3)
                     endif
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


