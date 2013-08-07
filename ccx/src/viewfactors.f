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
      subroutine viewfactors(textpart,iviewfile,istep,inpc,
     &  istat,n,key,iline,ipol,inl,ipoinp,inp,jobnamec,ipoinpc)
!
!     reading the input deck: *VIEWFACTOR 
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16),jobnamec(*)
!
      integer i,iviewfile,istep,n,istat,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),key,j,k,l,ipoinpc(0:*)
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in viscos: *VISCO can only be used'
         write(*,*) '  within a STEP'
         stop
      endif
!
      do i=2,n
         if(textpart(i)(1:4).eq.'READ') then
            if(iviewfile.eq.0) then
               iviewfile=-1
            else
               iviewfile=-abs(iviewfile)
            endif
         elseif(textpart(i)(1:9).eq.'WRITEONLY') then
            if(iviewfile.eq.0) then
               iviewfile=3
            else
               iviewfile=3*iviewfile/abs(iviewfile)
            endif
         elseif(textpart(i)(1:5).eq.'WRITE') then
            if(iviewfile.eq.0) then
               iviewfile=2
            else
               iviewfile=2*iviewfile/abs(iviewfile)
            endif
         elseif(textpart(i)(1:6).eq.'INPUT=') then
            jobnamec(2)(1:126)=textpart(i)(7:132)
            jobnamec(2)(127:132)='      '
            loop1: do j=1,126
               if(jobnamec(2)(j:j).eq.'"') then
                  do k=j+1,126
                     if(jobnamec(2)(k:k).eq.'"') then
                        do l=k-1,126
                           jobnamec(2)(l:l)=' '
                           exit loop1
                        enddo
                     endif
                     jobnamec(2)(k-1:k-1)=jobnamec(2)(k:k)
                  enddo
                  jobnamec(2)(126:126)=' '
               endif
            enddo loop1
         elseif(textpart(i)(1:7).eq.'OUTPUT=') then
            jobnamec(3)(1:125)=textpart(i)(8:132)
            jobnamec(3)(126:132)='      '
            loop2: do j=1,125
               if(jobnamec(3)(j:j).eq.'"') then
                  do k=j+1,125
                     if(jobnamec(3)(k:k).eq.'"') then
                        do l=k-1,125
                           jobnamec(3)(l:l)=' '
                           exit loop2
                        enddo
                     endif
                     jobnamec(3)(k-1:k-1)=jobnamec(3)(k:k)
                  enddo
                  jobnamec(3)(125:125)=' '
               endif
            enddo loop2
         else
            write(*,*) 
     &        '*WARNING in viewfactors: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

