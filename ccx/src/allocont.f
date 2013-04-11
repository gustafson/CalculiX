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
      subroutine allocont(ncont,ntie,tieset,nset,set,istartset,
     &  iendset,ialset,lakon,ncone)
!
!     counting the number of triangles needed for the 
!     triangulation of the contact master surfaces
!
      implicit none
!
      character*8 lakon(*)
      character*81 tieset(3,*),mastset,set(*),slavset
!
      integer ncont,ntie,i,j,k,nset,istartset(*),iendset(*),ialset(*),
     &  imast,nelem,jface,ncone,islav
!
!     number of master triangles
!
      ncont=0
!
!     number of slave nodes
!
      ncone=0
!
      do i=1,ntie
!
!        check for contact conditions
!
         if(tieset(1,i)(81:81).eq.'C') then
            mastset=tieset(3,i)
!
!           determining the master surface
!
            do j=1,nset
               if(set(j).eq.mastset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in triangucont: master surface'
               write(*,*) '       does not exist'
               stop
            endif
            imast=j
!
            do j=istartset(imast),iendset(imast)
               if(ialset(j).gt.0) then
c                  if(j.gt.istartset(imast)) then
c                     if(ialset(j).eq.ialset(j-1)) cycle
c                  endif
!
                  nelem=int(ialset(j)/10.d0)
                  jface=ialset(j)-10*nelem
!
                  if(lakon(nelem)(4:4).eq.'2') then
                     ncont=ncont+6
                  elseif(lakon(nelem)(4:4).eq.'8') then
                     ncont=ncont+2
                  elseif(lakon(nelem)(4:5).eq.'10') then
                     ncont=ncont+4
                  elseif(lakon(nelem)(4:4).eq.'4') then
                     ncont=ncont+1
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     if(jface.le.2) then
                        ncont=ncont+4
                     else
                        ncont=ncont+6
                     endif
                  elseif(lakon(nelem)(4:4).eq.'6') then
                     if(jface.le.2) then
                        ncont=ncont+1
                     else
                        ncont=ncont+2
                     endif
                  endif
!
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
!
                     nelem=int(k/10.d0)
                     jface=k-10*nelem
!
                     if(lakon(nelem)(4:4).eq.'2') then
                        ncont=ncont+6
                     elseif(lakon(nelem)(4:4).eq.'8') then
                        ncont=ncont+2
                     elseif(lakon(nelem)(4:5).eq.'10') then
                        ncont=ncont+4
                     elseif(lakon(nelem)(4:4).eq.'4') then
                        ncont=ncont+1
                     elseif(lakon(nelem)(4:5).eq.'15') then
                        if(jface.le.2) then
                           ncont=ncont+4
                        else
                           ncont=ncont+6
                        endif
                     elseif(lakon(nelem)(4:4).eq.'6') then
                        if(jface.le.2) then
                           ncont=ncont+1
                        else
                           ncont=ncont+2
                        endif
                     endif
!     
                  enddo
               endif
            enddo
!
!           counting the slave nodes
!
            slavset=tieset(2,i)
!
!           determining the master surface
!
            do j=1,nset
               if(set(j).eq.slavset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in triangucont: master surface'
               write(*,*) '       does not exist'
               stop
            endif
            islav=j
!
            do j=istartset(islav),iendset(islav)
               if(ialset(j).gt.0) then
                  ncone=ncone+1
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     ncone=ncone+1
                  enddo
               endif
            enddo
!
         endif
      enddo
!
      return
      end

