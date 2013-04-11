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
      subroutine readforce(zc,neq,nk,nev,nactdof,ikmpc,nmpc,ipompc,
     &  nodempc,mi,coefmpc)
!
!     reads a complex force (e.g. response of a fluid to a harmonic
!     structural excitation)
!     the force has to be stored in file 'dummy' in the form:
!
!     node,Fx-real,Fx-imag,Fy-real,Fy-imag,Fz-real,Fz-imag
!
!     only nodes in which a force is applied have to be stored. Modes
!     have to be separated by line starting with **
!
!     for cyclic symmetric structures the eigenmodes come in pairs
!     forces must be given for the first mode of each pair only
!
      implicit none
!
      integer mi(*),neq,nk,nev,i,j,k,nactdof(0:mi(2),*),ikmpc(*),nmpc,
     &  jdof,id,ist,ipompc(*),index,nodempc(3,*),node,istat
!
      real*8 coefmpc(*),comp(6)
!
      complex*16 zc(neq,*),force(3)
!
      open(27,file='dummy',status='unknown')
!
      do i=1,nev
c         do j=1,nk
         do
            read(27,*,iostat=istat) node,(comp(k),k=1,6)
            if(istat.ne.0) then
               exit
            endif
!
            do k=1,3
               force(k)=comp(2*k-1)*(1.d0,0.d0)+comp(2*k)*(0.d0,1.d0)
            enddo
!
            do k=1,3
               jdof=nactdof(k,node)
               if(jdof.ne.0) then
                  zc(jdof,i)=zc(jdof,i)-force(k)
               else
!     
!     node is a dependent node of a MPC: distribute
!     the forces among the independent nodes
!     (proportional to their coefficients)
!     
                  jdof=8*(node-1)+k
                  call nident(ikmpc,jdof,nmpc,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.jdof) then
                        ist=ipompc(id)
                        index=nodempc(3,ist)
                        if(index.eq.0) cycle
                        do
                           jdof=nactdof(nodempc(2,index),
     &                                  nodempc(1,index))
                           if(jdof.ne.0) then
                              zc(jdof,i)=zc(jdof,i)-
     &                             coefmpc(index)*force(k)/coefmpc(ist)
                           endif
                           index=nodempc(3,index)
                           if(index.eq.0) exit
                        enddo
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
!
      close(27)
!     
      return
      end
      
