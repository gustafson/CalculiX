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
!
!   Subroutine pre_extrapolate.f
!
!      Interpolates xstate values for the new integration points 
!      at the beginning of the new increment. 
!
!   by: Jaro Hokkanen
!
!
      subroutine interpolatestate(ne,ipkon,kon,lakon,ne0,mi,xstate,
     &  pslavsurf,nstate_,xstateini,islavsurf,islavsurfold,pslavsurfold)
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),kon(*),ne,iflag,i,n,mi(*),indexc,ne0,indexcj,
     &  nstate_,kk,nopespring,iface,nopespringj,ifacej,ielemslave,ll,
     &  numpts,islavsurf(2,*),islavsurfold(2,*)
!
      real*8 xstate(nstate_,mi(1),*),pslavsurf(3,*),pslavsurfold(3,*),
     &  xstateini(nstate_,mi(1),*)
!
      data iflag /2/
!     
!     Loop over the elements
!     
      i=ne0
      do
         i=i+1
         if(i.gt.ne) exit
!     
         if(ipkon(i).le.-1) then
            cycle
         else
            indexc=ipkon(i)
         endif
!     
!     Only contact elements
!     
         if((lakon(i)(1:1).ne.'E').or.(lakon(i)(7:7).ne.'C')) cycle
         nopespring=kon(indexc)
         kk=kon(indexc+nopespring+2)
         iface=islavsurf(1,kk)
!     
         ielemslave=int(iface/10.d0)
         lakonl=lakon(ielemslave)
!     
         n=i
!     
         do
            n=n+1
            if(n.gt.ne) exit
            indexcj=ipkon(n)
            nopespringj=kon(indexcj)
            ll=kon(indexcj+nopespringj+2)
            ifacej=islavsurf(1,ll)
            if(kk.ne.ll) exit
         enddo
         n=n-1
!     
         numpts=islavsurfold(2,kk+1)-islavsurfold(2,kk)
         if(numpts.gt.2) then
            call interpolateinface(kon,ipkon,i,n,kk,
     &           xstate,xstateini,numpts,nstate_,mi,pslavsurf,
     &           ne0,islavsurfold,pslavsurfold)
         endif
!     
         i=n  
!     
      enddo
!     
      return
      end
      
