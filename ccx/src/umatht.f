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
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     &  statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     &  cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     &  noel,npt,layer,kspt,kstep,kinc)
!
!     heat transfer material subroutine
!
!     INPUT:
!
!     statev(nstatv)      internal state variables at the start
!                         of the increment
!     temp                temperature at the start of the increment
!     dtemp               increment of temperature
!     dtemdx(ntgrd)       current values of the spatial gradients of the 
!                         temperature
!     time(1)             step time at the beginning of the increment
!     time(2)             total time at the beginning of the increment
!     dtime               time increment
!     predef              not used
!     dpred               not used
!     cmname              material name
!     ntgrd               number of spatial gradients of temperature
!     nstatv              number of internal state variables as defined
!                         on the *DEPVAR card
!     props(nprops)       user defined constants defined by the keyword
!                         card *USER MATERIAL,TYPE=THERMAL
!     nprops              number of user defined constants, as specified
!                         on the *USER MATERIAL,TYPE=THERMAL card
!     coords              global coordinates of the integration point
!     pnewd               not used
!     noel                element number
!     npt                 integration point number
!     layer               not used
!     kspt                not used
!     kstep               not used
!     kinc                not used
!
!     OUTPUT:
!     
!     u                   not used
!     dudt                not used
!     dudg(ntgrd)         not used
!     flux(ntgrd)         heat flux at the end of the increment
!     dfdt(ntgrd)         not used
!     dfdg(ntgrd,ntgrd)   variation of the heat flux with respect to the 
!                         spatial temperature gradient
!     statev(nstatv)      internal state variables at the end of the
!                         increment
!
      implicit none
!
      character*80 cmname
!
      integer ntgrd,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc
!
      real*8 u,dudt,dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     &  statev(nstatv),pnewdt,temp,dtemp,dtemdx(ntgrd),time(2),dtime,
     &  predef,dpred,props(nprops),coords(3),dfdg(ntgrd,ntgrd)
!
      return
      end
