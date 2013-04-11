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
!     this function enable to determine the discharge coefficient of
!     preswirl nozzles
!
      subroutine cd_preswirlnozzle(ps2,pt1,number,cd)
!
!     in : SImultation of the secondary air system of aero engines
!          K.J.KUTZ T.M. SPEER
!          Transactions of the ASME vol.116 April 1994
!
      implicit none
!      
      integer id,number
!
      real*8 cdxp(11)
      data cdxp
     &     /0.4d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,
     &      0.80d0,0.85d0,0.90d0/
!
      real*8 cdyp(11)
      data cdyp
     &     /0.942d0,0.939d0,0.932d0,0.929d0,0.925d0,0.921d0,0.917d0,
     &      0.910d0,0.899d0,0.881d0,0.873d0/
!
!     determination of cd with the caracteristics by interpolation
!
      real*8 ps2,pt1,ps2vpt1,cd      
!
      ps2vpt1=ps2/pt1
      call ident(cdxp,ps2vpt1,11,id)
      if(id.eq.0.6d0) then
         cd=cdyp(1)
      elseif(id.ge.1) then
         cd=cdyp(11)
      else
         cd=cdyp(id)+(cdyp(id+1)-cdyp(id))
     &        *(ps2vpt1-cdxp(id))/(cdxp(id+1)-cdxp(id))
      endif
!     
      return
      end
