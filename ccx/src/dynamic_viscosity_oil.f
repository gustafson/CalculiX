!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine dynamic_viscosity_oil(K_oil,Temp,rho_oil,nue_oil,
     &     dvi_oil)
!
!     this subroutine computes the different parameters for oil
!     source : Phd Thesis Roland Fischer Karlsruhe 1998
!     "Zweiphasenströmungen in Triebwerksleitungen - 
!     Theoretische und Experimentelle Untersuchubng von 
!     Luft/Öl Strömungen durch Blenden  "
!     Anhang A.3 Stoffwerte des verwendeten Triebwerköls 'MOBIL JET2'
!
      integer k_oil
!
      real*8 Temp,rho_oil,nue_oil,dvi_oil
!
      k_oil=k_oil
!
!     density [kg/m**3]
      rho_oil= 917D0-0.69d0*(Temp-273.15d0)
!
!     kinematic viscosity [m*2/s]
      nue_oil=0.00219d0*(Temp-273.15d0)**(-1.794d0)
!
!     dynamic viscosity [kg/(m*s)]
      dvi_oil=nue_oil/rho_oil

      return
      end
