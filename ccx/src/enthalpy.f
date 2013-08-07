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
!     This function computes the specific enthalpy for air as a function 
!     of the temperature
!
!     The values included in th is function as well as the methodology 
!     used to find them can be found in 
!
!     Y. Park and R.e. Sonntag
!     "Thermodynamic porperties of ideal gas air at 0.1MPa"
!      International Journal of Energy Research, vol 20,771-785(1996)

!     N.B.:
!     The temperatures range from 100K to 4000K
!
!     Enthalpy values in KJ/Kg in the document are given in J/Kg
!
!     Enthalpy values are given for a reference pressure P=1 bar
!     which explains why pressure has no influence on the result
!
      subroutine enthalpy(T,P,FARB,FARU,WAR,H) 
!
!     Input  :   T     TEMPERATURE (K)
!                P     PRESSURE (not used)
!                FARB  BURNT FUEL-TO-AIR RATIO (not used)
!                FARU  UNBURNT FUEL-TO-AIR RATIO (not used)
!                WAR   WATER-TO-AIR-RATIO (not used)
!     Result:    H     ENTHALPY (J/Kg)
!
!
      implicit none
!
      integer id,n40
!
      real*8 T,P,FARB,FARU,WAR, H
!
      real*8 tab_temperature(40)
      data tab_temperature/
     &     100 ,200 ,300 ,400 ,500 ,600 ,700 ,800 ,900 ,1000,
     &     1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     &     2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     &     3100,3200,3300,3400,3500,3600,3700,3800,3900,4000/
!
      real*8 tab_enthalpy(40)
      data tab_enthalpy/
     &     99972  ,200173 ,300473 ,401290 ,503345 ,607303 ,713551 ,
     &     822193 ,933140 ,1046202,1161158,1277784,1395877,1515258,
     &     1635776,1757301,1879726,2002959,2126924,2251557,2376802,
     &     2502611,2628944,2755764,2883040,3010742,3138846,3267327,
     &     3396165,3525340,3654835,3784631,3914716,4045073,4175691,
     &     4306558,4437661,4568993,4700543,4832303/
!
      data n40 /40/
!
      FARB=FARB
      FARU=FARU
      WAR=WAR
      P=P
!
!     linear interpolation
!
      call ident(tab_temperature,T,n40,id)
!
      if(id.le.1) then
         H=tab_enthalpy(1)
      elseif(id.ge.40) then
         H=tab_enthalpy(40)
      else
         H=tab_enthalpy(id)+
     &        (tab_enthalpy(id+1)-tab_enthalpy(id))
     &        *(T-tab_temperature(id))
     &        /(tab_temperature(id+1)-tab_temperature(id))
      endif
!
      write(*,*) 'enthalpy=',H
!
      return
      end

