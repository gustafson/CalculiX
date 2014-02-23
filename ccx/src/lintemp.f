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
      subroutine lintemp(t0,t1,konl,nope,jj,t0l,t1l)
!
!     calculates a trilinear approximation to the quadratic interpolation
!     of the temperatures in a C3D20 element (full integration). A
!     quadratic interpolation of the temperatures leads to quadratic
!     thermal stresses, which cannot be handled by the elements 
!     displacement functions (which lead to linear stresses). Thus,
!     the temperatures are approximated by a trilinear function.
!
      implicit none
!
      integer konl(20),nope,jj,i1,i
!
      real*8 t0(*),t1(*),t0l,t1l,a20l(20,27)
!                        
c      data ((a20l(i,j),i=1,20),j=1,10) /
      a20l=reshape((
     &/-0.088729832,-0.240369600,-0.059630393,-0.240369600,-0.240369600,
     & -0.059630393,-0.011270159,-0.059630393, 0.524865555, 0.066666663,
     &  0.066666663, 0.524865555, 0.066666663, 0.008467776, 0.008467776,
     &  0.066666663, 0.524865555, 0.066666663, 0.008467776, 0.066666663,
     & -0.164549715,-0.164549715,-0.149999995,-0.149999995,-0.149999995,
     & -0.149999995,-0.035450279,-0.035450279, 0.524865544, 0.295766106,
     &  0.066666668, 0.295766106, 0.066666668, 0.037567223, 0.008467777,
     &  0.037567223, 0.295766106, 0.295766106, 0.037567223, 0.037567223,
     & -0.240369600,-0.088729832,-0.240369600,-0.059630393,-0.059630393,
     & -0.240369600,-0.059630393,-0.011270159, 0.524865555, 0.524865555,
     &  0.066666663, 0.066666663, 0.066666663, 0.066666663, 0.008467776,
     &  0.008467776, 0.066666663, 0.524865555, 0.066666663, 0.008467776,
     & -0.164549715,-0.149999995,-0.149999995,-0.164549715,-0.149999995,
     & -0.035450279,-0.035450279,-0.149999995, 0.295766106, 0.066666668,
     &  0.295766106, 0.524865544, 0.037567223, 0.008467777, 0.037567223,
     &  0.066666668, 0.295766106, 0.037567223, 0.037567223, 0.295766106,
     & -0.157274855,-0.157274855,-0.157274855,-0.157274855,-0.092725137,
     & -0.092725137,-0.092725137,-0.092725137, 0.295766105, 0.295766105,
     &  0.295766105, 0.295766105, 0.037567224, 0.037567224, 0.037567224,
     &  0.037567224, 0.166666664, 0.166666664, 0.166666664, 0.166666664,
     & -0.149999995,-0.164549715,-0.164549715,-0.149999995,-0.035450279,
     & -0.149999995,-0.149999995,-0.035450279, 0.295766106, 0.524865544,
     &  0.295766106, 0.066666668, 0.037567223, 0.066666668, 0.037567223,
     &  0.008467777, 0.037567223, 0.295766106, 0.295766106, 0.037567223,
     & -0.240369600,-0.059630393,-0.240369600,-0.088729832,-0.059630393,
     & -0.011270159,-0.059630393,-0.240369600, 0.066666663, 0.066666663,
     &  0.524865555, 0.524865555, 0.008467776, 0.008467776, 0.066666663,
     &  0.066666663, 0.066666663, 0.008467776, 0.066666663, 0.524865555,
     & -0.149999995,-0.149999995,-0.164549715,-0.164549715,-0.035450279,
     & -0.035450279,-0.149999995,-0.149999995, 0.066666668, 0.295766106,
     &  0.524865544, 0.295766106, 0.008467777, 0.037567223, 0.066666668,
     &  0.037567223, 0.037567223, 0.037567223, 0.295766106, 0.295766106,
     & -0.059630393,-0.240369600,-0.088729832,-0.240369600,-0.011270159,
     & -0.059630393,-0.240369600,-0.059630393, 0.066666663, 0.524865555,
     &  0.524865555, 0.066666663, 0.008467776, 0.066666663, 0.066666663,
     &  0.008467776, 0.008467776, 0.066666663, 0.524865555, 0.066666663,
     & -0.164549715,-0.149999995,-0.035450279,-0.149999995,-0.164549715,
     & -0.149999995,-0.035450279,-0.149999995, 0.295766106, 0.037567223,
     &  0.037567223, 0.295766106, 0.295766106, 0.037567223, 0.037567223,
     &  0.295766106, 0.524865544, 0.066666668, 0.008467777, 0.066666668,
     & -0.157274855,-0.157274855,-0.092725137,-0.092725137,-0.157274855,
     & -0.157274855,-0.092725137,-0.092725137, 0.295766105, 0.166666664,
     &  0.037567224, 0.166666664, 0.295766105, 0.166666664, 0.037567224,
     &  0.166666664, 0.295766105, 0.295766105, 0.037567224, 0.037567224,
     & -0.149999995,-0.164549715,-0.149999995,-0.035450279,-0.149999995,
     & -0.164549715,-0.149999995,-0.035450279, 0.295766106, 0.295766106,
     &  0.037567223, 0.037567223, 0.295766106, 0.295766106, 0.037567223,
     &  0.037567223, 0.066666668, 0.524865544, 0.066666668, 0.008467777,
     & -0.157274855,-0.092725137,-0.092725137,-0.157274855,-0.157274855,
     & -0.092725137,-0.092725137,-0.157274855, 0.166666664, 0.037567224,
     &  0.166666664, 0.295766105, 0.166666664, 0.037567224, 0.166666664,
     &  0.295766105, 0.295766105, 0.037567224, 0.037567224, 0.295766105,
     & -0.124999996,-0.124999996,-0.124999996,-0.124999996,-0.124999996,
     & -0.124999996,-0.124999996,-0.124999996, 0.166666664, 0.166666664,
     &  0.166666664, 0.166666664, 0.166666664, 0.166666664, 0.166666664,
     &  0.166666664, 0.166666664, 0.166666664, 0.166666664, 0.166666664,
     & -0.092725137,-0.157274855,-0.157274855,-0.092725137,-0.092725137,
     & -0.157274855,-0.157274855,-0.092725137, 0.166666664, 0.295766105,
     &  0.166666664, 0.037567224, 0.166666664, 0.295766105, 0.166666664,
     &  0.037567224, 0.037567224, 0.295766105, 0.295766105, 0.037567224,
     & -0.149999995,-0.035450279,-0.149999995,-0.164549715,-0.149999995,
     & -0.035450279,-0.149999995,-0.164549715, 0.037567223, 0.037567223,
     &  0.295766106, 0.295766106, 0.037567223, 0.037567223, 0.295766106,
     &  0.295766106, 0.066666668, 0.008467777, 0.066666668, 0.524865544,
     & -0.092725137,-0.092725137,-0.157274855,-0.157274855,-0.092725137,
     & -0.092725137,-0.157274855,-0.157274855, 0.037567224, 0.166666664,
     &  0.295766105, 0.166666664, 0.037567224, 0.166666664, 0.295766105,
     &  0.166666664, 0.037567224, 0.037567224, 0.295766105, 0.295766105,
     & -0.035450279,-0.149999995,-0.164549715,-0.149999995,-0.035450279,
     & -0.149999995,-0.164549715,-0.149999995, 0.037567223, 0.295766106,
     &  0.295766106, 0.037567223, 0.037567223, 0.295766106, 0.295766106,
     &  0.037567223, 0.008467777, 0.066666668, 0.524865544, 0.066666668,
     & -0.240369600,-0.059630393,-0.011270159,-0.059630393,-0.088729832,
     & -0.240369600,-0.059630393,-0.240369600, 0.066666663, 0.008467776,
     &  0.008467776, 0.066666663, 0.524865555, 0.066666663, 0.066666663,
     &  0.524865555, 0.524865555, 0.066666663, 0.008467776, 0.066666663,
     & -0.149999995,-0.149999995,-0.035450279,-0.035450279,-0.164549715,
     & -0.164549715,-0.149999995,-0.149999995, 0.066666668, 0.037567223,
     &  0.008467777, 0.037567223, 0.524865544, 0.295766106, 0.066666668,
     &  0.295766106, 0.295766106, 0.295766106, 0.037567223, 0.037567223,
     & -0.059630393,-0.240369600,-0.059630393,-0.011270159,-0.240369600,
     & -0.088729832,-0.240369600,-0.059630393, 0.066666663, 0.066666663,
     &  0.008467776, 0.008467776, 0.524865555, 0.524865555, 0.066666663,
     &  0.066666663, 0.066666663, 0.524865555, 0.066666663, 0.008467776,
     & -0.149999995,-0.035450279,-0.035450279,-0.149999995,-0.164549715,
     & -0.149999995,-0.149999995,-0.164549715, 0.037567223, 0.008467777,
     &  0.037567223, 0.066666668, 0.295766106, 0.066666668, 0.295766106,
     &  0.524865544, 0.295766106, 0.037567223, 0.037567223, 0.295766106,
     & -0.092725137,-0.092725137,-0.092725137,-0.092725137,-0.157274855,
     & -0.157274855,-0.157274855,-0.157274855, 0.037567224, 0.037567224,
     &  0.037567224, 0.037567224, 0.295766105, 0.295766105, 0.295766105,
     &  0.295766105, 0.166666664, 0.166666664, 0.166666664, 0.166666664,
     & -0.035450279,-0.149999995,-0.149999995,-0.035450279,-0.149999995,
     & -0.164549715,-0.164549715,-0.149999995, 0.037567223, 0.066666668,
     &  0.037567223, 0.008467777, 0.295766106, 0.524865544, 0.295766106,
     &  0.066666668, 0.037567223, 0.295766106, 0.295766106, 0.037567223,
     & -0.059630393,-0.011270159,-0.059630393,-0.240369600,-0.240369600,
     & -0.059630393,-0.240369600,-0.088729832, 0.008467776, 0.008467776,
     &  0.066666663, 0.066666663, 0.066666663, 0.066666663, 0.524865555,
     &  0.524865555, 0.066666663, 0.008467776, 0.066666663, 0.524865555,
     & -0.035450279,-0.035450279,-0.149999995,-0.149999995,-0.149999995,
     & -0.149999995,-0.164549715,-0.164549715, 0.008467777, 0.037567223,
     &  0.066666668, 0.037567223, 0.066666668, 0.295766106, 0.524865544,
     &  0.295766106, 0.037567223, 0.037567223, 0.295766106, 0.295766106,
     & -0.011270159,-0.059630393,-0.240369600,-0.059630393,-0.059630393,
     & -0.240369600,-0.088729832,-0.240369600, 0.008467776, 0.066666663,
     &  0.066666663, 0.008467776, 0.066666663, 0.524865555, 0.524865555,
     &  0.066666663, 0.008467776, 0.066666663, 0.524865555, 0.066666663/
     &  ),(/20,27/))
!
      do i1=1,nope
         t0l=t0l+a20l(i1,jj)*t0(konl(i1))
         t1l=t1l+a20l(i1,jj)*t1(konl(i1))
      enddo
!
      return
      end
