!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2013 Guido Dhondt
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
c> soubroutine to get mu for current contact pair  
      subroutine getmu(mu,tietol,elcon,itie,ncmat_,ntmat_)
!     autor: Saskia Sitzmann 
       
       implicit none
!       
       integer itie,imat,ncmat_,ntmat_
!
       real*8 mu,tietol(2,*),elcon(0:ncmat_,ntmat_,*)
!       
       itie=itie+1
       imat=int(tietol(2,itie))
       if(ncmat_.lt.6)then
        mu=0.0
       else
        mu=elcon(6,1,imat)
       endif
       itie=itie-1
!
       return
      end
