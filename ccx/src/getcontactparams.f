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
      subroutine getcontactparams(mu,regmode,fkinv,p0,beta,tietol,elcon,
     &          itie,ncmat_,ntmat_)
!     autor: Saskia Sitzmann 
       
       implicit none
!       
       integer itie,imat,ncmat_,ntmat_,regmode
!
       real*8 mu,fkinv,p0,beta,tietol(2,*),elcon(0:ncmat_,ntmat_,*)
!       
       itie=itie+1
       imat=int(tietol(2,itie))
       if(ncmat_.lt.6)then
          mu=0.0
       else
          mu=elcon(6,1,imat)
       endif
c       write(*,*) 'el',elcon(3,1,imat),ncmat_
!
!      exponetial regularization
       if(ncmat_.gt.2)then
          if(elcon(3,1,imat).gt.1.4 .and.
     &         elcon(3,1,imat).lt.1.6 )then
             regmode=3
             p0=elcon(2,1,imat)
             beta=1.0/(elcon(1,1,imat))
             fkinv=0.0
             if(mu.gt.1.e-10)then
                write(*,*)'getcontactparams:'
                write(*,*)'exponential pressure overclosure',
     &               'with friction not yet supportes'
                stop
             endif
!     
!     linear regularization
          else if(elcon(3,1,imat).gt.2.4 .and. 
     &            elcon(3,1,imat).lt.2.6 )then
             regmode=1
             fkinv=1.0/elcon(2,1,imat)
             p0=0.0
             beta=0.0
!     
!     piecewiese linear regularization
          else if(elcon(3,1,imat).gt.3.4 .and.
     &            elcon(3,1,imat).lt.3.6 )then
             regmode=2
             p0=0.0
             beta=0.0
             fkinv=0.0
          else
             regmode=1
             fkinv=0.0
             p0=0.0
             beta=0.0
          endif
       else
          regmode=1
          fkinv=0.0
          p0=0.0
          beta=0.0
       endif 
       itie=itie-1
!
       return
      end
