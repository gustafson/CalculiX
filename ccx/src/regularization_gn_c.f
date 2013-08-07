!
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
c> regularization function for normal contactmortar
c>
c> @param [in]       lambdap   contact pressure in normal direction
c> @param [in]       divmode indicates whether funtion or derivate 
c>                             should be called
c>                    =0 function called
c>                    =1 derivative called    
c> @param [in]       regmode        selects regularization funtion
c>                    =1 perturbed Lagrange
c>                    =2 piece wise linear with given data points
c>                    =3 expontential contact law
c> @param [out]       gnc        result regularization function
c> @param [in]   aninvloc        stiffness constant for perturbed Lagrange
      subroutine regularization_gn_c(lambdap,divmode,regmode,
     &           gnc,aninvloc,p0,beta,elcon,nelcon,itie,ntmat_,
     &           plicon,nplicon,npmat_,ncmat_,tietol)
!
!       regularization function for normal contact
!       Author: Saskia Sitzmann
!     
      implicit none
!
      integer divmode,regmode,i,ndata,kode,npmat_,ncmat_,
     &  itie,ntmat_,nplicon(0:ntmat_,*),nelcon(2,*),
     &  imat
!
!
      real*8 lambdap,gnc,pn_d(40),gn_d(40),aninv1(40),t(40),
     &  beta,p0,aninvloc,elconloc(21),plconloc(802),t1l,
     &  elcon(0:ncmat_,ntmat_,*),plicon(0:2*npmat_,ntmat_,*),
     &  tietol(2,*)

      kode=-51
      t1l=0.0
      imat=int(tietol(2,itie+1))
!
c      beta=0.01
c      A=0.001
! hertzhueeber2d
c      A=1.0
c      beta=1.0/(log(100.0)/0.05)
! dbf 1
c       A=1.0
c       beta=0.0007
! dbf 2
c       A=1.0
c       beta=0.00009
! dbf 3
c       A=1.0
c       beta=0.00023
!
c      ndata=3
c      gn_d=(/0.0, 0.01, 0.0109, 0.0109, 0.0,0.0/)
c      pn_d=(/0.0, 100.0, 1000.0, 2000.0, 0.0,0.0/)
! contact law vera
c      ndata=8
c      gn_d=(/0.0, 0.0006, 0.0008, 0.001, 0.0012, 0.0014,
c     &      0.00155, 0.00155/)
c      pn_d=(/0.0, 0.2, 0.33, 0.56, 0.967, 1.73, 3.0, 10.0/)
! contact law vera*300 for dbf_coarse
c      ndata=8
c      do i=1, ndata
c        pn_d(i)=pn_d(i)*300
c      enddo
c      gn_d=(/0.0, 0.5, 1.0, 1.25, 1.5, 1.5/)
c      pn_d=(/0.0, 0.2, 0.5, 1.0, 2.5, 3.0/)
      gnc=0.0
!
!     perturbed Lagrange
!
      if(regmode.eq.1)then
         if(divmode.eq.0)then
            gnc=aninvloc*lambdap
         elseif(divmode.eq.1)then
            gnc=aninvloc
         else
            write(*,*)'error in regularzation_gn_c.f!'
            stop
         endif
!
!     multiple perturbed Lagrange
!
      else if(regmode.eq.2)then
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
         ndata=int(plconloc(801))
c          write(*,*)'ndata',ndata
         do i=1,ndata
          gn_d(i)=plconloc(2*i-1)
          pn_d(i)=plconloc(2*i)
c          write(*,*)'pn,gn', pn_d(i), gn_d(i)
         enddo 
         do i=1,ndata-1
            aninv1(i)=(gn_d(i+1)-gn_d(i))/(pn_d(i+1)-pn_d(i))
            t(i)=gn_d(i)-aninv1(i)*pn_d(i)
c         write(*,*)'mi',aninv1(i),'ti',t(i)
         enddo
         if(divmode.eq.0)then
            if(lambdap.lt.pn_d(1))then
               gnc=aninv1(1)*lambdap+t(1)
            endif
            do i=1,ndata-2
c               write(*,*) 'lp',lambdap,'pi',pn_d(i+1)
               if(pn_d(i+1).gt.lambdap.and.pn_d(i).le.lambdap)then
                  gnc=aninv1(i)*lambdap+t(i)
               endif
            enddo
            if(pn_d(ndata-1).le.lambdap)then
               gnc=aninv1(ndata-1)*lambdap+t(ndata-1)
            endif
         elseif(divmode.eq.1)then
            if(lambdap.lt.pn_d(1))then
               gnc=aninv1(1)
            endif
            do i=1,ndata-2
               if(pn_d(i+1).gt.lambdap.and.pn_d(i).le.lambdap)then
                  gnc=aninv1(i)
               endif
            enddo
            if(pn_d(ndata-1).le.lambdap)then
               gnc=aninv1(ndata-1)
            endif
         else
            write(*,*)'error in regularzation_gn_c.f!'
            stop
         endif
!
!     exponetial perturbed Lagrange
!
      else if (regmode.eq.3)then
         if(divmode.eq.0)then
            gnc=beta*log(lambdap+p0)-beta*log(p0)
         elseif(divmode.eq.1)then
            gnc=beta/(lambdap+p0)
         else
            write(*,*)'error in regularzation_gn_c.f!'
            stop
         endif
      endif
c      stop
      return
      end
