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
      subroutine applyboun(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon)
!
!     applies SPC and MPC conditions of the new increment to vold
!     and voldtu
!
      implicit none
!
      integer turbulent
!
      integer nrhcon(*),ielmat(*),ntmat_,nodeboun(*),isolidsurf(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:4,*),xbounact(*),shcon,
     &  voldtu(2,*),t1act(*),temp,r,dvi,xsolidsurf(*),reflength,
     &  refkin,reftuf,refvel,cp,voldaux(0:4,*),physcon(*)
!
!     inserting the boundary conditions (temperature, velocity
!     and pressure)
!
      do i=1,nboun
         node=nodeboun(i)
!
         vold(ndirboun(i),node)=xbounact(i)
!
!        update voldaux
!
         index=iponoel(node)
         if(index.le.0) cycle
         imat=ielmat(inoel(1,index))
         temp=vold(0,node)
         call materialdata_tg_sec(imat,ntmat_,temp,
     &        shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho,physcon)
!     
!     different treatment for gases and liquids
!     
         if(dabs(rho).lt.1.d-20) then
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            voldaux(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
c            write(*,*) 'iniiii2 ',node,vold(3,node),rho,voldaux(0,node)
         else
            voldaux(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
         endif
         voldaux(4,node)=rho
         do k=1,3
            voldaux(k,node)=rho*vold(k,node)
         enddo
c         write(*,*) 'bounnnn ',i,node
c         write(*,*) 'vold',(vold(k,node),k=0,4)
c         write(*,*) 'voldaux',(voldaux(k,node),k=0,4)
      enddo
c      write(*,*) 'vvvold ',(vold(i,307),i=0,4)
c      write(*,*) 'vvvoldaux ',(voldaux(i,307),i=0,4)
!
!     complete thermal field given
!
      if(ithermal.eq.1) then
         do i=1,nk
            if(iponoel(i).gt.0) then
               vold(0,i)=t1act(i)
            endif
         enddo
      endif
!
!     boundary conditions for turbulent kinetic energy and
!     turbulent frequency
!
      if(turbulent.eq.1) then
!
!        solid surface nodes
!
         do i=1,nsolidsurf
            node=isolidsurf(i)
            voldtu(1,i)=0.d0
            if(iponoel(node).eq.0) then
               write(*,*) '*ERROR in applyboun: no fluid element found'
               write(*,*) '       belonging to this node'
               stop
            endif
            imat=ielmat(inoel(1,iponoel(node)))
            temp=vold(0,node)
!
            call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,
     &          r,dvi,rhcon,nrhcon,rho)
!
!           separate treatment of gas and liquid
!            
            if(dabs(rho).lt.1.d-20) then
               voldtu(2,node)=800.d0*dvi/(voldaux(4,node)*
     &                           xsolidsurf(i)**2)
            else
               voldtu(2,node)=800.d0*dvi/(rho*xsolidsurf(i)**2)
            endif
!
         enddo
!
!        freestream nodes
!
!         check whether the freestream reference node belongs to
!        a fluid element
!
         if(iponoel(ifreenode).eq.0) then
            write(*,*) '*ERROR in applyboun: no fluid element found'
            write(*,*) '       belonging the freestream reference'
            write(*,*) '       node'
            stop
         endif
         imat=ielmat(inoel(1,iponoel(ifreenode)))
         temp=vold(0,ifreenode)
!     
         call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,
     &        r,dvi,rhcon,nrhcon,rho)
!     
         refvel=dsqrt(vold(1,ifreenode)**2+
     &        vold(2,ifreenode)**2+
     &        vold(3,ifreenode)**2)
!     
         reftuf=5.d0*refvel/reflength
         refkin=reftuf*dvi*1.d-3/voldaux(4,ifreenode)
!     
         do i=1,nfreestream
            node=ifreestream(i)
            voldtu(1,node)=refkin
            voldtu(2,node)=reftuf
         enddo
!     
      endif
c      write(*,*) 'after inserting the boundary conditions '
c      do i=1,nk
c         write(*,*) i,voldaux(1,i),voldaux(2,i),voldaux(3,i)
c      enddo
!     
      return
      end
