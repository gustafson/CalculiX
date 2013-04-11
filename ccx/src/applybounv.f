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
      subroutine applybounv(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &  compressible,ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,
     &  mi)
!
!     applies velocity boundary conditions
!
      implicit none
!
      integer turbulent,compressible
!
      integer nrhcon(*),ielmat(*),ntmat_,nodeboun(*),isolidsurf(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ismooth,indexi,nodei,nmpc,nodempc(3,*),ipompc(*),
     &  ist,ndir,ndiri,inomat(*),mi(2)
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:mi(2),*),xbounact(*),shcon,
     &  voldtu(2,*),t1act(*),temp,xsolidsurf(*),reflength,
     &  refkin,reftuf,refvel,voldaux(0:4,*),physcon(*),v(0:mi(2),*),
     &  rhoi,coefmpc(*),residu,size,correction
!
!     inserting the velocity boundary conditions
!
      do i=1,nboun
         if((ndirboun(i).lt.1).or.(ndirboun(i).gt.3)) cycle
!     
         node=nodeboun(i)
!
!        check whether fluid SPC
!
         imat=inomat(node)
         if(imat.eq.0) cycle
         if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
            temp=vold(0,node)
c            call materialdata_tg_sec(imat,ntmat_,temp,
c     &           shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho,physcon)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_)
         else
!     
!     determining rho from the solution field (for compressible fluids)
!     
            rho=voldaux(4,node)
         endif
         if(ismooth.eq.0) then
!     
!     in case of no smoothing (incompressible fluids or
!     pre-smoothing call of compressible fluids)
!     
            v(ndirboun(i),node)=xbounact(i)*rho
     &           -voldaux(ndirboun(i),node)
         else
!     
!     in case of smoothing: update voldaux (only for compressible
!     fluids)
!     
            voldaux(ndirboun(i),node)=xbounact(i)*rho
         endif
      enddo
!
!     taking velocity MPC's into account
!
      if(ismooth.eq.0) then
         do i=1,nmpc
            ist=ipompc(i)
            ndir=nodempc(2,ist)
            if((ndir.lt.1).or.(ndir.gt.3)) cycle
            node=nodempc(1,ist)
!
!           check whether fluid MPC
!
            imat=inomat(node)
            if(imat.eq.0) cycle
            if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
               temp=vold(0,node)
c               call materialdata_tg_sec(imat,ntmat_,temp,
c     &              shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho,physcon)
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_)
            else
!     
!     determining rho from the solution field (for compressible fluids)
!     
               rho=voldaux(4,node)+v(4,node)
            endif
!     
            index=nodempc(3,ist)
            residu=coefmpc(ist)*(voldaux(ndir,node)+v(ndir,node))/rho
            size=(coefmpc(ist)/rho)**2
            if(index.ne.0) then
               do
                  nodei=nodempc(1,index)
                  ndiri=nodempc(2,index)
!     
!     determining rho
!     
                  if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
                     imat=inomat(nodei)
                     if(imat.eq.0) cycle
                     temp=vold(0,nodei)
c                     call materialdata_tg_sec(imat,ntmat_,temp,
c     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rhoi,physcon)
                     call materialdata_rho(rhcon,nrhcon,imat,rhoi,
     &                    temp,ntmat_)
                  else
!     
!     determining rho from the solution field (for compressible fluids)
!     
                     rhoi=voldaux(4,nodei)+v(4,nodei)
                  endif
!     
                  residu=residu+coefmpc(index)*
     &                 (voldaux(ndiri,nodei)+v(ndiri,nodei))/rhoi
                  size=size+(coefmpc(index)/rhoi)**2
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
!
!           correcting all terms of the MPC
!
            residu=residu/size
!
            correction=-residu*coefmpc(ist)/rho
            v(ndir,node)=v(ndir,node)+correction
            index=nodempc(3,ist)
            if(index.ne.0) then
               do
                  nodei=nodempc(1,index)
                  ndiri=nodempc(2,index)
!     
!     determining rho
!     
                  if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
                     imat=inomat(nodei)
                     if(imat.eq.0) cycle
                     temp=vold(0,nodei)
c                     call materialdata_tg_sec(imat,ntmat_,temp,
c     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rhoi,physcon)
                     call materialdata_rho(rhcon,nrhcon,imat,rhoi,
     &                    temp,ntmat_)
                  else
!     
!     determining rho from the solution field (for compressible fluids)
!     
                     rhoi=voldaux(4,nodei)+v(4,nodei)
                  endif
!     
                  correction=-residu*coefmpc(index)/rhoi
                  v(ndiri,nodei)=v(ndiri,nodei)+correction
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
         enddo
      else
!
!        smoothing procedure: voldaux has already been updated
!
         do i=1,nmpc
            ist=ipompc(i)
            ndir=nodempc(2,ist)
            if((ndir.lt.1).or.(ndir.gt.3)) cycle
            node=nodempc(1,ist)
!
!           check whether fluid MPC
!
            imat=inomat(node)
            if(imat.eq.0) cycle
            if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
               temp=vold(0,node)
c               call materialdata_tg_sec(imat,ntmat_,temp,
c     &              shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho,physcon)
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_)
            else
!     
!     determining rho from the solution field (for compressible fluids)
!     
               rho=voldaux(4,node)+v(4,node)
            endif
!     
            index=nodempc(3,ist)
            residu=coefmpc(ist)*voldaux(ndir,node)/rho
            size=(coefmpc(ist)/rho)**2
            if(index.ne.0) then
               do
                  nodei=nodempc(1,index)
                  ndiri=nodempc(2,index)
!     
!     determining rho
!     
                  if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
                     imat=inomat(nodei)
                     if(imat.eq.0) cycle
                     temp=vold(0,nodei)
c                     call materialdata_tg_sec(imat,ntmat_,temp,
c     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rhoi,physcon)
                     call materialdata_rho(rhcon,nrhcon,imat,rhoi,
     &                    temp,ntmat_)
                  else
!     
!     determining rho from the solution field (for compressible fluids)
!     
                     rhoi=voldaux(4,nodei)
                  endif
!     
                  residu=residu+coefmpc(index)*
     &                 voldaux(ndiri,nodei)/rhoi
                  size=size+(coefmpc(index)/rhoi)**2
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
!
!           correcting all terms of the MPC
!
            residu=residu/size
!
            correction=-residu*coefmpc(ist)/rho
            voldaux(ndir,node)=voldaux(ndir,node)+correction
            index=nodempc(3,ist)
            if(index.ne.0) then
               do
                  nodei=nodempc(1,index)
                  ndiri=nodempc(2,index)
!     
!     determining rho
!     
                  if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
                     imat=inomat(nodei)
                     if(imat.eq.0) cycle
                     temp=vold(0,nodei)
c                     call materialdata_tg_sec(imat,ntmat_,temp,
c     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rhoi,physcon)
                     call materialdata_rho(rhcon,nrhcon,imat,rhoi,
     &                    temp,ntmat_)
                  else
!     
!     determining rho from the solution field (for compressible fluids)
!     
                     rhoi=voldaux(4,nodei)
                  endif
!     
                  correction=-residu*coefmpc(index)/rhoi
                  voldaux(ndiri,nodei)=voldaux(ndiri,nodei)+correction
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
         enddo
      endif
!     
      return
      end
      
