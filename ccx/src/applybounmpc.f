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
      subroutine applybounmpc(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &  compressible,ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,
     &  mi,matname)
!
!     applies temperature, pressure and velocity boundary conditions
!
      implicit none
!
      character*80 matname(*)
!
      integer turbulent,compressible,
     &  nrhcon(*),mi(*),ntmat_,nodeboun(*),isolidsurf(*),j,
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ismooth,indexi,nodei,nmpc,nodempc(3,*),ipompc(*),
     &  ist,ndir,ndiri,inomat(*),ielmat(mi(3),*)
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:mi(2),*),xbounact(*),
     &  vcontu(2,*),t1act(*),temp,xsolidsurf(*),reflength,r,cp,
     &  refkin,reftuf,refvel,vcon(0:4,*),physcon(*),v(0:mi(2),*),
     &  rhoi,coefmpc(*),residu,size,correction,shcon(0:3,ntmat_,*)
!
!     inserting SPC's
!
      do i=1,nboun
!     
         node=nodeboun(i)
!
!        check whether fluid SPC
!
         imat=inomat(node)
         if(imat.eq.0) cycle
!
!        storing the single point constraint
!
         vold(ndirboun(i),node)=xbounact(i)
!
!        deleting the conservative variables dependent on 
!        the SPC:
!        - if temperature: rho, rho*v_i, rho*eps_t
!        - if velocity v_j: rho*v_j
!        - if pressure and compressible: rho, rho*v_i, rho*eps_t
!
         if(ndirboun(i).eq.0) then
            do j=0,4
               vcon(j,node)=0.d0
            enddo
         elseif(ndirboun(i).lt.4) then
            vcon(ndirboun(i),node)=0.d0
         elseif(compressible.eq.0) then
            do j=0,4
               vcon(j,node)=0.d0
            enddo
         endif
!
      enddo
!
!     taking MPC's into account
!
      do i=1,nmpc
         ist=ipompc(i)
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
!     
         ndir=nodempc(2,ist)
         residu=0.d0
!     
         index=nodempc(3,ist)
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!     
               residu=residu+coefmpc(index)*vold(ndiri,nodei)
!     
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
         vold(ndir,node)=-residu/coefmpc(ist)
!     
!        deleting the conservative variables dependent on 
!        the SPC:
!        - if temperature: rho, rho*v_i, rho*eps_t
!        - if velocity v_j: rho*v_j
!        - if pressure and compressible: rho, rho*v_i, rho*eps_t
!
         if(ndir.eq.0) then
            do j=0,4
               vcon(j,node)=0.d0
            enddo
         elseif(ndir.lt.4) then
            vcon(ndir,node)=0.d0
         elseif(compressible.eq.0) then
            do j=0,4
               vcon(j,node)=0.d0
            enddo
         endif
!     
      enddo
!
!     update the conservative variables
!
!     SPC's
!
      do i=1,nboun
         node=nodeboun(i)
!
!        check whether fluid SPC
!
         imat=inomat(node)
         if(imat.eq.0) cycle
!
         if((vcon(0,node).eq.0.d0).or.
     &      ((vcon(4,node).eq.0.d0).and.(compressible.eq.1))) then
!
!           all conservative variables have to be updated
!
            temp=vold(0,node)
            call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &           nshcon,cp,physcon)
!     
!     different treatment for gases and liquids
!     
            if(compressible.eq.1) then
               r=shcon(3,1,imat)
               if(r.lt.1.d-10) then
                  write(*,*) '*ERROR in applybount: specific gas '
                  write(*,*) 'constant for material ',matname(imat)
                  write(*,*) 'is close to zero; maybe it has'
                  write(*,*) 'not been defined'
                  stop
               endif
               if(vold(0,node)-physcon(1).le.1.d-10) then
                  write(*,*)
     &               '*ERROR in applybount: absolute temperature '
                  write(*,*)
     &               '       is nearly zero; maybe absolute zero '
                  write(*,*) '       was wrongly defined or not defined'
                  write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
                  stop
               endif
               rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
               vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &              /2.d0)-vold(4,node)
            else
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_,ithermal)
               vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &              /2.d0)
            endif
            vcon(4,node)=rho
            do k=1,3
               vcon(k,node)=rho*vold(k,node)
            enddo
         elseif(vcon(1,node).eq.0.d0) then
            vcon(1,node)=vcon(4,node)*vold(1,node)
         elseif(vcon(2,node).eq.0.d0) then
            vcon(2,node)=vcon(4,node)*vold(2,node)
         elseif(vcon(3,node).eq.0.d0) then
            vcon(3,node)=vcon(4,node)*vold(3,node)
         endif
      enddo
!
!     MPC's
!
      do i=1,nmpc
         ist=ipompc(i)
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
!         
         if((vcon(0,node).eq.0.d0).or.
     &      ((vcon(4,node).eq.0.d0).and.(compressible.eq.1))) then
!
!           all conservative variables have to be updated
!
            temp=vold(0,node)
            call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &           nshcon,cp,physcon)
!     
!     different treatment for gases and liquids
!     
            if(compressible.eq.1) then
               r=shcon(3,1,imat)
               if(r.lt.1.d-10) then
                  write(*,*) '*ERROR in applybount: specific gas '
                  write(*,*) 'constant for material ',matname(imat)
                  write(*,*) 'is close to zero; maybe it has'
                  write(*,*) 'not been defined'
                  stop
               endif
               if(vold(0,node)-physcon(1).le.1.d-10) then
                  write(*,*) 
     &               '*ERROR in applybount: absolute temperature '
                  write(*,*) 
     &               '       is nearly zero; maybe absolute zero '
                  write(*,*) '       was wrongly defined or not defined'
                  write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
                  stop
               endif
               rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
               vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &              /2.d0)-vold(4,node)
            else
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_,ithermal)
               vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &              /2.d0)
            endif
            vcon(4,node)=rho
            do k=1,3
               vcon(k,node)=rho*vold(k,node)
            enddo
         elseif(vcon(1,node).eq.0.d0) then
            vcon(1,node)=vcon(4,node)*vold(1,node)
         elseif(vcon(2,node).eq.0.d0) then
            vcon(2,node)=vcon(4,node)*vold(2,node)
         elseif(vcon(3,node).eq.0.d0) then
            vcon(3,node)=vcon(4,node)*vold(3,node)
         endif
      enddo
!     
      return
      end
      
