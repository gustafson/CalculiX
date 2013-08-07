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
      subroutine applybount(nodeboun,ndirboun,nboun,xbounact,
     &  iponoel,vold,ipompc,nodempc,coefmpc,nmpc,inomat,matname,
     &  nshcon,shcon,nrhcon,rhcon,physcon,compressible,ntmat_,
     &  vcon,mi,ithermal)
!
!     applies temperature boundary conditions
!
      implicit none
!
      character*80 matname(*)
!
      integer nodeboun(*),ndirboun(*),i,nboun,node,iponoel(*),
     &  index,nodei,ndiri,ist,ipompc(*),nodempc(3,*),nmpc,
     &  ndir,inomat(*),imat,nshcon(*),nrhcon(*),k,compressible,
     &  ntmat_,mi(*),ithermal
!
      real*8 vold(0:mi(2),*),xbounact(*),residu,size,coefmpc(*),
     &  correction,
     &  temp,shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),physcon(*),
     &  cp,r,rho,vcon(0:4,*)
!
!     inserting the temperature boundary conditions
!
      do i=1,nboun
         if(ndirboun(i).ne.0) cycle
!
         node=nodeboun(i)
         if(inomat(node).eq.0) cycle
         vold(0,node)=xbounact(i)
!
!        update the conservative variables
!
         imat=inomat(node)
         temp=vold(0,node)
         call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &        nshcon,cp,physcon)
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
               write(*,*) '*ERROR in applybount: absolute temperature '
               write(*,*) '       is nearly zero; maybe absolute zero '
               write(*,*) '       was wrongly defined or not defined'
               write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
               stop
            endif
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
         endif
         vcon(4,node)=rho
         do k=1,3
            vcon(k,node)=rho*vold(k,node)
         enddo
      enddo
!
!     taking fluid temperature MPC's into account
!
      do i=1,nmpc
         ist=ipompc(i)
         ndir=nodempc(2,ist)
         if(ndir.ne.0) cycle
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
! 
         imat=inomat(node)
         if(imat.eq.0) cycle
c         index=iponoel(node)
c         if(index.le.0) cycle
!     
         index=nodempc(3,ist)
         residu=coefmpc(ist)*vold(ndir,node)
         size=coefmpc(ist)**2
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!     
               residu=residu+coefmpc(index)*vold(ndiri,nodei)
               size=size+coefmpc(index)**2
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
!     
!     correcting all terms of the MPC
!     
         residu=residu/size
!     
         correction=-residu*coefmpc(ist)
         vold(0,node)=vold(0,node)+correction
!
!        update the conservative variables
!
         imat=inomat(node)
         temp=vold(0,node)
         call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &        nshcon,cp,physcon)
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
               write(*,*) '*ERROR in applybount: absolute temperature '
               write(*,*) '       is nearly zero; maybe absolute zero '
               write(*,*) '       was wrongly defined or not defined'
               write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
               stop
            endif
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
         endif
         vcon(4,node)=rho
         do k=1,3
            vcon(k,node)=rho*vold(k,node)
         enddo
!
         index=nodempc(3,ist)
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!     
               correction=-residu*coefmpc(index)
               vold(ndiri,nodei)=vold(ndiri,nodei)+correction
!
!     update the conservative variables
!     
               imat=inomat(nodei)
               temp=vold(0,nodei)
               call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &              nshcon,cp,physcon)
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
                  if(vold(0,nodei)-physcon(1).le.1.d-10) then
                     write(*,*) 
     &                   '*ERROR in applybount: absolute temperature '
                     write(*,*) 
     &                   '       is nearly zero; maybe absolute zero '
                     write(*,*) 
     &                   '       was wrongly defined or not defined'
                     write(*,*) 
     &                   '       at all (*PHYSICAL CONSTANTS card)'
                     stop
                  endif
                  rho=vold(4,nodei)/(r*(vold(0,nodei)-physcon(1)))
                  vcon(0,nodei)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,nodei)**2+vold(2,nodei)**2+vold(3,nodei)**2)
     &                 /2.d0)-vold(4,nodei)
               else
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_,ithermal)
                  vcon(0,nodei)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,nodei)**2+vold(2,nodei)**2+vold(3,nodei)**2)
     &                 /2.d0)
               endif
               vcon(4,nodei)=rho
               do k=1,3
                  vcon(k,nodei)=rho*vold(k,nodei)
               enddo
!     
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
      enddo
!     
      return
      end
