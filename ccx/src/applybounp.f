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
      subroutine applybounp(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,voldtu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  voldcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &  ipompc,nodempc,coefmpc,nmpc,inomat,mi)
!
!     applies velocity boundary conditions
!
      implicit none
!
      integer turbulent
!
      integer nrhcon(*),ielmat(*),ntmat_,nodeboun(*),isolidsurf(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ipompc(*),nodempc(3,*),nmpc,ist,ndir,inomat(*),ndiri,
     &  nodei,mi(2)
!
      real*8 rhcon(0:1,ntmat_,*),vold(0:mi(2),*),xbounact(*),shcon,
     &  voldtu(2,*),t1act(*),temp,r,dvi,xsolidsurf(*),reflength,
     &  refkin,reftuf,refvel,cp,voldcon(0:4,*),physcon(*),v(0:mi(2),*),
     &  coefmpc(*),fixed_pres,size,correction,residu
!
!     inserting the pressure boundary conditions
!
      do i=1,nboun
         if(ndirboun(i).ne.4) cycle
!
         node=nodeboun(i)
         v(4,node)=xbounact(i)-vold(4,node)
      enddo
!
!     inserting the pressure MPC conditions
!
      do i=1,nmpc
         ist=ipompc(i)
         ndir=nodempc(2,ist)
         if(ndir.ne.4) cycle
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
!     
         index=nodempc(3,ist)
c         residu=coefmpc(ist)*(vold(ndir,node)+v(ndir,node))
         residu=0.d0
c         size=(coefmpc(ist))**2
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!
               residu=residu+coefmpc(index)*
     &              (vold(ndiri,nodei)+v(ndiri,nodei))
c               size=size+(coefmpc(index))**2
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
!
         v(ndir,node)=-residu/coefmpc(ist)-vold(ndir,node)
c!     
c!     correcting all terms of the MPC
c! 
c         residu=residu/size
c!     
c         correction=-residu*coefmpc(ist)
c         v(ndir,node)=v(ndir,node)+correction
c         index=nodempc(3,ist)
c         if(index.ne.0) then
c            do
c               nodei=nodempc(1,index)
c               ndiri=nodempc(2,index)
c!
c               correction=-residu*coefmpc(index)
c               v(ndiri,nodei)=v(ndiri,nodei)+correction
c               index=nodempc(3,index)
c               if(index.eq.0) exit
c            enddo
c         endif
      enddo
!     
      return
      end
