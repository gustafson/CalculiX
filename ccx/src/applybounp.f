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
     &  voldaux,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
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
     &  refkin,reftuf,refvel,cp,voldaux(0:4,*),physcon(*),v(0:mi(2),*),
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
         residu=coefmpc(ist)*(voldaux(ndir,node)+v(ndir,node))
         size=(coefmpc(ist))**2
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!
c               idof=8*(nodei-1)+ndiri
c               call nident(ikboun,idof,nboun,id)
c               if(id.ne.0) then
c                  if(ikboun(id).eq.idof) then
c                     index=nodempc(3,index)
c                     if(index.eq.0) exit
c                     cycle
c                  endif
c               endif
!
               residu=residu+coefmpc(index)*
     &              (voldaux(ndiri,nodei)+v(ndiri,nodei))
               size=size+(coefmpc(index))**2
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
         v(ndir,node)=v(ndir,node)+correction
         index=nodempc(3,ist)
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!
c               idof=8*(nodei-1)+ndiri
c               call nident(ikboun,idof,nboun,id)
c               if(id.ne.0) then
c                  if(ikboun(id).eq.idof) then
c                     index=nodempc(3,index)
c                     if(index.eq.0) exit
c                     cycle
c                  endif
c               endif
!
               correction=-residu*coefmpc(index)
               v(ndiri,nodei)=v(ndiri,nodei)+correction
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
      enddo
!     
      return
      end
