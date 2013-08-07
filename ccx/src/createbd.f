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
c> subroutine to calculate all combinations of p,q for current slave face l and store ist into contri.
c> goal is to calculate full  Bd[p,q]
c> \f$ B[p,q]=-<\Phi_q,\Psi_p> Id_d\f$.
c> @param [in] l                 current slave face
c> @param [in] ipkon             pointer for current element into field kon
c> @param [in] kon               Field containing the connectivity of the elements in succesive order
c> @param [in] lakon             element label 
c> @param [in] nodem             index of node on master side \f$ p \in S\f$
c> @param [in] nodesf            index of node on slave side \f$ q\in M \f$
c> @param [in] islavsurf         islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
c> @param [in] imastsurf         index of masterface corresponding to integration point i
c> @param [in] pmastsurf         field storing position and etal for integration points on master side
c> @param [in] itiefac           pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
c> @param [out] contr            field containing contributions for curretn face
c> @param [out] iscontr          (i) nodesf of contribution(i)
c> @param [out] imcontr          (i) nodem of contribution(i)
c> @param [in] iponoels          pointer for node i into field inoel...
c> @param [in] inoels            ...which stores 1D&2D elements belonging to node
c> @param [in] mi
c> @param [in] pslavsurf         field storing position xil, etal and weight for integration point on slave side
c> @param [in] pslavdual         (1:4,*)dual shape functions in nodesf 
      subroutine createbd(ict,l,ipkon,kon,lakon,co, vold, gapmints,
     &  islavsurf,imastsurf,pmastsurf,itiefac,contr,iscontr,imcontr,
     &  dcontr,idcontr1,idcontr2,gcontr,igcontr,iponoels,inoels,mi,
     &  pslavsurf,pslavdual,nslavnode,islavnode,nmastnode,imastnode,
     &  icounter,icounter2)
!         
!
!
!     Author: Saskia Sitzmann
!
      implicit none
!
      logical debug
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),konl(20),iflag,m,l,j,jj,
     &  indexe,islavsurf(2,*),iponoels(*),inoels(2,*),
     &  imastsurf(*),itiefac(2,*),ifaces,nelemens,jfaces,ifacem,
     &  mint2d,indexf,nopes1,nopes2,nodem,nodesf,nopes,
     &  locs,locm,mi(*),ns,mint2dloc1,mint2dloc2,
     &  ifs,ifm,nope1,nope2, iscontr(*),imcontr(*),getiface,
     &  jfacem,nelemenm,icounter,idummy,ifac,idcontr1(*),idcontr2(*),
     &  igcontr(*),icounter2,nslavnode(*),islavnode(*),ict,id,
     &  nmastnode(*),imastnode(*)
!
      real*8 pmastsurf(2,*),co(3,*),gapmints(*),
     &  vold(0:mi(2),*),weight,dx,help,
     &  ets,xis,etm,xim,xl2s(3,8),xsj2s(3),xsj2s2(3),
     &  shp2s(4,8),xs2s(3,2),xl2m(3,8),xsj2m(3),shp2m(7,8),xs2m(3,2),
     &  contribution,pslavsurf(3,*),pslavdual(16,*),contr(*), 
     &  dcontr(*),dcontribution,gcontr(*),gcontribution,
     &  shp2s2(7,8),xs2s2(3,2)
!
      include "gauss.f"
!            
      debug=.false.
      contribution = 0.d0
      dcontribution = 0.d0
      gcontribution = 0.d0
      icounter=0
      icounter2=0
      ifaces=islavsurf(1,l)
      nelemens = int(ifaces/10)
      jfaces = ifaces - nelemens*10
      indexe = ipkon(nelemens)
      ict=ict+1
      if(debug)write(*,*) 'createbd:l=',l,'tie',ict
      call getnumberofnodes(nelemens,jfaces,lakon,nope1,nopes1,idummy)
      mint2d=islavsurf(2,l+1)-islavsurf(2,l)
      if(mint2d==0) return
      indexf=islavsurf(2,l)
      if(debug)write(*,*) 'createbd:mint2d',mint2d
!     loop over all nodesf of current slave face      
      do j=1,nope1
            konl(j)=kon(ipkon(nelemens)+j)
      enddo
      do m=1,nopes1
         ifac=getiface(m,jfaces,nope1)
         do j=1,3
            xl2s(j,m)=co(j,konl(ifac))+
     &           vold(j,konl(ifac))       
         enddo
      enddo
!
      do j=1,nopes1
         do jj=1,nopes1
            dcontr(icounter2+nopes1*(j-1)+jj)=0.0
         enddo        
         gcontr(icounter2+j)=0.0
      enddo
!
      mint2dloc1=1
      mint2dloc2=1
      help=0.d0
!     lopp over all integration points created for current slave face  
      do 
         if (mint2dloc2>=mint2d) exit
!        find current master face
         ifacem=imastsurf(indexf+mint2dloc1)
         nelemenm= int (ifacem/10);
         jfacem=ifacem-nelemenm*10
         call getnumberofnodes(nelemenm,jfacem,lakon,
     &         nope2,nopes2,idummy)
!         find number of integration points belonging to master face     
         do
            if(ifacem==imastsurf(indexf+mint2dloc2+1))then
               mint2dloc2=mint2dloc2+1
               if(mint2dloc2==mint2d) exit     
            else
               exit
            endif  
         enddo
      if(debug)write(*,*) 'createbd:MF,loc1, loc2',ifacem,
     &    mint2dloc1,mint2dloc2       
         help=0.0
         do m=mint2dloc1,mint2dloc2
            xis=pslavsurf(1,indexf+m)
            ets=pslavsurf(2,indexf+m)
            weight=pslavsurf(3,indexf+m)
            ns=l
            iflag = 2
         if(debug)write(*,100) xis,ets,weight
 100     format('createbd:xs',3(1x,e15.8))    
            if(nopes1.eq.8) then
               call dualshape8q(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
            elseif(nopes1.eq.4) then
               call dualshape4q(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &              pslavdual,iflag)
            elseif(nopes1.eq.6) then
               call dualshape6tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,iflag)
            else
               call dualshape3tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &                    pslavdual,iflag)
            endif
            if(nopes1.eq.8) then
               call shape8q(xis,ets,xl2s,xsj2s2,xs2s2,shp2s2,iflag)
            elseif(nopes1.eq.4) then
               call shape4q(xis,ets,xl2s,xsj2s2,xs2s2,shp2s2,iflag)
            elseif(nopes1.eq.6) then
               call shape6tri(xis,ets,xl2s,xsj2s2,xs2s2,shp2s2,iflag)
            else
               call shape3tri(xis,ets,xl2s,xsj2s2,xs2s2,shp2s2,iflag)
            endif    
            xim = pmastsurf(1,indexf+m)
            etm = pmastsurf(2,indexf+m)
 101     format('createbd:xm',2(1x,e15.8))  
!     
            if(nopes2.eq.8) then
               call shape8q(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
            elseif(nopes2.eq.4) then
               call shape4q(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
            elseif(nopes2.eq.6) then
               call shape6tri(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
            else
               call shape3tri(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
            endif
            if(m==mint2dloc1)then
               do j=1,nopes1
                  do jj=1,nopes2
                     contr(icounter+nopes2*(j-1)+jj)=0.0
                  enddo
               enddo
            endif
            dx=dsqrt(xsj2s(1)**2+xsj2s(2)**2+xsj2s(3)**2)   
            do j=1,nopes1
               ifs=getiface(j,jfaces,nope1) 
               nodesf=kon(ipkon(nelemens)+ifs)
               locs=j
               gcontribution=shp2s(4,locs)
     &              *(gapmints(indexf+m))
     &              *pslavsurf(3,indexf+m) 
     &              *dx
               if(debug.and.j.eq.4)write(*,*)'ns',nodesf,'v',
     &          shp2s(4,locs),gapmints(indexf+m)
               gcontr(icounter2+j)=gcontr(icounter2+j)+gcontribution
            if(debug.and.j.eq.4)write(*,*) 'gap',gcontr(icounter2+j)
               if(m==1)then
                    call nident(islavnode(nslavnode(ict)+1),
     &                 nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
                    if(islavnode(nslavnode(ict)+id)==nodesf) then
                     igcontr(icounter2+j)=nslavnode(ict)+id
                    else
                       write(*,*)'createbd: node',nodesf
                       write(*,*)'was not catalogued properly in', 
     &                 'islavnode'
                       stop
                    endif
               endif
               do jj=1,nopes1
                  ifs=getiface(jj,jfaces,nope1) 
                  nodem=kon(ipkon(nelemens)+ifs)
                  dcontribution=shp2s(4,locs)
     &              *shp2s2(4,jj)  
     &              *pslavsurf(3,indexf+m)
     &              *dx
                  dcontr(icounter2+nopes1*(j-1)+jj)=
     &              dcontr(icounter2+nopes1*(j-1)+jj)+dcontribution                  
                  if(m==1)then
                    call nident(islavnode(nslavnode(ict)+1),
     &                 nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
                    if(islavnode(nslavnode(ict)+id)==nodesf) then
                     idcontr1(icounter2+nopes1*(j-1)+jj)=
     &                  nslavnode(ict)+id
                    else
                       write(*,*)'createbd: node',nodesf
                       write(*,*)'was not catalogued properly in', 
     &                 'islavnode'
                       stop
                    endif
!                    
                    call nident(islavnode(nslavnode(ict)+1),
     &                 nodem,(nslavnode(ict+1)-nslavnode(ict)),id)
                    if(islavnode(nslavnode(ict)+id)==nodem) then
                     idcontr2(icounter2+nopes1*(j-1)+jj)=
     &                 nslavnode(ict)+id
                    else
                       write(*,*)'createbd: node',nodem
                       write(*,*)'was not catalogued properly in', 
     &                 'islavnode'
                       stop
                    endif
!                    
                  endif
               enddo
               do jj=1,nopes2
                  ifm=getiface(jj,jfacem,nope2)
                  nodem=kon(ipkon(nelemenm)+ifm)
                  locm=jj
                  contribution=shp2s(4,locs)*shp2m(4,locm)
     &                 *pslavsurf(3,indexf+m)
     &                 *dx
                  contr(icounter+nopes2*(j-1)+jj)=
     &                 contr(icounter+nopes2*(j-1)+jj)+contribution
                  if(m==mint2dloc1)then
                    call nident(islavnode(nslavnode(ict)+1),
     &                 nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
                    if(islavnode(nslavnode(ict)+id)==nodesf) then
                     iscontr(icounter+nopes2*(j-1)+jj)=
     &                 nslavnode(ict)+id
                    else
                       write(*,*)'createbd: node',nodesf
                       write(*,*)'was not catalogued properly in', 
     &                 ' islavnode'
                       stop
                    endif                            
                    call nident(imastnode(nmastnode(ict)+1),
     &                 nodem,(nmastnode(ict+1)-nmastnode(ict)),id)
                    if(imastnode(nmastnode(ict)+id)==nodem) then
                     imcontr(icounter+nopes2*(j-1)+jj)=
     &                 nmastnode(ict)+id
                    else
                       write(*,*)'createbd: node',nodem
                       write(*,*)'was not catalogued properly in', 
     &                ' imastnode',nmastnode(ict)+id,
     &                  imastnode(nmastnode(ict)+id), nmastnode(ict)+1,
     &                  nmastnode(ict+1)
                       stop
                    endif
!                     
                  endif
                  contribution=0.d0
               enddo
            enddo
!            
         enddo
         mint2dloc1=mint2dloc2+1
         mint2dloc2=mint2dloc1
         icounter=icounter+nopes1*nopes2
      enddo
      icounter2=icounter2+nopes1*nopes1
!      
      if(debug)then
         write(*,*) 'createbd: contri,iscontr,imcontr',l
         do j=1, icounter
            write(*,*)contr(j),iscontr(j),imcontr(j)
         enddo 
         write(*,*) 'createbd: dcontri,idcontr',l
         do j=1, icounter2
            write(*,*) dcontr(j),idcontr1(j),idcontr2(j)
         enddo 
         write(*,*) 'createbd: gcontri,idcontr',l
         do j=1, nopes1
            write(*,*) gcontr(j),igcontr(j)
         enddo 
      endif
      ict=ict-1
      return 
      end
      
