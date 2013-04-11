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
      subroutine gencontrel(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,cg,straight,ifree,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,cs,
     &  elcon,istep,iinc,iit,ncmat_,ntmat_,
     &  vini,nmethod,islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,slavnor,slavtan,imastop,gap,
     &  islavact)
!
!     Calculating the normals and tangential vectors in the nodes of
!       the slave surface (slavnor, slavtan)
!     Determining the locations on the master surface opposite and
!       and locally orthogonal to the integration points on the slave
!       surface (imastsurf, pmastsurf)
!
!     Author: Li, Yang; Rakotonanahary, Samoela;
!
      implicit none
!
      character*1 c
      character*3 m11,m2,m3
      integer one,number_of_nodes
!
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,ifree,imastop(3,*),kmax(3),
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(10),nodeedge(2,10),iflag,kneigh,i,j,k,l,islav,isol,
     &  itri,kflag,n,ipos,nx(*),ny(*),ipointer(10),istep,iinc,
     &  nz(*),nstart,ifaceq(8,6),ifacet(6,4),index1,
     &  ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,iit,ncmat_,ntmat_,
     &  nnodelem,nface,nope,nodef(8),m1,km1,km2,km3,
     &  nmethod,islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,
     &  mint2d,m,jj,nopes,konl(20),id,islavact(*),indexnode(8)
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:4,*),p(3),ratio(9),
     &  totdist(10),dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),cs(17,*),
     &  beta,c0,elcon(0:ncmat_,ntmat_,*),vini(0:4,*),pmastsurf(2,*),
     &  pneigh(0:3,8),et,xi,weight,xl2(0:3,8),xsj2(3),shp2(7,8),
     &  xs2(3,2),slavnor(3,*),slavtan(6,*), xquad(2,8), xtri(2,6),dd,
     &  al,al1,al2,xn(3),s(3),xnabs(3),gap(*)
!
      include "gauss.f"
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             4,6,3,1/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      data iflag /2/
!     
!     new added data for the local coodinates for nodes
!
      data xquad /-1, -1,
     &           1, -1,
     &           1, 1,
     &           -1, 1,
     &           0, -1,
     &          1, 0,
     &           0, 1,
     &           -1, 0/
!
      data xtri /0, 0,
     &          1, 0,
     &          0, 1,
     &          0.5, 0,
     &          0.5, 0.5,
     &          0, 0.5/
!
!     maximum number of neighboring master triangles for a slave node
!
      kflag=2
      ifree = 0
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         kneigh=10
!
         slavset=tieset(2,i)
         ipos=index(slavset,' ')
         if(slavset(ipos-1:ipos-1).eq.'S') cycle
!     
!     determining the slave set
!     
         do j=1,nset
            if(set(j).eq.slavset) exit
         enddo
         if(j.gt.nset) then
            write(*,*) '*ERROR in gencontrel: contact slave set',
     &           slavset
            write(*,*) '       does not exist'
            stop
         endif
         islav=j
!     
         nstart=itietri(1,i)-1
         n=itietri(2,i)-nstart
         if(n.lt.kneigh) kneigh=n
         do j=1,n
            xo(j)=cg(1,nstart+j)
            x(j)=xo(j)
            nx(j)=j
            yo(j)=cg(2,nstart+j)
            y(j)=yo(j)
            ny(j)=j
            zo(j)=cg(3,nstart+j)
            z(j)=zo(j)
            nz(j)=j
         enddo
         call dsort(x,nx,n,kflag)
         call dsort(y,ny,n,kflag)
         call dsort(z,nz,n,kflag)
!
         do l = itiefac(1,i), itiefac(2,i)
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
!
! Decide the max integration points number, just consider 2D situation 
! 
            if(lakon(nelems)(4:5).eq.'8R') then
               mint2d=1
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:4).eq.'8') then
               mint2d=4
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:6).eq.'20R') then
               mint2d=4
               nopes=8
               nope=20
            elseif(lakon(nelems)(4:4).eq.'2') then
               mint2d=9
               nopes=8
               nope=20
            elseif(lakon(nelems)(4:5).eq.'10') then
               mint2d=3
               nopes=6
               nope=10
            elseif(lakon(nelems)(4:4).eq.'4') then
               mint2d=1
               nopes=3
               nope=4
!
!         treatment of wedge faces
!
            elseif(lakon(nelems)(4:4).eq.'6') then
               mint2d=1
               nope=6
               if(jfaces.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            elseif(lakon(nelems)(4:5).eq.'15') then
               nope=15
               if(jfaces.le.2) then
                  mint2d=3
                  nopes=6
               else
                  mint2d=4
                  nopes=8
               endif
            endif
!
!           actual position of the nodes belonging to the
!           slave surface
!
            do j=1,nope
               konl(j)=kon(ipkon(nelems)+j)
            enddo
!
            if((nope.eq.20).or.(nope.eq.8)) then
                do m=1,nopes
                   do j=1,3
                      xl2(j,m)=co(j,konl(ifaceq(m,jfaces)))+
     &                     vold(j,konl(ifaceq(m,jfaces)))
                   enddo
                enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
                do m=1,nopes
                   do j=1,3
                      xl2(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &                     vold(j,konl(ifacet(m,jfaces)))
                   enddo
                enddo
            else
                do m=1,nopes
                   do j=1,3
                      xl2(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &                     vold(j,konl(ifacew1(m,jfaces)))
                   enddo
                enddo
            endif

!	calculate the normal vector in the nodes belonging to the slave surface
!
            if(nopes.eq.8) then
               do m = 1, nopes
                 xi = xquad(1,m)
                 et = xquad(2,m)
                 call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                 dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
     &              + xsj2(3)*xsj2(3))
                 xsj2(1) = xsj2(1)/dd
                 xsj2(2) = xsj2(2)/dd
                 xsj2(3) = xsj2(3)/dd
!
                 if(nope.eq.20) then
                    node = konl(ifaceq(m,jfaces))
                 elseif(nope.eq.15) then
                    node=konl(ifacew2(m,jfaces))
                 endif
!
                 call nident(islavnode(nslavnode(i)+1), node, 
     &               nslavnode(i+1)-nslavnode(i), id)
                 index1=nslavnode(i)+id
                 indexnode(m)=index1
                 slavnor(1,index1) = slavnor(1,index1)
     &            +xsj2(1)
                 slavnor(2,index1) = slavnor(2,index1)
     &            +xsj2(2)
                 slavnor(3,index1) = slavnor(3,index1)
     &            +xsj2(3)
               enddo
            elseif(nopes.eq.4) then
               do m = 1, nopes
                 xi = xquad(1,m)
                 et = xquad(2,m)
                 call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                 dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &               + xsj2(3)*xsj2(3))
                 xsj2(1) = xsj2(1)/dd
                 xsj2(2) = xsj2(2)/dd
                 xsj2(3) = xsj2(3)/dd
!
                 if(nope.eq.8) then
                    node = konl(ifaceq(m,jfaces))
                 elseif(nope.eq.6) then
                    node=konl(ifacew1(m,jfaces))
                 endif
!
                 call nident(islavnode(nslavnode(i)+1), node, 
     &               nslavnode(i+1)-nslavnode(i), id)
!
                 index1=nslavnode(i)+id
                 indexnode(m)=index1
                 slavnor(1,index1) = slavnor(1,index1)
     &            +xsj2(1)
                 slavnor(2,index1) = slavnor(2,index1)
     &            +xsj2(2)
                 slavnor(3,index1) = slavnor(3,index1)
     &            +xsj2(3)
               enddo
            elseif(nopes.eq.6) then
              do m = 1, nopes
                 xi = xquad(1,m)
                 et = xquad(2,m)
                 call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                 dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &               + xsj2(3)*xsj2(3))
                 xsj2(1) = xsj2(1)/dd
                 xsj2(2) = xsj2(2)/dd
                 xsj2(3) = xsj2(3)/dd
!
                 if(nope.eq.10) then
                    node = konl(ifacet(m,jfaces))
                 elseif(nope.eq.15) then
                    node = konl(ifacew2(m,jfaces))
                 endif
!
                 call nident(islavnode(nslavnode(i)+1), node, 
     &               nslavnode(i+1)-nslavnode(i), id)
                 index1=nslavnode(i)+id
                 indexnode(m)=index1
                 slavnor(1,index1) = slavnor(1,index1)
     &            +xsj2(1)
                 slavnor(2,index1) = slavnor(2,index1)
     &            +xsj2(2)
                 slavnor(3,index1) = slavnor(3,index1)
     &            +xsj2(3)
               enddo
            else
              do m = 1, nopes
                 xi = xquad(1,m)
                 et = xquad(2,m)
                 call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                 dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &               + xsj2(3)*xsj2(3))
                 xsj2(1) = xsj2(1)/dd
                 xsj2(2) = xsj2(2)/dd
                 xsj2(3) = xsj2(3)/dd
!
                 if(nope.eq.6) then
                    node = konl(ifacew1(m,jfaces))
                 elseif(nope.eq.4) then
                    node = konl(ifacet(m,jfaces))
                 endif
!
                 call nident(islavnode(nslavnode(i)+1), node, 
     &               nslavnode(i+1)-nslavnode(i), id)
                 index1=nslavnode(i)+id
                 indexnode(m)=index1
                 slavnor(1,nslavnode(i)+id) = slavnor(1,index1)
     &            +xsj2(1)
                 slavnor(2,nslavnode(i)+id) = slavnor(2,index1)
     &            +xsj2(2)
                 slavnor(3,nslavnode(i)+id) = slavnor(3,index1)
     &            +xsj2(3)
               enddo
            endif
! 
!           determining the location on the master surface opposite and
!           locally orthogonal to the integration points on the slave
!           surface 
!
            islavsurf(2,l)=ifree
            do m = 1,mint2d
                ifree = ifree + 1
                if((lakon(nelems)(4:5).eq.'8R').or.
     &            ((lakon(nelems)(4:4).eq.'6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,m)
                  et=gauss2d1(2,m)
                  weight=weight2d1(m)
                elseif((lakon(nelems)(4:4).eq.'8').or.
     &               (lakon(nelems)(4:6).eq.'20R').or.
     &               ((lakon(nelems)(4:5).eq.'15').and.
     &               (nopes.eq.8))) then
                  xi=gauss2d2(1,m)
                  et=gauss2d2(2,m)
                  weight=weight2d2(m)
                elseif(lakon(nelems)(4:4).eq.'2') then
                  xi=gauss2d3(1,m)
                  et=gauss2d3(2,m)
                  weight=weight2d3(m)
                elseif((lakon(nelems)(4:5).eq.'10').or.
     &               ((lakon(nelems)(4:5).eq.'15').and.
     &               (nopes.eq.6))) then
                  xi=gauss2d5(1,m)
                  et=gauss2d5(2,m)
                  weight=weight2d5(m)
                elseif((lakon(nelems)(4:4).eq.'4').or.
     &               ((lakon(nelems)(4:4).eq.'6').and.
     &               (nopes.eq.3))) then
                  xi=gauss2d4(1,m)
                  et=gauss2d4(2,m)
                  weight=weight2d4(m)
                endif
!
                if(nopes.eq.8) then
                   call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                elseif(nopes.eq.4) then
                   call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                elseif(nopes.eq.6) then
                   call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                else
                   call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                endif
!
!               coordinates of the integration point
!
                do k=1,3
                    p(k)=0.d0
                    do j=1,nopes
                       p(k)=p(k)+xl2(k,j)*shp2(4,j)
                    enddo
                enddo
!     
!              determining the kneigh neighboring master contact
!              triangle centers of gravity
!   
               call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &             n,neigh,kneigh)
!     
!              normal vector in the integration point
!
               dd=dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
!
               do k=1,3
                  xn(k)=xsj2(k)/dd
               enddo
!     
               isol=0
!
               loop1: do k=1,kneigh
                  itri=neigh(k)+itietri(1,i)-1
                  loop2: do
                     al=-(straight(16,itri)+straight(13,itri)*p(1)
     &                 +straight(14,itri)*p(2)+straight(15,itri)*p(3))/
     &                 (straight(13,itri)*xn(1)+straight(14,itri)*xn(2)
     &                 +straight(15,itri)*xn(3))
!
                     do m1=1,3
                        al1=straight(4*m1-3,itri)*p(1)+
     &                      straight(4*m1-2,itri)*p(2)+
     &                      straight(4*m1-1,itri)*p(3)
                        al2=straight(4*m1-3,itri)*xn(1)+
     &                      straight(4*m1-2,itri)*xn(2)+
     &                      straight(4*m1-1,itri)*xn(3)
                        if(al1+al*al2+straight(4*m1,itri).gt.0.d0) then
                           if(al.lt.0.d0) cycle loop1
                           itri=imastop(m1,itri)
                           if(itri.eq.0) cycle loop1
                           cycle loop2
                        endif
                     enddo
!
                     isol=1
!
                     exit loop1
                  enddo loop2
               enddo loop1
!
               if(isol.eq.0) then
                  imastsurf(ifree)=0
                  pmastsurf(1,ifree)=0.d0
                  pmastsurf(2,ifree)=0.d0
!     
!                 no independent face was found: no corresponding
!                 master location
!
               else
!
!                 independent face found; all nodes belonging to 
!                 the face are active (only at the start of a new
!                 step)
!     
                  if((iinc.eq.1).and.(iit.eq.1)) then
                     do m1=1,nopes
                        islavact(indexnode(m1))=1
                     enddo
                  endif
!     
                  nelem=int(koncont(4,itri)/10.d0)
                  jface=koncont(4,itri)-10*nelem
!
                  indexe=ipkon(nelem)
                  if(lakon(nelem)(4:4).eq.'2') then
                     nnodelem=8
                     nface=6
                  elseif(lakon(nelem)(4:4).eq.'8') then
                     nnodelem=4
                     nface=6
                  elseif(lakon(nelem)(4:5).eq.'10') then
                     nnodelem=6
                     nface=4
                  elseif(lakon(nelem)(4:4).eq.'4') then
                     nnodelem=3
                     nface=4
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     if(jface.le.2) then
                        nnodelem=6
                     else
                        nnodelem=8
                     endif
                     nface=5
                     nope=15
                  elseif(lakon(nelem)(4:4).eq.'6') then
                     if(jface.le.2) then
                        nnodelem=3
                     else
                        nnodelem=4
                     endif
                     nface=5
                     nope=6
                  else
                     cycle
                  endif
!
!                 determining the nodes of the face
!
                  if(nface.eq.4) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacet(k,jface))
                     enddo
                  elseif(nface.eq.5) then
                     if(nope.eq.6) then
                        do k=1,nnodelem
                           nodef(k)=kon(indexe+ifacew1(k,jface))
                        enddo
                     elseif(nope.eq.15) then
                        do k=1,nnodelem
                           nodef(k)=kon(indexe+ifacew2(k,jface))
                        enddo
                     endif
                  elseif(nface.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifaceq(k,jface))
                     enddo
                  endif
!
                  do k=1,nnodelem
                     do j = 1,3
                        pneigh(j,k) = co(j,nodef(k))+vold(j,nodef(k))
                     enddo
                  enddo
                  call attach(pneigh,p,nnodelem,ratio,dist,xi,et)
                  imastsurf(ifree) = koncont(4,itri)
                  pmastsurf(1,ifree) = xi
                  pmastsurf(2,ifree) = et
               endif
            enddo
         enddo
!
!        normalizing the normals
!
         do l=nslavnode(i)+1,nslavnode(i+1)
            dd=dsqrt(slavnor(1,l)**2+slavnor(2,l)**2+
     &               slavnor(3,l)**2)
            do m=1,3
               slavnor(m,l)=slavnor(m,l)/dd
            enddo
!
!        determining the tangential directions
!
            do m=1,3
               xn(m)=slavnor(m,l)
               xnabs(m)=dabs(xn(m))
            enddo
            n=3
            kmax(1)=1
            kmax(2)=2
            kmax(3)=3
            kflag=2
!     
!     sorting the components of the normal
!     
            call dsort(xnabs,kmax,n,kflag)
!     
            km1=kmax(3)
            km2=km1+1
            if(km2.gt.3) km2=1
            km3=km2+1
            if(km3.gt.3) km3=1
!     
            slavtan(km1,l)=-slavnor(km2,l)
            slavtan(km2,l)=slavnor(km1,l)
            slavtan(km3,l)=0.d0
            dd=dsqrt(slavtan(km1,l)**2+slavtan(km2,l)**2)
            slavtan(km1,l)=slavtan(km1,l)/dd
            slavtan(km2,l)=slavtan(km2,l)/dd
!     
            slavtan(4,l)=xn(2)*slavtan(3,l)
     &           -xn(3)*slavtan(2,l)
            slavtan(5,l)=xn(3)*slavtan(1,l)
     &           -xn(1)*slavtan(3,l)
            slavtan(6,l)=xn(1)*slavtan(2,l)
     &           -xn(2)*slavtan(1,l)
!
!           determining the gap function (at the nodes)
!
            if(islavact(l).ne.0) then
               do m=1,3
!
!                 nodal coordinates
!
                  p(m)=co(m,islavnode(l))
                  xn(m)=slavnor(m,l)
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
                  call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &                 n,neigh,kneigh)
!     
                  isol=0
!     
                  loop3: do k=1,kneigh
                     itri=neigh(k)+itietri(1,i)-1
                     loop4: do
                      al=-(straight(16,itri)+straight(13,itri)*p(1)
     &                 +straight(14,itri)*p(2)+straight(15,itri)*p(3))/
     &                 (straight(13,itri)*xn(1)+straight(14,itri)*xn(2)
     &                 +straight(15,itri)*xn(3))
!     
!     intersection of normal to the slave surface with
!     the triangle plane of itri
!     
                        do m1=1,3
                           s(m1)=p(m1)+al*xn(m1)
                        enddo
!     
                        do m1=1,3
                           al1=straight(4*m1-3,itri)*p(1)+
     &                          straight(4*m1-2,itri)*p(2)+
     &                          straight(4*m1-1,itri)*p(3)
                           al2=straight(4*m1-3,itri)*xn(1)+
     &                          straight(4*m1-2,itri)*xn(2)+
     &                          straight(4*m1-1,itri)*xn(3)
                           if(al1+al*al2+straight(4*m1,itri).gt.0.d0) 
     &                        then
                              if(al.lt.0.d0) cycle loop3
                              itri=imastop(m1,itri)
                              if(itri.eq.0) cycle loop3
                              cycle loop4
                           endif
                        enddo
!     
                        isol=1
						gap(l)=al
!     
                        exit loop3
                     enddo loop4
                  enddo loop3
               enddo
            endif
         enddo
      enddo
!
            open(70,file='contact.frd',status='unknown')
            c='C'
            m11=' -1'
            m2=' -2'
            m3=' -3'
            p0='    0'
            p1='    1'
            p2='    2'
            p3='    3'
            p7='    7'
            p9999=' 9999'
            one=1
            write(70,'(a5,a1)') p1,c
            write(70,'(a5,a1,67x,i1)') p2,c,one
            number_of_nodes=0
            do i=1,itietri(2,ntie)
               number_of_nodes=max(number_of_nodes,koncont(1,i))
               number_of_nodes=max(number_of_nodes,koncont(2,i))
               number_of_nodes=max(number_of_nodes,koncont(3,i))
            enddo
            do i=1,number_of_nodes
               write(70,'(a3,i10,1p,3e12.5)') m11,i,(co(j,i),j=1,3)
            enddo
            write(70,'(a3)') m3
            write(70,'(a5,a1,67x,i1)') p3,c,one
            do i=1,itietri(2,ntie)
               write(70,'(a3,i10,2a5)')m11,i,p7,p0
               write(70,'(a3,3i10)') m2,(koncont(j,i),j=1,3)
            enddo
            write(70,'(a3)') m3
            write(70,'(a5)') p9999
            close(70)
!     
      return
      end
