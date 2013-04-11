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
     &  lakon,set,cg,straight,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,
     &  iinc,iit,
     &  islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,slavnor,slavtan,imastop,
     &  islavact,mi,ncont,ipe,ime,pslavsurf,pslavdual)
!
!     Calculating the normals and tangential vectors in the nodes of
!       the slave surface (slavnor, slavtan)
!     Determining the gap at the regular integration points of the
!       slave surface
!     Determining the coefficients of the dual shape functions on
!       the slave surface
!
!     Author: Li, Yang; Rakotonanahary, Samoela;
!
      implicit none
!
      character*1 c
      character*3 m11,m2,m3
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,ifree,imastop(3,*),kmax(3),ncont,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(10),iflag,kneigh,i,j,k,l,islav,isol,
     &  itri,kflag,ntri,ipos,nx(*),ny(*),iinc,
     &  nz(*),nstart,ifaceq(8,6),ifacet(6,4),index1,ifreeintersec,
     &  ifacew1(4,5),ifacew2(8,5),nelemm,jfacem,indexe,iit,
     &  nnodelem,nface,nope,nodef(8),m1,km1,km2,km3,number,
     &  islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,mi(2),
     &  mint2d,m,nopes,konl(20),id,islavact(*),indexnode(8),
     &  itria(4),ntria,itriacorner(4,4),inodesin(3*ncont),line,
     &  nnodesin,inodesout(3*ncont),nnodesout,iactiveline(3,3*ncont),
     &  nactiveline,intersec(2,6*ncont),ipe(*),ime(4,*),nintpoint,k1,j1,
     &  ipiv(4),info,ipnt,one,number_of_nodes,itel
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  xntersec(3,6*ncont),xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &  pmastsurf(2,*),xl2m(3,8),et,xi,weight,xl2s(3,8),xsj2(3),
     &  shp2(7,8),
     &  xs2(3,2),slavnor(3,*),slavtan(6,*), xquad(2,8), xtri(2,6),dd,
     &  al,al1,al2,xn(3),xnabs(3),slavstraight(20),
     &  pslavdual(16,*),diag_els(4),m_els(10),contribution,work(4)
!
      real*4 rand
      real*8 pslavsurf(3,*),err,pnodesin(3,3*ncont)
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
      open(40,file='contact.fbd',status='unknown')
!
!     maximum number of neighboring master triangles for a slave node
!
      kflag=2
      ifree = 0
!     
      err=1.d-6
!
!     storing the triangulation of the master surfaces
!     
      open(70,file='TriMasterContactPair.frd',status='unknown')
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
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         kneigh=1
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
!        ntri: number of triangles in the triangulation of the master
!        surface corresponding to tie i
!
         nstart=itietri(1,i)-1
         ntri=itietri(2,i)-nstart
         if(ntri.lt.kneigh) kneigh=ntri
         do j=1,ntri
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
         call dsort(x,nx,ntri,kflag)
         call dsort(y,ny,ntri,kflag)
         call dsort(z,nz,ntri,kflag)
c      do j=1,ntri
c         write(*,*) j,x(j),y(j),z(j),nx(j),ny(j),nz(j)
c      enddo
!     
         do l = itiefac(1,i), itiefac(2,i)
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
!     
!     initialization for Dualshape Coefficient matrix
!     
            ipnt=0
            do k=1,4
               diag_els(k)=0.0
c               do j=k,4
               do j=1,k
                  ipnt=ipnt+1
                  m_els(ipnt)=0.0
               enddo
            enddo
!     
!     Decide on the max integration point number, just consider 2D situation 
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
!     treatment of wedge faces
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
!     actual position of the nodes belonging to the
!     slave surface
!     
            do j=1,nope
               konl(j)=kon(ipkon(nelems)+j)
            enddo
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopes
                  do j=1,3
                     xl2s(j,m)=co(j,konl(ifaceq(m,jfaces)))+
     &                    vold(j,konl(ifaceq(m,jfaces)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopes
                  do j=1,3
                     xl2s(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &                    vold(j,konl(ifacet(m,jfaces)))
                  enddo
               enddo
            else
               do m=1,nopes
                  do j=1,3
                     xl2s(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &                    vold(j,konl(ifacew1(m,jfaces)))
                  enddo
               enddo
            endif
            
!     calculate the normal vector in the nodes belonging to the slave surface
!     
            if(nopes.eq.8) then
               do m = 1, nopes
                  xi = xquad(1,m)
                  et = xquad(2,m)
                  call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
                  dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
     &                 + xsj2(3)*xsj2(3))
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
     &                 nslavnode(i+1)-nslavnode(i), id)
                  index1=nslavnode(i)+id
                  indexnode(m)=index1
                  slavnor(1,index1) = slavnor(1,index1)
     &                 +xsj2(1)
                  slavnor(2,index1) = slavnor(2,index1)
     &                 +xsj2(2)
                  slavnor(3,index1) = slavnor(3,index1)
     &                 +xsj2(3)
               enddo
            elseif(nopes.eq.4) then
               do m = 1, nopes
                  xi = xquad(1,m)
                  et = xquad(2,m)
                  call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
                  dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                 + xsj2(3)*xsj2(3))
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
     &                 nslavnode(i+1)-nslavnode(i), id)
!     
                  index1=nslavnode(i)+id
                  indexnode(m)=index1
                  slavnor(1,index1) = slavnor(1,index1)
     &                 +xsj2(1)
                  slavnor(2,index1) = slavnor(2,index1)
     &                 +xsj2(2)
                  slavnor(3,index1) = slavnor(3,index1)
     &                 +xsj2(3)
               enddo
            elseif(nopes.eq.6) then
               do m = 1, nopes
                  xi = xquad(1,m)
                  et = xquad(2,m)
                  call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
                  dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                 + xsj2(3)*xsj2(3))
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
     &                 nslavnode(i+1)-nslavnode(i), id)
                  index1=nslavnode(i)+id
                  indexnode(m)=index1
                  slavnor(1,index1) = slavnor(1,index1)
     &                 +xsj2(1)
                  slavnor(2,index1) = slavnor(2,index1)
     &                 +xsj2(2)
                  slavnor(3,index1) = slavnor(3,index1)
     &                 +xsj2(3)
               enddo
            else
               do m = 1, nopes
                  xi = xquad(1,m)
                  et = xquad(2,m)
                  call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
                  dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                 + xsj2(3)*xsj2(3))
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
     &                 nslavnode(i+1)-nslavnode(i), id)
                  index1=nslavnode(i)+id
                  indexnode(m)=index1
                  slavnor(1,nslavnode(i)+id) = slavnor(1,index1)
     &                 +xsj2(1)
                  slavnor(2,nslavnode(i)+id) = slavnor(2,index1)
     &                 +xsj2(2)
                  slavnor(3,nslavnode(i)+id) = slavnor(3,index1)
     &                 +xsj2(3)
               enddo
            endif
!     
!     determining the gap contribution of the integration points
!     and the coefficient for the slave dualshape functions
!     
            do m = 1,mint2d
               ifree = ifree + 1
               if((lakon(nelems)(4:5).eq.'8R').or.
     &              ((lakon(nelems)(4:4).eq.'6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,m)
                  et=gauss2d1(2,m)
                  weight=weight2d1(m)
               elseif((lakon(nelems)(4:4).eq.'8').or.
     &                 (lakon(nelems)(4:6).eq.'20R').or.
     &                 ((lakon(nelems)(4:5).eq.'15').and.
     &                 (nopes.eq.8))) then
                  xi=gauss2d2(1,m)
                  et=gauss2d2(2,m)
                  weight=weight2d2(m)
               elseif(lakon(nelems)(4:4).eq.'2') then
                  xi=gauss2d3(1,m)
                  et=gauss2d3(2,m)
                  weight=weight2d3(m)
               elseif((lakon(nelems)(4:5).eq.'10').or.
     &                 ((lakon(nelems)(4:5).eq.'15').and.
     &                 (nopes.eq.6))) then
                  xi=gauss2d5(1,m)
                  et=gauss2d5(2,m)
                  weight=weight2d5(m)
               elseif((lakon(nelems)(4:4).eq.'4').or.
     &                 ((lakon(nelems)(4:4).eq.'6').and.
     &                 (nopes.eq.3))) then
                  xi=gauss2d4(1,m)
                  et=gauss2d4(2,m)
                  weight=weight2d4(m)
               endif
!     
               if(nopes.eq.8) then
                  call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
                  contribution=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+
     &                 xsj2(3)**2)
                  ipnt=0
                  do k=1,4
                     diag_els(k)=diag_els(k)+shp2(4,k)*contribution
                     do j=1,k
                        ipnt=ipnt+1
                        m_els(ipnt)=m_els(ipnt)+shp2(4,k)*shp2(4,j)*
     &                       contribution
                     enddo
                  enddo
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               else
                  call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               endif
!     
!     Calculate the Mass matrix for compilation of the dualshapefunction 
!     pslavdual(16,*)
!     
               do k=1,3
                  p(k)=0.d0
                  do j=1,nopes
                     p(k)=p(k)+xl2s(k,j)*shp2(4,j)
                  enddo
               enddo
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
               call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &              ntri,neigh,kneigh)
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
                  itel=0
                  loop2: do
                     itel=itel+1
                     al=-(straight(16,itri)+straight(13,itri)*p(1)
     &                 +straight(14,itri)*p(2)+straight(15,itri)*p(3))/
     &                 (straight(13,itri)*xn(1)+straight(14,itri)*xn(2)
     &                 +straight(15,itri)*xn(3))
!     
                     do m1=1,3
                        al1=straight(4*m1-3,itri)*p(1)+
     &                       straight(4*m1-2,itri)*p(2)+
     &                       straight(4*m1-1,itri)*p(3)
                        al2=straight(4*m1-3,itri)*xn(1)+
     &                       straight(4*m1-2,itri)*xn(2)+
     &                       straight(4*m1-1,itri)*xn(3)
                        if(al1+al*al2+straight(4*m1,itri).gt.1.d-10)then
c                           if(al.lt.1.d-10) cycle loop1
                           itri=imastop(m1,itri)
                           if(itel.gt.ntri) then
                              write(*,*) '*INFO in gencontrel: circular'
                              write(*,*) 
     &                           '      hopping; no contact established'
                              exit loop1
                           endif
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
               if(isol.ne.0) then
!     
!     independent face found; all nodes belonging to 
!     the face are active (only at the start of a new
!     step)
!     
c               WRITE(*,*) "gencontrel",al
                  if((iinc.eq.1).and.(iit.eq.1)) then
                     do m1=1,nopes
                        islavact(indexnode(m1))=1
                     enddo
                  endif
               endif
            enddo
!
!     computation of psladual            
!     
!     compute inverse of me_ls
!     factorisation
!     
            call dsptrf('U',4,m_els,ipiv,info)
!     inverse
            call dsptri('U',4,m_els,ipiv,work,info)
!     
!     stack of pslavdual multiplication with diag_els
!     
            pslavdual(1,l)=diag_els(1)*m_els(1)
            pslavdual(2,l)=diag_els(1)*m_els(2)
            pslavdual(3,l)=diag_els(1)*m_els(4)
            pslavdual(4,l)=diag_els(1)*m_els(7)
            pslavdual(5,l)=diag_els(2)*m_els(2)
            pslavdual(6,l)=diag_els(2)*m_els(3)
            pslavdual(7,l)=diag_els(2)*m_els(5)
            pslavdual(8,l)=diag_els(2)*m_els(8)
            pslavdual(9,l)=diag_els(3)*m_els(4)
            pslavdual(10,l)=diag_els(3)*m_els(5)
            pslavdual(11,l)=diag_els(3)*m_els(6)
            pslavdual(12,l)=diag_els(3)*m_els(9)
            pslavdual(13,l)=diag_els(4)*m_els(7)
            pslavdual(14,l)=diag_els(4)*m_els(8)
            pslavdual(15,l)=diag_els(4)*m_els(9)
            pslavdual(16,l)=diag_els(4)*m_els(10)
         enddo
!     
!     FIRST SLAVE SURFACE LOOP DONE
!     
!     normalizing the normals
!     
         do l=nslavnode(i)+1,nslavnode(i+1)
            dd=dsqrt(slavnor(1,l)**2+slavnor(2,l)**2+
     &           slavnor(3,l)**2)
            do m=1,3
               slavnor(m,l)=slavnor(m,l)/dd
            enddo
!     
!     determining the tangential directions
!     
            do m=1,3
               xn(m)=slavnor(m,l)
               xnabs(m)=dabs(xn(m))
            enddo
            number=3
            kmax(1)=1
            kmax(2)=2
            kmax(3)=3
            kflag=2
!     
!     sorting the components of the normal
!     
            call dsort(xnabs,kmax,number,kflag)
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
         enddo
!     
      enddo 
!     
      return
      end
