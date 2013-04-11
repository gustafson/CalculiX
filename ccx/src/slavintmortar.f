!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine slavintmortar(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,cg,straight,nintpoint,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,
     &  iinc,iit,
     &  islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,slavnor,slavtan,imastop,gapmints,
     &islavact,mi,ncont,ipe,ime,pslavsurf,pslavdual,i,l,ntri)
!
!     Determining the location of the integration points in slave
!     surface ifaces. This location depends on the triangulation of 
!     the opposite master surface. For the slave surface the local
!     coordinates and the integration weight is stored in pslavsurf,
!     the label of the opposite master face and the local coordinates
!     of the opposite point on the master faces are stored in 
!     islavsurf and pmastsurf, respectively
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
      integer ntie,nset,nintpoint,imastop(3,*),kmax(3),ncont,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(10),iflag,kneigh,i,j,k,l,islav,isol,
     &  itri,kflag,n,ipos,nx(*),ny(*),iinc,
     &  nz(*),nstart,ifaceq(8,6),ifacet(6,4),index1,ifreeintersec,
     &  ifacew1(4,5),ifacew2(8,5),nelemm,jfacem,indexe,iit,
     &  nnodelem,nface,nope,nodef(8),m1,km1,km2,km3,number,
     &  islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,mi(*),
     &  mint2d,m,nopes,konl(20),id,islavact(*),indexnode(8),
     &  itria(4),ntria,itriacorner(4,4),inodesin(3*ncont),line,
     &  nnodesin,inodesout(3*ncont),nnodesout,iactiveline(3,3*ncont),
     &  nactiveline,intersec(2,6*ncont),ipe(*),ime(4,*),k1,j1,
     &  ipiv(4),info,ipnt,ntri,nintpfirst,
     &  compt,il
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  xntersec(3,6*ncont),xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &  pmastsurf(2,*),xl2m(3,8),et,xi,weight,xl2s(3,8),xsj2(3),
     &  shp2(7,8),pmiddle(3),
     &  xs2(3,2),slavnor(3,*),slavtan(6,*),dd,
     &  al,al1,al2,xn(3),xnabs(3),gapmints(*),slavstraight(20),
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
      kneigh=1
      err=1.d-6
c      err=1.d-2
      nintpfirst=nintpoint
      compt=0
c      WRITE(*,*) "SLAVINTMORTAR"
!     
!     Research of the contact integration points
!     
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
!     
!     Decide the max integration points number, just consider 2D situation 
!     
            if(lakon(nelems)(4:5).eq.'8R') then
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:4).eq.'8') then
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:6).eq.'20R') then
               nopes=8
               nope=20
            elseif(lakon(nelems)(4:4).eq.'2') then
c               nopes=8
c               nope=20
               nopes=4
               nope=20
            elseif(lakon(nelems)(4:5).eq.'10') then
               nopes=6
               nope=10
            elseif(lakon(nelems)(4:4).eq.'4') then
               nopes=3
               nope=4
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelems)(4:4).eq.'6') then
               nope=6
               if(jfaces.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            elseif(lakon(nelems)(4:5).eq.'15') then
               nope=15
               if(jfaces.le.2) then
                  nopes=6
               else
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
     &                     vold(j,konl(ifaceq(m,jfaces)))
c     &                    vold(j,konl(ifaceq(m,jfaces)))+err*rand(0)
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopes
                  do j=1,3
                     xl2s(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &                     vold(j,konl(ifacet(m,jfaces)))
c     &                    vold(j,konl(ifacet(m,jfaces)))+err*rand(0)
                  enddo
               enddo
            else
               do m=1,nopes
                  do j=1,3
                     xl2s(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &                     vold(j,konl(ifacew1(m,jfaces)))
c     &                    vold(j,konl(ifacew1(m,jfaces)))+err*rand(0)
                  enddo
               enddo
            endif
!
!           slightly reducing the size of the slave surface in
!           an aleatoric way
!
            do j=1,3
               pmiddle(j)=0.d0
               do m=1,nopes
                  pmiddle(j)=pmiddle(j)+xl2s(j,m)
               enddo
               pmiddle(j)=pmiddle(j)/nopes
            enddo
            do j=1,3
               do m=1,nopes
                  xl2s(j,m)=xl2s(j,m)-err*rand(0)*(xl2s(j,m)-pmiddle(j))
               enddo
            enddo
!
!     calculate the mean normal vector on the Slave Surface
!     
            do k=1,3
               xn(k)=0.d0
            enddo
            if(nopes.eq.8) then
               do m = 1, nopes
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
                  do k=1,3
                     xn(k)=slavnor(k,index1)+xn(k)
                  enddo
               enddo
            elseif(nopes.eq.4) then
               do m = 1, nopes
!     
                  if((nope.eq.8).or.(nope.eq.20)) then
                     node = konl(ifaceq(m,jfaces))
                  elseif(nope.eq.6) then
                     node=konl(ifacew1(m,jfaces))
                  endif
!     
                  call nident(islavnode(nslavnode(i)+1), node, 
     &                 nslavnode(i+1)-nslavnode(i), id)
!     
                  index1=nslavnode(i)+id
                  do k=1,3
                     xn(k)=slavnor(k,index1)+xn(k)
                  enddo
               enddo
            elseif(nopes.eq.6) then
               do m = 1, nopes
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
                  do k=1,3
                     xn(k)=slavnor(k,index1)+xn(k)
                  enddo
               enddo
            else
               do m = 1, nopes
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
                  do k=1,3
                     xn(k)=slavnor(k,index1)+xn(k)
                  enddo
               enddo
            endif
!     
!     normalizing the mean normal on the Slave surface
!     
            dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
            do k=1,3
               xn(k)=xn(k)/dd
            enddo
!     
!     determine the equations of the triangle/quadrilateral
!     (mean)plane and of the planes boardering the 
!     triangle/quadrilateral
!     
            if(nopes.eq.3) then
               call straighteq3d(xl2s,slavstraight)
            else
               call approxplane(xl2s,slavstraight,xn)
            endif
!     
!     determine the triangles corresponding to the corner
!     nodes
!     
            ntria=0
            do j=1,4
               itria(j)=0
               do k=1,4
                  itriacorner(j,k)=0
               enddo
            enddo
!     
            do j=1,nopes
               call neartriangle(xl2s(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &           ntri,neigh,kneigh,itietri,ntie,straight,imastop,itri,i)
               if(itri.eq.0) cycle
!      
!
               node = konl(ifaceq(j,jfaces))
               call nident(islavnode(nslavnode(i)+1), node, 
     &           nslavnode(i+1)-nslavnode(i), id)
               if (islavact(nslavnode(i)+id).eq.-1) then
                      islavact(nslavnode(i)+id)=0
               endif
!
!
               call nident(itria,itri,ntria,id)
               if(id.gt.0) then
                  if(itria(id).eq.itri) then
                     itriacorner(j,id)=1
                     cycle
                  endif
               endif
!     
!     triangle was not covered yet: add to stack
!     
               ntria=ntria+1
               do k=ntria,id+2,-1
                  itria(k)=itria(k-1)
                  do m=1,j-1
                     itriacorner(m,k)=itriacorner(m,k-1)
                  enddo
               enddo
               itria(id+1)=itri
               itriacorner(j,id+1)=1
               do m=1,j-1
                  itriacorner(m,id+1)=0
               enddo               
            enddo
!  
            nnodesin=0
            nnodesout=0
            nactiveline=0
            ifreeintersec=0
!     
!     treating the corner triangles first
!     
            do j=1,ntria
               itri=itria(j)
               nelemm=int(koncont(4,itri)/10.d0)
               jfacem=koncont(4,itri)-10*nelemm
!     
               indexe=ipkon(nelemm)
               if(lakon(nelemm)(4:4).eq.'2') then
                  nnodelem=8
                  nface=6
               elseif(lakon(nelemm)(4:4).eq.'8') then
                  nnodelem=4
                  nface=6
               elseif(lakon(nelemm)(4:5).eq.'10') then
                  nnodelem=6
                  nface=4
               elseif(lakon(nelemm)(4:4).eq.'4') then
                  nnodelem=3
                  nface=4
               elseif(lakon(nelemm)(4:5).eq.'15') then
                  if(jfacem.le.2) then
                     nnodelem=6
                  else
                     nnodelem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelemm)(4:4).eq.'6') then
                  if(jfacem.le.2) then
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
!     determining the nodes of the face
!     
               if(nface.eq.4) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacet(k,jfacem))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew1(k,jfacem))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew2(k,jfacem))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifaceq(k,jfacem))
                  enddo
               endif
!     
               do k1=1,nnodelem
                  do j1 = 1,3
                     xl2m(j1,k1) = co(j1,nodef(k1))+vold(j1,nodef(k1))
                  enddo
               enddo
!     
               call treattriangle(inodesin,nnodesin,inodesout,
     &              nnodesout,nopes,slavstraight,xn,co,xl2s,ipe,ime,
     &              iactiveline,nactiveline,intersec,xntersec,
     &              ifreeintersec,itri,koncont,itriacorner(1,j),
     &              nintpoint,pslavsurf,ncont,imastsurf,pmastsurf,
     &              xl2m,nnodelem,vold,mi,pnodesin,straight,gapmints,l)
            enddo
!     
!     retrieving all triangles by neighborhood search
!     
            do
               line=iactiveline(1,1)
               if(nactiveline.eq.0) exit
               if(ime(2,line).eq.iactiveline(2,1)) then
                  itri=imastop(ime(3,line),ime(2,line))
               else
                  itri=ime(2,line)
               endif
!     
!     corners of the Slave surface have already been treated
!     
               if(itri.eq.0) then
                nactiveline=nactiveline-1
                  do il=1,nactiveline
                   do k=1,3
                      iactiveline(k,il)=iactiveline(k,il+1)
                   enddo
                  enddo
                cycle
               endif
               do j=1,4
                  itriacorner(j,1)=0
               enddo
!     
               nelemm=int(koncont(4,itri)/10.d0)
               jfacem=koncont(4,itri)-10*nelemm
!     
               indexe=ipkon(nelemm)
               if(lakon(nelemm)(4:4).eq.'2') then
                  nnodelem=8
                  nface=6
               elseif(lakon(nelemm)(4:4).eq.'8') then
                  nnodelem=4
                  nface=6
               elseif(lakon(nelemm)(4:5).eq.'10') then
                  nnodelem=6
                  nface=4
               elseif(lakon(nelemm)(4:4).eq.'4') then
                  nnodelem=3
                  nface=4
               elseif(lakon(nelemm)(4:5).eq.'15') then
                  if(jfacem.le.2) then
                     nnodelem=6
                  else
                     nnodelem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelemm)(4:4).eq.'6') then
                  if(jfacem.le.2) then
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
!     determining the nodes of the face
!     
               if(nface.eq.4) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacet(k,jfacem))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew1(k,jfacem))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew2(k,jfacem))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifaceq(k,jfacem))
                  enddo
               endif
!     
               do k1=1,nnodelem
                  do j1 = 1,3
                     xl2m(j1,k1) = co(j1,nodef(k1))+vold(j1,nodef(k1))
                  enddo
               enddo
!     
               compt=compt+1
               call treattriangle(inodesin,nnodesin,inodesout,
     &              nnodesout,nopes,slavstraight,xn,co,xl2s,ipe,ime,
     &              iactiveline,nactiveline,intersec,xntersec,
     &              ifreeintersec,itri,koncont,itriacorner,nintpoint,
     &              pslavsurf,ncont,imastsurf,pmastsurf,
     &              xl2m,nnodelem,vold,mi,pnodesin,straight,gapmints,l)
            enddo
!     
         islavsurf(2,l+1)=nintpoint
!     
!
!       
      return
      end
