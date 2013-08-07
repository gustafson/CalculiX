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
!
c>     Determining the location of the integration points in slave
c>     surface ifaces. This location depends on the triangulation of 
c>     the opposite master surface. For the slave surface the local
c>     coordinates and the integration weight is stored in pslavsurf,
c>     the label of the opposite master face and the local coordinates
c>     of the opposite point on the master faces are stored in 
c>     islavsurf and pmastsurf, respectively
c> @param   [in]     tieset      name and dependent surface of tie set
c> @param   [in]     ntie        number of contraints
c> @param   [in]     itietri     (1,i)pointer to node where trangulation starts for i (2,i) pointer to end
c> @param   [in]     ipkon       pointer into field kon
c> @param   [in]     kon         Field containing the connectivity of the elements in succesive order
c> @param   [in]     set         (i)name of set_i
c> @param   [in]     cg          field containing centers of gravity
c> @param   [in]     straight    (1:4 5:8 9:13,i)coeffs of plane equation for edges of triagle_i (13:16,i) coeffs of plane containing triagle
c> @param   [in]     koncont     (1:3,i) nodes of triagle_i (4,i) element face
c> @param   [in]     co          field containing the coordinates of all nodes
c> @param   [in]     vold        field containing the displacements
c> @param   [in,out] x,y,z       ONLY HELP FIELD
c> @param   [in,out] xo,yo,zo    ONLY HELP FIELD
c> @param   [in,out] nx,ny,nz    ONLY HELP FIELD
c> @param   [out]    nintpoint
c> @param   [in]     nset
c> @param   [in]     lakon       element label 
c> @param   [in]     iinc        index increment
c> @param   [in]     iit         index iteration
c> @param   [in]     itiefac     pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
c> @param   [in,out]     islavsurf   islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i
c> @param   [in]     nslavnode   (i) for contraint i pointer into field islavnode
c> @param   [in]     imastop     (l,i) for edge l in triagle i neightbouring triangle
c> @param   [out]     imastsurf   pointer into pmastsurf    
c> @param   [out]     pmastsurf   field storing position and etal for integration points on master side 
C> @param   [in]     islavnode   fields containing nodes of slace surfaces
C> @param   [in]     slavnor     normal vektors in the nods of slave surface
C> @param   [in]     slavtan      tangetial vektors in the nodes of slave surface NOT USED!
C> @param   [in]     mi
C> @param   [in]     ncont
C> @param   [in]     ipe
C> @param   [in]     ime
C> @param   [in]     pslavsurf
C> @param   [in]     pslavdual   (1:4,i)dual shape functions for face i 
c> @param   [in]     islavact    active set
c> @param   [out]    pslavsurf   field storing position xil, etal and weight for integration point on slave side
c> @param   [in]     i           current tie
c> @param   [in]     l           current face
c> @param   [in]     ntri        # triangles
c> @param   [out]    gapmints    stores gaps between master and slave side
c> @todo l.353 check if this can be improved. is it possible to reduce complexity working with islavact???    
      subroutine slavintmortar(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,cg,straight,nintpoint,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,
     &  iinc,iit,
     &  islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,slavnor,slavtan,imastop,gapmints,
     &  islavact,mi,ncont,ipe,ime,pslavsurf,pslavdual,i,l,ntri)

!     Author: Li, Yang; Rakotonanahary, Samoela; Sitzmann,Saskia
!
      implicit none
!     
      logical debug,nogap
!     
      character*8 lakon(*)
      character*81 tieset(3,*),set(*)
!     
      integer ntie,nset,nintpoint,imastop(3,*),ncont,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(10),iflag,kneigh,i,j,k,l,
     &  itri,ipos,nx(*),ny(*),
     &  nz(*),index1,ifreeintersec,
     &  nelemm,jfacem,indexe,iit,iinc,
     &  nnodelem,nope,m1,
     &  islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,mi(*),
     &  m,nopes,konl(20),id,islavact(*),indexnode(8),
     &  itria(4,2),ntria,itriacorner(4,4),line,
     &  iactiveline(3,3*ncont),
     &  nactiveline,ipe(*),ime(4,*),k1,j1,
     &  info,ntri,nintpfirst,nodem(8),
     &  compt,il,ifac,getiface,ifacem,idummy,nopemm
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),
     &  xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &  pmastsurf(2,*),xl2m(3,8),xl2s(3,8),
     &  pmiddle(3),xl2sr(3,8),xl2sp(3,8),
     &  slavnor(3,*),slavtan(6,*),dd,xns(3,8),areaslav,
     &  al,xn(3),gapmints(*),slavstraight(20),
     &  pslavdual(16,*),err2,dist,distmin,
     &  pslavsurf(3,*),err,pnodesin(3,3*ncont),xquad(2,8), 
     &  xtri(2,6),xi,et,xsj2(3),xs2(3,2),shp2(7,8),anglesm
!     
      include "gauss.f"
!
      debug=.false.
      data iflag /2/
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
      kneigh=1
      err=0.1
      areaslav=0.0
      err2=1.d-4
      nintpfirst=nintpoint
      compt=0
      islavsurf(2,l)=nintpoint
      if(debug)WRITE(30,*) '#SLAVINTMORTAR iinc',iinc, 'face',l
c     WRITE(20,*) '#SLAVINTMORTAR iinc',iinc, 'face',l      
c     WRITE(*,*) '#SLAVINTMORTAR iit',iit, 'face',l
!     
!     Research of the contact integration points
!     
      ifaces = islavsurf(1,l)
      nelems = int(ifaces/10)
      jfaces = ifaces - nelems*10
!     
!     get nope,nopes
!     
      call getnumberofnodes(nelems,jfaces,lakon,nope,nopes,idummy)
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
      do j=1,nope
         konl(j)=kon(ipkon(nelems)+j)
      enddo
!     
      do m=1,nopes
         do j=1,3
            ifac=getiface(m,jfaces,nope)
            xl2s(j,m)=co(j,konl(ifac))+
     &           vold(j,konl(ifac))     
         enddo
      enddo  
!     
!     slightly reducing the size of the slave surface in
!     an aleatoric way
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
            xl2sr(j,m)=xl2s(j,m)-0.5*err*(xl2s(j,m)-pmiddle(j))
c     xl2sr(j,m)=xl2s(j,m)
         enddo
      enddo
!     
!     calculate the mean normal vector on the Slave Surface
!     
      do k=1,3
         xn(k)=0.d0
      enddo
!     
      do m = 1, nopes
         if(nopes.eq.4 .or. nopes.eq.8)then
            xi = xquad(1,m)
            et = xquad(2,m)
         else
            xi = xtri(1,m)
            et = xtri(2,m)
         endif
         if(nopes.eq.8)then
            call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
         elseif(nopes.eq.4)then
            call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
         elseif(nopes.eq.6)then
            call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
         else
            call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
         endif   
         dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
     &        + xsj2(3)*xsj2(3))
         xsj2(1) = xsj2(1)/dd
         xsj2(2) = xsj2(2)/dd
         xsj2(3) = xsj2(3)/dd
!     
         do k=1,3
            xn(k) = xn(k)
     &           +xsj2(k)
         enddo
      enddo            
!     
!     normalizing the mean normal on the Slave surface
!     
      dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
      do k=1,3
         xn(k)=xn(k)/dd
      enddo
c     write(*,*)'slavm xn',(xn(k),k=1,3)
!     
!     determine the equations of the triangle/quadrilateral
!     (mean)plane and of the planes boardering the 
!     triangle/quadrilateral
!     
      if(nopes.eq.3) then
         call straighteq3d(xl2s,slavstraight)
         do k=1,3
            xn(k)= slavstraight(4*nopes+k)
         enddo               
      else
         call approxplane(xl2s,slavstraight,xn)
      endif
!     
!     Project slave nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1, nopes
         al=-xn(1)*xl2s(1,j)-xn(2)*
     &        xl2s(2,j)-xn(3)*xl2s(3,j)-
     &        slavstraight(nopes*4+4)
         if(nopes.eq.4)then
            do k=1,3
               xl2sp(k,j)= xl2s(k,j)+al*xn(k)
            enddo
         else
            do k=1,3
               xl2sp(k,j)= xl2s(k,j)
            enddo
         endif
      enddo 
!     
!     determine the triangles corresponding to the corner
!     nodes
!     
      ntria=0
      do j=1,4
         itria(j,1)=0
         itria(j,2)=0
         do k=1,4
            itriacorner(j,k)=0
         enddo
      enddo
!     check for nogap-nodes
      do j=1,nopes
         call neartriangle(xl2sr(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &        ntri,neigh,kneigh,itietri,ntie,straight,imastop,itri,i,
     &        debug)
         ifac= getiface(j,jfaces,nope)
         node= konl(ifac)
         if(debug) then
            write(20,*) neigh(1),neigh(1)+itietri(1,i)-1
            write(20,*) 'itri',itri,'node',node
         endif
         anglesm=xn(1)*straight(13,itri)
     &        +xn(2)*straight(14,itri)
     &        +xn(3)*straight(15,itri)
         if(anglesm.lt.-0.2)then
            call nident(islavnode(nslavnode(i)+1), node, 
     &           nslavnode(i+1)-nslavnode(i), id)
c     write(20,*) ' node',node,islavnode(nslavnode(i)+id),
c     &            id,islavact(nslavnode(i)+id) 
            if(itri.ne.0.and.islavact(nslavnode(i)+id).eq.-3) then  
               islavact(nslavnode(i)+id)=0
            endif
            if(itri.eq.0.and.islavact(nslavnode(i)+id).gt.-1) then  
               islavact(nslavnode(i)+id)=-3
            endif
         else
            call nident(islavnode(nslavnode(i)+1), node, 
     &           nslavnode(i+1)-nslavnode(i), id)
c     write(20,*) ' node',node,islavnode(nslavnode(i)+id),
c     &            id,islavact(nslavnode(i)+id) 
            if(itri.ne.0.and.islavact(nslavnode(i)+id).gt.-1) then  
               islavact(nslavnode(i)+id)=-3
            endif
         endif
      enddo
c     
      do j=1,3
         do m=1,nopes
            xl2sr(j,m)=xl2s(j,m)-2*err*(xl2s(j,m)-pmiddle(j))
         enddo
      enddo
      distmin=1.1   
      do j=1,nopes
         call neartriangle(xl2sr(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &        ntri,neigh,kneigh,itietri,ntie,straight,imastop,itri,i,
     &        debug)
         ifac= getiface(j,jfaces,nope)
         node= konl(ifac) 
         if(itri.eq.0) then  
            cycle
         endif
         dist= -(straight(13,itri)*xl2sr(1,j)+
     &        straight(14,itri)*xl2sr(2,j)+
     &        straight(15,itri)*xl2sr(3,j)+
     &        straight(16,itri))/
     &        (straight(13,itri)*xn(1)+
     &        straight(14,itri)*xn(2)+
     &        straight(15,itri)*xn(3))
         if(dist.lt.distmin)distmin=dist
         if(debug)write(20,*) 'j',j,'dist',dist,distmin
         ifacem=koncont(4,itri)
         if(debug)write(20,*)'noder ',node, 'itri',itri,
     &        'ifacem',ifacem
!     
!     
         call nident(itria(1:4,1),itri,ntria,id)
         if(id.gt.0) then
            if(itria(id,1).eq.itri) then
               itriacorner(j,id)=1
               cycle
            endif
         endif
         call nident(itria(1:4,2),ifacem,ntria,id)
         if(id.gt.0) then
            if(itria(id,2).eq.ifacem) then
               itriacorner(j,id)=1
               cycle
            endif
         endif
!     
!     triangle was not covered yet: add to stack
!     
!     angle criteria 
!     
         anglesm=xn(1)*straight(13,itri)
     &        +xn(2)*straight(14,itri)
     &        +xn(3)*straight(15,itri)
         if(debug)write(20,*)'cos alpa',anglesm
!     
c     if(dist.lt.1.0)then
c     if(anglesm.lt.0.0)then
         if(anglesm.lt.-0.1)then
            ntria=ntria+1
            do k=ntria,id+2,-1
               itria(k,1)=itria(k-1,1)
               itria(k,2)=itria(k-1,2)
               do m=1,j-1
                  itriacorner(m,k)=itriacorner(m,k-1)
               enddo
            enddo
            itria(id+1,1)=itri
            itria(id+1,2)=ifacem
            itriacorner(j,id+1)=1
            do m=1,j-1
               itriacorner(m,id+1)=0
            enddo 
         endif              
      enddo
      if(debug)then 
         write(20,*)'itria ifacem n1 n2 n3 n4 face', l
         do k=1,ntria
            write(20,*) itria (k,1:2),itriacorner(1:nopes,k)
         enddo
      endif
!     
c     if(distmin.gt.1.0) then
c     write(20,*) 'face',l,'distmin',distmin
c     write(20,*) 'no integrationpoints generated,',
c     &      'too much dist!'
c     islavsurf(2,l+1)=nintpoint
c     return
c     endif           
      nactiveline=0
      ifreeintersec=0
!     
!     treating the corner elements first
!     
      do j=1,ntria
         if(debug)write(20,*) 'corner triangle j',j
         itri=itria(j,1)
         ifacem=koncont(4,itri)
         nelemm=int(ifacem/10.d0)
         jfacem=ifacem-10*nelemm
         if(debug)write(20,*)itri,itria(j,2), ifacem,nelemm,jfacem
         call getnumberofnodes(nelemm,jfacem,lakon,nopemm,
     &        nnodelem,idummy)     
!     
!     determining the nodes of the face
!     
         do j1=1,nopemm
            konl(j1)=kon(ipkon(nelemm)+j1)
         enddo
         do k1=1,nnodelem
            ifac=getiface(k1,jfacem,nopemm)
            nodem(k1)=konl(ifac)
            do j1=1,3
               xl2m(j1,k1)=co(j1,konl(ifac))+
     &              vold(j1,konl(ifac))
            enddo
         enddo 
         dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         if(debug)then
            write(20,*) 'dd',dd    
            write(20,100)(xn(k),k=1,3)
            write(20,100)(slavstraight(nopes*4+k),k=1,3)
            write(20,*) 'SIM: xl2s'
            do j1=1,nopes
               write(20,*)(xl2s(k,j1),k=1,3)
            enddo    
            write(20,*) 'SIM: xl2m'
            do j1=1,nnodelem
               write(20,*)(xl2m(k,j1),k=1,3)
            enddo
         endif                                        
 100     format('SIM: xns',3(3x,e15.8))
 101     format(3(e15.8))
!     
         if(debug) write(20,*) 'TT: itri',nelemm
!     
!     Project master nodes to meanplane, needed for Sutherland-Hodgman
!     
         call treattriangleS(
     &        nopes,slavstraight,xn,xns,co,xl2s,xl2sp,
     &        ipe,ime,iactiveline,nactiveline,
     &        ifreeintersec,ifacem,itriacorner(1,j),
     &        nintpoint,pslavsurf,ncont,imastsurf,pmastsurf,
     &        xl2m,nnodelem,nodem,mi,pnodesin,straight,gapmints,l,
     &        areaslav,debug)
      enddo
!     
!     retrieving all triangles by neighborhood search
!     
      do
         line=iactiveline(1,1)
         if(nactiveline.eq.0) exit
         if(koncont(4,ime(2,line)).eq.iactiveline(2,1)) then
            itri=imastop(ime(3,line),ime(2,line))
         else
            itri=ime(2,line)
         endif
!     check whether still in contact tie
         if(itri.gt.itietri(2,i) .or. itri.lt.itietri(1,i))then
            if(itri.ne.0)then
c     write(*,*)'sim:tie',i,'face',l,'itiri',itri
c     write(*,*)' mintri',itietri(1,i),'maxtri',itietri(2,i)
            endif
            itri=0
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
         ifacem=koncont(4,itri)
         nelemm=int(koncont(4,itri)/10.d0)
         jfacem=koncont(4,itri)-10*nelemm
!     
         indexe=ipkon(nelemm)
         call getnumberofnodes(nelemm,jfacem,lakon,nopemm,
     &        nnodelem,idummy)
         
!     
!     determining the nodes of the face
!     
         do j1=1,nopemm
            konl(j1)=kon(ipkon(nelemm)+j1)
         enddo
         do k1=1,nnodelem
            ifac=getiface(k1,jfacem,nopemm)
            nodem(k1)=konl(ifac)
            do j1=1,3
               xl2m(j1,k1)=co(j1,konl(ifac))+
     &              vold(j1,konl(ifac))
            enddo
         enddo
         compt=compt+1
         if(debug)then
            write(20,*) 'dd',dd ,'ifacem', ifacem   
            write(20,100)(xn(k),k=1,3)
            write(20,*) 'SIM: xl2s'
            do j1=1,nopes
               write(20,*)(xl2s(k,j1),k=1,3)
            enddo    
            write(20,*) 'SIM: xl2m'
            do j1=1,nnodelem
               write(20,*)(xl2m(k,j1),k=1,3)
            enddo
         endif  
         call treattriangleS(
     &        nopes,slavstraight,xn,xns,co,xl2s,xl2sp,
     &        ipe,ime,iactiveline,nactiveline,
     &        ifreeintersec,ifacem,itriacorner,nintpoint,
     &        pslavsurf,ncont,imastsurf,pmastsurf,
     &        xl2m,nnodelem,nodem,mi,pnodesin,straight,gapmints,l,
     &        areaslav,debug)
      enddo
      islavsurf(2,l+1)=nintpoint
      if(debug) then
         if(areaslav.lt.1.e-12)write(*,*)'areaslav(',l,')=',
     &        areaslav
      endif
c     write(20,*)'nintp',nintpoint-nintpfirst
!     check for N-nodes
c     if(nintpoint-nintpfirst.gt.0)then
      do j=1,nope
         konl(j)=kon(ipkon(nelems)+j)
      enddo
c     nogap=.false.
c     do j=1,nopes
c     ifac= getiface(j,jfaces,nope)
c     node= konl(ifac)
c     call nident(islavnode(nslavnode(i)+1), node, 
c     &              nslavnode(i+1)-nslavnode(i), id)
c     if(islavact(nslavnode(i)+id).gt.-1) then  
c     nogap=.true.
c     endif
c     enddo
      if(nintpoint-nintpfirst.gt.0)then
         do j=1,nopes
            ifac= getiface(j,jfaces,nope)
            node= konl(ifac)
            call nident(islavnode(nslavnode(i)+1), node, 
     &           nslavnode(i+1)-nslavnode(i), id)
c     write(20,*) ' node',node,islavnode(nslavnode(i)+id),
c     &            id,islavact(nslavnode(i)+id)
            if(islavact(nslavnode(i)+id).eq.-3) then  
               islavact(nslavnode(i)+id)=-1
            endif
         enddo
      endif
      if(debug)then
         write(20,*) 'mint2d', (nintpoint-islavsurf(2,l))
         write(20,*)'areaslav(',l,')=',areaslav
      endif
!     
      return
      end
      
