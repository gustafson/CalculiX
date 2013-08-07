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
c>     Calculating the normals and tangential vectors in the nodes of
c>       the slave surface (slavnor, slavtan)
c>
c>
c> @param   [in]     tieset      name and dependent surface of tie set
c> @param   [in]     ntie        number of contraints
c> @param   [in]     itietri     (1,i)pointer to node where trangulation starts for i (2,i) pointer to end
c> @param   [in]     ipkon       pointer into field kon
c> @param   [in]     set         (i)name of set_i
c> @param   [in]     cg          field containing centers of gravity
c> @param   [in]     straight    (1:4 5:8 9:13,i)coeffs of plane equation for edges of triagle_i (13:16,i) coeffs of plane containing triagle
c> @param   [in]     koncont     (1:3,i) nodes of triagle_i (4,i) element face
c> @param   [in]     co          field containing the coordinates of all nodes
c> @param   [in]     vold        field containing the displacements
c> @param   [in,out] x,y,z       ONLY HELP FIELD
c> @param   [in,out] xo,yo,zo    ONLY HELP FIELD
c> @param   [in,out] nx,ny,nz    ONLY HELP FIELD
c> @param   [in]     kon         Field containing the connectivity of the elements in succesive order
c> @param   [in]     lakon       element label 
c> @param   [in]     iinc        index increment
c> @param   [in]     iit         index iteration
c> @param   [in]     itiefac     pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
c> @param   [in]     islavsurf   islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
c> @param   [in]     nslavnode   (i) for contraint i pointer into field islavnode
c> @param   [in]     imastop     (l,i) for edge l in triagle i neightbouring triangle
c> @param   [in]     imastsurf   pointer into pmastsurf    NOT USED
c> @param   [in]     pmastsurf   field storing position and etal for integration points on master side NOT USED
C> @param   [in]     islavnode   fields containing nodes of slace surfaces
C> @param   [out]    slavnor     normal vektors in the nods of slave surface
C> @param   [out]    slavtan     tangetial vektors in the nodes of slave surface 
C> @param   [in]       mi        NOT USED
C> @param   [in]     ncont
C> @param   [in]       ipe       NOT USED
C> @param   [in]       ime       NOT USED
C> @param   [in]       pslavsurf NOT USED
C> @param   [out]    pslavdual	 (1:4,i)dual shape functions for face i 

      subroutine gencontrel(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,cg,straight,
     &  koncont,co,vold,nset,
     &  iinc,iit,
     &  islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,slavnor,slavtan,imastop,
     &  mi,ncont,ipe,ime,pslavsurf,pslavdual)
!
!     Calculating the normals and tangential vectors in the nodes of
!       the slave surface (slavnor, slavtan)
!     Determining the coefficients of the dual shape functions on
!       the slave surface
!
!     Author: Li, Yang; Rakotonanahary, Samoela;
!
      implicit none
!
      logical checkbiorthogonality
!
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,ifree,imastop(3,*),kmax(3),ncont,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  iflag,kneigh,i,j,l,islav,
     &  itri,kflag,ipos,iinc,iit,ijk,
     &  index1,
     &  nface,nope,m1,km1,km2,km3,number,
     &  islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,mi(*),
     &  mint2d,m,nopes,konl(20),id,indexnode(8),
     &  line,
     &  ipe(*),ime(4,*),nintpoint,
     &  ipiv(4),ifac,getiface,nodesf,ifs,
     &   flagtan
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),
     &  pmastsurf(2,*),et,xi,xl2s(3,8),xsj2(3),
     &  shp2(7,8),t1(6),t2(6),xlnode(3),
     &  xs2(3,2),slavnor(3,*),slavtan(6,*), xquad(2,8), xtri(2,6),dd,
     &  al2,xn(3),xnabs(3),e(3,3),
     &  pslavdual(16,*),pslavsurf(3,*),err,xs2m(3,2),xsj2m(3)
!     
      include "gauss.f"
!
      data iflag /2/
      ijk=0
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
      data e /1.d0 , 0.d0 , 0.d0,
     &        0.d0 , 1.d0 , 0.d0,
     &        0.d0 , 0.d0 , 1.d0/
!
      checkbiorthogonality=.false.
      flagtan=7
!
      open(40,file='contact.fbd',status='unknown')
      open(50,file='slavtan.fbd',status='unknown')
      open(20,file='slavintmortar.out',status='unknown')
      open(30,file='intpoints.out',status='unknown')
!
!     maximum number of neighboring master triangles for a slave node
!
      kflag=2
      ifree = 0     
      err=1.d-6
      do i=1,ntie  
         do l=nslavnode(i)+1,nslavnode(i+1)
            do m=1,3
               slavnor(m,l)=0.0
            enddo
         enddo
      enddo
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         kneigh=1     
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

         do l = itiefac(1,i), itiefac(2,i)
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
            call getnumberofnodes(nelems,jfaces,lakon,nope,
     &           nopes,mint2d) 
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            do j=1,nope
               konl(j)=kon(ipkon(nelems)+j)
            enddo
!     
            do m=1,nopes
               ifac=getiface(m,jfaces,nope)
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifac))+
     &                 vold(j,konl(ifac))   
               enddo
            enddo          
!     calculate the normal vector in the nodes belonging to the slave surface
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
     &                 + xsj2(3)*xsj2(3))
                  xsj2(1) = xsj2(1)/dd
                  xsj2(2) = xsj2(2)/dd
                  xsj2(3) = xsj2(3)/dd                 
                  ifac=getiface(m,jfaces,nope)
                  node = konl(ifac)
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
!
 105     format(4(1x,e15.8))
         enddo
!     
!     FIRST SLAVE SURFACE LOOP DONE
!     
!     normalizing the normals
!     
         do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
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
               xlnode(m)=co(m,node)+vold(m,node)
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
! tan5    
!
      if(flagtan==5.or.flagtan==4.or.flagtan==1)then
            km1=kmax(3)
            km2=km1+1
            if(km2.gt.3) km2=1
            km3=km2+1
            if(km3.gt.3) km3=1     
            t1(km1)=-slavnor(km3,l)
            t1(km3)=slavnor(km1,l)
            t1(km2)=0.d0
            dd=dsqrt(t1(km1)**2+t1(km3)**2)
            t1(km1)=t1(km1)/dd
            t1(km3)=t1(km3)/dd
!     
            t1(4)=xn(2)*t1(3)
     &           -xn(3)*t1(2)
            t1(5)=xn(3)*t1(1)
     &           -xn(1)*t1(3)
            t1(6)=xn(1)*t1(2)
     &           -xn(2)*t1(1)   
            t2(km1)=-slavnor(km2,l)
            t2(km2)=slavnor(km1,l)
            t2(km3)=0.d0
            dd=dsqrt(t2(km1)**2+t2(km2)**2)
            t2(km1)=t2(km1)/dd
            t2(km2)=t2(km2)/dd
     
            t2(4)=xn(2)*t2(3)
     &           -xn(3)*t2(2)
            t2(5)=xn(3)*t2(1)
     &           -xn(1)*t2(3)
            t2(6)=xn(1)*t2(2)
     &           -xn(2)*t2(1)     
            do m=1,3
               if(flagtan==5)then
                  slavtan(m,l)=(t1(m)+t2(m))
               elseif(flagtan==1) then
                  slavtan(m,l)=(t1(m))
               else
                  slavtan(m,l)=(t2(m))
               endif
            enddo
            do m=4,6
               if(flagtan==5)then
                  slavtan(m,l)=(t1(m)+t2(m))
               elseif(flagtan==1) then
                  slavtan(m,l)=(t1(m))
               else
                  slavtan(m,l)=(t2(m))
               endif
            enddo
            dd=dsqrt(slavtan(1,l)**2+ slavtan(2,l)**2
     &        +slavtan(3,l)**2)
            do m=1,3
               slavtan(m,l)=slavtan(m,l)/dd
            enddo 
            dd=dsqrt(slavtan(4,l)**2+ slavtan(5,l)**2
     &         +slavtan(6,l)**2)
            do m=4,6
               slavtan(m,l)=slavtan(m,l)/dd
            enddo
          elseif(flagtan==6)then
! A.Popp
             if(abs(xn(2))>1.d-6)then
                slavtan(1,l)=1.d0
                slavtan(3,l)=1.d0
                slavtan(2,l)=(-xn(1)-xn(3))/xn(2)
             elseif(abs(xn(3))>1.d-6)then
                slavtan(1,l)=1.d0
                slavtan(2,l)=1.d0
                slavtan(3,l)=(-xn(1)-xn(2))/xn(3)
             elseif(abs(xn(1))>1.d-6)then
                slavtan(2,l)=1.d0
                slavtan(3,l)=1.d0
                slavtan(1,l)=(-xn(2)-xn(3))/xn(1)       
             else
                write(*,*)'gencontrel: something wrong with slavnor'
             endif
             dd=dsqrt(slavtan(1,l)**2+ slavtan(2,l)**2
     &        +slavtan(3,l)**2)
             do m=1,3
                slavtan(m,l)=slavtan(m,l)/dd
             enddo 
             slavtan(4,l)=xn(2)*slavtan(3,l)
     &            -xn(3)*slavtan(2,l)
             slavtan(5,l)=xn(3)*slavtan(1,l)
     &            -xn(1)*slavtan(3,l)
             slavtan(6,l)=xn(1)*slavtan(2,l)
     &            -xn(2)*slavtan(1,l)
             dd=dsqrt(slavtan(4,l)**2+ slavtan(5,l)**2
     &            +slavtan(6,l)**2)
             do m=4,6
                slavtan(m,l)=slavtan(m,l)/dd
             enddo
          elseif(flagtan==7)then
               if(1.d0 - dabs(xn(1)).lt.1.5231d-6) then       
!           
!     calculating the local directions on master surface
!
                  slavtan(1,l)=-xn(3)*xn(1)
                  slavtan(2,l)=-xn(3)*xn(2)
                  slavtan(3,l)=1.d0-xn(3)*xn(3)
               else
                  slavtan(1,l)=1.d0-xn(1)*xn(1)
                  slavtan(2,l)=-xn(1)*xn(2)
                  slavtan(3,l)=-xn(1)*xn(3)
               endif
             dd=dsqrt(slavtan(1,l)**2+ slavtan(2,l)**2
     &        +slavtan(3,l)**2)
             do m=1,3
                slavtan(m,l)=slavtan(m,l)/dd
             enddo 
               slavtan(4,l)=-(xn(2)*slavtan(3,l)-xn(3)*slavtan(2,l))
               slavtan(5,l)=-(xn(3)*slavtan(1,l)-xn(1)*slavtan(3,l))
               slavtan(6,l)=-(xn(1)*slavtan(2,l)-xn(2)*slavtan(1,l)) 
          endif
!
         ijk=ijk+1
         write(50,100) ijk,(xlnode(m),m=1,3)
         ijk=ijk+1
         write(50,100) ,ijk,(xlnode(m)+0.25*slavnor(m,l),m=1,3)
         ijk=ijk+1
         write(50,100) ,ijk,(xlnode(m)+0.25*slavtan(m,l),m=1,3)
         ijk=ijk+1
         write(50,100) ,ijk,(xlnode(m)+0.25*slavtan(3+m,l),m=1,3)
         write(50,101) ijk-3,ijk-3,ijk-2
         write(50,101) ijk-2,ijk-3,ijk-1
         write(50,101) ijk-1,ijk-3,ijk
 100     format('PNT ',i10,'P',3(1x,e15.8))
 101     format('LINE ',i10,'L',i10,'P ',i10,'P')
         enddo
!     
      enddo
!
      close(50)
      return
      end
