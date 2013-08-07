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
c>     Determining the coefficients of the dual shape functions on
c>       the slave surface
c>
c>
c> @param   [in]     tieset      name and dependent surface of tie set
c> @param   [in]     ntie        number of contraints
c> @param   [in]     itietri     (1,i)pointer to node where trangulation starts for i (2,i) pointer to end
c> @param   [in]     ipkon       pointer into field kon
c> @param   [in]     se          (i)name of set_i
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
C> @param   [in]       mi        NOT USED
C> @param   [in]     ncont
C> @param   [in]       ipe       NOT USED
C> @param   [in]       ime       NOT USED
C> @param   [in]       pslavsurf NOT USED
C> @param   [out]    pslavdual	 (1:4,i)dual shape functions for face i 

      subroutine gendualcoeffs(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,cg,straight,
     &  koncont,co,vold,nset,
     &  iinc,iit,islavact,
     &  islavsurf,imastsurf,pmastsurf,itiefac,
     &  islavnode,nslavnode,imastop,
     &  mi,ncont,ipe,ime,pslavsurf,pslavdual)
!
!     Determining the coefficients of the dual shape functions on
!       the slave surface
!
!     Author: Sitzmann,Saskia ;
!
      implicit none
!      
      logical checkbiorthogonality, checknorm
!
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,ifree,imastop(3,*),kmax(3),ncont,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(10),iflag,kneigh,i,ii,j,k,l,islav,isol,
     &  itri,kflag,ntri,ipos,iinc,islavact(*),
     &  nstart,index1,ifreeintersec,
     &  nelemm,jfacem,indexe,iit,
     &  nnodelem,nface,nope,nodef(8),m1,km1,km2,km3,number,
     &  islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &  imastsurf(*),itiefac(2,*),ifaces,nelems,jfaces,mi(*),
     &  mint2d,m,nopes,konl(20),id,indexnode(8),
     &  itria(4),ntria,itriacorner(4,4),inodesin(3*ncont),line,
     &  nnodesin,inodesout(3*ncont),nnodesout,iactiveline(3,3*ncont),
     &  nactiveline,intersec(2,6*ncont),ipe(*),ime(4,*),nintpoint,k1,j1,
     &  ipiv(4),info,ipnt,one,number_of_nodes,itel,ifac,getiface,
     &  lnode(2,4),nnogap,n1,n2,n3,g(3),iscontr(20),
     &  imcontr(4*4),jj,locm,locs,nodesf,nodem,ifs,
     &  ifm,ns, indexf,icounter,idummy,modf
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  xntersec(3,6*ncont),ph1(3),ph2(3),areax,areay,areaz,area,
     &  pmastsurf(2,*),xl2m(3,8),ets,xis,weight,xl2s(3,8),xsj2(3),
     &  shp2(7,8),shp2s(4,8),dx,help,
     &  xs2(3,2), xquad(2,8), xtri(2,6),
     &  xn(3),xnabs(3),slavstraight(20),
     &  pslavdual(16,*),diag_els(4),m_els(10),contribution,work(4),
     &  contr(20),xs2m(3,2),xsj2m(3),shp2m(7,8),etm,xim,
     &  pslavsurf(3,*)
!     
      include "gauss.f"
!
      data iflag /2/
!     
      checkbiorthogonality=.false.
      checknorm=.false.
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
!     
!     ntri: number of triangles in the triangulation of the master
!     surface corresponding to tie i
!     
         do l = itiefac(1,i), itiefac(2,i)
!     
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
            call getnumberofnodes(nelems,jfaces,lakon,nope,
     &           nopes,idummy) 
!     
!     initialization for Dualshape Coefficient matrix
!     
            ipnt=0
            do k=1,4
               diag_els(k)=0.0
               do j=1,k
                  ipnt=ipnt+1
                  m_els(ipnt)=0.0
               enddo
            enddo
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            checknorm=.false.            
            do j=1,nope
               konl(j)=kon(ipkon(nelems)+j)
            enddo
            nnogap=0
            do m=1,nopes
               ifac=getiface(m,jfaces,nope)
               lnode(1,m)=konl(ifac)
               call nident(islavnode(nslavnode(i)+1),
     &            konl(ifac),(nslavnode(i+1)-nslavnode(i)),id)
               if(islavnode(nslavnode(i)+id)==konl(ifac)) then
                lnode(2,m)=islavact(nslavnode(i)+id)
                if(lnode(2,m).lt.0) nnogap=nnogap+1
                if(checknorm)then
                write(*,*) 'node',lnode(1,m),lnode(2,m)
                endif
               else
                write(*,*)'createbd: node',konl(ifac)
                write(*,*)'was not catalogued properly in islavnode'
                stop
               endif  
               
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifac))+
     &                 vold(j,konl(ifac))
               enddo
            enddo
            ph1(1)=xl2s(1,2)-xl2s(1,1)
            ph1(2)=xl2s(2,2)-xl2s(2,1)
            ph1(3)=xl2s(3,2)-xl2s(3,1)
            ph2(1)=xl2s(1,3)-xl2s(1,1)
            ph2(2)=xl2s(2,3)-xl2s(2,1)
            ph2(3)=xl2s(3,3)-xl2s(3,1)
            areax=((ph1(2)*ph2(3))-(ph2(2)*ph1(3)))**2
            areay=(-(ph1(1)*ph2(3))+(ph2(1)*ph1(3)))**2
            areaz=((ph1(1)*ph2(2))-(ph2(1)*ph1(2)))**2
            area=dsqrt(areax+areay+areaz)/2. 
            mint2d=islavsurf(2,l+1)-islavsurf(2,l)
            if(mint2d==0) cycle
            indexf=islavsurf(2,l)
            do m=1,mint2d
               xis=pslavsurf(1,indexf+m)
               ets=pslavsurf(2,indexf+m)
               weight=pslavsurf(3,indexf+m)
               ns=l
               iflag = 2
!     
!     determining the gap contribution of the integration points
!     and the coefficient for the slave dualshape functions
!     
               if(nopes.eq.8) then
C     > @todo calculation of coeffs of dual shape funktion for quadratic elements not jet implemented
c     >     case nopes.eq.8
                  call shape8q(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
                  dx=dsqrt(xsj2(1)**2+xsj2(2)**2+
     &                 xsj2(3)**2)
                  ipnt=0
                  do k=1,4
                     diag_els(k)=diag_els(k)+shp2(4,k)*weight
     &                    *dx
                     do j=1,k
                        ipnt=ipnt+1
                        m_els(ipnt)=m_els(ipnt)+shp2(4,k)*shp2(4,j)
     &                       *weight
     &                       *dx  
                     enddo
                  enddo
!                  
               elseif(nopes.eq.6) then
C     > @todo calculation of coeffs of dual shape funktion for quadratic elements not yet implemented
c     >     case nopes.eq.6
                  call shape6tri(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
               else
c     > @todo simple calculation of coeffs of dual shape funktion for 3tri not jet implemented
                  call shape3tri(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
                  dx=dsqrt(xsj2(1)**2+xsj2(2)**2+
     &                 xsj2(3)**2)
                  help=3.d0*shp2(4,1)-shp2(4,2)-shp2(4,3)
                  ipnt=0
                  do k=1,3
                     diag_els(k)=diag_els(k)+shp2(4,k)*weight
     &                    *dx
                     do j=1,k
                        ipnt=ipnt+1
                        m_els(ipnt)=m_els(ipnt)+shp2(4,k)*shp2(4,j)
     &                       *weight
     &                       *dx  
                     enddo
                  enddo
               endif 
            enddo
!     
!     Calculate the Mass matrix for compilation of the dualshapefunction 
!     pslavdual(16,*)
!     
!     compute inverse of me_ls
!     factorisation
!     
            if(checknorm .or.checkbiorthogonality)then
               write(*,*)'diag_els',l
               write(*,105)(diag_els(j),j=1,nopes)
            endif
            call dsptrf('U',nopes,m_els,ipiv,info)
!     inverse
            call dsptri('U',nopes,m_els,ipiv,work,info)
            if(checknorm .or. checkbiorthogonality)then
            endif
!     
!     stack of pslavdual multiplication with diag_els
!     
            if(nopes.eq.4)then  
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
            elseif(nopes.eq.3)then
               pslavdual(1,l)=diag_els(1)*m_els(1)
               pslavdual(2,l)=diag_els(1)*m_els(2)
               pslavdual(3,l)=diag_els(1)*m_els(4)
               
               pslavdual(5,l)=diag_els(2)*m_els(2)
               pslavdual(6,l)=diag_els(2)*m_els(3)
               pslavdual(7,l)=diag_els(2)*m_els(5)
               
               pslavdual(9,l)=diag_els(3)*m_els(4)
               pslavdual(10,l)=diag_els(3)*m_els(5)
               pslavdual(11,l)=diag_els(3)*m_els(6)    
            endif
           if(checknorm)then
               write(*,*)'pslavdual1, face',l
               do ii=1,nopes
                     help=0.0
                     do j=1,nopes
                       help=help+pslavdual((ii-1)*4+j,l)
                     enddo
                  write(*,105) (pslavdual((ii-1)*4+j,l),j=1,nopes),help
               enddo
            endif         
            if(nnogap.gt.0 .and. nopes.eq.4) then
             if(nnogap==1) then
               do ii=1,nopes
                if(lnode(2,ii).lt.0) exit
               enddo
                 n1=modf(nopes,ii-2)
                 n2=modf(nopes,ii)
                do jj=1,nopes
                 pslavdual((n1-1)*4+jj,l)=pslavdual((n1-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                 pslavdual((n2-1)*4+jj,l)=pslavdual((n2-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                enddo
                do jj=1,nopes
                 pslavdual((ii-1)*4+jj,l)=0.0
                enddo
             elseif(nnogap==2)then
               do ii=1,nopes
                if(lnode(2,ii).lt.0) then
                  n1=modf(nopes,ii-2)
                  n2=modf(nopes,ii)
c                 write(*,*) 'ii',ii,'n',n1,n2
                  if(lnode(2,n1).lt.0) n1=n2
                  if(lnode(2,n2).lt.0) n2=n1
                 do jj=1,nopes
                  pslavdual((n1-1)*4+jj,l)=pslavdual((n1-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                  pslavdual((n2-1)*4+jj,l)=pslavdual((n2-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                 enddo
                 do jj=1,nopes
                  pslavdual((ii-1)*4+jj,l)=0.0
                 enddo
                endif
               enddo
             elseif(nnogap==3)then
               do ii=1,nopes
               if(lnode(2,ii).gt.-1) n1=ii
               enddo
               do ii=1,nopes
                if(lnode(2,ii).lt.0) then
                 do jj=1,nopes
                  pslavdual((n1-1)*4+jj,l)=pslavdual((n1-1)*4+jj,l)
     &             +pslavdual((ii-1)*4+jj,l)
                 enddo
                 do jj=1,nopes
                  pslavdual((ii-1)*4+jj,l)=0.0
                 enddo
                endif
               enddo
             endif
            endif
            if(nnogap.gt.0 .and. nopes.eq.3) then
             if(nnogap==1) then
               do ii=1,nopes
                if(lnode(2,ii).lt.0) exit
               enddo
                 n1=modf(nopes,ii-2)
                 n2=modf(nopes,ii)
                do jj=1,nopes
                 pslavdual((n1-1)*4+jj,l)=pslavdual((n1-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                 pslavdual((n2-1)*4+jj,l)=pslavdual((n2-1)*4+jj,l)
     &            +0.5*pslavdual((ii-1)*4+jj,l)
                enddo
                do jj=1,nopes
                 pslavdual((ii-1)*4+jj,l)=0.0
                enddo
             elseif(nnogap==2)then
               do ii=1,nopes
               if(lnode(2,ii).gt.-1) n1=ii
               enddo
               do ii=1,nopes
                if(lnode(2,ii).lt.0) then
                 do jj=1,nopes
                  pslavdual((n1-1)*4+jj,l)=pslavdual((n1-1)*4+jj,l)
     &             +pslavdual((ii-1)*4+jj,l)
                 enddo
                 do jj=1,nopes
                  pslavdual((ii-1)*4+jj,l)=0.0
                 enddo
                endif
               enddo

             endif
            endif
            if(checknorm)then
               write(*,*)'pslavdual2, face',l
               do ii=1,nopes
                     help=0.0
                     do j=1,nopes
                       help=help+pslavdual((ii-1)*4+j,l)
                     enddo
                  write(*,105) (pslavdual((ii-1)*4+j,l),j=1,nopes),help
               enddo
            endif
 105        format(4(1x,e15.8))
         enddo
!     
!     FIRST SLAVE SURFACE LOOP DONE
!     
         if(checkbiorthogonality)then
            do l = itiefac(1,i), itiefac(2,i)
               ifaces = islavsurf(1,l)
               nelems = int(ifaces/10)
               jfaces = ifaces - nelems*10
               call getnumberofnodes(nelems,jfaces,lakon,nope,
     &              nopes,idummy)
               mint2d=islavsurf(2,l+1)-islavsurf(2,l)
               if(mint2d==0) cycle
               indexf=islavsurf(2,l)
!     
               do j=1,nope
                  konl(j)=kon(ipkon(nelems)+j)
               enddo
               do m=1,nopes
                  do j=1,3
                     ifac=getiface(m,jfaces,nope)
                     xl2s(j,m)=co(j,konl(ifac))+
     &                    vold(j,konl(ifac))    
                  enddo
               enddo
               ph1(1)=xl2s(1,2)-xl2s(1,1)
               ph1(2)=xl2s(2,2)-xl2s(2,1)
               ph1(3)=xl2s(3,2)-xl2s(3,1)
               ph2(1)=xl2s(1,3)-xl2s(1,1)
               ph2(2)=xl2s(2,3)-xl2s(2,1)
               ph2(3)=xl2s(3,3)-xl2s(3,1)
               areax=((ph1(2)*ph2(3))-(ph2(2)*ph1(3)))**2
               areay=(-(ph1(1)*ph2(3))+(ph2(1)*ph1(3)))**2
               areaz=((ph1(1)*ph2(2))-(ph2(1)*ph1(2)))**2
               area=dsqrt(areax+areay+areaz)/2. 
               do m=1,mint2d
                  xis=pslavsurf(1,indexf+m)
                  ets=pslavsurf(2,indexf+m)
                  weight=pslavsurf(3,indexf+m)
                  ns=l
                  iflag = 2
 100              format('createbd:xs',3(1x,e15.8))
                  if(nopes.eq.8) then
                     call dualshape8q(xis,ets,xl2s,xsj2,xs2,shp2s,iflag)
                  elseif(nopes.eq.4) then
                     call dualshape4q(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                    pslavdual,iflag)
                  elseif(nopes.eq.6) then
                     call dualshape6tri
     &                     (xis,ets,xl2s,xsj2,xs2,shp2s,iflag)
                  else
                     call dualshape3tri(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                    pslavdual,iflag)
                  endif
                  xim=pslavsurf(1,indexf+m)
                  etm=pslavsurf(2,indexf+m)
 101              format('createbd:xm',2(1x,e15.8))
!     
                  if(nopes.eq.8) then
                     call shape8q(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
                  elseif(nopes.eq.4) then
                     call shape4q(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
                  elseif(nopes.eq.6) then
                     call shape6tri(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
                  else
                     call shape3tri(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
                  endif
                  if(m==1)then
                     do j=1,nopes
                        do jj=1,nopes
                           contr(nopes*(j-1)+jj)=0.0
                        enddo
                     enddo
                  endif 
                  do j=1,nopes
                     ifs=getiface(j,jfaces,nope)
                     nodesf=kon(ipkon(nelems)+ifs)
                     locs=j 
                     do jj=1,nopes
                        ifm=getiface(jj,jfaces,nope)
                        nodem=kon(ipkon(nelems)+ifm)
                        locm=jj
                        contribution=shp2s(4,locs)*shp2m(4,locm)
     &                       *pslavsurf(3,indexf+m)
     &                       *dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2) 
                        contr(nopes*(j-1)+jj)=contr(nopes*(j-1)+jj)
     &                       +contribution
                        iscontr(nopes*(j-1)+jj)=nodesf
                        imcontr(nopes*(j-1)+jj)=nodem
                        contribution=0.d0
                     enddo
                  enddo
               enddo 
               if(l.eq.13 .or. l.eq.36 .or. l.eq.59) then
                  write(*,*) 'checkbiorth.: contri,iscontr,imcontr',l
                  do j=1, nopes
                     write(*,105)(contr((j-1)*nopes+jj),jj=1,nopes)
                  enddo 
               endif
            enddo
         endif
!         
         checknorm=.false.
         if(checknorm)then
            do l = itiefac(1,i), itiefac(2,i)
               ifaces = islavsurf(1,l)
               nelems = int(ifaces/10)
               jfaces = ifaces - nelems*10
               call getnumberofnodes(nelems,jfaces,lakon,nope,
     &              nopes,idummy)
               mint2d=islavsurf(2,l+1)-islavsurf(2,l)
               if(mint2d==0) cycle
               indexf=islavsurf(2,l)
!     
               do j=1,nope
                  konl(j)=kon(ipkon(nelems)+j)
               enddo
               do m=1,nopes
                  do j=1,3
                     ifac=getiface(m,jfaces,nope)
                     xl2s(j,m)=co(j,konl(ifac))+
     &                    vold(j,konl(ifac))
                  enddo
               enddo 
               ph1(1)=xl2s(1,2)-xl2s(1,1)
               ph1(2)=xl2s(2,2)-xl2s(2,1)
               ph1(3)=xl2s(3,2)-xl2s(3,1)
               ph2(1)=xl2s(1,3)-xl2s(1,1)
               ph2(2)=xl2s(2,3)-xl2s(2,1)
               ph2(3)=xl2s(3,3)-xl2s(3,1)
               areax=((ph1(2)*ph2(3))-(ph2(2)*ph1(3)))**2
               areay=(-(ph1(1)*ph2(3))+(ph2(1)*ph1(3)))**2
               areaz=((ph1(1)*ph2(2))-(ph2(1)*ph1(2)))**2
               area=dsqrt(areax+areay+areaz)/2. 
               do m=1,mint2d
                  xis=pslavsurf(1,indexf+m)
                  ets=pslavsurf(2,indexf+m)
                  weight=pslavsurf(3,indexf+m)
                  ns=l
                  iflag = 2
                  if(nopes.eq.8) then
                     call dualshape8q(xis,ets,xl2s,xsj2,xs2,shp2s,
     &                    iflag)
                  elseif(nopes.eq.4) then
                     call dualshape4q(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                    pslavdual,iflag)
                  elseif(nopes.eq.6) then
                     call dualshape6tri(xis,ets,xl2s,xsj2,xs2,shp2s,
     &                    iflag)
                  else
                     call dualshape3tri(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                    pslavdual,iflag)
                  endif
                  if(m==1)then
                     do j=1,nopes
                        contr(j)=0.0
                     enddo
                  endif 
                  dx= dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
!     
                  do j=1,nopes
                     ifs=getiface(j,jfaces,nope)    
                     nodesf=kon(ipkon(nelems)+ifs)
                     locs=j
                     contribution=shp2s(4,locs)*
     &                    weight  
     &                    *dx
                     contr(j)=contr(j)
     &                    +contribution
                     iscontr(j)=nodesf
                     contribution=0.d0
                  enddo   
               enddo  
                  write(*,*) 'checknorm: contri,iscontr,imcontr',l
                  do j=1, nopes
                     write(*,*)contr(j),iscontr(j)
                  enddo 
!     
            enddo
         endif
      enddo
!     
      return
      end
