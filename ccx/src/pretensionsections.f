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
      subroutine pretensionsections(inpc,textpart,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &  labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,lakon,
     &  kon,ipkon,set,nset,istartset,iendset,ialset,co,ics)
!
!     reading the input deck: *PRE-TENSION SECTION
!
      implicit none
!
      logical twod
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 surface,set(*)
      character*132 textpart(16)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &  n,i,j,key,nk,node,ifacequad(3,4),ifacetria(3,3),npt,
     &  mpcfreeold,ikmpc(*),ilmpc(*),id,idof,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),irefnode,lathyp(3,6),inum,
     &  jn,jt,jd,iside,nelem,jface,nnodelem,nface,nodef(8),nodel(8),
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),indexpret,
     &  k,ipos,nkold,nope,m,kon(*),ipkon(*),indexe,iset,nset,idir,
     &  istartset(*),iendset(*),ialset(*),index1,ics(2,*)
!
      real*8 coefmpc(*),xn(3),xt(3),xd(3),dd,co(3,*)
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
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
!     nodes per face for quad elements
!
      data ifacequad /1,2,5,
     &                2,3,6,
     &                3,4,7,
     &                4,1,8/
!
!     nodes per face for tria elements
!
      data ifacetria /1,2,4,
     &                2,3,5,
     &                3,1,6/
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in pretensionsections.f: *EQUATION should'
         write(*,*) '       be placed before all step definitions'
         stop
      endif
!
      do i=2,n
         if(textpart(i)(1:8).eq.'SURFACE=') then
            surface=textpart(i)(9:88)
            ipos=index(surface,' ')
            surface(ipos:ipos)='T'
         elseif(textpart(i)(1:5).eq.'NODE=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) irefnode
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            if((irefnode.gt.nk).or.(irefnode.le.0)) then
               write(*,*) '*ERROR in pretensionsections.f:'
               write(*,*) '       node ',irefnode,' is not defined'
               stop
            endif
         endif
      enddo
!
!     checking whether the surface exists and is an element face
!     surface
!
      iset=0
      do i=1,nset
         if(set(i).eq.surface) then
            iset=i
            exit
         endif
      enddo
      if(iset.eq.0) then
         write(*,*) '*ERROR in pretensionsections: nonexistent surface'
         write(*,*) '       or surface consists of nodes'
         call inputerror(inpc,ipoinpc,iline)
      endif
!         
!     reading the normal vector and normalizing
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      do i=1,3
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) xn(i)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      enddo
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!
!     finding a unit vector xt perpendicular to the normal vector
!     using a unit vector in x or in y 
!
      if(dabs(xn(1)).lt.0.95d0) then
         xt(1)=1.d0-xn(1)*xn(1)
         xt(2)=-xn(1)*xn(2)
         xt(3)=-xn(1)*xn(3)
      else
         xt(1)=-xn(2)*xn(1)
         xt(2)=1.d0-xn(2)*xn(2)
         xt(3)=-xn(2)*xn(3)
      endif
      dd=dsqrt(xt(1)*xt(1)+xt(2)*xt(2)+xt(3)*xt(3))
      do i=1,3
         xt(i)=xt(i)/dd
      enddo
!
!     xd=xn x xt
!         
      xd(1)=xn(2)*xt(3)-xn(3)*xt(2)
      xd(2)=xn(3)*xt(1)-xn(1)*xt(3)
      xd(3)=xn(1)*xt(2)-xn(2)*xt(1)
!
!     generating a Latin hypercube
!     checking which DOF's of xn, xt and xd are nonzero
!
      do inum=1,6
         if((dabs(xn(lathyp(1,inum))).gt.1.d-3).and.
     &      (dabs(xt(lathyp(2,inum))).gt.1.d-3).and.
     &      (dabs(xd(lathyp(3,inum))).gt.1.d-3)) exit
      enddo
      jn=lathyp(1,inum)
      jt=lathyp(2,inum)
      jd=lathyp(3,inum)
!
!     generating the MPCs
!         
      indexpret=0
      nkold=nk
      m=iendset(iset)-istartset(iset)+1
c!
c!     check whether any MPC was defined in the nodes belonging to
c!     the pre-tension surface
c!
c!     loop over all element faces belonging to the surface
c!      
c      do k=1,m
c         iside=ialset(istartset(iset)+k-1)
c         nelem=int(iside/10.d0)
c         indexe=ipkon(nelem)
c         jface=iside-10*nelem
c!
c!        nnodelem: #nodes in the face
c!        the nodes are stored in nodef(*)
c!
c         if(lakon(nelem)(4:4).eq.'2') then
c            nnodelem=8
c            nface=6
c         elseif(lakon(nelem)(3:4).eq.'D8') then
c            nnodelem=4
c            nface=6
c         elseif(lakon(nelem)(4:5).eq.'10') then
c            nnodelem=6
c            nface=4
c            nope=10
c         elseif(lakon(nelem)(4:4).eq.'4') then
c            nnodelem=3
c            nface=4
c            nope=4
c         elseif(lakon(nelem)(4:5).eq.'15') then
c            if(jface.le.2) then
c               nnodelem=6
c            else
c               nnodelem=8
c            endif
c            nface=5
c            nope=15
c         elseif(lakon(nelem)(3:4).eq.'D6') then
c            if(jface.le.2) then
c               nnodelem=3
c            else
c               nnodelem=4
c            endif
c            nface=5
c            nope=6
c         elseif((lakon(nelem)(2:2).eq.'8').or.
c     &          (lakon(nelem)(4:4).eq.'8')) then
c            nnodelem=3
c            nface=4
c            nope=8
c            if(lakon(nelem)(4:4).eq.'8') then
c               jface=jface-2
c            endif
c         elseif((lakon(nelem)(2:2).eq.'6').or.
c     &          (lakon(nelem)(4:4).eq.'6')) then
c            nnodelem=3
c            nface=3
c            if(lakon(nelem)(4:4).eq.'6') then
c               jface=jface-2
c            endif
c         else
c            cycle
c         endif
c!     
c!     determining the nodes of the face
c!     
c         if(nface.eq.3) then
c            do i=1,nnodelem
c               nodef(i)=kon(indexe+ifacetria(i,jface))
c            enddo
c         elseif(nface.eq.4) then
c            if(nope.eq.8) then
c               do i=1,nnodelem
c                  nodef(i)=kon(indexe+ifacequad(i,jface))
c               enddo
c            else
c               do i=1,nnodelem
c                  nodef(i)=kon(indexe+ifacet(i,jface))
c               enddo
c            endif
c         elseif(nface.eq.5) then
c            if(nope.eq.6) then
c               do i=1,nnodelem
c                  nodef(i)=kon(indexe+ifacew1(i,jface))
c               enddo
c            elseif(nope.eq.15) then
c               do i=1,nnodelem
c                  nodef(i)=kon(indexe+ifacew2(i,jface))
c               enddo
c            endif
c         elseif(nface.eq.6) then
c            do i=1,nnodelem
c               nodef(i)=kon(indexe+ifaceq(i,jface))
c            enddo
c         endif
c!
c!        loop over the nodes belonging to the face      
c!         
c         do i=1,nnodelem
c            node=nodef(i)
c!
c            idof=8*(node-1)+jt
c            call nident(ikmpc,idof,nmpc,id)
c            if(id.gt.0) then
c               if(ikmpc(id).eq.idof) then
c!
c!                 MPC was defined in node: error
c!
c                  write(*,*) '*ERROR in pretensionsections:'
c                  write(*,*) '       a non-pretension MPC is defined'
c                  write(*,*) '       in node ',node,'.'
c                  write(*,*) '       This is not allowed'
c                  stop
c               endif
c            endif
c         enddo
c      enddo
!
!     number of distinct pre-strain nodes for the present keyword
!
      npt=0
!
!     loop over all element faces belonging to the surface
!      
      do k=1,m
         twod=.false.
         iside=ialset(istartset(iset)+k-1)
         nelem=int(iside/10.d0)
         indexe=ipkon(nelem)
         jface=iside-10*nelem
!
!        nnodelem: #nodes in the face
!        the nodes are stored in nodef(*)
!
         if(lakon(nelem)(4:4).eq.'2') then
            nnodelem=8
            nface=6
         elseif(lakon(nelem)(3:4).eq.'D8') then
            nnodelem=4
            nface=6
         elseif(lakon(nelem)(4:5).eq.'10') then
            nnodelem=6
            nface=4
            nope=10
         elseif(lakon(nelem)(4:4).eq.'4') then
            nnodelem=3
            nface=4
            nope=4
         elseif(lakon(nelem)(4:5).eq.'15') then
            if(jface.le.2) then
               nnodelem=6
            else
               nnodelem=8
            endif
            nface=5
            nope=15
         elseif(lakon(nelem)(3:4).eq.'D6') then
            if(jface.le.2) then
               nnodelem=3
            else
               nnodelem=4
            endif
            nface=5
            nope=6
         elseif((lakon(nelem)(2:2).eq.'8').or.
     &          (lakon(nelem)(4:4).eq.'8')) then
            nnodelem=3
            nface=4
            nope=8
            if(lakon(nelem)(4:4).eq.'8') then
               twod=.true.
               jface=jface-2
            endif
         elseif((lakon(nelem)(2:2).eq.'6').or.
     &          (lakon(nelem)(4:4).eq.'6')) then
            nnodelem=3
            nface=3
            if(lakon(nelem)(4:4).eq.'6') then
               twod=.true.
               jface=jface-2
            endif
         else
            cycle
         endif
!     
!     determining the nodes of the face
!     
         if(nface.eq.3) then
            do i=1,nnodelem
               nodef(i)=kon(indexe+ifacetria(i,jface))
               nodel(i)=ifacetria(i,jface)
            enddo
         elseif(nface.eq.4) then
            if(nope.eq.8) then
               do i=1,nnodelem
                  nodef(i)=kon(indexe+ifacequad(i,jface))
                  nodel(i)=ifacequad(i,jface)
               enddo
            else
               do i=1,nnodelem
                  nodef(i)=kon(indexe+ifacet(i,jface))
                  nodel(i)=ifacet(i,jface)
               enddo
            endif
         elseif(nface.eq.5) then
            if(nope.eq.6) then
               do i=1,nnodelem
                  nodef(i)=kon(indexe+ifacew1(i,jface))
                  nodel(i)=ifacew1(i,jface)
               enddo
            elseif(nope.eq.15) then
               do i=1,nnodelem
                  nodef(i)=kon(indexe+ifacew2(i,jface))
                  nodel(i)=ifacew2(i,jface)
               enddo
            endif
         elseif(nface.eq.6) then
            do i=1,nnodelem
               nodef(i)=kon(indexe+ifaceq(i,jface))
               nodel(i)=ifaceq(i,jface)
            enddo
         endif
!
!        loop over the nodes belonging to the face      
!         
         do i=1,nnodelem
            node=nodef(i)
            call nident2(ics,node,npt,id)
            if(id.gt.0) then
               if(ics(1,id).eq.node) then
!
!                 node was already treated: replacing the node
!                 by the partner node
!
                  kon(indexe+nodel(i))=ics(2,id)
                  cycle
               endif
            endif
c
c            call nident(ikmpc,idof,nmpc,id)
c            if(id.gt.0) then
c               if(ikmpc(id).eq.idof) then
c!
c!                 node was already treated: replacing the node
c!                 by the partner node
c!
cc                  kon(indexe+nodel(i))=nodempc(1,nodempc(3,
cc     &                 nodempc(3,nodempc(3,ipompc(ilmpc(id))))))
c                  index1=ipompc(ilmpc(id))
c                  do
c                     if(nodempc(1,index1).ne.node) then
c                        kon(indexe+nodel(i))=nodempc(1,index1)
c                        exit
c                     else
c                        index1=nodempc(3,index1)
c                     endif
c                  enddo
c!
c                  cycle
c               endif
c            endif
!
!           generating a partner node
!
            nk=nk+1
!
!           coordinates for the new node
!
            do j=1,3
               co(j,nk)=co(j,node)
            enddo
!
!           updating the topology
!
            kon(indexe+nodel(i))=nk
!
!           updating ics
!
            npt=npt+1
            do j=npt,id+2,-1
               ics(1,j)=ics(1,j-1)
               ics(2,j)=ics(2,j-1)
            enddo
            ics(1,id+1)=node
            ics(2,id+1)=nk
!
!           first MPC perpendicular to the normal direction
!
c            idof=8*(node-1)+jt
            idof=8*(nk-1)+jt
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in equations: increase nmpc_'
               stop
            endif
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)='                    '
!
!           updating ikmpc and ilmpc
!
            do j=nmpc,id+2,-1
               ikmpc(j)=ikmpc(j-1)
               ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!
            idir=jt
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=jt
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=jt
               coefmpc(mpcfree)=xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xt(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=xt(idir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
            nodempc(3,mpcfreeold)=0
!
!           second MPC perpendicular to the normal direction
!
            if(.not.twod) then
c               idof=8*(node-1)+jd
               idof=8*(nk-1)+jd
               call nident(ikmpc,idof,nmpc,id)
!     
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) '*ERROR in equations: increase nmpc_'
                  stop
               endif
               labmpc(nmpc)='                    '
               ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
!     
               idir=jd
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=node
                  nodempc(1,mpcfree)=nk
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
!     
               idir=idir+1
               if(idir.eq.4) idir=1
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=node
                  nodempc(1,mpcfree)=nk
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
!     
               idir=idir+1
               if(idir.eq.4) idir=1
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=node
                  nodempc(1,mpcfree)=nk
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
!     
               idir=jd
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=nk
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
!     
               idir=idir+1
               if(idir.eq.4) idir=1
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=nk
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
!     
               idir=idir+1
               if(idir.eq.4) idir=1
               if(dabs(xd(idir)).gt.1.d-10) then
c                  nodempc(1,mpcfree)=nk
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=xd(idir)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
               endif
               nodempc(3,mpcfreeold)=0
            endif
!     
!     MPC in normal direction
!     
!           check whether initialized
!
            if(indexpret.eq.0) then
c               idof=8*(node-1)+jn
               idof=8*(nk-1)+jn
               call nident(ikmpc,idof,nmpc,id)
!
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) '*ERROR in equations: increase nmpc_'
                  stop
               endif
               labmpc(nmpc)='PRETENSION          '
               ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
            else
               nodempc(3,indexpret)=mpcfree
            endif
!
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=node
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=-xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
c               nodempc(1,mpcfree)=nk
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=xn(idir)
               indexpret=mpcfree
               mpcfree=nodempc(3,mpcfree)
            endif
!
         enddo
      enddo
!
      nodempc(3,indexpret)=mpcfree
      nodempc(1,mpcfree)=irefnode
      nodempc(2,mpcfree)=1
      coefmpc(mpcfree)=1.d0*(nk-nkold)
      mpcfreeold=mpcfree
      mpcfree=nodempc(3,mpcfree)
      nodempc(3,mpcfreeold)=0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
c      do i=1,nmpc
c         call writempc(ipompc,nodempc,coefmpc,labmpc,i)
c      enddo
c      do i=1,nmpc
c         write(*,*) i,ikmpc(i),ilmpc(i)
c      enddo
!
      return
      end

