!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine gen3dboun(ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &  xboun,iamboun,typeboun,iponoel,inoel,iponoelmax,kon,ipkon,
     &  lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,co,
     &  nmethod,iperturb,istep,vold,mi)
!
!     connects nodes of 1-D and 2-D elements, for which SPC's were
!     defined, to the nodes of their expanded counterparts
!
      implicit none
!
      logical fixed
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
!
      integer ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),nodeact,
     &  ndirboun(*),idim,ier,matz,ncgnodes,lstart,lend,linc,nnodes,
     &  iamboun(*),iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),ne,
     &  iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikmpc(*),ilmpc(*),rig(*),ntrans,inotr(2,*),nbounold,i,node,
     &  index,ielem,j,indexe,indexk,idir,iamplitude,irotnode,nk,nk_,
     &  newnode,idof,id,mpcfreenew,k,nam,nmethod,iperturb,ndepnodes,
     &  idepnodes(80),l,iexpnode,indexx,irefnode,imax,isol,mpcfreeold,
     &  nod,impc,istep,nrhs,ipiv(3),info,m,mi(*),itr,idirref,id1
!
      real*8 xboun(*),xnor(*),coefmpc(*),trab(7,*),val,co(3,*),
     &  xnoref(3),dmax,d(3,3),e(3,3,3),alpha,q(3),w(3),xn(3),
     &  a1(3),a2(3),dd,c1,c2,c3,ww,c(3,3),vold(0:mi(2),*),a(3,3),
     &  e1(3),e2(3),t1(3),b(3,3),x(3),y(3),fv1(3),dot,
     &  fv2(3),z(3,3),xi1,xi2,xi3,u(3,3),r(3,3)
!
      data d /1.,0.,0.,0.,1.,0.,0.,0.,1./
      data e /0.,0.,0.,0.,0.,-1.,0.,1.,0.,
     &        0.,0.,1.,0.,0.,0.,-1.,0.,0.,
     &        0.,-1.,0.,1.,0.,0.,0.,0.,0./
!
      label='                    '
      fixed=.false.
!
      nbounold=nboun
      do i=1,nbounold
         node=nodeboun(i)
         if(node.gt.iponoelmax) then
            if(ndirboun(i).gt.4) then
               write(*,*) '*ERROR: in gen3dboun: node ',node,
     &              ' does not'
               write(*,*) '        belong to a beam nor shell'
               write(*,*) '        element and consequently has no'
               write(*,*) '        rotational degrees of freedom.'
               write(*,*) '        This may be caused by applying'
               write(*,*) '        a local coordinate system to a'
               write(*,*) '        node in which a rotational boundary'
               write(*,*) '        conditions is defined. In CalculiX'
               write(*,*) '        this is not allowed.'
               call exit(201)
            endif
            cycle
         endif
         index=iponoel(node)
         if(index.eq.0) then
            if(ndirboun(i).gt.4) then
               write(*,*) '*ERROR: in gen3dboun: node ',node,
     &              ' does not'
               write(*,*) '        belong to a beam nor shell'
               write(*,*) '        element and consequently has no'
               write(*,*) '        rotational degrees of freedom.'
               write(*,*) '        This may be caused by applying'
               write(*,*) '        a local coordinate system to a'
               write(*,*) '        node in which a rotational boundary'
               write(*,*) '        conditions is defined. In CalculiX'
               write(*,*) '        this is not allowed.'
               call exit(201)
            endif
            cycle
         endif
         ielem=inoel(1,index)
         j=inoel(2,index)
         indexe=ipkon(ielem)
         indexk=iponor(2,indexe+j)
         idir=ndirboun(i)
         val=xboun(i)
         if(nam.gt.0) iamplitude=iamboun(i)
!
         if(rig(node).ne.0) then
!
!           existing knot
!
            if(idir.gt.4) then
               if(rig(node).lt.0) then
                  write(*,*) '*ERROR in gen3dboun: in node ',node
                  write(*,*) '       a rotational DOF is constrained'
                  write(*,*) '       by a SPC; however, the elements'
                  write(*,*) '       to which this node belongs do not'
                  write(*,*) '       have rotational DOFs'
                  call exit(201)
               endif
               j=idir-4
               irotnode=rig(node)
               type='B'
               call bounadd(irotnode,j,j,val,nodeboun,
     &              ndirboun,xboun,nboun,nboun_,iamboun,
     &              iamplitude,nam,ipompc,nodempc,coefmpc,
     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &              type,typeboun,nmethod,iperturb,fixed,vold,
     &              irotnode,mi,label)
            endif
         else
c!
c!           check for rotational DOFs defined in any but the first step
c!
c            if(idir.gt.4) then
c!
c!              create a knot: determine the knot
c!
c               ndepnodes=0
c               if(lakon(ielem)(7:7).eq.'L') then
c                  do k=1,3
c                     ndepnodes=ndepnodes+1
c                     idepnodes(ndepnodes)=knor(indexk+k)
c                  enddo
c                  idim=1
c               elseif(lakon(ielem)(7:7).eq.'B') then
c                  do k=1,8
c                     ndepnodes=ndepnodes+1
c                     idepnodes(ndepnodes)=knor(indexk+k)
c                  enddo
c                  idim=3
c               else
c                  write(*,*) 
c     &           '*ERROR in gen3dboun: a rotational DOF was applied'
c                  write(*,*) 
c     &           '*      to node',node,' without rotational DOFs'
c                  call exit(201)
c               endif
c!
c!              remove all MPC's in which the knot nodes are
c!              dependent nodes
c!              
c               do k=1,ndepnodes
c                  nod=idepnodes(k)
c                  do l=1,3
c                     idof=8*(nod-1)+l
c                     call nident(ikmpc,idof,nmpc,id)
c                     if(id.gt.0) then
c                        if(ikmpc(id).eq.idof) then
c                           impc=ilmpc(id)
c                           call mpcrem(impc,mpcfree,nodempc,nmpc,
c     &                       ikmpc,ilmpc,labmpc,coefmpc,ipompc)
c                        endif
c                     endif
c                  enddo
c               enddo
c!
c!              generate a rigid body knot
c!
c               irefnode=node
c               nk=nk+1
c               if(nk.gt.nk_) then
c                  write(*,*) '*ERROR in rigidbodies: increase nk_'
c                  call exit(201)
c               endif
c               irotnode=nk
c               rig(node)=irotnode
c               nk=nk+1
c               if(nk.gt.nk_) then
c                  write(*,*) '*ERROR in rigidbodies: increase nk_'
c                  call exit(201)
c               endif
c               iexpnode=nk
c               do k=1,ndepnodes
c                  call knotmpc(ipompc,nodempc,coefmpc,irefnode,
c     &                 irotnode,iexpnode,
c     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
c     &                 nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
c     &                 idepnodes,typeboun,co,xboun,istep,k,ndepnodes,
c     &                 idim,e1,e2,t1)
c               enddo
c!
c!              determine the location of the center of gravity of
c!              the section and its displacements
c!
c               do l=1,3
c                  q(l)=0.d0
c                  w(l)=0.d0
c               enddo
c               if(ndepnodes.eq.3) then
c                  do k=1,ndepnodes,2
c                     nod=idepnodes(k)
c                     do l=1,3
c                        q(l)=q(l)+co(l,nod)
c                        w(l)=w(l)+vold(l,nod)
c                     enddo
c                  enddo
c                  do l=1,3
c                     q(l)=q(l)/2.d0
c                     w(l)=w(l)/2.d0
c                  enddo
c               else
c                  do k=1,ndepnodes
c                     nod=idepnodes(k)
c                     do l=1,3
c                        q(l)=q(l)+co(l,nod)
c                        w(l)=w(l)+vold(l,nod)
c                     enddo
c                  enddo
c                  do l=1,3
c                     q(l)=q(l)/ndepnodes
c                     w(l)=w(l)/ndepnodes
c                  enddo
c               endif
c!
c!              check whether the displacements are zero
c!
c               dd=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
c               if(dd.lt.1.d-20) then
c                  do l=1,3
c                     vold(l,irefnode)=0.d0
c                     vold(l,irotnode)=0.d0
c                     vold(l,iexpnode)=0.d0
c                  enddo
c               else
c!
c!              determine the first displacements of iexpnode
c!
c                  alpha=0.d0
c                  ncgnodes=0
c                  do k=1,ndepnodes
c                     nod=idepnodes(k)
c                     dd=(co(1,nod)-q(1))**2
c     &                    +(co(2,nod)-q(2))**2
c     &                    +(co(3,nod)-q(3))**2
c                     if(dd.lt.1.d-20) then
c                        ncgnodes=ncgnodes+1
c                        cycle
c                     endif
c                     alpha=alpha+dsqrt(
c     &                    ((co(1,nod)+vold(1,nod)-q(1)-w(1))**2
c     &                    +(co(2,nod)+vold(2,nod)-q(2)-w(2))**2
c     &                    +(co(3,nod)+vold(3,nod)-q(3)-w(3))**2)/dd)
c                  enddo
c                  if(ndepnodes-ncgnodes.gt.0) then
c                     alpha=alpha/(ndepnodes-ncgnodes)
c                  endif
c!     
c!     determine the displacements of irotnodes
c!     
c                  do l=1,3
cc                     do m=1,3
cc                        a(l,m)=0.d0
cc                     enddo
c                     xn(l)=0.d0
c                  enddo
c!
c                  ncgnodes=0
c                  do k=1,ndepnodes
c                     nod=idepnodes(k)
c                     dd=0.d0
c                     do l=1,3
c                        a1(l)=co(l,nod)-q(l)
c                        a2(l)=vold(l,nod)-w(l)
c                        dd=dd+a1(l)*a1(l)
c                     enddo
c                     dd=dsqrt(dd)
c                     if(dd.lt.1.d-10) then
c                        ncgnodes=ncgnodes+1
c                        cycle
c                     endif
c                     do l=1,3
c                        a1(l)=a1(l)/dd
c                        a2(l)=a2(l)/dd
c                     enddo
c                     xn(1)=xn(1)+(a1(2)*a2(3)-a1(3)*a2(2))
c                     xn(2)=xn(2)+(a1(3)*a2(1)-a1(1)*a2(3))
c                     xn(3)=xn(3)+(a1(1)*a2(2)-a1(2)*a2(1))
cc                     do l=1,3
cc                        do m=1,3
cc                           a(l,m)=a(l,m)+a1(l)*a1(m)
cc                        enddo
cc                     enddo
c                  enddo
c!     
c                  if(ndepnodes-ncgnodes.gt.0) then
c                     do l=1,3
cc                        do m=1,3
cc                           a(l,m)=a(l,m)/(ndepnodes-ncgnodes)
cc                        enddo
c                        xn(l)=xn(l)/(ndepnodes-ncgnodes)
c                     enddo
c                  endif
c!     
c                  dd=0.d0
c                  do l=1,3
c                     dd=dd+xn(l)*xn(l)
c                  enddo
c                  dd=dsqrt(dd)
c                  do l=1,3
c                     xn(l)=dasin(dd/alpha)*xn(l)/dd
c                  enddo
c!     
c!     determine the displacements of irefnode
c!     
c                  ww=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
c!     
c                  c1=dcos(ww)
c                  if(ww.gt.1.d-10) then
c                     c2=dsin(ww)/ww
c                  else
c                     c2=1.d0
c                  endif
c                  if(ww.gt.1.d-5) then
c                     c3=(1.d0-c1)/ww**2
c                  else
c                     c3=0.5d0
c                  endif
c!     
c!     rotation matrix c
c!     
c                  do k=1,3
c                     do l=1,3
c                        r(k,l)=c1*d(k,l)+
c     &                       c2*(e(k,1,l)*xn(1)+e(k,2,l)*xn(2)+
c     &                       e(k,3,l)*xn(3))+c3*xn(k)*xn(l)
c                     enddo
c                  enddo
c!     
c!     copying the displacements
c!     
c                  do l=1,3
c                     vold(l,irefnode)=w(l)
c                     vold(l,irotnode)=xn(l)
c                  enddo
c                  vold(1,iexpnode)=alpha-1.d0
c!     
c!     correction of the expansion values for beam sections
c!     
c                  if(idim.eq.2) then
c!     
c!     initializing matrices b and c
c!     
c                     do l=1,3
c                        do m=1,3
c                           b(l,m)=0.d0
c                           c(l,m)=0.d0
c                        enddo
c                     enddo
c!     
c!     solving a least squares problem to determine 
c!
c!              start meanrotationmpc
c!              change: mean rotation MPC instead of KNOT
c!
c               idirref=idir-4
c!
c               if(lakon(ielem)(7:7).eq.'L') then
c                  lstart=3
c                  lend=1
c                  linc=-2
c               elseif(lakon(ielem)(7:7).eq.'B') then
c                  lstart=4
c                  lend=1
c                  linc=-1
c               endif
c!
c!              check for transformations
c!
c               if(ntrans.le.0) then
c                  itr=0
c               elseif(inotr(1,node).eq.0) then
c                  itr=0
c               else
c                  itr=inotr(1,node)
c               endif
c!
c!              determine a unit vector on the rotation axis
c!
c!     the transpose of the deformation gradient:
c!     c.F^T=b
c!     
c                     do k=1,ndepnodes
c                        nod=idepnodes(k)
c                        do l=1,3
c                           x(l)=co(l,nod)-q(l)
c                           y(l)=x(l)+vold(l,nod)-w(l)
c                        enddo
c                        do l=1,3
c                           do m=1,3
c                              c(l,m)=c(l,m)+x(l)*x(m)
c                              b(l,m)=b(l,m)+x(l)*y(m)
c                           enddo
c                        enddo
c                     enddo
c!     
c!     solving the linear equation system
c!     
c                     m=3
c                     nrhs=3
c                     call dgesv(m,nrhs,c,m,ipiv,b,m,info)
c                     if(info.ne.0) then
c                        write(*,*) '*ERROR in gen3dforc:'
c                        write(*,*) '       singular system of equations'
c                        call exit(201)
c                     endif
c!     
c!     now b=F^T
c!     
c!     constructing the right stretch tensor
c!     U=F^T.R
c!     
c                     do l=1,3
c                        do m=l,3
c                           u(l,m)=b(l,1)*r(1,m)+b(l,2)*r(2,m)+
c     &                          b(l,3)*r(3,m)
c                        enddo
c                     enddo
c                     u(2,1)=u(1,2)
c                     u(3,1)=u(1,3)
c                     u(3,2)=u(2,3)
c!     
c!     determining the eigenvalues and eigenvectors of U
c!     
c                     m=3
c                     matz=1
c                     ier=0
c                     call rs(m,m,u,w,matz,z,fv1,fv2,ier)
c                     if(ier.ne.0) then
c                        write(*,*) 
c     &                       '*ERROR in knotmpc while calculating the'
c                        write(*,*) '       eigenvalues/eigenvectors'
c                        call exit(201)
c                     endif
c!     
c                     if((dabs(w(1)-1.d0).lt.dabs(w(2)-1.d0)).and.
c     &                    (dabs(w(1)-1.d0).lt.dabs(w(3)-1.d0))) then
c                        l=2
c                        m=3
c                     elseif((dabs(w(2)-1.d0).lt.dabs(w(1)-1.d0)).and.
c     &                       (dabs(w(2)-1.d0).lt.dabs(w(3)-1.d0))) then
c                        l=1
c                        m=3
c                     else
c                        l=1
c                        m=2
c                     endif
c                     xi1=datan2
c     &                   ((z(1,l)*e2(1)+z(2,l)*e2(2)+z(3,l)*e2(2)),
c     &                    (z(1,l)*e1(1)+z(2,l)*e1(2)+z(3,l)*e1(2)))
c                     xi2=w(l)-1.d0
c                     xi3=w(m)-1.d0
c!     
c                     vold(1,iexpnode)=xi1
c                     vold(2,iexpnode)=xi2
c                     vold(3,iexpnode)=xi3
c                  endif
c               endif
c!     
c!     apply the boundary condition
c!     
c               idir=idir-4
c               type='B'
c               call bounadd(irotnode,idir,idir,val,nodeboun,
c     &              ndirboun,xboun,nboun,nboun_,iamboun,
c     &              iamplitude,nam,ipompc,nodempc,coefmpc,
c     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
c     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
c     &              type,typeboun,nmethod,iperturb,fixed,vold,
c     &              irotnode,mi,label)
c!
c!              check for shells whether the rotation about the normal
c!              on the shell has been eliminated
c!               
c               if(lakon(ielem)(7:7).eq.'L') then
c                  indexx=iponor(1,indexe+j)
c                  do j=1,3
c                     xnoref(j)=xnor(indexx+j)
c                  enddo
c                  dmax=0.d0
c                  imax=0
c                  do j=1,3
c                     if(dabs(xnoref(j)).gt.dmax) then
c                        dmax=dabs(xnoref(j))
c                        imax=j
c                     endif
c                  enddo
c!     
c!                 check whether a SPC suffices
c!
c                  if(dabs(1.d0-dmax).lt.1.d-3) then
c                     val=0.d0
c                     if(nam.gt.0) iamplitude=0
c                     type='R'
c                     call bounadd(irotnode,imax,imax,val,nodeboun,
c     &                    ndirboun,xboun,nboun,nboun_,iamboun,
c     &                    iamplitude,nam,ipompc,nodempc,coefmpc,
c     &                    nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
c     &                    ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
c     &                    type,typeboun,nmethod,iperturb,fixed,vold,
c     &                    irotnode,mi,label)
c                  else
c!     
c!                    check for an unused rotational DOF
c!     
c                     isol=0
c                     do l=1,3
c                        idof=8*(node-1)+4+imax
c                        call nident(ikboun,idof,nboun,id)
c                        if((id.gt.0).and.(ikboun(id).eq.idof)) then
c                           imax=imax+1
c                           if(imax.gt.3) imax=imax-3
c                           cycle
c                        endif
c                        isol=1
c                        exit
c                     enddo
c!     
c!     if one of the rotational dofs was not used so far,
c!     it can be taken as dependent side for fixing the
c!     rotation about the normal. If all dofs were used,
c!     no additional equation is needed.
c!     
c                     if(isol.eq.1) then
c                        idof=8*(irotnode-1)+imax
c                        call nident(ikmpc,idof,nmpc,id)
c                        nmpc=nmpc+1
c                        if(nmpc.gt.nmpc_) then
c                           write(*,*) 
c     &                          '*ERROR in gen3dboun: increase nmpc_'
c                           call exit(201)
c                        endif
c!     
c                        ipompc(nmpc)=mpcfree
c                        labmpc(nmpc)='                    '
c!     
c                        do l=nmpc,id+2,-1
c                           ikmpc(l)=ikmpc(l-1)
c                           ilmpc(l)=ilmpc(l-1)
c                        enddo
c                        ikmpc(id+1)=idof
c                        ilmpc(id+1)=nmpc
c!     
c                        nodempc(1,mpcfree)=irotnode
c                        nodempc(2,mpcfree)=imax
c                        coefmpc(mpcfree)=xnoref(imax)
c                        mpcfree=nodempc(3,mpcfree)
c                        imax=imax+1
c                        if(imax.gt.3) imax=imax-3
c                        nodempc(1,mpcfree)=irotnode
c                        nodempc(2,mpcfree)=imax
c                        coefmpc(mpcfree)=xnoref(imax)
c                        mpcfree=nodempc(3,mpcfree)
c                        imax=imax+1
c                        if(imax.gt.3) imax=imax-3
c                        nodempc(1,mpcfree)=irotnode
c                        nodempc(2,mpcfree)=imax
c                        coefmpc(mpcfree)=xnoref(imax)
c                        mpcfreeold=mpcfree
c                        mpcfree=nodempc(3,mpcfree)
c                        nodempc(3,mpcfreeold)=0
c                     endif
c                  endif
c               endif
c               cycle
!
            if(idir.gt.4) then
!
!              start meanrotationmpc
!              change: mean rotation MPC instead of KNOT
!
!              if a mean rotation MPC has already been created 
!              for idof, ilboun(id) contains the index of the
!              SPC created for the angle value
!
               idof=8*(node-1)+idir
               call nident(ikboun,idof,nboun,id)
               if(ilboun(id).ne.i) cycle

               idirref=idir-4
!
               if(lakon(ielem)(7:7).eq.'L') then
                  lstart=3
                  lend=1
                  linc=-2
!
!                 calculating the normal vector =
!                 vector along the drilling direction
!
                  do j=1,3
                     xn(j)=co(j,knor(indexk+3))-co(j,knor(indexk+1))
                  enddo
                  dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
                  do j=1,3
                     xn(j)=xn(j)/dd
                  enddo
!                  
               elseif(lakon(ielem)(7:7).eq.'B') then
                  lstart=4
                  lend=1
                  linc=-1
               endif
!
!              check for transformations
!
               if(ntrans.le.0) then
                  itr=0
               elseif(inotr(1,node).eq.0) then
                  itr=0
               else
                  itr=inotr(1,node)
               endif
!
!              determine a unit vector on the rotation axis
!
               if(itr.eq.0) then
                  do j=1,3
                     do k=1,3
                        a(j,k)=0.d0
                     enddo
                     a(j,j)=1.d0
                  enddo
               else
                  call transformatrix(trab(1,itr),co(1,node),a)
               endif
!
!              check whether the rotation vector does not
!              have a component along the drilling direction
!
               if(lakon(ielem)(7:7).eq.'L') then
                  dot=a(1,idirref)*xn(1)+a(2,idirref)*xn(2)+
     &                a(3,idirref)*xn(3)
                  if(dot.gt.0.05) then
                     do k=5,7
                        if(k.eq.idir) then
                           if(dabs(xboun(i)).lt.1.d-10) cycle
                           write(*,*) '*ERROR in gen3dboun: rotation'
                           write(*,*) '       vector in node ',node
                           write(*,*) '       and direction ',idir-1
                           write(*,*) '       has a significant'
                           write(*,*) 
     &                       '       component along the drilling'
                           write(*,*) '       direction; this is not'
                           write(*,*) '       allowed'
                           call exit(201)
                        endif
                        idof=8*(node-1)+k
                        call nident(ikboun,idof,nboun,id1)
                        if(id1.gt.0) then
                           if(ikboun(id1).eq.idof) then
                              if(dabs(xboun(ilboun(id1))).lt.1.d-10) 
     &                             cycle
                           endif
                        endif
                        write(*,*) '*ERROR in gen3dboun: rotation'
                        write(*,*) '       vector in node ',node
                        write(*,*) '       and direction ',idir-1
                        write(*,*) '       has a significant'
                        write(*,*) 
     &                       '       component along the drilling'
                        write(*,*) '       direction; this is not'
                        write(*,*) '       allowed'
                        call exit(201)
                     enddo
                  endif
               endif
!
!              specific label for mean rotations for beams and
!              shells
!
               label='MEANROTBS           '
               nnodes=0
               do j=lstart,lend,linc
                  nodeact=knor(indexk+j)
                  do k=1,3
                     nnodes=nnodes+1
                     call usermpc(ipompc,nodempc,coefmpc,
     &                    labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                    nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                    nboun,nboun_,nnodes,nodeact,co,label,
     &                    typeboun,iperturb,node,idirref,xboun)
                  enddo
               enddo
!
!              rotation value term
!
               nodeact=nk+1
               do k=1,3
                  co(k,nodeact)=a(k,idirref)
               enddo
               nnodes=nnodes+1
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,nnodes,nodeact,co,label,
     &              typeboun,iperturb,node,idirref,xboun)
!
!              inhomogeneous term
!
               nodeact=0
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,nnodes,nodeact,co,label,
     &              typeboun,iperturb,node,idirref,xboun)
!
!              end meanrotationmpc
!
!              SPC angle term
!
               if(nodeact.ne.-1) then
                  idir=1
                  type='B'
                  call bounadd(nk,idir,idir,val,nodeboun,
     &              ndirboun,xboun,nboun,nboun_,iamboun,
     &              iamplitude,nam,ipompc,nodempc,coefmpc,
     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &              type,typeboun,nmethod,iperturb,fixed,vold,
     &              nk,mi,label)
!
!                 storing the index of the SPC with the angle
!                 value in ilboun(id)
!
                  ilboun(id)=nboun
               endif
!
               cycle
            endif
!     
!           2d element shell element: generate MPC's
!
!           u(n_1)+u(n_3)=2*u(n)
!
            if(lakon(ielem)(7:7).eq.'L') then
               newnode=knor(indexk+1)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=knor(indexk+3)
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
!
!              fixing the temperature degrees of freedom
!
               if(idir.eq.0) then
!
!                 t(n_3)=t(n)
!
                  newnode=knor(indexk+3)
                  idof=8*(newnode-1)+idir
                  call nident(ikmpc,idof,nmpc,id)
                  if((id.le.0).or.(ikmpc(id).ne.idof)) then
                     nmpc=nmpc+1
                     if(nmpc.gt.nmpc_) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
                        call exit(201)
                     endif
                     labmpc(nmpc)='                    '
                     ipompc(nmpc)=mpcfree
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                     nodempc(1,mpcfree)=newnode
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
                        call exit(201)
                     endif
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=-1.d0
                     mpcfreenew=nodempc(3,mpcfree)
                     if(mpcfreenew.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
                        call exit(201)
                     endif
                     nodempc(3,mpcfree)=0
                     mpcfree=mpcfreenew
                  endif
               endif
            elseif(lakon(ielem)(7:7).eq.'B') then
!
!                       1d beam element: generate MPC's
!
!              u(n_1)+u(n_2)+u(n_3)+u(n_4)=4*u(n)
!
               newnode=knor(indexk+1)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  do k=2,4
                     nodempc(1,mpcfree)=knor(indexk+k)
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
                        call exit(201)
                     endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
!
!              fixing the temperature degrees of freedom
!
               if(idir.eq.0) then
                  do k=2,4
!
!                    t(n_k)=t(n), k=2,4
!
                     newnode=knor(indexk+k)
                     idof=8*(newnode-1)+idir
                     call nident(ikmpc,idof,nmpc,id)
                     if((id.le.0).or.(ikmpc(id).ne.idof)) then
                        nmpc=nmpc+1
                        if(nmpc.gt.nmpc_) then
                           write(*,*) 
     &                          '*ERROR in gen3dboun: increase nmpc_'
                           call exit(201)
                        endif
                        labmpc(nmpc)='                    '
                        ipompc(nmpc)=mpcfree
                        do j=nmpc,id+2,-1
                           ikmpc(j)=ikmpc(j-1)
                           ilmpc(j)=ilmpc(j-1)
                        enddo
                        ikmpc(id+1)=idof
                        ilmpc(id+1)=nmpc
                        nodempc(1,mpcfree)=newnode
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=1.d0
                        mpcfree=nodempc(3,mpcfree)
                        if(mpcfree.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dboun: increase nmpc_'
                           call exit(201)
                        endif
                        nodempc(1,mpcfree)=node
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=-1.d0
                        mpcfreenew=nodempc(3,mpcfree)
                        if(mpcfreenew.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dboun: increase nmpc_'
                           call exit(201)
                        endif
                        nodempc(3,mpcfree)=0
                        mpcfree=mpcfreenew
                     endif
                  enddo
               endif
            else
!     
!                       2d plane stress, plane strain or axisymmetric
!                       element: MPC in all but z-direction
!
               newnode=knor(indexk+2)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if(((id.le.0).or.(ikmpc(id).ne.idof)).and.
     &              (idir.ne.3)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-1.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
            endif
         endif
      enddo
!
      return
      end


