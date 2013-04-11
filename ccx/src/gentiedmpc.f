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
      subroutine gentiedmpc(tieset,ntie,itietri,ipkon,kon,
     &  lakon,set,istartset,iendset,ialset,cg,straight,
     &  koncont,co,xo,yo,zo,x,y,z,nx,ny,nz,nset,
     &  ifaceslave,istartfield,iendfield,ifield,
     &  ipompc,nodempc,coefmpc,nmpc,nmpctied,mpcfree,ikmpc,ilmpc,
     &  labmpc,ithermal,tietol,cfd,ncont,imastop)
!
!     generates MPC's for the slave tied contact nodes
!
      implicit none
!
      character*1 c
      character*3 m11,m2,m3
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,istartset(*),iendset(*),ialset(*),
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &  neigh(1),iflag,kneigh,i,j,k,l,isol,itri,ll,kflag,n,nx(*),
     &  ny(*),ipointer(1),nz(*),nstart,ifaceq(8,6),ifacet(6,4),
     &  ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,imastop(3,*),
     &  nnodelem,nface,nope,nodef(8),idof,kstart,kend,jstart,id,
     &  jend,ifield(*),istartfield(*),iendfield(*),ifaceslave(*),
     &  ipompc(*),nodempc(3,*),nmpc,nmpctied,mpcfree,ikmpc(*),
     &  ilmpc(*),ithermal(2),cfd,ncont,mpcfreeold,m,one,number_of_nodes,
     &  itriold,itrinew,ntriangle,ntriangle_,itriangle(100)
!
      real*8 cg(3,*),straight(16,*),co(3,*),p(3),
     &  dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),pl(3,8),
     &  ratio(8),xi,et,coefmpc(*),tietol(2,*),tolloc
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
!     opening a file to store the nodes which are not connected
!
c      open(9,file='nodes_not_connected.fbd',status='unknown',err=51)
c      close(9,status='delete',err=52)
c      open(9,file='nodes_not_connected.fbd',status='unknown',err=51)
      open(40,file='WarnNodeMissMasterIntersect.nam',status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeMissMasterIntersect'
!
!     storing the triangulation of the master surfaces
!     
      open(70,file='TriMasterContactTie.frd',status='unknown')
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
      nmpctied=nmpc
!
!     calculating a typical element size
!
      tolloc=0.d0
      do i=1,ncont
         tolloc=tolloc+dabs(straight(1,i)*cg(1,i)+
     &               straight(2,i)*cg(2,i)+
     &               straight(3,i)*cg(3,i)+
     &               straight(4,i))
      enddo
      tolloc=0.025*tolloc/ncont
!
!     determining for which dofs MPC's have to be generated
!
      if(cfd.eq.1) then
         if(ithermal(2).le.1) then
            kstart=1
            kend=4
         else
            kstart=0
            kend=4
         endif
      else
         if(ithermal(2).le.1) then
            kstart=1
            kend=3
         elseif(ithermal(2).eq.2) then
            kstart=0
            kend=0
         else
            kstart=0
            kend=3
         endif
      endif
!
!     maximum number of neighboring master triangles for a slave node
!
      kflag=2
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'T') cycle
         iflag=0
         kneigh=1
         slavset=tieset(2,i)
!
!        default tolerance if none is specified 
!
         if(tietol(1,i).lt.1.d-10) tietol(1,i)=tolloc
!     
!     determining the slave set
!     
         if(ifaceslave(i).eq.0) then
c            ipos=index(slavset,' ')
c            slavset(ipos:ipos)='S'
            do j=1,nset
               if(set(j).eq.slavset) then
                  exit
               endif
            enddo
c            if(j.gt.nset) then
c               write(*,*) 
c     &              '*ERROR in gentiedmpc: tied contact slave set',
c     &              slavset
c               write(*,*) '       does not exist'
c               stop
c            endif
            jstart=istartset(j)
            jend=iendset(j)
         else
            jstart=istartfield(i)
            jend=iendfield(i)
         endif
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
         do j=jstart,jend
            if(((ifaceslave(i).eq.0).and.(ialset(j).gt.0)).or.
     &          (ifaceslave(i).eq.1)) then
!     
               if(ifaceslave(i).eq.0) then
                  node=ialset(j)
               else
                  node=ifield(j)
               endif
!     
c               write(*,*) 'gentiedmpc ',j,node
               do k=1,3
                  p(k)=co(k,node)
               enddo
!     
!              determining the kneigh neighboring master contact
!              triangle centers of gravity
!   
               call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &             n,neigh,kneigh)
!     
               isol=0
!     
c               do k=1,kneigh
c                  itri=neigh(k)+itietri(1,i)-1
c!     
c                  totdist(k)=0.d0
c!     
c                  do l=1,3
c                     ll=4*l-3
c                     dist=straight(ll,itri)*p(1)+
c     &                    straight(ll+1,itri)*p(2)+
c     &                    straight(ll+2,itri)*p(3)+
c     &                    straight(ll+3,itri)
c                     if(dist.gt.0.d0) then
c                        totdist(k)=totdist(k)+dist
c                     endif
c                  enddo
cc                  write(*,*) 'gentiedmpc ',k,itri,koncont(4,itri),
cc     &                       totdist(k)
c                  totdist(k)=dsqrt(totdist(k)**2+
c     &                (straight(13,itri)*p(1)+
c     &                 straight(14,itri)*p(2)+
c     &                 straight(15,itri)*p(3)+
c     &                 straight(16,itri))**2)
cc                  cgdist=dsqrt((p(1)-cg(1,itri))**2+
cc     &                         (p(2)-cg(2,itri))**2+
cc     &                         (p(3)-cg(3,itri))**2)
cc                  write(*,*) 'gentiedmpc ',k,itri,koncont(4,itri),
cc     &                       totdist(k),cgdist
c!     
c                  if(totdist(k).le.tietol(1,i)) then
c                     isol=k
c                     exit
c                  endif
c               enddo
!     
               isol=0
!
               itriold=0
               itri=neigh(1)+itietri(1,i)-1
               ntriangle=0
               ntriangle_=100
!
               loop1: do
                  do l=1,3
                     ll=4*l-3
                     dist=straight(ll,itri)*p(1)+
     &                    straight(ll+1,itri)*p(2)+
     &                    straight(ll+2,itri)*p(3)+
     &                    straight(ll+3,itri)
c                     if(dist.gt.1.d-6) then
                     if(dist.gt.tietol(1,i)) then
                        itrinew=imastop(l,itri)
                        if(itrinew.eq.0) then
c                           write(*,*) '**border reached'
                           isol=-1
                           exit loop1
                        elseif(itrinew.eq.itriold) then
c                           write(*,*) '**solution in between triangles'
                           isol=itri
                           exit loop1
                        else
                           call nident(itriangle,itrinew,ntriangle,id)
                           if(id.gt.0) then
                              if(itriangle(id).eq.itrinew) then
c                                 write(*,*) '**circular path; no solution'
                                 isol=-2
                                 exit loop1
                              endif
                           endif
                           ntriangle=ntriangle+1
                           if(ntriangle.gt.ntriangle_) then
c                              write(*,*) '**too many iterations'
                              isol=-3
                              exit loop1
                           endif
                           do k=ntriangle,id+2,-1
                              itriangle(k)=itriangle(k-1)
                           enddo
                           itriangle(id+1)=itrinew
                           itriold=itri
                           itri=itrinew
                           cycle loop1
                        endif
                     elseif(l.eq.3) then
c                              write(*,*) '**regular solution'
                        isol=itri
                        exit loop1
                     endif
                  enddo
               enddo loop1
!
!              if an opposite triangle is found: check the distance
!              perpendicular to the triangle
!
               if(isol.gt.0) then
                  dist=dsqrt(straight(13,itri)*p(1)+
     &                 straight(14,itri)*p(2)+
     &                 straight(15,itri)*p(3)+
     &                 straight(16,itri))**2
                  if(dist.gt.tietol(1,i)) isol=0
               endif
!
               if(isol.le.0) then
!
!                 no MPC is generated
!
                  write(*,*) '*WARNING in gentiedmpc: no tied MPC'
                  write(*,*) '         generated for node ',node
                  if(isol.eq.0) then
                     write(*,*) '         master face too far away'
                     write(*,*) '         distance: ',dist
                     write(*,*) '         tolerance: ',tietol(1,i)
                  else
                     write(*,*) '         no corresponding master face'
                     write(*,*) '         found; tolerance: ',
     &                       tietol(1,i)
                  endif
                  write(40,*) node
                else
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
!                 attaching the node with coordinates in p
!                 to the face
!
                  do k=1,nnodelem
                     do l=1,3
                        pl(l,k)=co(l,nodef(k))
                     enddo
                  enddo
                  call attach(pl,p,nnodelem,ratio,dist,xi,et)
                  do k=1,3
                     co(k,node)=p(k)
                  enddo
!
!                 generating MPC's
!
                  do l=kstart,kend
                     idof=8*(node-1)+l
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,*) '*WARNING in gentiedmpc:'
                           write(*,*) '         DOF ',l,' of node ',
     &                          node,' is not active;'
                           write(*,*) '         no tied constraint ',
     &                                'is generated'
                           write(40,*) node
                           cycle
                        endif
                     endif
!
                     nmpc=nmpc+1
                     labmpc(nmpc)='                    '
                     ipompc(nmpc)=mpcfree
!     
!                    updating ikmpc and ilmpc
!     
                     do m=nmpc,id+2,-1
                        ikmpc(m)=ikmpc(m-1)
                        ilmpc(m)=ilmpc(m-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
!     
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=l
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*)
     &                    '*ERROR in gentiedmpc: increase memmpc_'
                        stop
                     endif
                     do k=1,nnodelem
                        nodempc(1,mpcfree)=nodef(k)
                        nodempc(2,mpcfree)=l
                        coefmpc(mpcfree)=-ratio(k)
                        mpcfreeold=mpcfree
                        mpcfree=nodempc(3,mpcfree)
                        if(mpcfree.eq.0) then
                           write(*,*)
     &                      '*ERROR in gentiedmpc: increase memmpc_'
                           stop
                        endif
                     enddo
                     nodempc(3,mpcfreeold)=0

c                     call writempc(ipompc,nodempc,coefmpc,labmpc,nmpc)
                  enddo
!
               endif
!     
            else
               node=ialset(j-2)
               do
                  node=node-ialset(j)
                  if(node.ge.ialset(j-1)) exit
!     
                  do k=1,3
                     p(k)=co(k,node)
                  enddo
!     
!                 determining the kneigh neighboring master contact
!                 triangle centers of gravity
!     
                  call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &                 n,neigh,kneigh)
!     
                  isol=0
!     
c                  do k=1,kneigh
c                     itri=neigh(k)+itietri(1,i)-1
c!     
c                     totdist(k)=0.d0
c!     
c                     do l=1,3
c                        ll=4*l-3
c                        dist=straight(ll,itri)*p(1)+
c     &                       straight(ll+1,itri)*p(2)+
c     &                       straight(ll+2,itri)*p(3)+
c     &                       straight(ll+3,itri)
c                        if(dist.gt.0.d0) then
c                           totdist(k)=totdist(k)+dist
c                        endif
c                     enddo
cc                  write(*,*) 'gentiedmpc ',k,itri,koncont(4,itri),
cc     &                       totdist(k)
c                     totdist(k)=dsqrt(totdist(k)**2+
c     &                    (straight(13,itri)*p(1)+
c     &                    straight(14,itri)*p(2)+
c     &                    straight(15,itri)*p(3)+
c     &                    straight(16,itri))**2)
cc                     cgdist=dsqrt((p(1)-cg(1,itri))**2+
cc     &                    (p(2)-cg(2,itri))**2+
cc     &                    (p(3)-cg(3,itri))**2)
cc                  write(*,*) 'gentiedmpc ',k,itri,koncont(4,itri),
cc     &                       totdist(k),cgdist
c!     
c                     if(totdist(k).le.tietol(1,i)) then
c                        isol=k
c                        exit
c                     endif
c                  enddo
!     
               isol=0
!
               itriold=0
               itri=neigh(1)+itietri(1,i)-1
               ntriangle=0
               ntriangle_=100
!
               loop2: do
                  do l=1,3
                     ll=4*l-3
                     dist=straight(ll,itri)*p(1)+
     &                    straight(ll+1,itri)*p(2)+
     &                    straight(ll+2,itri)*p(3)+
     &                    straight(ll+3,itri)
c                     if(dist.gt.1.d-6) then
                     if(dist.gt.tietol(1,i)) then
                        itrinew=imastop(l,itri)
                        if(itrinew.eq.0) then
c                           write(*,*) '**border reached'
                           isol=-1
                           exit loop2
                        elseif(itrinew.eq.itriold) then
c                           write(*,*) '**solution in between triangles'
                           isol=itri
                           exit loop2
                        else
                           call nident(itriangle,itrinew,ntriangle,id)
                           if(id.gt.0) then
                              if(itriangle(id).eq.itrinew) then
c                                 write(*,*) '**circular path; no solution'
                                 isol=-2
                                 exit loop2
                              endif
                           endif
                           ntriangle=ntriangle+1
                           if(ntriangle.gt.ntriangle_) then
c                              write(*,*) '**too many iterations'
                              isol=-3
                              exit loop2
                           endif
                           do k=ntriangle,id+2,-1
                              itriangle(k)=itriangle(k-1)
                           enddo
                           itriangle(id+1)=itrinew
                           itriold=itri
                           itri=itrinew
                           cycle loop2
                        endif
                     elseif(l.eq.3) then
c                              write(*,*) '**regular solution'
                        isol=itri
                        exit loop2
                     endif
                  enddo
               enddo loop2
!
!              if an opposite triangle is found: check the distance
!              perpendicular to the triangle
!
               if(isol.gt.0) then
                  dist=dsqrt(straight(13,itri)*p(1)+
     &                 straight(14,itri)*p(2)+
     &                 straight(15,itri)*p(3)+
     &                 straight(16,itri))**2
                  if(dist.gt.tietol(1,i)) isol=0
               endif
!     
!     check whether distance is larger than tietol(1,i):
!     no element is generated
!     
                  if(isol.eq.0) then
!
!                 no MPC is generated
!
                     write(*,*) '*WARNING in gentiedmpc: no tied MPC'
                     write(*,*) '         generated for node ',node
                     if(isol.eq.0) then
                        write(*,*) '         master face too far away'
                        write(*,*) '         distance: ',dist
                        write(*,*) '         tolerance: ',tietol(1,i)
                     else
                      write(*,*) '         no corresponding master face'
                      write(*,*) '         found; tolerance: ',
     &                                tietol(1,i)
                     endif
                     write(40,*) node
                  else
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
!     determining the nodes of the face
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
!                 attaching the node with coordinates in p
!                 to the face
!
                     do k=1,nnodelem
                        do l=1,3
                           pl(l,k)=co(l,nodef(k))
                        enddo
                     enddo
                     call attach(pl,p,nnodelem,ratio,dist,xi,et)
                     do k=1,3
                        co(k,node)=p(k)
                     enddo
!     
!                    generating MPC's
!     
                     do l=kstart,kend
                        idof=8*(node-1)+l
                        call nident(ikmpc,idof,nmpc,id)
                        if(id.gt.0) then
                           if(ikmpc(id).eq.idof) then
                              write(*,*) '*WARNING in gentiedmpc:'
                              write(*,*) '         DOF ',l,' of node ',
     &                             node,' is not active;'
                              write(*,*) '         no tied constraint ',
     &                             'is generated'
                              write(40,*) node
                              cycle
                           endif
                        endif
!     
                        nmpc=nmpc+1
                        labmpc(nmpc)='                    '
                        ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
                        do m=nmpc,id+2,-1
                           ikmpc(m)=ikmpc(m-1)
                           ilmpc(m)=ilmpc(m-1)
                        enddo
                        ikmpc(id+1)=idof
                        ilmpc(id+1)=nmpc
!     
                        nodempc(1,mpcfree)=node
                        nodempc(2,mpcfree)=l
                        coefmpc(mpcfree)=1.d0
                        mpcfree=nodempc(3,mpcfree)
                        if(mpcfree.eq.0) then
                           write(*,*)
     &                      '*ERROR in gentiedmpc: increase memmpc_'
                           stop
                        endif
                        do k=1,nnodelem
                           nodempc(1,mpcfree)=nodef(k)
                           nodempc(2,mpcfree)=l
                           coefmpc(mpcfree)=-ratio(k)
                           mpcfreeold=mpcfree
                           mpcfree=nodempc(3,mpcfree)
                           if(mpcfree.eq.0) then
                              write(*,*)
     &                    '*ERROR in gentiedmpc: increase memmpc_'
                              stop
                           endif
                        enddo
                        nodempc(3,mpcfreeold)=0
                     enddo
                  endif
!     
               enddo
            endif
         enddo
      enddo
!
!     number of tied MPC's
!
      nmpctied=nmpc-nmpctied
!
      close(40)
!     
      return
!
c 51   write(*,*) '*ERROR in openfile: could not open file ',
c     &  'nodes_not_connected.fbd'
c      stop
c 52   write(*,*) '*ERROR in openfile: could not delete file ',
c     &  'nodes_not_connected.fbd'
c      stop
!
      end

