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
      subroutine gencontelem(tieset,ntie,itietri,ne,ipkon,kon,
     &  lakon,set,istartset,iendset,ialset,cg,straight,ifree,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,ielmat,cs,
     &  elcon,istep,iinc,iit,ncmat_,ntmat_,ifcont1,ifcont2,ne0,
     &  vini,nmethod,mi,imastop)
!
!     generate contact elements for the slave contact nodes
!
      implicit none
!
      character*8 lakon(*)
c      character*18 cfile
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,istartset(*),iendset(*),ialset(*),ifree,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),ne,node,
     &  neigh(1),nodeedge(2,10),iflag,kneigh,i,j,k,l,islav,isol,
     &  itri,ll,kflag,n,ipos,nx(*),ny(*),ipointer(10),istep,iinc,
     &  nz(*),nstart,ielmat(*),material,ifaceq(8,6),ifacet(6,4),
     &  ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,iit,
     &  nnodelem,nface,nope,nodef(8),ncmat_,ntmat_,ifcont1(*),
     &  ifcont2(*),ne0,ifaceref,isum,nmethod,mi(2),iteller,
     &  imastop(3,*), itriangle(100),ntriangle,ntriangle_,itriold,
     &  itrinew,id
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  totdist(10),dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),cs(17,*),
     &  beta,c0,elcon(0:ncmat_,ntmat_,*),vini(0:mi(2),*)
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
      data iteller /0/
      save iteller
!
!     maximum number of neighboring master triangles for a slave node
!
      kflag=2
!
!     opening a file to store the contact spring elements
!    
c      iteller=iteller+1
c      cfile(1:18)='contactelem       '
c      if(iteller.lt.10) then
c         write(cfile(12:12),'(i1)') iteller
c         cfile(13:16)='.inp'
c      elseif(iteller.lt.100) then
c         write(cfile(12:13),'(i2)') iteller
c         cfile(14:17)='.inp'
c      elseif(iteller.lt.1000) then
c         write(cfile(12:14),'(i3)') iteller
c         cfile(15:18)='.inp'
c      else
c         write(*,*) '*ERROR in gencontelem: more than 1000'
c         write(*,*) '       contact element files'
c         stop
c      endif
c      open(27,file=cfile,status='unknown')
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         iflag=0
         kneigh=1
         slavset=tieset(2,i)
         material=int(cs(1,i))
!     
!     determining the slave set
!     
         do j=1,nset
            if(set(j).eq.slavset) exit
         enddo
         if(j.gt.nset) then
            write(*,*) '*ERROR in gencontelem: contact slave set',
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
         do j=istartset(islav),iendset(islav)
            if(ialset(j).gt.0) then
!     
               node=ialset(j)
!     
               do k=1,3
                  p(k)=co(k,node)+vold(k,node)
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
c                     if(dist.gt.0.d0) then
                     if(dist.gt.1.d-6) then
                        itrinew=imastop(l,itri)
                        if(itrinew.eq.0) then
c                           write(*,*) '**border reached'
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
                                 exit loop1
                              endif
                           endif
                           ntriangle=ntriangle+1
                           if(ntriangle.gt.ntriangle_) then
c                              write(*,*) '**too many iterations'
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
c               do k=1,kneigh
c                  itri=neigh(k)+itietri(1,i)-1
c!     
c                  ipos=0
c                  totdist(k)=0.d0
c                  nodeedge(1,k)=0
c                  nodeedge(2,k)=0
c!     
c                  do l=1,3
c                     ll=4*l-3
c                     dist=straight(ll,itri)*p(1)+
c     &                    straight(ll+1,itri)*p(2)+
c     &                    straight(ll+2,itri)*p(3)+
c     &                    straight(ll+3,itri)
c                     if(dist.gt.0.d0) then
c                        totdist(k)=totdist(k)+dist
c                        if(ipos.eq.0) then
c                           nodeedge(1,k)=koncont(l,itri)
c                           if(l.ne.3) then
c                              nodeedge(2,k)=koncont(l+1,itri)
c                           else
c                              nodeedge(2,k)=koncont(1,itri)
c                           endif
c                        else
c                           if((nodeedge(1,k).eq.koncont(l,itri)).or.
c     &                      (nodeedge(2,k).eq.koncont(l,itri)))then
c                              nodeedge(1,k)=koncont(l,itri)
c                              nodeedge(2,k)=0
c                           else
c                              if(l.ne.3) then
c                                 nodeedge(1,k)=koncont(l+1,itri)
c                              else
c                                 nodeedge(1,k)=koncont(1,itri)
c                              endif
c                           endif
c                        endif
c                        ipos=ipos+1
c                     endif
c                  enddo
c!     
c                  if(totdist(k).le.0.d0) then
c                     isol=k
c                     exit
c                  endif
c               enddo
c!
c!              if no independent face was found, a small
c!              tolerance is applied
c!
c               if(isol.eq.0) then
c                  do k=1,kneigh
c                     ipointer(k)=neigh(k)+itietri(1,i)-1
c                  enddo
c                  call dsort(totdist,ipointer,kneigh,kflag)
c                  do k=1,kneigh
c                     itri=ipointer(k)
c                     dist=dabs(straight(1,itri)*cg(1,itri)+
c     &                         straight(2,itri)*cg(2,itri)+
c     &                         straight(3,itri)*cg(3,itri)+
c     &                         straight(4,itri))
c                     if(totdist(k).lt.1.d-3*dist) then
c                        isol=k
c                        exit
c                     endif
c                  enddo
c               endif
!
!              check whether distance is larger than c0:
!              no element is generated
!
               if(isol.ne.0) then
                  dist=straight(13,itri)*p(1)+
     &                 straight(14,itri)*p(2)+
     &                 straight(15,itri)*p(3)+
     &                 straight(16,itri)
                  beta=elcon(1,1,material)
                  if(beta.gt.0.d0) then
                     c0=dlog(100.d0)/beta
                  else
                     c0=0.d0
                  endif
                  if(dist.gt.c0) then
                     isol=0
!
!                    adjusting the bodies at the start of the
!                    calculation such that they touch
!
                  elseif((istep.eq.1).and.(iinc.eq.1).and.
     &                   (iit.le.0).and.(dist.lt.0.d0).and.
     &                   (nmethod.eq.1)) then
                     do k=1,3
                        vold(k,node)=vold(k,node)-
     &                        dist*straight(12+k,itri)
                        vini(k,node)=vold(k,node)
                    enddo
                  endif
               endif
!     
               if(isol.eq.0) then
!
!                 no independent face was found: no spring
!                 element is generated
!
               else
!     
!                 plane spring
!     
                  ne=ne+1
                  ipkon(ne)=ifree
                  lakon(ne)='ESPRNGC '
                  ielmat(ne)=material
                  nelem=int(koncont(4,itri)/10.d0)
                  jface=koncont(4,itri)-10*nelem
!
!                 storing the face in ifcont1 and the
!                 element number in ifcont2
!
                  ifcont1(ne-ne0)=koncont(4,itri)
                  ifcont2(ne-ne0)=ne
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
                     kon(ifree+k)=nodef(k)
                  enddo
                  ifree=ifree+nnodelem+1
                  kon(ifree)=node
                  ifree=ifree+1
!
                  write(lakon(ne)(8:8),'(i1)') nnodelem+1
c                  write(*,*) 'new elem',ne,(nodef(k),k=1,nnodelem),node
                  if((nnodelem.eq.3).or.(nnodelem.eq.6)) then
c                     write(27,100)
c 100                 format('*ELEMENT,TYPE=C3D4')
c                     write(27,*) ne,',',nodef(1),',',nodef(2),',',
c     &                           nodef(3),',',node
                  else
c                     write(27,101)
c 101                 format('*ELEMENT,TYPE=C3D6')
c                     write(27,*) ne,',',nodef(2),',',node,',',nodef(3),
c     &                    ',',nodef(1),',',node,',',nodef(4)
                  endif
               endif
!     
            else
               node=ialset(j-2)
               do
                  node=node-ialset(j)
                  if(node.ge.ialset(j-1)) exit
!     
                  do k=1,3
                     p(k)=co(k,node)+vold(k,node)
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
                  do k=1,kneigh
                     itri=neigh(k)+itietri(1,i)-1
!     
                     ipos=0
                     totdist(k)=0.d0
                     nodeedge(1,k)=0
                     nodeedge(2,k)=0
!     
                     do l=1,3
                        ll=4*l-3
                        dist=straight(ll,itri)*p(1)+
     &                       straight(ll+1,itri)*p(2)+
     &                       straight(ll+2,itri)*p(3)+
     &                       straight(ll+3,itri)
                        if(dist.gt.0.d0) then
                           totdist(k)=totdist(k)+dist
                           if(ipos.eq.0) then
                              nodeedge(1,k)=koncont(l,itri)
                              if(l.ne.3) then
                                 nodeedge(2,k)=koncont(l+1,itri)
                              else
                                 nodeedge(2,k)=koncont(1,itri)
                              endif
                           else
                              if((nodeedge(1,k).eq.koncont(l,itri)).or.
     &                         (nodeedge(2,k).eq.koncont(l,itri)))then
                                 nodeedge(1,k)=koncont(l,itri)
                                 nodeedge(2,k)=0
                              else
                                 if(l.ne.3) then
                                    nodeedge(1,k)=koncont(l+1,itri)
                                 else
                                    nodeedge(1,k)=koncont(1,itri)
                                 endif
                              endif
                           endif
                           ipos=ipos+1
                        endif
                     enddo
!     
                     if(totdist(k).le.0.d0) then
                        isol=k
                        exit
                     endif
                  enddo
!
!              if no independent face was found, a small
!              tolerance is applied
!
                  if(isol.eq.0) then
                     do k=1,kneigh
                        ipointer(k)=neigh(k)+itietri(1,i)-1
                     enddo
                     call dsort(totdist,ipointer,kneigh,kflag)
                     do k=1,kneigh
                        itri=ipointer(k)
                        dist=straight(1,itri)*cg(1,itri)+
     &                       straight(2,itri)*cg(2,itri)+
     &                       straight(3,itri)*cg(3,itri)+
     &                       straight(4,itri)
                        if(totdist(k).lt.1.d-3*dist) then
                           isol=k
                           exit
                        endif
                     enddo
                  endif
!     
                  if(isol.eq.0) then
                  else
!     
!                   plane spring
!     
                     ne=ne+1
                     ipkon(ne)=ifree
                     lakon(ne)='ESPRNGC '
                     ielmat(ne)=material
                     nelem=int(koncont(4,itri)/10.d0)
                     jface=koncont(4,itri)-10*nelem
!
!                    storing the face in ifcont1 and the
!                    element number in ifcont2
!
                     ifcont1(ne-ne0)=koncont(4,itri)
                     ifcont2(ne-ne0)=ne
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
                     do k=1,nnodelem
                        kon(ifree+k)=nodef(k)
                     enddo
                     ifree=ifree+nnodelem+1
                     kon(ifree)=node
                     ifree=ifree+1
                     write(lakon(ne)(8:8),'(i1)') nnodelem+1
c                   write(*,*) 'new elem',ne,(nodef(k),k=1,nnodelem),node
                  endif
!     
               enddo
            endif
         enddo
      enddo
!
!     sorting all used independent faces
!
      n=ne-ne0
      call isortii(ifcont1,ifcont2,n,kflag)
!
!     replace the faces by the number of times they were used in 
!     contact spring elements
!
      i=1
      loop: do
         ifaceref=ifcont1(i)
         isum=1
         j=i+1
         if(j.gt.ne-ne0) exit loop
         do
            if(ifcont1(j).eq.ifaceref) then
               isum=isum+1
               j=j+1
               if(j.gt.ne-ne0) exit loop
               cycle
            else
               do k=i,j-1
                  ifcont1(k)=isum
               enddo
               i=j
               exit
            endif
         enddo
      enddo loop
      do k=i,j-1
         ifcont1(k)=isum
      enddo
!
!     sorting in the original order
!
      call isortii(ifcont2,ifcont1,n,kflag)
!
!     storing the number of dependent nodes as last entry
!     in the topology
!
      do i=ne0+1,ne
         read(lakon(i)(8:8),'(i1)') nope
         kon(ipkon(i)+nope+1)=ifcont1(i-ne0)
c         write(*,*) 'gencontelem',i,ifcont1(i-ne0)
      enddo
!     
c      close(27)
!
      return
      end

