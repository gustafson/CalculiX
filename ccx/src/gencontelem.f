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
     &  lakon,cg,straight,ifree,
     &  koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,cs,
     &  elcon,istep,iinc,iit,ncmat_,ntmat_,ne0,
     &  vini,nmethod,mi,imastop,nslavnode,islavnode,islavsurf,
     &  itiefac,areaslav,iponoels,inoels,springarea,ikmpc,
     &  ilmpc,nmpc,ipompc,nodempc,coefmpc,set,nset,istartset,
     &  iendset,ialset,tietol,reltime)
!
!     generate contact elements for the slave contact nodes
!
      implicit none
!
      character*8 lakon(*)
c      character*18 cfile
      character*81 tieset(3,*),slavset,set(*),noset
!
      integer ntie,ifree,
     &  itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),ne,node,
     &  neigh(1),iflag,kneigh,i,j,k,l,isol,iset,idummy,
     &  itri,ll,kflag,n,nx(*),ny(*),istep,iinc,
     &  nz(*),nstart,ielmat(*),material,ifaceq(8,6),ifacet(6,4),
     &  ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,iit,
     &  nnodelem,nface,nope,nodef(8),ncmat_,ntmat_,index1,
     &  ne0,nmethod,mi(2),iteller,ifaces,jfaces,
     &  imastop(3,*), itriangle(100),ntriangle,ntriangle_,itriold,
     &  itrinew,id,nslavnode(*),islavnode(*),islavsurf(2,*),
     &  itiefac(2,*),iponoels(*),inoels(3,*),konl(20),nelems,m,
     &  mint2d,nopes,idof,index2,ikmpc(*),ilmpc(*),nmpc,
     &  nodempc(3,*),ipompc(*),ipos,nset,istartset(*),iendset(*),
     &  ialset(*)
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),cs(17,*),
     &  beta,c0,elcon(0:ncmat_,ntmat_,*),vini(0:mi(2),*),weight,
     &  areaslav(*),springarea(2,*),xl2(3,8),area,xi,et,shp2(7,8),
     &  xs2(3,2),xsj2(3),coefmpc(*),adjust,tietol(2,*),reltime
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
!     flag for shape functions
!
      data iflag /2/
!
      data iteller /0/
      save iteller
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
         kneigh=1
         slavset=tieset(2,i)
         material=int(tietol(2,i))
!
!        check whether an adjust node set has been defined
!        only checked at the start of the first step
!
c         if((istep.eq.1).and.(iinc.eq.1).and.(iit.le.0)) then
         if((istep.eq.1).and.(iit.lt.0)) then
            iset=0
            if(tieset(1,i)(1:1).ne.' ') then
               noset(1:80)=tieset(1,i)(1:80)
               noset(81:81)=' '
               ipos=index(noset,' ')
               noset(ipos:ipos)='N'
               do iset=1,nset
                  if(set(iset).eq.noset) exit
               enddo
               kflag=1
               call isortii(ialset(istartset(iset)),idummy,
     &            iendset(iset)-istartset(iset)+1,kflag)
            endif
         endif
!     
!        determine the area of the slave surfaces
!     
         do l = itiefac(1,i), itiefac(2,i)
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
!     
!     Decide on the max integration points number, just consider 2D situation 
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
                     xl2(j,m)=co(j,konl(ifaceq(m,jfaces)))+
     &                    vold(j,konl(ifaceq(m,jfaces)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopes
                  do j=1,3
                     xl2(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &                    vold(j,konl(ifacet(m,jfaces)))
                  enddo
               enddo
            else
               do m=1,nopes
                  do j=1,3
                     xl2(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &                    vold(j,konl(ifacew1(m,jfaces)))
                  enddo
               enddo
            endif
!     
!           calculating the area of the slave face
!
            area=0.d0
            do m=1,mint2d
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
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               else
                  call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               endif
               area=area+weight*dsqrt(xsj2(1)**2+xsj2(2)**2+
     &              xsj2(3)**2)
            enddo
            areaslav(l)=area
         enddo
!
!        search a master face for each slave node and generate a contact
!        spring element if successful
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
         kflag=2
         call dsort(x,nx,n,kflag)
         call dsort(y,ny,n,kflag)
         call dsort(z,nz,n,kflag)
!
         do j=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(j)
!
!                 calculating the area corresponding to the
!                 slave node; is made up of the area
!                 of the neighboring slave faces
!
            area=0.d0
            index1=iponoels(node)
            do
               if(index1.eq.0) exit
               area=area+areaslav(inoels(1,index1))/
     &              inoels(2,index1)
               index1=inoels(3,index1)
            enddo
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
!              check whether distance is larger than c0:
!              no element is generated
!
               if(isol.ne.0) then
                  dist=straight(13,itri)*p(1)+
     &                 straight(14,itri)*p(2)+
     &                 straight(15,itri)*p(3)+
     &                 straight(16,itri)
!
!                 check for an adjust parameter (only at the start
!                 of the first step)
!
c                  if((istep.eq.1).and.(iinc.eq.1).and.(iit.le.0)) then
                  if((istep.eq.1).and.(iit.lt.0)) then
                     if(iset.ne.0) then
!
!                       check whether node belongs to the adjust node
!                       set
!
                        call nident(ialset(istartset(iset)),node,
     &                      iendset(iset)-istartset(iset)+1,id)
                        if(id.gt.0) then
                           if(ialset(istartset(iset)+id-1).eq.node) then
                              do k=1,3
                                 co(k,node)=co(k,node)-
     &                                      dist*straight(12+k,itri)
                              enddo
                              dist=0.d0
                           endif
                        endif
                     elseif(dabs(tietol(1,i)).ge.2.d0) then
!
!                       adjust parameter
!
                        adjust=dabs(tietol(1,i))-2.d0
                        if(dist.le.adjust) then
                           do k=1,3
                              co(k,node)=co(k,node)-
     &                             dist*straight(12+k,itri)
                           enddo
                           dist=0.d0
                        endif
                     endif
                  endif
!                           
                  beta=elcon(1,1,material)
                  if(beta.gt.0.d0) then
                     c0=dlog(100.d0)/beta
                  else
                     if(dabs(area).gt.0.d0) then
                        c0=1.d-6*dsqrt(area)
                     else
                        c0=1.d-10
                     endif
                  endif
                  if(dist.gt.c0) then
                     isol=0
!
!                    adjusting the bodies at the start of the
!                    calculation such that they touch
!
c                  elseif((istep.eq.1).and.(iinc.eq.1).and.
c     &                   (iit.le.0).and.(dist.lt.0.d0).and.
c     &                   (nmethod.eq.1)) then
c                     do k=1,3
c                        vold(k,node)=vold(k,node)-
c     &                        dist*straight(12+k,itri)
c                        vini(k,node)=vold(k,node)
c                    enddo
                  endif
               endif
!     
               if(isol.ne.0) then
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
!                 storing the area corresponding to the slave node
!                 and the clearance if penetration takes place,
!                 i.e. dist <0 at the start of every step
!                
                  springarea(1,j)=area
                  if(iit.lt.0.d0) then
                     if(dist.lt.0.d0) then
                        springarea(2,j)=dist
                     else
                        springarea(2,j)=0.d0
                     endif
                  endif
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
                  kon(ifree)=j
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
         enddo
      enddo
!
!     
c      close(27)
!
      return
      end

