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
      subroutine rhs(co,nk,kon,ipkon,lakon,ne,
     &  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &  fext,nactdof,neq,nmethod,
     &  ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,alcon,
     &  nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,
     &  iprestr,vold,iperturb,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,npmat_,ttime,time,istep,iinc,dtime,
     &  physcon,ibody)
!
!     filling the right hand side load vector b
!
!     b contains the contributions due to mechanical forces only
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer kon(*),ipompc(*),nodempc(3,*),ipobody(2,*),nbody,
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),ikmpc(*),
     &  ilmpc(*),nactdof(0:3,*),konl(20),nelcon(2,*),ibody(3,*),
     &  nrhcon(*),nalcon(2,*),ielmat(*),ielorien(*),ipkon(*)
!
      integer nk,ne,nmpc,nforc,nload,neq,nmethod,nom,
     &  ithermal,iprestr,iperturb,i,j,k,idist,jj,
     &  id,ist,index,jdof1,jdof,node1,ntmat_,indexe,nope,norien,
     &  iexpl,idof1,iinc,istep,icalccg
!
      integer nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),coefmpc(*),xforc(*),xload(2,*),p1(3,2),
     &  p2(3,2),fext(*),bodyf(3),elcon(0:21,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &  t0(*),t1(*),vold(0:3,*),ff(60),time,ttime,dtime
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*)
!
      real*8 om(2),physcon(3)
!
      icalccg=0
!
      do i=1,neq
         fext(i)=0.d0
      enddo
!
!        distributed forces (body forces or thermal loads or
!        residual stresses or distributed face loads)
!
      if((nbody.ne.0).or.(ithermal.ne.0).or.
     &     (iprestr.ne.0).or.(nload.ne.0)) then
         idist=1
      else
         idist=0
      endif
!
      if((ithermal.le.1).or.(ithermal.eq.3)) then
!
!     mechanical analysis: loop over all elements
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j) 
         enddo
!
!       assigning centrifugal forces
!
         if(nbody.gt.0) then
c            om=0.d0
            nom=0
            om(1)=0.d0
            om(2)=0.d0
            bodyf(1)=0.d0
            bodyf(2)=0.d0
            bodyf(3)=0.d0
!     
            index=i
            do
               j=ipobody(1,index)
               if(j.eq.0) exit
               if(ibody(1,j).eq.1) then
                  nom=nom+1
                  if(nom.gt.2) then
                     write(*,*)'*ERROR in rhs: no more than two centri-'
                     write(*,*)'       fugal loading cards allowed'
                     stop
                  endif
                  om(nom)=xbody(1,j)
                  p1(1,nom)=xbody(2,j)
                  p1(2,nom)=xbody(3,j)
                  p1(3,nom)=xbody(4,j)
                  p2(1,nom)=xbody(5,j)
                  p2(2,nom)=xbody(6,j)
                  p2(3,nom)=xbody(7,j)
!     
!     assigning gravity forces
!     
               elseif(ibody(1,j).eq.2) then
                  bodyf(1)=bodyf(1)+xbody(1,j)*xbody(2,j)
                  bodyf(2)=bodyf(2)+xbody(1,j)*xbody(3,j)
                  bodyf(3)=bodyf(3)+xbody(1,j)*xbody(4,j)
!     
!     assigning newton gravity forces
!     
               elseif(ibody(1,j).eq.3) then
                  call newton(icalccg,ne,ipkon,lakon,kon,t0,co,rhcon,
     &               nrhcon,ntmat_,physcon,i,cgr,bodyf,ielmat,ithermal,
     &               vold)
               endif
               index=ipobody(2,index)
               if(index.eq.0) exit
            enddo
         endif
!     
         call e_c3d_rhs(co,nk,konl,lakon(i),p1,p2,om,bodyf,nbody,
     &        ff,i,nmethod,rhcon,ielmat,ntmat_,vold,iperturb,
     &        nelemload,sideload,xload,nload,idist,ttime,time,istep,
     &        iinc,dtime)
!
         do jj=1,3*nope
!
            j=(jj-1)/3+1
            k=jj-3*(j-1)
!
            node1=kon(indexe+j)
            jdof1=nactdof(k,node1)
!
!            distributed forces
!
            if(idist.ne.0) then
               if(jdof1.eq.0) then
                  if(nmpc.ne.0) then
                     idof1=(node1-1)*7+k
                     call nident(ikmpc,idof1,nmpc,id)
                     if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                        id=ilmpc(id)
                        ist=ipompc(id)
                        index=nodempc(3,ist)
                        do
                           jdof1=nactdof(nodempc(2,index),
     &                          nodempc(1,index))
                           if(jdof1.ne.0) then
                              fext(jdof1)=fext(jdof1)
     &                             -coefmpc(index)*ff(jj)/coefmpc(ist)
                           endif
                           index=nodempc(3,index)
                           if(index.eq.0) exit
                        enddo
                     endif
                  endif
                  cycle
               endif
               fext(jdof1)=fext(jdof1)+ff(jj)
            endif
!
         enddo
      enddo
!
      else
!
!     thermal analysis: loop over all elements
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j) 
         enddo
!     
         call e_c3d_rhs_th(co,nk,konl,lakon(i),
     &        ff,i,nmethod,t0,t1,vold,nelemload,
     &        sideload,xload,nload,idist,dtime,
     &        ttime,time,istep,iinc)
!
         do jj=1,nope
!
            j=jj
!
            node1=kon(indexe+j)
            jdof1=nactdof(0,node1)
!
!            distributed forces
!
            if(idist.ne.0) then
               if(jdof1.eq.0) then
                  if(nmpc.ne.0) then
                     idof1=(node1-1)*7
                     call nident(ikmpc,idof1,nmpc,id)
                     if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                        id=ilmpc(id)
                        ist=ipompc(id)
                        index=nodempc(3,ist)
                        do
                           jdof1=nactdof(nodempc(2,index),
     &                          nodempc(1,index))
                           if(jdof1.ne.0) then
                              fext(jdof1)=fext(jdof1)
     &                             -coefmpc(index)*ff(jj)/coefmpc(ist)
                           endif
                           index=nodempc(3,index)
                           if(index.eq.0) exit
                        enddo
                     endif
                  endif
                  cycle
               endif
               fext(jdof1)=fext(jdof1)+ff(jj)
            endif
!
         enddo
      enddo
!
      endif
!
!        point forces
!      
      do i=1,nforc
         if(ndirforc(i).gt.3) cycle
         jdof=nactdof(ndirforc(i),nodeforc(1,i))
         if(jdof.ne.0) then
            fext(jdof)=fext(jdof)+xforc(i)
         endif
      enddo
c      write(*,*) 'rhs '
c      write(*,'(6(1x,e11.4))') (fext(i),i=1,neq)
!
      return
      end
