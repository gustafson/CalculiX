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
      subroutine gen3delem(kon,ipkon,lakon,ne,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,
     &  nboun_,nodeboun,ndirboun,xboun,iamboun,nam,
     &  inotr,trab,nk,nk_,iponoel,inoel,iponor,xnor,thicke,thickn,
     &  knor,istep,offset,t0,t1,ikforc,ilforc,rig,nforc,
     &  nforc_,nodeforc,ndirforc,xforc,iamforc,nelemload,sideload,
     &  nload,ithermal,ntrans,co,ixfree,ikfree,inoelfree,iponoelmax,
     &  iperturb,tinc,tper,tmin,tmax,ctrl,typeboun,nmethod,nset,set,
     &  istartset,iendset,ialset,prop,ielprop)
!
!     generates three-dimensional elements:
!         for isochoric elements
!         for plane stress
!         for plane strain
!         for plate and shell elements
!         for beam elements
!
      implicit none
!
      character*1 typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*81 set(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,ikmpc(*),
     &  ilmpc(*),kon(*),ipkon(*),ne,mpc,indexe,i,j,k,node,idof,
     &  id,mpcfreeold,ikboun(*),ilboun(*),nboun,nboun_,kflag,idummy,
     &  iterm(500),nterm,neigh(7,8),l,m,nodeboun(*),ndirboun(*),nk,
     &  nk_,index,iponoel(*),inoel(3,*),inoelfree,istep,nmpcold,
     &  ikforc(*),ilforc(*),nodeforc(2,*),ndirforc(*),iamforc(*),
     &  nelemload(*),nforc,nforc_,ithermal(2),nload,iamboun(*),
     &  ntrans,inotr(2,*),nam,iponoelmax,iperturb,numnod,itransaxial,
     &  rig(*),nmethod,nset,istartset(*),iendset(*),ialset(*),
     &  ielprop(*)
!
      integer iponor(2,*),knor(*),
     & ixfree,ikfree
!
      real*8 coefmpc(*),thicke(2,*),xnor(*),thickn(2,*),tinc,tper,tmin,
     &  tmax,offset(2,*),t0(*),t1(*),xforc(*),trab(7,*),co(3,*),
     &  xboun(*),pi,ctrl(26),prop(*)
!
      data neigh /1,9,2,17,5,12,4,2,9,1,18,6,10,3,
     &            3,11,4,19,7,10,2,4,11,3,20,8,12,1,
     &            5,13,6,17,1,16,8,6,13,5,18,2,14,7,
     &            7,15,8,19,3,14,6,8,15,7,20,4,16,5/
!
      pi=4.d0*datan(1.d0)
!
!     catalogueing the element per node relationship for shell/beam
!     elements and transferring the nodal thickness to the elements
!
!     inoelfree=1 means that there is at least one 1D or 2D element
!     in the structure. Otherwise inoelfree=0.
!
      if((istep.eq.1).and.(inoelfree.eq.1)) then
!
         itransaxial=0
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            if((lakon(i)(1:2).ne.'C3').and.(lakon(i)(1:1).ne.'D').and.
     &         (lakon(i)(1:1).ne.'G').and.(lakon(i)(1:1).ne.'E')) then
               if(lakon(i)(1:1).eq.'B') then
                  numnod=3
               elseif((lakon(i)(2:2).eq.'6').or.
     &                (lakon(i)(4:4).eq.'6')) then
                  numnod=6
               else
                  numnod=8
               endif
               indexe=ipkon(i)
               do j=1,numnod
                  node=kon(indexe+j)
                  iponoelmax=max(iponoelmax,node)
                  inoel(1,inoelfree)=i
                  inoel(2,inoelfree)=j
                  inoel(3,inoelfree)=iponoel(node)
                  iponoel(node)=inoelfree
                  inoelfree=inoelfree+1
                  if(lakon(i)(1:2).ne.'CA') then
                     if(thickn(1,node).gt.0.d0)
     &                    thicke(1,indexe+j)=thickn(1,node)
                     if(thickn(2,node).gt.0.d0)
     &                    thicke(2,indexe+j)=thickn(2,node)
                  endif
                  if(thicke(1,indexe+j).le.0.d0) then
                     if(lakon(i)(1:1).eq.'C') then
                        thicke(1,indexe+j)=1.d0
                     else
                        write(*,*)'*ERROR in gen3delem: first thickness'
                        write(*,*)'       in node ',j,' of element ',i
                        write(*,*)'       is zero'
                        stop
                     endif
                  endif
                  if((lakon(i)(1:1).eq.'B').and.
     &                 (thicke(2,indexe+j).le.0.d0)) then
                     write(*,*) '*ERROR in gen3delem: second thickness'
                     write(*,*)'       in node ',j,' of beam element ',i
                     write(*,*)'       is zero'
                     stop
                  endif
               enddo
            endif
         enddo
!
!        checking whether any rotational degrees of freedom are fixed
!        by SPC's, MPC's or loaded by bending moments or torques
!        in the end, rig(i)=0 if no rigid knot is defined in node i,
!        else rig(i)=the rotational node of the knot. The value -1 is
!        a dummy. 
!
         do i=1,nboun
            if(ndirboun(i).gt.3) rig(nodeboun(i))=-1
         enddo
         do i=1,nforc
            if(ndirforc(i).gt.3) rig(nodeforc(1,i))=-1
         enddo
         do i=1,nmpc
            index=ipompc(i)
            do
               if(index.eq.0) exit
               if(nodempc(2,index).gt.3) then
                  rig(nodempc(1,index))=-1
               endif
               index=nodempc(3,index)
            enddo
         enddo
!
!     calculating the normals in nodes belonging to shells/beams
!
         nmpcold=nmpc
!
         call gen3dnor(nk,nk_,co,iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,thicke,offset,iponor,xnor,knor,rig,iperturb,tinc,
     &     tper,tmin,tmax,ctrl,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,nboun_,
     &     nodeboun,ndirboun,xboun,iamboun,typeboun,nam,ntrans,inotr,
     &     trab,ikfree,ixfree,nmethod,ithermal)
!
      endif
!
      if(istep.eq.1) then
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
!
!        incompressible elements
!
            if(lakon(i)(1:7).eq.'C3D20RH') then
               indexe=ipkon(i)
               do j=1,8
                  node=kon(indexe+j)
                  mpc=0
                  do k=1,3
                     idof=7*(node-1)+k
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           if(labmpc(ilmpc(i))(1:9).eq.'ISOCHORIC') then
!
!                          existing mpc
!
                              mpc=ilmpc(id)
                              exit
                           else
                              cycle
                           endif
                        endif
                     endif
!
!                 new mpc
!
                     nmpc=nmpc+1
                     if(nmpc.gt.nmpc_) then
                        write(*,*) '*ERROR in gen3delem: increase nmpc_'
                        stop
                     endif
                     mpc=nmpc
                     do l=nmpc,id+2,-1
                        ikmpc(l)=ikmpc(l-1)
                        ilmpc(l)=ilmpc(l-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                  enddo
!
                  if(mpc.eq.0) then
                     write(*,*) '*WARNING in gen3delem: no free DOF in'
                     write(*,*) '         node ',node,' for isochoric'
                     write(*,*) '         MPC application'
                     cycle
                  endif
!
                  if(ipompc(mpc).eq.0) then
!
!                 new mpc
!
                     labmpc(mpc)(1:9)='ISOCHORIC'
                     ipompc(mpc)=mpcfree
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=k
                     mpcfree=nodempc(3,mpcfree)
!
                     do k=1,7
                        do l=1,3
                           nodempc(1,mpcfree)=kon(indexe+neigh(k,j))
                           nodempc(2,mpcfree)=l
                           mpcfree=nodempc(3,mpcfree)
                        enddo
                     enddo
!
!                 add nonhomogeneous term
!                  
                     nk=nk+1
                     if(nk.gt.nk_) then
                        write(*,*) '*ERROR in gen3delem: increase nk_'
                        stop
                     endif
                     mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                     nodempc(3,mpcfreeold)=0
                     idof=7*(nk-1)+k
                     call nident(ikboun,idof,nboun,id)
                     nboun=nboun+1
                     if(nboun.gt.nboun_) then
                        write(*,*)'*ERROR in gen3delem: increase nboun_'
                        stop
                     endif
                     nodeboun(nboun)=nk
                     ndirboun(nboun)=j
                     typeboun(nboun)='I'
                     do l=nboun,id+2,-1
                        ikboun(l)=ikboun(l-1)
                        ilboun(l)=ilboun(l-1)
                     enddo
                     ikboun(id+1)=idof
                     ilboun(id+1)=nboun
!
                  else
!
                     index=nodempc(3,ipompc(mpc))
                     nterm=0
                     do
                        if(index.eq.0) exit
                        nterm=nterm+1
                        if(nterm.gt.500) then
                           write(*,*) '*ERROR in gen3delem:'
                           write(*,*) '       increase nterm_'
                           stop
                        endif
                        iterm(nterm)=
     &                       7*(nodempc(1,index)-1)+nodempc(2,index)
                        index=nodempc(3,nodempc(3,nodempc(3,index)))
                     enddo
                     kflag=1
                     call isortii(iterm,idummy,nterm,kflag)
!
                     do k=1,7
                        do l=1,3
                           m=7*(kon(indexe+neigh(k,j))-1)+l
                           call nident(iterm,m,nterm,id)
                           if(id.ne.0) then
                              if(iterm(id).eq.m) then
                                 cycle
                              endif
                           endif
                           mpcfreeold=mpcfree
                           mpcfree=nodempc(3,mpcfree)
                           nodempc(3,mpcfree)=ipompc(mpc)
                           ipompc(mpc)=mpcfreeold
                           nodempc(1,mpcfreeold)=kon(indexe+neigh(k,j))
                           nodempc(2,mpcfreeold)=l
                        enddo
                     enddo
!
                  endif
               enddo
            elseif((lakon(i)(1:2).eq.'CP').or.
     &             (lakon(i)(1:1).eq.'S').or.
     &             (lakon(i)(1:2).eq.'CA')) then
!
               call gen3dfrom2d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &           thicke,offset,ntrans,inotr,trab,ikboun,ilboun,nboun,
     &           nboun_,nodeboun,ndirboun,xboun,iamboun,typeboun,ipompc,
     &           nodempc,coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,
     &           nk,nk_,co,rig,nmethod,iperturb,ithermal)
!            
            elseif(lakon(i)(1:1).eq.'B') then
               call gen3dfrom1d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &              thicke,ntrans,inotr,trab,nk,nk_,co,offset)
            endif
!
            if(lakon(i)(1:4).eq.'CPE6') then
               lakon(i)(1:7)='C3D15 E'
            elseif(lakon(i)(1:5).eq.'CPE8R') then
               lakon(i)(1:7)='C3D20RE'
            elseif(lakon(i)(1:4).eq.'CPE8') then
               lakon(i)(1:7)='C3D20 E'
            elseif(lakon(i)(1:4).eq.'CPS6') then
               lakon(i)(1:7)='C3D15 S'
            elseif(lakon(i)(1:5).eq.'CPS8R') then
               lakon(i)(1:7)='C3D20RS'
            elseif(lakon(i)(1:4).eq.'CPS8') then
               lakon(i)(1:7)='C3D20 S'
            elseif(lakon(i)(1:4).eq.'CAX6') then
               lakon(i)(1:7)='C3D15 A'
            elseif(lakon(i)(1:5).eq.'CAX8R') then
               lakon(i)(1:7)='C3D20RA'
            elseif(lakon(i)(1:4).eq.'CAX8') then
               lakon(i)(1:7)='C3D20 A'
            elseif(lakon(i)(1:2).eq.'S6') then
               lakon(i)(1:7)='C3D15 L'
            elseif(lakon(i)(1:3).eq.'S8R') then
               lakon(i)(1:7)='C3D20RL'
            elseif(lakon(i)(1:2).eq.'S8') then
               lakon(i)(1:7)='C3D20 L'
            elseif(lakon(i)(1:4).eq.'B32R') then
               lakon(i)(1:7)='C3D20RB'
            elseif(lakon(i)(1:1).eq.'B') then
               lakon(i)(1:7)='C3D20 B'
            endif
         enddo
!
!        filling the new KNOT MPC's (needs the coordinates
!        of the expanded nodes)
!
         if(inoelfree.ne.0) then
            call fillknotmpc(co,ipompc,nodempc,coefmpc,labmpc,
     &           nmpc,nmpcold)
            call gen3dprop(prop,ielprop,iponoel,inoel,iponoelmax,kon,
     &           ipkon,lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,
     &           nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,
     &           trab,nam,nk,nk_,co,nmethod,iperturb)
         endif
!     
      endif
!
!        generating MPC's to connect shells and beams with solid
!        elements
!
      if((inoelfree.ne.0).and.(istep.eq.1)) then
         call gen3dconnect(kon,ipkon,lakon,ne,iponoel,inoel,
     &     iponoelmax,rig,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,
     &     nmpc_,mpcfree,ikmpc,ilmpc,labmpc)
      endif
!
      if(inoelfree.ne.0) then
!
!           multiplying existing boundary conditions
!
         call gen3dboun(ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &     xboun,iamboun,typeboun,iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,
     &     co,nmethod,iperturb)
!
!        updating the nodal surfaces: establishing links between the user
!        defined nodes and the newly generated nodes (mid-nodes
!        for 2d elements, mean of corner nodes for 1d elements)
!
         if(istep.eq.1) then
            call gen3dsurf(iponoel,inoel,iponoelmax,kon,ipkon,
     &        lakon,ne,iponor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &        mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,
     &        nk_,co,nmethod,iperturb,nset,set,istartset,iendset,ialset)
         endif
!
!        updating the MPCs: establishing links between the user
!        defined nodes and the newly generated nodes (mid-nodes
!        for 2d elements, mean of corner nodes for 1d elements)
!
         if(istep.eq.1) then
            call gen3dmpc(ipompc,nodempc,coefmpc,nmpc,nmpc_,mpcfree,
     &        ikmpc,ilmpc,labmpc,iponoel,inoel,iponoelmax,kon,ipkon,
     &        lakon,ne,iponor,xnor,knor,rig)
         endif
!
!        updating the temperatures
!
         if(ithermal(1).gt.0) then
            call gen3dtemp(iponoel,inoel,iponoelmax,kon,ipkon,lakon,ne,
     &           iponor,xnor,knor,t0,t1,thicke,offset,rig,nk,nk_,co,
     &           istep,ithermal)
         endif
!
!        updating the concentrated loading
!
         call gen3dforc(ikforc,ilforc,nforc,nforc_,nodeforc,
     &     ndirforc,xforc,iamforc,ntrans,inotr,trab,rig,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,iponoel,inoel,
     &     iponoelmax,kon,ipkon,lakon,ne,iponor,xnor,knor,nam,nk,nk_,
     &     co,thicke)
      endif
!
      return
      end


