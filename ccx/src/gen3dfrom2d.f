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
      subroutine gen3dfrom2d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &  thicke,offset,ntrans,inotr,trab,ikboun,ilboun,nboun,nboun_,
     &  nodeboun,ndirboun,xboun,iamboun,typeboun,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,nk,nk_,co,rig,nmethod,
     &  iperturb,ithermal)
!
!     expands 2d element i into a 3d element
!
!     generates additional MPC's for plane stress, plane strain and
!     axisymmetric elements
!
      implicit none
!
      logical axial
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*)
!
      integer kon(*),ipkon(*),ne,iponor(2,*),knor(*),ntrans,inotr(2,*),
     &  ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),ndirboun(*),
     &  iamboun(*),nam,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikmpc(*),ilmpc(*),nk,nk_,i,rig(*),nmethod,iperturb,ishift
!
      integer indexe,j,nodel(8),indexx,indexk,k,nedge,nodes(3,8),
     &  iamplitude,l,newnode,idir,idof,id,m,mpcfreenew,node,ithermal(2),
     &  jmin,jmax
! 
      real*8 xnor(*),thicke(2,*),offset(2,*),trab(7,*),xboun(*),
     &  coefmpc(*),co(3,*)
!
      real*8 thicks(8),xnors(3,8),dc,ds,val,x,y
!
!              check for axial elements  
!
      if(lakon(i)(1:2).eq.'CA') then
         axial=.true.
      else
         axial=.false.
      endif
!
      indexe=ipkon(i)
!
!           localizing the nodes, thicknesses and normals for the
!           2-D element
!            
      if((lakon(i)(2:2).eq.'6').or.
     &     (lakon(i)(4:4).eq.'6')) then
         nedge=3
         ishift=15
      else
         nedge=4
         ishift=20
      endif
!
      do j=1,2*nedge
         nodel(j)=kon(indexe+j)
         kon(indexe+ishift+j)=nodel(j)
         indexk=iponor(2,indexe+j)
         thicks(j)=thicke(1,indexe+j)
         do k=1,3
            nodes(k,j)=knor(indexk+k)
         enddo
      enddo
!
!           generating the 3-D element topology for shell and plane
!           stress/strain elements
!
      if(lakon(i)(1:2).ne.'CA') then
         do j=1,2*nedge
            indexx=iponor(1,indexe+j)
            do k=1,3
               xnors(k,j)=xnor(indexx+k)
            enddo
            if(ntrans.gt.0) then
               do k=1,3
                  inotr(1,nodes(k,j))=inotr(1,nodel(j))
               enddo
            endif
         enddo
!
         do k=1,nedge
            kon(indexe+k)=nodes(1,k)
!
            do j=1,3
               co(j,nodes(1,k))=co(j,nodel(k))
     &              -thicks(k)*xnors(j,k)*(.5d0+offset(1,i))
            enddo
         enddo
         do k=1,nedge
            kon(indexe+nedge+k)=nodes(3,k)
            do j=1,3
               co(j,nodes(3,k))=co(j,nodel(k))
     &              +thicks(k)*xnors(j,k)*(.5d0-offset(1,i))
            enddo
         enddo
         do k=nedge+1,2*nedge
            kon(indexe+nedge+k)=nodes(1,k)
            do j=1,3
               co(j,nodes(1,k))=co(j,nodel(k))
     &              -thicks(k)*xnors(j,k)*(.5d0+offset(1,i))
            enddo
         enddo
         do k=nedge+1,2*nedge
            kon(indexe+2*nedge+k)=nodes(3,k)
            do j=1,3
               co(j,nodes(3,k))=co(j,nodel(k))
     &              +thicks(k)*xnors(j,k)*(.5d0-offset(1,i))
            enddo
         enddo
         do k=1,nedge
            kon(indexe+4*nedge+k)=nodes(2,k)
            do j=1,3
               co(j,nodes(2,k))=co(j,nodel(k))
     &              -thicks(k)*xnors(j,k)*offset(1,i)
            enddo
         enddo
      else
!
!           generating the 3-D element topology for axisymmetric elements
!
         dc=dcos(thicks(1)/2.d0)
         ds=dsin(thicks(1)/2.d0)
         do j=1,nedge
            indexk=iponor(2,indexe+j)
            x=co(1,nodel(j))
            y=co(2,nodel(j))
!
            node=knor(indexk+1)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=-x*ds
            kon(indexe+j)=node
!
            node=knor(indexk+2)
            co(1,node)=x
            co(2,node)=y
            co(3,node)=0.d0
            kon(indexe+4*nedge+j)=node
!
            node=knor(indexk+3)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=x*ds
            kon(indexe+nedge+j)=node
         enddo
!
         do j=nedge+1,2*nedge
            indexk=iponor(2,indexe+j)
            x=co(1,nodel(j))
            y=co(2,nodel(j))
!
            node=knor(indexk+1)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=-x*ds
            kon(indexe+nedge+j)=node
!
            node=knor(indexk+3)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=x*ds
            kon(indexe+2*nedge+j)=node
         enddo
      endif
!
!           additional SPC's due to plane strain/plane stress/axisymmetric
!           conditions
!
      do j=1,2*nedge
         if(lakon(i)(1:1).ne.'S') then
!
!                    fixing the middle plane
!
            if(rig(nodel(j)).gt.0) cycle
!
            if(ithermal(2).ne.2) then
               val=0.d0
               k=3
               if(nam.gt.0) iamplitude=0
               type='M'
               call bounadd(nodes(2,j),k,k,val,nodeboun,
     &              ndirboun,xboun,nboun,nboun_,iamboun,
     &              iamplitude,nam,ipompc,nodempc,coefmpc,
     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &              labmpc,type,typeboun,nmethod,iperturb)
            endif
!     
!           specifying that the side planes do the same
!           as the middle plane (in all directions for
!           plane strain and axisymmetric elements, in the
!           plane for plane stress elements)
!
            if(ithermal(2).le.1) then
               jmin=1
               jmax=3
            elseif(ithermal(2).eq.2) then
               jmin=0
               jmax=0
            else
               jmin=0
               jmax=3
            endif
!
            do l=1,3,2
               newnode=nodes(l,j)
               do idir=jmin,jmax
                  if((idir.eq.3).and.(lakon(i)(1:3).eq.'CPS'))
     &                 cycle
                  idof=7*(newnode-1)+idir
                  call nident(ikmpc,idof,nmpc,id)
                  if((id.le.0).or.(ikmpc(id).ne.idof)) then
                     nmpc=nmpc+1
                     if(nmpc.gt.nmpc_) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase nmpc_'
                        stop
                     endif
                     labmpc(nmpc)='                    '
                     ipompc(nmpc)=mpcfree
                     do m=nmpc,id+2,-1
                        ikmpc(m)=ikmpc(m-1)
                        ilmpc(m)=ilmpc(m-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                     nodempc(1,mpcfree)=newnode
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase nmpc_'
                        stop
                     endif
                     nodempc(1,mpcfree)=nodes(2,j)
                     if((lakon(i)(2:2).eq.'A').and.(idir.eq.3))
     &                    then
                        nodempc(2,mpcfree)=1
                     else
                        nodempc(2,mpcfree)=idir
                     endif
                     if(lakon(i)(2:2).eq.'A') then
                        if(idir.eq.1) then
                           coefmpc(mpcfree)=-dc
                        elseif(idir.eq.3) then
                           if(l.eq.1) then
                              coefmpc(mpcfree)=ds
                           else
                              coefmpc(mpcfree)=-ds
                           endif
                        else
                           coefmpc(mpcfree)=-1.d0
                        endif
                     else
                        coefmpc(mpcfree)=-1.d0
                     endif
                     mpcfreenew=nodempc(3,mpcfree)
                     if(mpcfreenew.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase nmpc_'
                        stop
                     endif
                     nodempc(3,mpcfree)=0
                     mpcfree=mpcfreenew
                  endif
               enddo
            enddo
         endif
      enddo
!            
      return
      end


