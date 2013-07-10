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
      subroutine checkspcmpc(lakon,ipkon,kon,ntie,tieset,nset,set,
     &  itiefac,islavsurf,islavnode,
     &  imastnode,nslavnode,nmastnode,
     &  slavnor,slavtan,islavact,
     &  nboun,ndirboun,nodeboun,xboun,
     &  nmpc,ipompc,nodempc,coefmpc,
     &  ikboun,ilboun,ikmpc,ilmpc,
     &  nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
     &  nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
     &  islavborder)
c>     check whether SPC's and MPC's in salve nodes are compatible
c>     with mortar contact    
c>     
c>
c>     author: Sitzmann,Saskia
c>
c>     islavmpc(2,j)=1  directional blocking
c>     islavmpc(2,j)=2  cyclic symmetry
c>
c>     imastmpc(2,j)=1  directional blocking
c>     imastmpc(2,j)=2  cyclic symmetry
c>     imastmpc(2,j)=3  spc with displacement
!
        implicit none
!
      logical nodeslavsurf,debug, incompatible,nogap,twod
!
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,mastset,set(*)
!
      integer ntie,i,j,k,l,nset,dir,dirind,dirdep,
     &  ifaces,nelems,jfaces,ifacem,nelemm,
     &  jfacem,indexe,nopes,nopem,ipkon(*),kon(*),id,
     &  ifaceq(8,6),ifacet(7,4),ifacew1(4,5),ifacew2(8,5),node,
     &  itiefac(2,*),islavsurf(2,*),islavnode(*),imastnode(*),
     &  nslavnode(ntie+1),nmastnode(ntie+1),ifacecount,islav,imast,
     &  ipos,index1,ifreenoelold,
     &  numbern,numberf,iface,kflag,nk,
     &  islavact(*),
     &  nboun,ndirboun(*),nodeboun(*),
     &  nmpc,ipompc(*),nodempc(3,*),index,
     &  ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &  nslavspc(2,*),islavspc(2,*),nsspc,nslavmpc(2,*),islavmpc(2,*),
     &  nsmpc,nmastspc(2,*),imastspc(2,*),nmspc,nmastmpc(2,*),
     &  imastmpc(2,*),nmmpc,islavborder(*),isspc,imspc,ismpc,immpc,
     &  nodem,nodes,ist,zs(3),dof,node2,nsl,nc
!      
      real*8  xboun(*),coefmpc(*),nn,n(3),fixed_disp,coefdep,
     &  slavnor(3,*),slavtan(6,*),v(3),sp
!
       debug=.false.
!
       if(nsspc.gt.0 .or. nsmpc.gt.0)then
          do i=1,ntie
             twod=.false.
             if(tieset(1,i)(81:81).ne.'C') cycle
             nsl= nslavnode(i+1)-nslavnode(i)+1
             nc=0
             do l=nslavnode(i)+1,nslavnode(i+1)
                nc=nc+(nslavspc(2,l)-nslavspc(1,l))
                nc=nc+(nslavmpc(2,l)-nslavmpc(1,l))
             enddo
             if(nc.gt.0.88*nsl)then
                twod=.true.
             endif
             write(*,*)'csm:tie',i,'nc',nc,'nsl',nsl                         
             do l=nslavnode(i)+1,nslavnode(i+1)
                node=islavnode(l)
                debug=.false.
                if(node.eq.18)debug=.true.
                do k=1,3
                   n(k)=slavnor(k,l)
                enddo
                nn=n(1)*n(1)+n(2)*n(2)+n(3)*n(3)
                incompatible=.false.
                if(islavact(l).eq.-1) then
                  nogap=.true.
                else
                  nogap=.false.
                endif
                if(((nslavspc(2,l)-(nslavspc(1,l))).gt.0 .or.
     &               (nslavmpc(2,l)-(nslavmpc(1,l))).gt.0).and.
     &               islavborder(l).lt.0)then
                  if(debug)then
                   write(*,*) 'checkspcmpc: only SPCs and MPCs on',
     &                  '  border slave nodes allowed'
                   write(*,*) 'node',node, 'has no LM'
                  endif
                   incompatible=.true.
                endif
                do j=nslavspc(1,l)+1,nslavspc(2,l)
                   v(1)=0.0
                   v(2)=0.0
                   v(3)=0.0
                   ist=islavspc(1,j)
                   dir=ndirboun(ist)
                   v(dir)=1.0
                   sp=v(1)*n(1)+v(2)*n(2)+v(3)*n(3)
                   islavspc(2,j)=1
!     tolerance of around 1 degrees
                   if(.not.(sp.lt.-0.08 .or. sp.gt.0.08) )then
                    do k=1,3
                     n(k)=n(k)-sp*v(k)
                    enddo
                   endif
                   if(sp.lt.-0.08 .or. sp.gt.0.08 )then
                    if(debug)then
                      write(*,*) 'checkspcmpc: normal direction cannot', 
     &                     'be blocked'
                      write(*,*) 'node',node, 'has no LM'
                    endif
                      incompatible=.true.
                      islavspc(2,j)=-2
                   endif           
                   if(xboun(ist).gt.1.e-16 .or.
     &                  xboun(ist).lt.-1.e-16)then
                     if(debug)then
                      write(*,*) 'checkspcmpc: no displacements allowed'
                      write(*,*) 'so far:node',node, 'has no LM'
                     endif
                      incompatible=.true.
                      islavspc(2,j)=-2
                   endif
                   if((islavspc(2,j).eq.1).and.(.not.twod))then
                      incompatible=.true.
                   endif
                enddo                
! check for directional blocking                
                do j=nslavmpc(1,l)+1,nslavmpc(2,l)
                   v(1)=0.0
                   v(2)=0.0
                   v(3)=0.0
                   islavmpc(2,j)=1
                   ist=islavmpc(1,j)
                   dirdep=nodempc(2,ist)
                   coefdep=coefmpc(ist)
                   v(dirdep)=coefdep
                   index=nodempc(3,ist)
                   fixed_disp=0.d0
                   if(index.ne.0) then
                      do
                         if(nodempc(1,index).eq.node)then
                            dirind=nodempc(2,index)
                            v(dirind)=coefmpc(index)
                         else
                           dof=8*(nodempc(1,index)-1)+nodempc(2,index)
                           call nident(ikboun, dof, 
     &                        nboun, id)
                           if(id.gt.0 .and. ikboun(id).eq.dof)then
                            if(xboun(ilboun(id)).gt.1.e-16 .or.
     &                       xboun(ilboun(id)).lt.-1.e-16)then
                       if(debug)then
                      write(*,*) 'checkspcmpc: no displacements allowed'
                      write(*,*) 'so far:node',node, 'has no LM'
                       endif
                             incompatible=.true.
                             islavmpc(2,j)=-2
                             v(dirdep)=0.0
                            endif
                           else
                            v(dirdep)=0.0
                           endif                           
                         endif
                         index=nodempc(3,index)
                         if(index.eq.0) exit              
                      enddo
                      sp=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
                      if(sp>1.e-16)then
                       do k=1,3
                         v(k)=v(k)/sp 
                       enddo
                       sp=v(1)*n(1)+v(2)*n(2)+v(3)*n(3)
                       if(.not.(sp.lt.-0.08 .or. sp.gt.0.08) )then
                        do k=1,3
                         n(k)=n(k)-sp*v(k)
                        enddo
!     tolerance of around 2 degrees
                       else
                         if(debug)then
                         write(*,*) 'checkspcmpc: normal direction can', 
     &                        'not be blocked'
                         write(*,*) 'node',node, 'has no LM'
                         endif
                         incompatible=.true.
                           islavmpc(2,j)=-2
                       endif
                      else
                       islavmpc(2,j)=-2

                      endif 
                   endif
                   if((islavmpc(2,j).eq.1).and.(.not.twod))then
                   incompatible=.true.
                   endif                   
                enddo
! check for zyclic symmetry 
                zs(1)=0
                zs(2)=0
                zs(3)=0              
                do j=nslavmpc(1,l)+1,nslavmpc(2,l)
                 if(islavmpc(2,j).ne.1)then
                   ist=islavmpc(1,j)
                   dirdep=nodempc(2,ist)
                   coefdep=coefmpc(ist)
                   index=nodempc(3,ist)
                   fixed_disp=0.d0
                   zs(dirdep)=nodempc(1,index)
                   call nident(islavnode(nslavnode(i)+1), zs(dirdep), 
     &              nslavnode(i+1)-nslavnode(i), id)
                   if(id>0)then
                    if(islavnode(nslavnode(i)+id).ne.zs(dirdep))then
                        if(debug)then
                         write(*,*) 'checkspcmpc: only zyclic ', 
     &                        'symmetry on slave nodes supported'
                         write(*,*) 'node',node, 'has no LM'
                         endif
                           incompatible=.true.
                           islavmpc(2,j)=-2
                    endif
                   else
                         if(debug)then
                          write(*,*) 'checkspcmpc: only zyclic ', 
     &                        'symmetry on slave nodes supported'
                         write(*,*) 'node',node, 'has no LM'
                         endif
                           incompatible=.true.
                           islavmpc(2,j)=-2     
                   endif              
c check whether ind node is slave node at border
                   if(index.ne.0) then
                      do
                         if(nodempc(1,index).ne.zs(dirdep))then
                  if(debug)then
                         write(*,*) 'checkspcmpc: only zyclic symmetry', 
     &                        'and directional blocking supported'
                         write(*,*) 'node',node, 'has no LM'
                  endif
                           incompatible=.true.
                           islavmpc(2,j)=-2
                          if(debug)then
                          write(*,*) nodempc(1,index),zs(dirdep)
                          endif
                         endif
                         index=nodempc(3,index)
                         if(index.eq.0) exit              
                      enddo
                      if(.not.incompatible.and.islavmpc(2,j).ne.1)then
                       islavmpc(2,j)=2
                      endif
                   endif
                 endif 
                 if((islavmpc(2,j).eq.2) .and. (.not.twod))then
!                   make cyclic symmetry master node NoLM-node, too
                    islavact(l)=-2
                    islavact(nslavnode(i)+id)=-2                
                  endif
                enddo
!
                if(debug)write(*,*) zs(1),zs(2),zs(3)
!                
                if(.not.(zs(1).eq.zs(2) .or. zs(2).eq.zs(3)
     &                   .or. zs(1).eq.zs(3)) )then
                  if(debug)then
                    write(*,*) 'checkspcmpc: only zyclic symmetry', 
     &                        'and directional blocking supported'
                    write(*,*) 'node',node, 'has no LM'
                    write(*,*) zs(1),zs(2),zs(3)
                  endif
                    incompatible=.true.
                endif 
!                                
                if(incompatible)then
                   islavact(l)=-2
                else
!                adapt slavnormal due to spc's
                 do k=1,3
                   slavnor(k,l)=n(k)
                 enddo                  
                endif
                if(debug)then
                write(*,*)'checkspcmpc:node',node,'act',islavact(l)
                do j=nslavspc(1,l)+1,nslavspc(2,l)
                 write(*,*) 'spc',islavspc(1,j),'type',islavspc(2,j)
                enddo
                do j=nslavmpc(1,l)+1,nslavmpc(2,l)
                 write(*,*) 'mpc',islavmpc(1,j),'type',islavmpc(2,j)
                enddo
                write(*,*) '***************************************'
                endif
             enddo
          enddo
!          
       endif
      if(nmspc.gt.0 .or. nmmpc.gt.0)then
          do i=1,ntie
             if(tieset(1,i)(81:81).ne.'C') cycle
             do l=nmastnode(i)+1,nmastnode(i+1)   
                node=imastnode(l)             
! check for directional blocking and "hidden" SPCs            
                do j=nmastmpc(1,l)+1,nmastmpc(2,l)
                   v(1)=0.0
                   v(2)=0.0
                   v(3)=0.0
                   imastmpc(2,j)=1
                   ist=imastmpc(1,j)
                   dirdep=nodempc(2,ist)
                   coefdep=coefmpc(ist)
                   v(dirdep)=coefdep
                   index=nodempc(3,ist)
                   fixed_disp=0.d0
                   k=1
                   if(index.ne.0) then
                      do
                         k=k+1
                         if(nodempc(1,index).eq.node)then
                            dirind=nodempc(2,index)
                            v(dirind)=coefmpc(index)
                         else
                           dof=8*(nodempc(1,index)-1)+nodempc(2,index)
                           call nident(ikboun, dof, 
     &                        nboun, id)
                           if(id.gt.0 .and. ikboun(id).eq.dof)then
                            if(xboun(ilboun(id)).gt.1.e-16 .or.
     &                       xboun(ilboun(id)).lt.-1.e-16)then
                      write(*,*) 'checkspcmpc:nodem',node
                      write(*,*) 'displacement',xboun(ilboun(id))
                             incompatible=.true.
                             if(imastmpc(2,j).ne.-2)imastmpc(2,j)=3
                            endif
                           else
                            imastmpc(2,j)=-2
                           endif                           
                         endif
                         index=nodempc(3,index)
                         if(index.eq.0) exit              
                      enddo
                   endif
!                   
                enddo
! check for zyclic symmetry 
                   zs(1)=0
                   zs(2)=0
                   zs(3)=0              
                do j=nmastmpc(1,l)+1,nmastmpc(2,l)
                 if(imastmpc(2,j).ne.1.and.imastmpc(2,j).ne.3)then
                   ist=imastmpc(1,j)
                   dirdep=nodempc(2,ist)
                   coefdep=coefmpc(ist)
                   index=nodempc(3,ist)
                   fixed_disp=0.d0
                   zs(dirdep)=nodempc(1,index)
                   call nident(imastnode(nmastnode(i)+1), zs(dirdep), 
     &              nmastnode(i+1)-nmastnode(i), id)
                   if(id>0)then
                    if(imastnode(nmastnode(i)+id).ne.zs(dirdep))then
c                         write(*,*) 'checkspcmpc: only zyclic ', 
c     &                        'symmetry on master nodes supported'
c                         write(*,*) 'nodem',node
                    endif
                   else
c                          write(*,*) 'checkspcmpc: only zyclic ', 
c     &                        'symmetry on master nodes supported'
c                         write(*,*) 'nodem',node   
                   endif              
                   if(index.ne.0) then
                      do
                         if(nodempc(1,index).ne.zs(dirdep).and. 
     &                    nodempc(1,index).ne.node)then
c                         write(*,*) 'checkspcmpc: only zyclic symmetry', 
c     &                        'and directional blocking supported'
c                         write(*,*) 'nodem',node
c                          write(*,*) nodempc(1,index),zs(dirdep)
                         endif
                         index=nodempc(3,index)
                         if(index.eq.0) exit              
                      enddo
                      if(.not.incompatible.and.imastmpc(2,j).ne.1)then
                       imastmpc(2,j)=2
                      endif
                   endif
                 endif  
                enddo
!
                if(.not.(zs(1).eq.zs(2) .or. zs(2).eq.zs(3)
     &                   .or. zs(1).eq.zs(3)) )then
c                    write(*,*) 'checkspcmpc: only zyclic symmetry', 
c     &                        'and directional blocking supported'
c                    write(*,*) 'nodem',node
c                    write(*,*) zs(1),zs(2),zs(3)
                endif 
!                
                if(debug)then
                write(*,*)'checkspcmpc:nodem',node
c                do j=nmastspc(1,l)+1,nmastspc(2,l)
c                 write(*,*) 'spc',imastspc(1,j),'type',imastspc(2,j)
c                enddo
                do j=nmastmpc(1,l)+1,nmastmpc(2,l)
                 write(*,*) 'mpc',imastmpc(1,j),'type',imastmpc(2,j)
                enddo
                write(*,*) '***************************************'
                endif
             enddo
          enddo
!          
       endif 
!
       if(ntie.gt.1)then
          do i=1,ntie
             if(tieset(1,i)(81:81).ne.'C') cycle
             do l=nslavnode(i)+1,nslavnode(i+1)
                node=islavnode(l)
                if(islavact(l).gt.-1)then
                 do j=1,ntie
                  if(j.ne.i)then
                   if(tieset(1,j)(81:81).ne.'C') cycle
                   call nident(islavnode(nslavnode(j)+1), node, 
     &              nslavnode(j+1)-nslavnode(j), id)
                   if(id>0)then
                    if(islavnode(nslavnode(j)+id).eq.node)then
                     islavact(l)=-2
                     if(debug)then
                  write(*,*)'checkspcmpc: node',node,'tie1s',i,'tie2s',j
                  write(*,*)'in more than one contact tie and set NoLM!'
                     endif
                    endif
                   endif                   
                   call nident(imastnode(nmastnode(j)+1), node, 
     &              nmastnode(j+1)-nmastnode(j), id)
                   if(id>0)then
                    if(imastnode(nmastnode(j)+id).eq.node)then
                     islavact(l)=-2
                     if(debug)then
                  write(*,*)'checkspcmpc: node',node,'tie1s',i,'tie2m',j
                  write(*,*)'in more than one contact tie and set NoLM!'
                     endif
                    endif
                   endif                   
                  endif
                 enddo
                endif
             enddo
          enddo
!
       endif
!             
       return
       end
