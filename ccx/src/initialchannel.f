!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2022 Guido Dhondt
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
!     Solve the Bresse equation for the turbulent stationary flow
!     in channels with a non-erosive bottom
!     
      subroutine initialchannel(itg,ieg,ntg,lakon,v,ipkon,kon,nflow,
     &     ikboun,nboun,prop,ielprop,ndirboun,nodeboun,xbounact,
     &     ielmat,ntmat_,shcon,nshcon,physcon,rhcon,nrhcon,ipobody,
     &     ibody,xbodyact,co,nbody,network,vold,set,istep,iit,mi,
     &     ineighe,ilboun,ttime,time,itreated,iponoel,inoel,istack,
     &     sfr,hfr,sba,hba,ndata,jumpup,jumpdo)
!
      implicit none
!
      character*1 mode
      character*8 lakon(*)
      character*81 set(*)
!      
      integer mi(*),itg(*),ieg(*),ntg,nflow,ipkon(*),kon(*),ikboun(*),
     &     nboun,ielprop(*),ndirboun(*),nodeboun(*),ielmat(mi(3),*),
     &     ntmat_,nshcon(*),nrhcon(*),ipobody(2,*),ibody(3,*),nbody,
     &     network,istep,iit,ineighe(*),ilboun(*),i,j,nelem,indexe,
     &     node1,node2,id,itreated(*),id1,id2,nup,index,iponoel(*),
     &     inoel(2,*),nmid,ndo,inv,nelemio,nelup,node,imat,neldo,
     &     istack(2,*),nstack,nel,ndata,jumpup(*),jumpdo(*)
!
      real*8 v(0:mi(2),*),prop(*),xbounact(*),shcon(0:3,ntmat_,*),
     &     physcon(*),rhcon(0:1,ntmat_,*),xbodyact(7,*),co(3,*),
     &     vold(0:mi(2),*),ttime,time,xflow,g(3),dg,temp,cp,dvi,r,
     &     rho,sfr(*),hfr(*),sba(*),hba(*),epsilon
!
      if(network.le.2) then
        write(*,*) '*ERROR: a network channel canot be used for'
        write(*,*) '        temperature calculations only'
        write(*,*)
        call exit(201)
      endif
!
      epsilon=1.d-5
!
!     identify the number of "real" channel neighbors (i.e. without
!     CHANNEL INOUT elements) per node: stored in field ineighe(*)
!
      do i=1,ntg
        ineighe(i)=0
      enddo
!
      do i=1,nflow
        nelem=ieg(i)
        if(lakon(nelem)(2:5).ne.'LICH') cycle
        if(lakon(nelem)(6:7).eq.'IO') cycle
        indexe=ipkon(nelem)
!
        node1=kon(indexe+1)
        call nident(itg,node1,ntg,id)
        ineighe(id)=ineighe(id)+1
!
        node2=kon(indexe+3)
        call nident(itg,node2,ntg,id)
        ineighe(id)=ineighe(id)+1
!
      enddo
!
!     setting all values to zero
!
      do i=1,ntg
        node=itg(i)
        do j=0,mi(2)
          v(j,i)=0.d0
        enddo
      enddo
!
!     applying the boundary conditions
!
      do j=1,nboun
        v(ndirboun(j),nodeboun(j))=xbounact(j)
      enddo
!     
!     determine the gravity vector
!     
      do j=1,3
        g(j)=0.d0
      enddo
      if(nbody.gt.0) then
        index=nelem
        do
          j=ipobody(1,index)
          if(j.eq.0) exit
          if(ibody(1,j).eq.2) then
            g(1)=g(1)+xbodyact(1,j)*xbodyact(2,j)
            g(2)=g(2)+xbodyact(1,j)*xbodyact(3,j)
            g(3)=g(3)+xbodyact(1,j)*xbodyact(4,j)
          endif
          index=ipobody(2,index)
          if(index.eq.0) exit
        enddo
      endif
      dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!
      nstack=0
!
!     major loop: looking for SLUICE GATE elements
!
      loop1: do i=1,nflow
        nelem=ieg(i)
        if(((lakon(nelem)(6:7).ne.'SG').and.
     &     (lakon(nelem)(6:7).ne.'WE')).or.
     &       (itreated(i).eq.1)) cycle
!
!       untreated SLUICE GATE element found
!
        indexe=ipkon(nelem)
        node1=kon(indexe+1)
        call nident(itg,node1,ntg,id1)
        node2=kon(indexe+3)
        call nident(itg,node2,ntg,id2)
!
!       looking for a SLUICE GATE element connected on one side to
!       a CHANNEL INOUT element (as only element)
!
        if((ineighe(id1).gt.1).and.(ineighe(id2).gt.1)) cycle
!
!       new branch found
!
!       determine the upstream node nup of the element
!
        if(ineighe(id1).eq.1) then
          nup=node1
          inv=1
        else
          nup=node2
          inv=-1
        endif
!
!       determine the upstream element nelup
!
        index=iponoel(nup)
        do
          if(index.eq.0) then
            write(*,*) '*ERROR: node',nup
            write(*,*) '        is only connected to one element'
            write(*,*) '        with number',nelem
            write(*,*)
            call exit(201)
          endif
          if(inoel(1,index).ne.nelem) then
            nelup=inoel(1,index)
            exit
          else
            index=inoel(2,index)
          endif
        enddo
!
!       define mode to be "forward"
!
        mode='F'
!
!       mass flow is taken from the IO element upstream of the sluice gate
!
        xflow=dabs(v(1,kon(ipkon(nelup)+2)))
!
!       loop over all elements in present branch
!
        loop2: do
!
!         if F (forward): nelup and nup known
!
        if(mode.eq.'F') then
          if(nelem.eq.0) then
              call nident(itg,nup,ntg,id)
              if(ineighe(id).eq.1) then
                write(*,*) '*INFO: branch finished'
!
!               IO-element: determine the mass flow
!
                index=iponoel(nup)
                do
                  if(inoel(1,index).ne.nelup) then
                    nelem=inoel(1,index)
                    if(nup.eq.kon(ipkon(nelem)+1)) then
                      v(1,kon(ipkon(nelem)+2))=xflow
                    else
                      v(1,kon(ipkon(nelem)+2))=-xflow
                    endif
                    cycle loop1
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) exit
                enddo
!
              elseif(ineighe(id).gt.2) then
                write(*,*) '*ERROR in initialchannel: branch split'
                write(*,*)
                call exit(201)
              endif
!     
!     one "true" element connected downstream
!     loop over all elements connected to nup
!     
              index=iponoel(nup)
              do
                if(inoel(1,index).ne.nelup) then
                  if(lakon(inoel(1,index))(6:7).ne.'IO') then
                    nelem=inoel(1,index)
                  else
!     
!     add flow
!     
                    nelemio=inoel(1,index)
                    if((lakon(nelup)(6:7).eq.'SG').or.
     &                 (lakon(nelup)(6:7).eq.'WE')) then
                      write(*,*)
     &                     '*ERROR in initialchannel: no IO element'
                      write(*,*) '       allowed immediately downstream'
                      write(*,*) '       of a SLUICE GATE element.'
                      write(*,*) '       faulty element:',nelemio
                      write(*,*) 
                      call exit(201)
                    endif
                    if(nup.eq.kon(ipkon(nelemio)+3)) then
                      xflow=xflow+v(1,kon(ipkon(nelemio)+2))
                    else
                      xflow=xflow-v(1,kon(ipkon(nelemio)+2))
                    endif
                  endif
                endif
                index=inoel(2,index)
                if(index.eq.0) exit
              enddo
            endif
!     
!     actual element = nelem
!     determining the middle and downstream node nmid and ndo
!     
!     if the actual flow is from kon(indexe+1) to kon(indexe+3) then
!     inv=1
!     if the actual flow is from kon(indexe+3) to kon(indexe+1) then
!     inv=-1
!     
            indexe=ipkon(nelem)
            nmid=kon(indexe+2)
            if(kon(indexe+1).eq.nup) then
              ndo=kon(indexe+3)
              inv=1
            else
              ndo=kon(indexe+1)
              inv=-1
            endif
          else
!     
!     B (backward): either nelem and ndo known or
!                   neldo and ndo known
!
            if(nelem.eq.0) then
              call nident(itg,ndo,ntg,id)
              if(ineighe(id).eq.1) then
!
!              IO-element: determine the mass flow
!     
                index=iponoel(ndo)
                do
                  if(inoel(1,index).ne.neldo) then
                    nelem=inoel(1,index)
                    if(ndo.eq.kon(ipkon(nelem)+3)) then
                      v(1,kon(ipkon(nelem)+2))=xflow
                    else
                      v(1,kon(ipkon(nelem)+2))=-xflow
                    endif
                    exit
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) exit
                enddo
!
                if(nstack.gt.0) then
                  mode='F'
                  nelup=istack(1,nstack)
                  nup=istack(2,nstack)
                  nelem=0
                  xflow=dabs(v(1,kon(ipkon(nelup)+2)))
                  nstack=nstack-1
                  cycle loop2
                else
                  write(*,*) '*INFO: branch finished'
                  cycle loop1
                endif
              elseif(ineighe(id).gt.2) then
                write(*,*) '*ERROR in initialchannel: branch split'
                write(*,*)
                call exit(201)
              endif
!     
!     if nelem is zero, neldo is known and nelem has to be
!     determined
!     
              index=iponoel(ndo)
              do
                if(inoel(1,index).ne.neldo) then
                  if(lakon(inoel(1,index))(6:7).ne.'IO') then
                    nelem=inoel(1,index)
                  else
!     
!     add flow
!     
                    nelemio=inoel(1,index)
                    if(lakon(nelup)(6:7).eq.'RE') then
                      write(*,*)
     &                     '*ERROR in initialchannel: no IO element'
                      write(*,*) '       allowed immediately upstream'
                      write(*,*) '       of a RESERVOIR element.'
                      write(*,*) '       faulty element:',nelemio
                      write(*,*) 
                      call exit(201)
                    endif
                    if(kon(ipkon(nelemio)+3).eq.ndo) then
                      xflow=xflow-v(1,kon(ipkon(nelemio)+2))
                    else
                      xflow=xflow+v(1,kon(ipkon(nelemio)+2))
                    endif
                  endif
                endif
                index=inoel(2,index)
                if(index.eq.0) exit
              enddo
            endif
!
!           determine nmid and nup
!
            indexe=ipkon(nelem)
            nmid=kon(indexe+2)
            if(kon(indexe+3).eq.ndo) then
              nup=kon(indexe+1)
              inv=1
            else
              nup=kon(indexe+3)
              inv=-1
            endif
!     
!           determining the upstream element (needed in istack
!           in certain cases
!     
            nelup=0
            index=iponoel(nup)
            do
              if(index.eq.0) exit
              if(inoel(1,index).ne.nelem) then
                if(lakon(inoel(1,index))(6:7).ne.'IO') then
                  nelup=inoel(1,index)
                  exit
                endif
              endif
              index=inoel(2,index)
            enddo
          endif
!
!         taking the temperature of the upstream node for the
!         material properties
!
          temp=v(0,nup)
          imat=ielmat(1,nelem)
!     
          call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &         dvi,rhcon,nrhcon,rho)
!     
!     treating the specific element type
!     
          call nident(ieg,nelem,nflow,nel)
          if((lakon(nelem)(6:7).eq.'SG').or.
     &       (lakon(nelem)(6:7).eq.'WE')) then
            call sluicegate(nelem,ielprop,prop,nup,nmid,ndo,co,g,dg,
     &           mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,
     &           nboun,mi,v,ipkon,kon,inv,epsilon,lakon)
          elseif(lakon(nelem)(6:7).eq.'  ') then
            call straightchannel(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &           dg,mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,
     &           nboun,mi,v,ipkon,kon,ndata,nel,sfr((nel-1)*ndata+1),
     &           hfr((nel-1)*ndata+1),sba((nel-1)*ndata+1),
     &           hba((nel-1)*ndata+1),jumpup,jumpdo,inv,epsilon)
          elseif(lakon(nelem)(6:7).eq.'RE') then
            call reservoir(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &           dg,mode,xflow,rho,dvi,nelup,mi,v,inv,epsilon,istack,
     &           nstack)
          elseif((lakon(nelem)(6:7).eq.'CO').or.
     &           (lakon(nelem)(6:7).eq.'EL').or.
     &           (lakon(nelem)(6:7).eq.'ST').or.
     &           (lakon(nelem)(6:7).eq.'DR')) then
            call contraction(nelem,ielprop,prop,nup,nmid,ndo,dg,
     &           mode,xflow,rho,nelup,neldo,istack,nstack,
     &           mi,v,inv,epsilon,co)
          else
            write(*,*) '*ERROR in initialchannel:'
            write(*,*) '       element of type ',lakon(nelem)
            write(*,*) '       is not known'
            write(*,*)
            call exit(201)
          endif
          itreated(nel)=1
        enddo loop2
      enddo loop1
!
      return
      end
      
