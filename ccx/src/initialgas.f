!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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

!     calculate the initial conditions for the gas!     
!     
      subroutine initialgas(itg,ieg,ntm,ntg,ac,bc,lakon,v,
     &     ipkon,kon,nflow,ikboun,nboun,prop,ielprop,
     &     nactdog,ndirboun,nodeboun,xbounact,
     &     ielmat,ntmat_,shcon,nshcon,physcon,ipiv,nteq,
     &     rhcon,nrhcon,ipobody,ibody,xbodyact,co,nbody,network,
     &     iin_abs,vold,set,istep,iit,mi)
!     
      implicit none
!     
      logical identity,calcinitialpressure,gravity,gaspipe
!
      character*8 lakon(*)
      character*81 set(*)
!           
      integer ieg(*),nflow,i,j,ntg,ielmat(*),ntmat_,id,ntm,node1,node2,
     &     nelem,index,nshcon(*),ipkon(*),kon(*),ikboun(*),nboun,idof,
     &     nodem,idirf(5),nactdog(0:3,*),imat,ielprop(*),
     &     nodef(5),ndirboun(*),nodeboun(*),itg(*),node,kflag,ipiv(*),
     &     nrhs,info,idof1,idof2,nteq,nrhcon(*),ipobody(2,*),ibody(3,*),
     &     nbody,numf,network,iin_abs,icase,index2,index1,nelem1,nelem2,
     &     node11,node21,node12,node22,istep,iit,mi(2)
!     
      real*8 ac(ntm,ntm), bc(ntm),prop(*),shcon(0:3,ntmat_,*),
     &     f,df(5),xflow,xbounact(*),v(0:mi(2),*),cp,r,tg1,
     &     tg2,gastemp,physcon(*),pressmin,dvi,rho,g(3),
     &     rhcon(0:1,ntmat_,*),co(3,*),xbodyact(7,*),kappa,
     &     a,Tt,Pt,Ts,pressmax,constant,vold(0:mi(2),*)
!
      kflag=1
!
!     initialization of the pressure to zero except in the first 
!     iteration of the first step
!
c      if(network.ne.1) then
c         if((istep.eq.1).and.(iit.eq.-1).and. (iin_abs.ne.0)) then
c            do i=1,ntg
c               node=itg(i)
c                  v(2,node)=0.d0
c            enddo
c         endif
c      endif
!
      do j=1,nboun
         v(ndirboun(j),nodeboun(j))=xbounact(j)
      enddo
!     
!     no initial pressures for purely thermal calculations
!
      if(network.ne.0) then
!   
!        determining whether pressure initial conditions 
!        are provided for all nodes
!     
         pressmin=-1.d0
         pressmax=0.d0
         constant=1.55d0

         do i=1,ntg
            node=itg(i)
            if(v(2,node).lt.1.d-10) then
               v(2,node)=0.d0
            else
               if(pressmin.lt.0.d0) then
                  pressmin=v(2,node)
               elseif(v(2,node).lt.pressmin) then
                  pressmin=v(2,node)
               endif
!
               if(v(2,node).gt.pressmax)then
                  pressmax=v(2,node)
               endif
!
            endif
         enddo
!         
         if(pressmin.lt.0.d0) then
            write(*,*) '*ERROR in initialgas: minimum initial pressure'
            write(*,*) '       is smaller than zero'
            stop
         endif
!
!        in nodes in which no initial pressure is given v(2,*)
!        is replaced by -n, where n is the number of elements the
!        node belongs to: allows to find boundary nodes of the 
!        network
!
         gaspipe=.false.
         calcinitialpressure=.false.
!
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
!     
            if ((((lakon(nelem)(1:5).eq.'DGAPF')
     &           .or.(lakon(nelem)(1:5).eq.'DGAPI')).and.(iin_abs.eq.0))
     &           .or.((lakon(nelem)(1:3).eq.'DRE')
     &           .and.(lakon(nelem)(1:7).ne.'DREWAOR')
     &           .and.(iin_abs.eq.0))) then 
!     
!     In the case of a element of type GASPIPE or RESTRICTOR 
!     (except TYPE= RESTRICTOR WALL ORIFICE)
!     the number of pipes connected to node 1 and 2
!     are computed and stored in nactdog(3,node1)
!     respectively nactdog(3,node2)
!     
               gaspipe=.true.
               if(node1.ne.0) then
                  if (nactdog(3,node1).ge.0) then
!
                     if(node2.ne.0)then
                        nactdog(3,node1)=nactdog(3,node1)+1
                     endif
                  endif
               endif
               if(node2.ne.0) then
                  if (nactdog(3,node2).ge.0) then
                     if(node1.ne.0) then
                        nactdog(3,node2)=nactdog(3,node2)+1
                     endif
                  endif
               endif
            else
               if(iin_abs.eq.0) then
!     
!     for all other elements (different from GASPIPE or 
!     RESTRICTOR), including RESTRICTOR WALL ORIFICE 
!     nactdog(3,nodei)=-1
!     which means that they are connected to chambers
!     i.e. static and total values are equal
!     
                  if (node1.ne.0) then
                     nactdog(3,node1)=-1
                  endif
                  if(node2.ne.0) then
                     nactdog(3,node2)=-1
                  endif
               endif
            endif
!     
            if((node1.eq.0).or.(node2.eq.0)) cycle
            if(v(2,node1).lt.1.d-10) then
               v(2,node1)=v(2,node1)-1.d0
               calcinitialpressure=.true.
            endif
            if(v(2,node2).lt.1.d-10) then
               v(2,node2)=v(2,node2)-1.d0
               calcinitialpressure=.true.
            endif
         enddo
!     
         if(calcinitialpressure) then
!     
!           assigning values to the boundary nodes of the network
!     
            do i=1,ntg
               node=itg(i)
!     boundary condition nodes
               if(nactdog(2,node).eq.0) then 
                  v(2,node)=dtan(v(2,node)
     &                 *constant/pressmax) 
                  cycle 
               endif
!     initial condition nodes
               if((nactdog(2,node).ne.0)
     &              .and.(v(2,node).gt.0d0)) then
                  v(2,node) = dtan(v(2,node)
     &                 *constant/pressmax) 
                  cycle
               endif
!     nodes neither defined as *BOUNDARY nor as *INITIAL CONDITIONS
               if(abs(v(2,node)+1.d0).lt.1.d-10) then
                  v(2,node)=0.95d0*pressmin
                  pressmin=0.95d0*pressmin
                  v(2,node)=dtan(v(2,node)
     &                 *constant/pressmax)
               endif              
            enddo
!     
!           initializing ac and bc
!     
            do i=1,nteq
               do j=1,nteq
                  ac(i,j)=0.d0
               enddo
               bc(i)=0.d0
            enddo
!     
!           mass flow and temperature dofs: 1 on the diagonal
!     
            do i=1,ntg
               node=itg(i)
               do j=0,1
                  if(nactdog(j,node).ne.0) then
                     idof=nactdog(j,node)
                     ac(idof,idof)=1.d0
                  endif
               enddo
            enddo
!     
!           pressure dofs
!     
            do i=1,nflow
               nelem=ieg(i)
               index=ipkon(nelem)
               node1=kon(index+1)
               node2=kon(index+3)
               if((node1.eq.0).or.(node2.eq.0)) cycle
               idof1=nactdog(2,node1)
               idof2=nactdog(2,node2)
               if(idof1.ne.0) then
                  if(v(2,node1).gt.0.d0) then
                     ac(idof1,idof1)=1.d0
                     bc(idof1)=v(2,node1)
                  else
                     ac(idof1,idof1)=ac(idof1,idof1)+1.d0
                     if(v(2,node2).gt.0.d0) then
                        bc(idof1)=bc(idof1)+v(2,node2)
                     else
                        ac(idof1,idof2)=ac(idof1,idof2)-1.d0
                     endif
                  endif
               endif
               if(idof2.ne.0) then
                  if(v(2,node2).gt.0.d0) then
                     ac(idof2,idof2)=1.d0
                     bc(idof2)=v(2,node2)
                  else
                     ac(idof2,idof2)=ac(idof2,idof2)+1.d0
                     if(v(2,node1).gt.0.d0) then
                        bc(idof2)=bc(idof2)+v(2,node1)
                     else
                        ac(idof2,idof1)=ac(idof2,idof1)-1.d0
                     endif
                  endif
               endif
            enddo
!
!           solving the system
!     
            nrhs=1
            call dgesv(nteq,nrhs,ac,ntm,ipiv,bc,ntm,info)
            if(info.ne.0) then
               write(*,*) '*ERROR in initialgas: singular matrix'
               stop
            endif
!     
!           storing the initial conditions in v
!     
            do i=1,ntg
               node=itg(i)
               if(nactdog(2,node).eq.0) then
                  v(2,node)=pressmax
     &                 *datan(v(2,node))/constant
                  cycle
               endif
               v(2,node)=pressmax
     &              *datan(bc(nactdog(2,node)))/constant
            enddo
         endif
!
!     for ELEMENT TYPE BRANCH
!     initial pressures in branch 1 and 2 are equal and set to the 
!     maximum of the two pressures
!
         do i=1,nflow
            nelem=ieg(i)
            if(lakon(i)(2:5).ne.'REBR') cycle
            index=ielprop(nelem)
!     
            nelem1=prop(index+2)
            index1=ipkon(nelem1)
            node11=kon(index1+1)
            node21=kon(index1+3)
!     
            nelem2=prop(index+3)
            index2=ipkon(nelem2)
            node12=kon(index2+1)
            node22=kon(index2+3)
!     
            if(node11.eq.node12) then
               node1=node21
               node2=node22
            elseif(node11.eq.node22) then
               node1=node21
               node2=node12
            elseif(node21.eq.node12) then
               node1=node11
               node2=node22
            elseif(node21.eq.node22) then
               node1=node11
               node2=node12
            endif
!     
            if(v(2,node1).ge.v(2,node2)) then
               v(2,node2)=v(2,node1)
            else
               v(2,node1)=v(2,node2)
            endif
         enddo
      endif
!
!     temperature initial conditions
!
      do i=1,ntg
         node=itg(i)
         if (nactdog(0,node).eq.0) cycle
         if (v(0,node)+physcon(1).lt.1.d-10) then
            write(*,*)
     &           '*WARNING in initialgas : the initial temperature for n
     &ode',node
            write(*,*) 
     &           'is O Kelvin or less; the default is taken (293 K)'
            write(*,*)
            v(0,node)=293.d0+physcon(1)
         endif
      enddo
!  
!     determining the initial mass flow in those nodes for which no
!     flux boundary conditions are defined
!   
      if(network.ne.0) then
      do i=1,nflow
         nelem=ieg(i)
         index=ipkon(nelem)
         node1=kon(index+1)
         node2=kon(index+3)
         if((node1.eq.0).or.(node2.eq.0)) cycle
         nodem=kon(index+2)
         idof=8*(nodem-1)+1
         call nident(ikboun,idof,nboun,id)
         if (id.gt.0) then
            if (ikboun(id).eq.idof) then 
               cycle
            endif
         endif
!
!        for liquids: determine the gravity vector
!
         if(lakon(nelem)(2:3).eq.'LI') then
            gravity=.false.
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
                     gravity=.true.
                  endif
                  index=ipobody(2,index)
                  if(index.eq.0) exit
               enddo
            endif
            if(.not.gravity) then
               write(*,*) '*ERROR in initialgas: no gravity vector'
               write(*,*) '       was defined for liquid element',nelem
               stop
            endif
         endif
!
         tg1=v(0,node1)
         tg2=v(0,node2)
         gastemp=(tg1+tg2)/2.d0
         imat=ielmat(nelem)
         call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,dvi,
     &        rhcon,nrhcon,rho)
!
!     If inlet and outlet pressure are equal this leads to a massflow rate 
!     equal to 0 and in turn possibly to a singular matrix configuration
!
         if(v(2,node1).eq.v(2,node2))then
!
!     if neither inlet nor outlet pressure are active D.O.F: error-> stop
!
            if((nactdog(2,node1).eq.0)
     &           .and.(nactdog(2,node2).eq.0)) then
               WRITE(*,*) '*****************************************'
               write(*,*) '*ERROR:in subroutine initialgas.f'
               write(*,*) '       in element', nelem
               write(*,*) '       Inlet and outlet pressures are '
               write(*,*) '       boundary conditions '
               write(*,*) '       node1',node1,' pressure',v(2,node1)
               write(*,*) '       node2',node2,' pressure',v(2,node2)
               stop
!
!     if inlet pressure is an active degree of freedom
!
            else if((nactdog(2,node1).ne.0)
     &              .and.(nactdog(2,node2).eq.0))then
               WRITE(*,*) '*****************************************'
               write(*,*) '*WARNING:in subroutine initialgas.f'
               write(*,*) '       in element', nelem
               write(*,*) '       Inlet pressure initial condition '
               write(*,*) '       is changed '
               write(*,*) '       node1',node1,' given initial pressure'
     &              ,v(2,node1)
               v(2,node1)=1.1*v(2,node1)
               write(*,*) '       node1',node1,' new initial pressure',
     &              v(2,node1)               
               write(*,*) '       node2',node2,' pressure',v(2,node2)
!
!     if outlet pressure is an active D.O.F.
!
            else if((nactdog(2,node1).eq.0)
     &              .and.(nactdog(2,node2).ne.0))then
               WRITE(*,*) '*****************************************'
               write(*,*) '*WARNING:in subroutine initialgas.f'
               write(*,*) '       in element', nelem
               write(*,*) '       Outlet pressure initial condition '
               write(*,*) '       is changed '
               write(*,*) '       node1',node1,' pressure'
     &             ,v(2,node1)    
               write(*,*) '       node2',node2,'given intial pressure',
     &              v(2,node2)
               v(2,node2)=0.9*v(2,node2)
               write(*,*) '       node2',node2,' new initial pressure',
     &              v(2,node2)
!
!     if both inlet and outlet pressures are active D.O.F.
!
            else if((nactdog(2,node1).ne.0)
     &              .and.(nactdog(2,node2).ne.0))then
               WRITE(*,*) '*****************************************'
               write(*,*) '*WARNING:in subroutine initialgas.f'
               write(*,*) '       in element', nelem
               write(*,*) '       Inlet and outlet pressure '
               write(*,*) '       initial condition are changed '
               write(*,*) '       node1',node1,' given initial pressure'
     &             ,v(2,node1)  
               v(2,node1)=1.05*v(2,node2)
               write(*,*) '       node1',node1,' new intial pressure',
     &              v(2,node1)
               write(*,*) '       node2',node2,' given initial pressure'
     &             ,v(2,node2)
               v(2,node2)=0.95*v(2,node2)
               write(*,*) '       node2',node2,' new intial pressure',
     &              v(2,node2)
            endif
         endif

         call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,rho,physcon,g,co,dvi,numf,
     &        vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi)
!
         v(1,nodem)=xflow
!
         if(lakon(nelem)(2:4).ne.'LIP') then
            if(v(1,nodem).eq.0d0) then
               WRITE(*,*) '*****************************************'
               write(*,*) '*ERROR:in subroutine initialgas.f'
               write(*,*) '       in element', nelem,lakon(nelem)(1:6)
               write(*,*) '       mass flow rate value = 0 !'
               write(*,*) '       node1',node1,' pressure',v(2,node1)
               write(*,*) '       node2',node2,' pressure',v(2,node2)
               stop
            endif
            if (v(1,nodem).lt.0) then
               WRITE(*,*) '*****************************************'
               write(*,*) '*WARNING: in subroutine initialgas.f'
               write(*,*) '        in element', nelem
               write(*,*) '        mass flow rate value .le. 0 !'
               write(*,*) '        node1',node1,'pressure',v(2,node1)
               write(*,*) '        node2',node2,'pressure',v(2,node2)
               write(*,*) '        check element definition'
            endif
         endif
      enddo
!
!     calculating the static temperature for nodes belonging to gas pipes
!     and restrictors (except RESTRICTOR WALL ORIFICE)
!  
      if(gaspipe) then
         node=nelem
      endif
      if(iin_abs.eq.0) then
         node=nelem
      endif
      if (gaspipe.and.(iin_abs.eq.0)) then
!
!        nactdog(3,node) is set to zero for chamber nodes
!
         do i=1,ntg
            node=itg(i)
            if(nactdog(3,node).lt.0) nactdog(3,node)=0
            if(nactdog(3,node).gt.2) nactdog(3,node)=-nactdog(3,node)
         enddo
!
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
!
!     if exactly one or two pipes are connected to a node then
!     the number of the element the node belongs to is stored 
!     in nactdog(3,nodei)
!     
            if(node1.gt.0) then
               if((nactdog(3,node1).eq.1).or.
     &            (nactdog(3,node1).eq.2)) then
                  if(node2.ne.0)then
                     nactdog(3,node1)=nelem
                  endif
               elseif(nactdog(3,node1).lt.-2) then  
                  nactdog(3,node1)=0
                  write(*,*) '*WARNING :in subroutine initialgas.f'
                  write(*,*) '          more than 2 elements GASPIPE'
                  write(*,*) '          or RESTRICTOR are connected '
                  write(*,*) '          to the same node',node1
                  write(*,*) '          The common node is a chamber'
                  write(*,*) '          Total and static parameters are'
                  write(*,*) '          equal'
               endif
            endif
!     
            if(node2.gt.0) then
               if((nactdog(3,node2).eq.1).or.
     &            (nactdog(3,node2).eq.2)) then
                  if (node1.ne.0) then
                     nactdog(3,node2)=nelem
                  endif
               elseif(nactdog(3,node2).lt.-2) then  
                  nactdog(3,node2)=0
                  write(*,*) '*WARNING :in subroutine initialgas.f'
                  write(*,*) '          more than 2 elements GASPIPE'
                  write(*,*) '          or RESTRICTOR are connected '
                  write(*,*) '          to the same node',node2
                  write(*,*) '          The common node is a chamber'
                  write(*,*) '          Total and static parameters are'
                  write(*,*) '          equal'
               endif
            endif
         enddo
!
!        for a non-chamber node, the gaspipe element this node belongs
!        to is stored in nactdog(3,node)
!
         do i=1,ntg
            node=itg(i)
            if(nactdog(3,node).le.0) then           
               nactdog(3,node)=-1
            endif
         enddo
!
!     The static temperature is calculated and stored in v(3,node)
!     total temperatures are supposed equal (adiabatic pipe)
!     
         do i=1,ntg
            node=itg(i)
            if(nactdog(3,node).gt.0) then 
!     
               nelem=nactdog(3,node)
               index=ielprop(nelem)
               nodem=kon(ipkon(nelem)+2)
!     
               imat=ielmat(nelem)
               call materialdata_tg(imat,ntmat_,v(0,node),
     &              shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
               kappa=cp/(cp-R)
               xflow=v(1,nodem)
               Tt=v(0,node)
               Pt=v(2,node)
!      
               if((lakon(nelem)(2:5).eq.'GAPF')
     &            .or.(lakon(nelem)(2:5).eq.'GAPI')) then
                  A=prop(index+1)
                  if((lakon(nelem)(2:6).eq.'GAPFA') 
     &            .or.(lakon(nelem)(2:6).eq.'GAPIA')) then
                     icase=0
                  elseif((lakon(nelem)(2:6).eq.'GAPFI')
     &                .or.(lakon(nelem)(2:6).eq.'GAPII')) then
                     icase=1
                  endif     
!               
               elseif(lakon(nelem)(2:3).eq.'RE') then
                  index2=ipkon(nelem)
                  node1=kon(index2+1)
                  node2=kon(index2+3)
                  if(lakon(nelem)(4:5).eq.'EX') then
                     if((lakon(int(prop(index+4)))(2:6).eq.'GAPFA')
     &             .or.(lakon(int(prop(index+4)))(2:6).eq.'GAPIA')) then
                        icase=0
                     elseif((lakon(int(prop(index+4)))(2:6).eq.'GAPFI')
     &                 .or.(lakon(int(prop(index+4)))(2:6).eq.'GAPII')) 
     &                  then
                        icase=1
                     endif
                  else
                     icase=0
                  endif
!
                  if(lakon(nelem)(4:5).eq.'BE') then
                     a=prop(index+1)
!
                  elseif(lakon(nelem)(4:5).eq.'BR') then
                     if(lakon(nelem)(4:6).eq.'BRJ') then
                        if(nelem.eq.nint(prop(index+2)))then
                              A=prop(index+5)
                        elseif(nelem.eq.nint(prop(index+3))) then
                              A=prop(index+6)
                        endif
                     elseif(lakon(nelem)(4:6).eq.'BRS') then
                        if(nelem.eq.nint(prop(index+2)))then
                              A=prop(index+5)
                        elseif(nelem.eq.nint(prop(index+3))) then
                              A=prop(index+6)
                        endif
                     endif
!                     
                  else
                     if(node.eq.node1) then
                        a=prop(index+1)
                     elseif(node.eq.node2) then
                        a=prop(index+2)
                     endif
                  endif
               endif
!     
               if(v(3,node).eq.0) then
                  call ts_calc(xflow,Tt,Pt,kappa,r,a,Ts,icase)
                  v(3,node)=Ts
               endif
!     
!     if the element is not of gaspipe or branch type,
!     total and static temperatures are equal
!
            elseif(node.ne.0) then
               v(3,node)=v(0,node)
            endif
         enddo
      endif
      endif
!    
!     initialisation of bc
!     
      do i=1,ntm
         bc(i)=0.d0
      enddo
!
      return
      end
