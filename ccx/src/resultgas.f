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
!     
!     construction of the B matrix      
!     
      subroutine resultgas(itg,ieg,ntg,ntm,
     &     bc,nload,sideload,nelemload,xloadact,lakon,ntmat_,
     &     voldgas,shcon,nshcon,ipkon,kon,co,nflow,ikboun,xbounact,
     &     nboun,iinc,istep,dtime,ttime,time,ilboun,ndirboun,nodeboun,
     &     ikforc,ilforc,xforcact,nforc,ielmat,nteq,prop,ielprop,
     &     nactdog,nacteq,iin,physcon,camt,camf,camp,rhcon,nrhcon,
     &     ipobody,ibody,xbodyact,nbody,dtheta,vold)
!     
      implicit none
!     
      logical identity
      character*8 lakonl,lakon(*)
      character*20 sideload(*)
!     
      integer itg(*),ieg(*),ntg,nteq,nflow,nload,ielmat(*),iflag,
     &     nelemload(2,*),nope,nopes,mint2d,i,j,k,l,nrhcon(*),
     &     node,imat,ntmat_,id,ntm,ifaceq(8,6),ifacet(6,4),numf,
     &     ifacew(8,5),node1,node2,nshcon(*),nelem,ig,index,konl(20),
     &     ipkon(*),kon(*),ikboun(*),nboun,ndirboun(*),nodeboun(*),idof,
     &     iinc,istep,jltyp,nfield,ilboun(*),ikforc(*),ipobody(2,*),
     &     ilforc(*),nforc,nodem,idirf(5),ieq,nactdog(0:3,*),nbody,
     &     nacteq(0:3,*),ielprop(*),nodef(5),iin,kflag,ibody(3,*),case,
     &     inv, index2
!     
      real*8 bc(ntm),xloadact(2,*),cp,h(2),physcon(3),r,dvi,rho,
     &     xl2(0:3,8),coords(3),dxsj2,temp,xi,et,weight,xsj2(3),
     &     gastemp,voldgas(0:3,*),shcon(0:3,ntmat_,*),co(3,*),shp2(4,8),
     &     xbounact(*),field,prop(*),tg1,tg2,dtime,ttime,time,g(3),
     &     xforcact(*),areaj,xflow,tvar(2),f,df(5),camt,camf,camp,
     &     rhcon(0:1,ntmat_,*),xbodyact(7,*),sinktemp,kappa,a,T,Tt,Pt,
     &     M,dtheta,ts1,ts2,xs2(3,2),xk1,xk2,xdenom1,xdenom2,expon,pt1,
     &     pt2,dt1,dt2,xcst,xnum1,xnum2,Qred_crit,xflow_crit,
     &     nelem0,nodem0,xflow0,nelem1,nodem1,xflow1,
     &     nelem2,nodem2,xflow2,R1,R2,Rout,Rin,Uout,Uin,heat,pi,
     &     Cp_cor,U,Ct,nelemswirl,vold(0:3,*)
!     
      real*8 gauss2d1(2,1),gauss2d2(2,4),gauss2d3(2,9),gauss2d4(2,1),
     &     gauss2d5(2,3),gauss3d1(3,1),gauss3d2(3,8),gauss3d3(3,27),
     &     gauss3d4(3,1),gauss3d5(3,4),gauss3d6(3,15),gauss3d7(3,2),
     &     gauss3d8(3,9),gauss3d9(3,18),weight2d1(1),weight2d2(4),
     &     weight2d3(9),weight2d4(1),weight2d5(3),weight3d1(1),
     &     weight3d2(8),weight3d3(27),weight3d4(1),weight3d5(4),
     &     weight3d6(15),weight3d7(2),weight3d8(9),weight3d9(18)
!     
      include "gauss.f"
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data iflag /2/
!     
      kflag=2
!     
      tvar(1)=time
      tvar(2)=ttime+dtime
!     
      pi=4.d0*datan(1.d0)
!     
!     calculating the maximum change in the solution
!     
      camt=0.d0
      camf=0.d0
      camp=0.d0
!     
      do i=1,ntg
         node=itg(i)
         do j=0,2
            if(nactdog(j,node).eq.0) cycle
            idof=nactdog(j,node)
            if(j.eq.0) then
               if(dabs(bc(idof)).gt.camt) camt=dabs(bc(idof))
            elseif(j.eq.1) then
               if(dabs(bc(idof)).gt.camf) camf=dabs(bc(idof))
            else
               if(dabs(bc(idof)).gt.camp) camp=dabs(bc(idof))
            endif
         enddo
      enddo
!     
!     updating voldgas
!     
      do i=1,ntg
         node=itg(i)
         do j=0,2
            if(nactdog(j,node).eq.0) cycle
            voldgas(j,node)=voldgas(j,node)+bc(nactdog(j,node))*dtheta
         enddo
      enddo
!     
!     testing the validity of the pressures
!     
      do i=1,ntg
         node=itg(i)
         if(voldgas(2,node).lt.0) then
            write(*,*) 'wrong pressure'
            iin=0
            return
         endif
      enddo
!     
!     testing validity of temperatures
!     
      do i=1,ntg
         node=itg(i)
         if(voldgas(0,node).lt.0) then
            iin=0
            return
         endif
      enddo
!     
!     testing the validity of the solution for branches elements
!     and restrictor. Since the element properties is dependent on
!     
!     
      do i=1, nflow
         nelem=ieg(i)
         if ((lakon(nelem)(4:5).eq.'ATR').or. 
     &        (lakon(nelem)(4:5).eq.'RTA')) then
            xflow=voldgas(1,kon(ipkon(nelem)+2))
            if(xflow.lt.0d0)then
               Write(*,*)'*WARNING in resultgas.f'
               write(*,*)'Element',nelem,'of TYPE ABSOLUTE TO RELATIVE'
               write(*,*)'The flow direction is no more conform '
               write(*,*)'to element definition'
               write(*,*)'Check the pertinence of the results'
            endif
         elseif(lakon(nelem)(2:3).eq.'RE') then
!     
            if(lakon(nelem)(4:5).ne.'BR') then
               nodem=kon(ipkon(nelem)+2)
               xflow=voldgas(1,nodem)
               if (xflow.lt.0) then
                  Write(*,*)'*WARNING in resultgas.f'
                  write(*,*)'Element',nelem,'of TYPE RESTRICTOR'
                  write(*,*)'The flow direction is no more conform '
                  write(*,*)'to element definition'
                  write(*,*)'Check the pertinence of the results'
               endif
!     
            elseif(lakon(nelem)(4:5).eq.'BR') then
               index=ielprop(nelem)
!     
               nelem0=int(prop(index+1))
               nodem0=kon(ipkon(nelem0)+2)
               xflow0=voldgas(1,nodem0)
!     
               nelem1=int(prop(index+2))
               nodem1=kon(ipkon(nelem1)+2)
               xflow1=voldgas(1,nodem1)
!     
               nelem2=int(prop(index+3))
               nodem2=kon(ipkon(nelem2)+2)
               xflow2=voldgas(1,nodem2)
!     
               if((xflow0.lt.0).or.(xflow1.lt.0).or.(xflow2.lt.0)) then
                  Write(*,*)'*WARNING in resultgas.f'
                  write(*,*)'Element',nelem,'of TYPE BRANCH'
                  write(*,*)'The flow direction is no more conform '
                  write(*,*)'to element definition'
                  write(*,*)'Check the pertinence of the results'
               endif
            endif
         endif
      enddo
!     
!     node1 or node2 do not belong to GASPIPE or RESTRICTOR element
!     
      do i=1,ntg
         node=itg(i)
         nelem=nactdog(3,node)
         if(nelem.le.0) then
            voldgas(3,node)=voldgas(0,node)
         endif
      enddo
!     
!     iteratively solving Tt=T+0.5*v**2/(2*Cp) to obtain T static
!     
      do i=1,ntg
         node=itg(i)
         nelem=nactdog(3,node)
!     
         if (nelem.gt.0) then 
!     
            nodem=kon(ipkon(nelem)+2)
            T=voldgas(3,node)
            Tt=voldgas(0,node)
            Pt=voldgas(2,node)
            xflow=voldgas(1,nodem)
!     
            imat=ielmat(nelem)
            call materialdata_tg(imat,ntmat_,voldgas(3,node),
     &           shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!     
            index=ielprop(nelem)
            kappa=(cp/(cp-R))
!     
            if(lakon(nelem)(2:5).eq.'GAPI') then
               A=prop(index+1)
               if(lakon(nelem)(2:6).eq.'GAPIA') then
                  case=0
               elseif(lakon(nelem)(2:6).eq.'GAPII') then
                  case=1
               endif  
            elseif(lakon(nelem)(2:3).eq.'OR') then
               A=prop(index+1)
               case=0
!     
            elseif(lakon(nelem)(2:3).eq.'RE') then
               index2=ipkon(nelem)
               node1=kon(index2+1)
               node2=kon(index2+3)
!     
               if(lakon(nelem)(4:5).eq.'EX') then
                  if(lakon(int(prop(index+4)))(2:6).eq.'GAPIA') then
                     case=0
                  elseif(lakon(int(prop(index+4)))(2:6).eq.'GAPII')then
                     case=1
                  endif
               else
                  case=0
               endif
!     
!     defining the sections
!     
               if(lakon(nelem)(4:5).eq.'BE') then
                  a=prop(index+1)
!     
               elseif(lakon(nelem)(4:5).eq.'BR') then
                  if(lakon(nelem)(4:6).eq.'BRJ') then
                     if(nelem.eq.int(prop(index+2)))then
                        A=prop(index+5)
                     elseif(nelem.eq.int(prop(index+3))) then
                        A=prop(index+6)
                     endif
                  elseif(lakon(nelem)(4:6).eq.'BRS') then
                     if(nelem.eq.int(prop(index+2)))then
                        A=prop(index+5)
                     elseif(nelem.eq.int(prop(index+3))) then
                        A=prop(index+6)
                     endif
                  endif
!     
               else
!     
                  if(node.eq.node1) then
                     a=prop(index+1)
                  elseif(node.eq.node2) then
                     a=prop(index+2)
                  endif
               endif
            endif

            if(xflow.lt.0d0) then
               inv=-1
            else
               inv=1
            endif
!     
            if(case.eq.0) then
               Qred_crit=dsqrt(kappa/R)*(1.+0.5*(kappa-1.))
     &              **(-0.5d0*(kappa+1.)/(kappa-1.))
            else
               Qred_crit=dsqrt(1/R)*(1.+0.5*(kappa-1.)/kappa)
     &              **(-0.5d0*(kappa+1.)/(kappa-1.))
            endif
            xflow_crit=inv*Qred_crit*Pt*A/dsqrt(Tt)             
!     
            if(dabs(voldgas(1,nodem)).ge.dabs(xflow_crit)) then
               if(case.eq.1) then
                  voldgas(1,nodem)=xflow_crit
               endif
            endif
!     
            call ts_calc(xflow,Tt,Pt,kappa,r,a,T,case)
!     
            voldgas(3,node)=T
!     
!     Mach number
!     
            M=(2.d0/(kappa-1.d0)*(Tt/T-1.d0))**0.5d0
!     
         endif
!     
      enddo
!     
!     reinitialisation of the Bc matrix
!     
      do i=1,nteq
         bc(i)=0.d0
      enddo
!     
!     determining the residual
!     
      do i=1,nflow
         nelem=ieg(i)
         index=ipkon(nelem)
         node1=kon(index+1)
         nodem=kon(index+2)
         node2=kon(index+3)
         xflow=voldgas(1,nodem)
!     
         if(node1.eq.0) then
            tg1=voldgas(0,node2)
            tg2=tg1
            ts1=voldgas(3,node2)
            ts2=Ts1
         elseif(node2.eq.0) then
            tg1=voldgas(0,node1)
            tg2=tg1
            ts1=voldgas(3,node1)
            ts2=ts1
         else
            tg1=voldgas(0,node1)
            tg2=voldgas(0,node2)
            ts1=voldgas(3,node1)
            ts2=voldgas(3,node2)
         endif
         gastemp=(ts1+ts2)/2.d0
!     
         imat=ielmat(nelem)
!     
         call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,dvi,
     &        rhcon,nrhcon,rho)
!     
!     Definitions of the constant for isothermal flow elements
!     
         if(lakon(nelem)(2:6).eq.'GAPII') then
            if((node1.ne.0).and.(node2.ne.0)) then
               A=prop(ielprop(nelem)+1)
               pt1=voldgas(2,node1)
               pt2=voldgas(2,node2)
!     
               if(pt1.ge.pt2)then
                  if((voldgas(3,nodem)).ge.(pt2/pt1))then
                     pt2=voldgas(3,nodem)*pt1
                  endif
                  tg1=voldgas(0,node1)
                  ts1=voldgas(3,node1)
                  tg2=voldgas(0,node2)
                  ts2=voldgas(3,node2)
               else
                  pt1=voldgas(2,node2)
                  pt2=voldgas(2,node1)
                  if(voldgas(3,nodem).ge.(pt2/pt1))then
                     pt2=voldgas(3,nodem)*pt1
                  endif
                  tg1=voldgas(0,node2)
                  ts1=voldgas(3,node2)
                  tg2=voldgas(0,node1)
                  ts2=voldgas(3,node1)
               endif
!     
               dt1=tg1/ts1-1d0
               dt2=tg2/ts2-1d0
               xcst=2.d0*Cp*A**2/(R**2)
               expon=2.d0*kappa/(kappa-1.d0)
               xk1=pt1**2*(ts1/tg1)**expon
               xk2=pt2**2*(ts2/tg2)**expon
!     
               xnum1=xcst*dt1*xk1-xflow**2*ts1
               xdenom1=xcst*xk1*(1.d0-expon*dt1)/ts1+2.d0*xflow**2
               xnum2=xcst*dt2*xk2-xflow**2*ts2
               xdenom2=xcst*xk2*(1.d0-expon*dt2)/ts2+2.d0*xflow**2
!     
            endif
         endif
!     
         if(node1.ne.0) then
!     
!     energy equation contribution node1
!     
            if (nacteq(0,node1).ne.0) then
               ieq=nacteq(0,node1)
!     
               if(nacteq(3,node1).eq.0) then
                  if (xflow.lt.0d0)then
                     bc(ieq)=bc(ieq)+cp*(tg1-tg2)*xflow
                  endif
!     
               elseif((nacteq(3,node1).eq.node2)) then
!     
                  bc(ieq)=(ts2+xnum2/xdenom2-ts1-xnum1/xdenom1)
!     
               endif
            endif
!     
!     mass equation contribution node1
!     
            if (nacteq(1,node1).ne.0) then
               ieq=nacteq(1,node1)
               bc(ieq)=bc(ieq)-xflow
            endif
         endif
!     
         if(node2.ne.0) then
!     
!     energy equation contribution node2
!     
            if (nacteq(0,node2).ne.0) then
               ieq=nacteq(0,node2)
!     
               if(nacteq(3,node2).eq.0) then
                  if (xflow.gt.0d0)then
                     bc(ieq)=bc(ieq)-cp*(tg2-tg1)*xflow
                  endif
!     
               elseif(nacteq(3,node2).eq.node1) then
!     
                  bc(ieq)=(ts2+xnum2/xdenom2-ts1-xnum1/xdenom1) 
!     
               endif
            endif
!     
!     mass equation contribution node2
!     
            if (nacteq(1,node2).ne.0) then
               ieq=nacteq(1,node2)
               bc(ieq)=bc(ieq)+xflow
            endif
         endif
!     
!     element equation
!     
         if (nacteq(2,nodem).ne.0) then
            ieq=nacteq (2,nodem)
!     
!     for liquids: determine the gravity vector
!     
            if(lakon(nelem)(2:3).eq.'LI') then
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
            endif
!     
            call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &           nactdog,identity,
     &           ielprop,prop,kflag,voldgas,xflow,f,nodef,idirf,df,
     &           cp,r,rho,physcon,g,co,dvi,numf,vold)
            bc(ieq)=-f
         endif
      enddo
!     
!     convection with the walls
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'FC') then
            nelem=nelemload(1,i)
            lakonl=lakon(nelem)
            node=nelemload(2,i)
            ieq=nacteq(0,node)
            if(ieq.eq.0) then 
               cycle
            endif
!     
            call nident(itg,node,ntg,id)
!     
!     calculate the area
!     
            read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
!     
            if(lakonl(4:4).eq.'2') then
               nope=20
               nopes=8
            elseif(lakonl(4:4).eq.'8') then
               nope=8
               nopes=4
            elseif(lakonl(4:5).eq.'10') then
               nope=10
               nopes=6
            elseif(lakonl(4:4).eq.'4') then
               nope=4
               nopes=3
            elseif(lakonl(4:5).eq.'15') then
               nope=15
            else
               nope=6
            endif
!     
            if(lakonl(4:5).eq.'8R') then
               mint2d=1
            elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &              then
               if(lakonl(7:7).eq.'A') then
                  mint2d=2
               else
                  mint2d=4
               endif
            elseif(lakonl(4:4).eq.'2') then
               mint2d=9
            elseif(lakonl(4:5).eq.'10') then
               mint2d=3
            elseif(lakonl(4:4).eq.'4') then
               mint2d=1
            endif
!     
            if(lakonl(4:4).eq.'6') then
               mint2d=1
               if(ig.le.2) then
                  nopes=3

               else
                  nopes=4
               endif
            endif
            if(lakonl(4:5).eq.'15') then
               if(ig.le.2) then
                  mint2d=3
                  nopes=6
               else
                  mint2d=4
                  nopes=8
               endif
            endif
!     
!     connectivity of the element
!     
            index=ipkon(nelem)
            if(index.lt.0) then
               write(*,*) '*ERROR in radflowload: element ',nelem
               write(*,*) '       is not defined'
               stop
            endif
            do k=1,nope
               konl(k)=kon(index+k)
            enddo
!     
!     coordinates of the nodes belonging to the face
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do k=1,nopes
                  xl2(0,k)=voldgas(0,konl(ifaceq(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &                    voldgas(j,konl(ifaceq(k,ig)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do k=1,nopes
                  xl2(0,k)=voldgas(0,konl(ifacet(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &                    voldgas(j,konl(ifacet(k,ig)))
                  enddo
               enddo
            else
               do k=1,nopes
                  xl2(0,k)=voldgas(0,konl(ifacew(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &                    voldgas(j,konl(ifacew(k,ig)))
                  enddo
               enddo
            endif
!     
!     integration to obtain the area and the mean
!     temperature
!     
            do l=1,mint2d
               if((lakonl(4:5).eq.'8R').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,l)
                  et=gauss2d1(2,l)
                  weight=weight2d1(l)
               elseif((lakonl(4:4).eq.'8').or.
     &                 (lakonl(4:6).eq.'20R').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                  xi=gauss2d2(1,l)
                  et=gauss2d2(2,l)
                  weight=weight2d2(l)
               elseif(lakonl(4:4).eq.'2') then
                  xi=gauss2d3(1,l)
                  et=gauss2d3(2,l)
                  weight=weight2d3(l)
               elseif((lakonl(4:5).eq.'10').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                  xi=gauss2d5(1,l)
                  et=gauss2d5(2,l)
                  weight=weight2d5(l)
               elseif((lakonl(4:4).eq.'4').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                  xi=gauss2d4(1,l)
                  et=gauss2d4(2,l)
                  weight=weight2d4(l)
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
!     
               dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &              xsj2(3)*xsj2(3))
               areaj=dxsj2*weight
!     
               temp=0.d0
               do k=1,3
                  coords(k)=0.d0
               enddo
               do j=1,nopes
                  temp=temp+xl2(0,j)*shp2(4,j)
                  do k=1,3
                     coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
               enddo
!     
               sinktemp=voldgas(0,node)
               if(sideload(i)(5:6).ne.'NU') then
                  h(1)=xloadact(1,i)
               else
                  read(sideload(i)(2:2),'(i1)') jltyp
                  jltyp=jltyp+10
                  call film(h,sinktemp,temp,istep,
     &                 iinc,tvar,nelem,l,coords,jltyp,field,nfield,
     &                 sideload(i),node,areaj,voldgas)
               endif
               if(lakonl(5:7).eq.'0RA') then
                  bc(ieq)=bc(ieq)+
     &                 2.d0*(temp-sinktemp)*h(1)*dxsj2*weight
               else
                  bc(ieq)=bc(ieq)+
     &                 (temp-sinktemp)*h(1)*dxsj2*weight
               endif
            enddo
         endif
      enddo
!     
!     prescribed heat generation        
!     
      do i=1,ntg
         node=itg(i)
         idof=7*(node-1)
         call nident(ikforc,idof,nforc,id)
         if(id.gt.0) then
            if(ikforc(id).eq.idof) then
               ieq=nacteq(0,node)
               bc(ieq)=bc(ieq)+xforcact(ilforc(id))
               cycle
            endif
         endif
      enddo
!     
!     in the case of forced vortices when, when temperature change 
!     is required , an additionnal heat input is added in the energy equation for the 
!     downstream node
!     
      do i=1,nflow
         nelem=ieg(i)
         if(lakon(nelem)(2:3).ne.'VO') cycle
!     
!     free vortex and no temperature change
!     
         if((lakon(nelem)(2:5).eq.'VOFR').and.
     &        (prop(ielprop(nelem)+8).eq.0)) cycle
!     
!     free vortex and temperature change in the absolute system
!     
         if((lakon(nelem)(2:5).eq.'VOFR').and.
     &        (prop(ielprop(nelem)+8).eq.1)) cycle
!     
!     forced vortex and no temperature change
!     
         if((lakon(nelem)(2:5).eq.'VOFO').and.
     &        (prop(ielprop(nelem)+6).eq.0)) cycle
!     
         nodem=kon(ipkon(nelem)+2)
         xflow=voldgas(1,nodem)
         if(xflow.gt.0d0) then
            node1=kon(ipkon(nelem)+1)
            node2=kon(ipkon(nelem)+3)
         else
            node1=kon(ipkon(nelem)+1)
            node2=kon(ipkon(nelem)+3)
         endif
!     
         if(xflow.gt.0d0) then
            R1=prop(ielprop(nelem)+2)
            R2=prop(ielprop(nelem)+1)
            if(R1.gt.R2) then
               Rout=R2
               Rin=R1
            else
               Rout=R2
               Rin=R1
            endif
         else
            R1=prop(ielprop(nelem)+2)
            R2=prop(ielprop(nelem)+1)
            if(R1.gt.R2) then
               Rout=R1
               Rin=R2
            else
               Rout=R1
               Rin=R2
            endif
         endif
!     
!     computing temperature corrected Cp=Cp(T) coefficient
!     
!     
         Tg1=voldgas(0,node1)
         Tg2=voldgas(0,node2)
         gastemp=(Tg1+Tg2)/2.d0
!
         imat=ielmat(nelem)
         call materialdata_tg(imat,ntmat_,gastemp,
     &        shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
         
       
!
         call cp_corrected(cp,Tg1,Tg2,cp_cor)
!     
         Uout=Pi/30*prop(ielprop(nelem)+5)*Rout
         Uin=Pi/30*prop(ielprop(nelem)+5)*Rin
!     
!     free and forced vortices with temperature 
!     change in the relative system of coordinates 
!     
         if((lakon(nelem)(2:5).eq.'VOFR') .and.
     &        (prop(ielprop(nelem)+8).eq.(-1))) then
!     
            Uout=Pi/30*prop(ielprop(nelem)+7)*Rout
            Uin=Pi/30*prop(ielprop(nelem)+7)*Rin
!     
            heat=0.5d0*Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
         elseif (((lakon(nelem)(2:5).eq.'VOFO')
     &           .and.(prop(ielprop(nelem)+6).eq.(-1)))) then
!     
            Uout=Pi/30*prop(ielprop(nelem)+5)*Rout
            Uin=Pi/30*prop(ielprop(nelem)+5)*Rin
!     
            heat=0.5d0*Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
!     forced vortices with temperature change in the absolute system
!     
         elseif((lakon(nelem)(2:5).eq.'VOFO')
     &           .and.((prop(ielprop(nelem)+6).eq.1))) then
!     
            heat=Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
         endif
!     
!     including the resulting additional heat flux in the energy equation
!     
         if(xflow.gt.0d0)then
            ieq=nacteq(0,node2)
            if(nacteq(0,node2).ne.0)then
               bc(ieq)=bc(ieq)+heat
!     
            endif
         else
            ieq=nacteq(0,node1)
            if(nacteq(0,node1).ne.0)then
               bc(ieq)=bc(ieq)+heat
            endif
         endif
      enddo
!     
!     transfer element ABSOLUTE TO RELATIVE / RELATIVE TO ABSOLUTE
!     
      do i= 1, nflow
         nelem=ieg(i)
!     
         if((lakon(nelem)(2:4).eq.'ATR').or.
     &        (lakon(nelem)(2:4).eq.'RTA'))  then
!     
            nodem=kon(ipkon(nelem)+2)
            xflow=voldgas(1,nodem)
            if(xflow.gt.0d0) then
               node1=kon(ipkon(nelem)+1)
               node2=kon(ipkon(nelem)+3)

            else
               node1=kon(ipkon(nelem)+1)
               node2=kon(ipkon(nelem)+3)

            endif
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
!     
            Tg1=voldgas(0,node1)
            Tg2=voldgas(0,node2)
            gastemp=(Tg1+Tg2)/2.d0
!
            imat=ielmat(nelem)
            call materialdata_tg(imat,ntmat_,gastemp,
     &           shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!
            call cp_corrected(cp,Tg1,Tg2,cp_cor)
!     
            index=ielprop(nelem)
            U=prop(index+1)
            ct=prop(index+2)
!     
            if(ct.eq.0) then
               nelemswirl=prop(index+3)
               index2=ielprop(nelemswirl)
!     
!     previous element is a preswirl nozzle
!     
               if(lakon(nelemswirl)(2:5).eq.'ORPN') then
                  ct=prop(index2+4)
!     
!     previous element is a forced vortex
!     
               elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
                  ct=prop(index2+7)
!     
!     previous element is a free vortex
!     
               elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
                  ct=prop(index2+9)
               endif
            endif                
!     
            if(lakon(nelem)(2:4).eq.'ATR') then
               heat=Cp/Cp_cor*(0.5d0*(U**2-2d0*U*Ct)*xflow)
!     
            elseif(lakon(nelem)(2:4).eq.'RTA') then
               heat=Cp/Cp_cor*(-0.5d0*(U**2-2d0*U*Ct)*xflow)
            endif
!     
!     including the resulting additional heat flux in the energy equation
!     
            if(xflow.gt.0d0)then
               ieq=nacteq(0,node2)
               if(nacteq(0,node2).ne.0)then
                  bc(ieq)=bc(ieq)+heat
!     
               endif
            else
               ieq=nacteq(0,node1)
               if(nacteq(0,node1).ne.0)then
                  bc(ieq)=bc(ieq)+heat
               endif
            endif
         endif
      enddo         
!     
      return
      end
