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
      subroutine generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,nr,nz,
     &  rcs0,zcs0,ncs_,cs,labmpc,istep,istat,n,
     &  mcs,ithermal,triangulation,csab,xn,yn,zn,phi,nodel,ncsnodes,
     &  nodesonaxis,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,lcs,
     &  kontri,straight,ne,ipkon,kon,lakon,ifacetet,inodface)
!
!     generate cyclic mpc's
!
      implicit none
!
      logical triangulation,nodesonaxis,interpolation
!
      character*8 lakon(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),nneigh,ne,ipkon(*),kon(*),
     &  istep,istat,n,j,k,nk,nmpc,nmpc_,mpcfree,ics(*),nterms,
     &  nr(*),nz(*),nodel,noder,ikmpc(*),ilmpc(*),kontri(3,*),
     &  number,idof,ndir,node,ncsnodes,id,mpcfreeold,ncs_,
     &  mcs,ithermal,nrcg(*),nzcg(*),jcs(*),lcs(*),nodef(8),
     &  netri,ifacetet(*),inodface(*)
!
      real*8 tolloc,co(3,*),coefmpc(*),rcs(*),zcs(*),rcs0(*),zcs0(*),
     &  csab(7),xn,yn,zn,xap,yap,zap,rp,zp,al(3,3),ar(3,3),cs(17,*),phi,
     &  x2,y2,z2,x3,y3,z3,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),
     &  straight(9,*),ratio(8)
!
      nneigh=1
!
      xap=co(1,nodel)-csab(1)
      yap=co(2,nodel)-csab(2)
      zap=co(3,nodel)-csab(3)
!     
      zp=xap*xn+yap*yn+zap*zn
      rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
      call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,node,nneigh)
      noder=abs(ics(node))
!
!     check whether node is on axis
!
      if(noder.eq.nodel) then
         return
      endif
!
      interpolation=.false.
!
      if(rp.gt.1.d-10) then
         x2=(xap-zp*xn)/rp
         y2=(yap-zp*yn)/rp
         z2=(zap-zp*zn)/rp
         x3=yn*z2-y2*zn
         y3=x2*zn-xn*z2
         z3=xn*y2-x2*yn
      endif
!
         if((tolloc.ge.0.d0).and.
     &        (tolloc.le.dsqrt((rp-rcs0(node))**2+(zp-zcs0(node))**2)))
     &        then
!
!           the nodal positions on the dependent and independent
!           sides of the mpc's do no agree: interpolation is
!           necessary. 
!
            write(*,*) '*WARNING in generatecycmpcs: no cyclic'
            write(*,*) '         symmetric partner found for'
            write(*,*) '         node ',nodel,'. A new partner'
            write(*,*) '         is created and connected to an'
            write(*,*) '         existing element side'
            write(*,*)
!
            interpolation=.true.
!
            if(.not.triangulation) then
               call triangulate(co,nk,ics,rcs0,zcs0,ncsnodes,
     &              rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,
     &              straight,ne,ipkon,kon,lakon,lcs,netri,ifacetet,
     &              inodface)
               triangulation=.true.
            endif
!
            label='CYCLIC              '
            if(mcs.lt.10) then
               write(label(7:7),'(i1)') mcs
            elseif(mcs.lt.100) then
               write(label(7:8),'(i2)') mcs
            else
               write(*,*)'*ERROR in generatecycmpcs: no more than 99'
               write(*,*)'       cyclic symmetry definitions allowed'
               stop
            endif
!
            noder=nk+1
            co(1,noder)=csab(1)+zp*xn+rp*(x2*dcos(phi)+x3*dsin(phi))
            co(2,noder)=csab(2)+zp*yn+rp*(y2*dcos(phi)+y3*dsin(phi))
            co(3,noder)=csab(3)+zp*zn+rp*(z2*dcos(phi)+z3*dsin(phi))
!
            call linkdissimilar(co,nk,ics,csab,xn,yn,zn,ncsnodes,
     &         rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,straight,
     &         lcs,nodef,ratio,nterms,rp,zp,netri,
     &         nodesonaxis,nodel,noder,ifacetet,inodface)
c            write(*,*) 'before: ',co(1,nodel),co(2,nodel),co(3,nodel)
           co(1,nodel)=co(1,noder)+rp*(x2*(1.d0-dcos(phi))-x3*dsin(phi))
           co(2,nodel)=co(2,noder)+rp*(y2*(1.d0-dcos(phi))-y3*dsin(phi))
           co(3,nodel)=co(3,noder)+rp*(z2*(1.d0-dcos(phi))-z3*dsin(phi))
c            write(*,*) 'after: ',co(1,nodel),co(2,nodel),co(3,nodel)
         else
            if(ics(node).lt.0) return
         endif
!     
!     generating the mechanical MPC's; the generated MPC's are for
!     nodal diameter 0. For other nodal diameters the MPC's are
!     changed implicitly in mastructcs and mafillsmcs
!     
      call transformatrix(csab,co(1,nodel),al)
      call transformatrix(csab,co(1,noder),ar)
!     
      do ndir=1,3
         nmpc=nmpc+1
         labmpc(nmpc)='CYCLIC              '
         if(mcs.lt.10) then
            write(labmpc(nmpc)(7:7),'(i1)') mcs
         elseif(mcs.lt.100) then
            write(labmpc(nmpc)(7:8),'(i2)') mcs
         else
            write(*,*)'*ERROR in generatecycmpcs: no more than 99'
            write(*,*)'       cyclic symmetry definitions allowed'
            stop
         endif
         ipompc(nmpc)=mpcfree 
         number=ndir-1 
!     
!     determining which direction to use for the
!     dependent side: should not occur on the dependent
!     side in another MPC and should have a nonzero
!     coefficient
!     
         do j=1,3
            number=number+1
            if(number.gt.3) number=1
            idof=7*(nodel-1)+number
            call nident(ikmpc,idof,nmpc-1,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  cycle
               endif
            endif
            if(dabs(al(number,ndir)).lt.1.d-10) cycle
            exit
         enddo
         if(j.gt.3) then
            write(*,*) '*ERROR in generatecycmpcs: cyclic MPC in node'
            write(*,*) '       ',nodel,' cannot be created: all'
            write(*,*) '       DOFs in the node are used as'
            write(*,*) '       dependent nodes in other MPCs'
            stop
         endif
         number=number-1
!     
!     updating ikmpc and ilmpc
!     
         do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
         enddo
         ikmpc(id+1)=idof
         ilmpc(id+1)=nmpc
!     
         do j=1,3
            number=number+1
            if(number.gt.3) number=1
            if(dabs(al(number,ndir)).lt.1.d-10) cycle
            nodempc(1,mpcfree)=nodel
            nodempc(2,mpcfree)=number
            coefmpc(mpcfree)=al(number,ndir)
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
               stop
            endif
         enddo
         do j=1,3
            number=number+1
            if(number.gt.3) number=1
            if(dabs(ar(number,ndir)).lt.1.d-10) cycle
            if(.not.interpolation) then
               nodempc(1,mpcfree)=noder
               nodempc(2,mpcfree)=number
               coefmpc(mpcfree)=-ar(number,ndir)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
                  stop
               endif
            else
               do k=1,nterms
                  nodempc(1,mpcfree)=nodef(k)
                  nodempc(2,mpcfree)=number
                  coefmpc(mpcfree)=-ar(number,ndir)*ratio(k)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) '*ERROR in generatecycmpcs: increase nmp
     &c_'
                     stop
                  endif
               enddo
            endif
         enddo
         nodempc(3,mpcfreeold)=0
      enddo
!     
!     generating the thermal MPC's; the generated MPC's are for nodal
!     diameter 0. 
!     
      nmpc=nmpc+1
      labmpc(nmpc)='CYCLIC              '
      if(mcs.lt.10) then
         write(labmpc(nmpc)(7:7),'(i1)') mcs
      elseif(mcs.lt.100) then
         write(labmpc(nmpc)(7:8),'(i2)') mcs
      else
         write(*,*)'*ERROR in generatecycmpcs: no more than 99'
         write(*,*)'       cyclic symmetry definitions allowed'
         stop
      endif
      ipompc(nmpc)=mpcfree 
      idof=7*(nodel-1)
      call nident(ikmpc,idof,nmpc-1,id)
      if(id.gt.0) then
         if(ikmpc(id).eq.idof) then
            write(*,*) '*ERROR in generatecycmpcs: temperature'
            write(*,*) '       in node',nodel,'is already used'
            stop
         endif
      endif
!     
!     updating ikmpc and ilmpc
!     
      do j=nmpc,id+2,-1
         ikmpc(j)=ikmpc(j-1)
         ilmpc(j)=ilmpc(j-1)
      enddo
      ikmpc(id+1)=idof
      ilmpc(id+1)=nmpc
!     
      nodempc(1,mpcfree)=nodel
      nodempc(2,mpcfree)=0
      coefmpc(mpcfree)=1.d0
      mpcfree=nodempc(3,mpcfree)
      if(mpcfree.eq.0) then
         write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
         stop
      endif
      if(.not.interpolation) then
         nodempc(1,mpcfree)=noder
         nodempc(2,mpcfree)=0
         coefmpc(mpcfree)=-1.d0
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         if(mpcfree.eq.0) then
            write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
            stop
         endif
      else
         do k=1,nterms
            nodempc(1,mpcfree)=nodef(k)
            nodempc(2,mpcfree)=0
            coefmpc(mpcfree)=-ratio(k)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
               stop
            endif
         enddo
      endif
      nodempc(3,mpcfreeold)=0
!     
      return
      end

