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
     &  mcs,ithermal,triangulation,csab,xn,yn,zn,phi,noded,ncsnodes,
     &  nodesonaxis,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,lcs,
     &  kontri,straight,ne,ipkon,kon,lakon,ifacetet,inodface)
!
!     generate cyclic mpc's
!
      implicit none
!     
      logical triangulation,nodesonaxis,interpolation
!     
      character*1 c
      character*3 m1,m2,m3
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*20 labmpc(*),label
!     
      integer ipompc(*),nodempc(3,*),nneigh,ne,ipkon(*),kon(*),
     &     istep,istat,n,j,k,nk,nmpc,nmpc_,mpcfree,ics(*),nterms,
     &     nr(*),nz(*),noded,nodei,ikmpc(*),ilmpc(*),kontri(3,*),
     &     number,idof,ndir,node,ncsnodes,id,mpcfreeold,ncs_,
     &     mcs,ithermal,nrcg(*),nzcg(*),jcs(*),lcs(*),nodef(8),
     &     netri,ifacetet(*),inodface(*),lathyp(3,6),inum,one,i,
     &     noden(10)
!     
      real*8 tolloc,co(3,*),coefmpc(*),rcs(*),zcs(*),rcs0(*),zcs0(*),
     &  csab(7),xn,yn,zn,xap,yap,zap,rp,zp,al(3,3),ar(3,3),cs(17,*),phi,
     &  x2,y2,z2,x3,y3,z3,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),
     &  straight(9,*),ratio(8)
!
      save netri
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
!     
c      nneigh=1
      nneigh=10
!     
      xap=co(1,noded)-csab(1)
      yap=co(2,noded)-csab(2)
      zap=co(3,noded)-csab(3)
!     
      zp=xap*xn+yap*yn+zap*zn
      rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
c      call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,node,nneigh)
      call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,noden,nneigh)
c      nodei=abs(ics(node))
      node=noden(1)
      nodei=abs(ics(noden(1)))
c      do i=1,10
c      write(*,*) 'dep node: ',i,noded,rp,zp,abs(ics(noden(i))),
c     &     rcs0(noden(i)),zcs0(noden(i))
c      enddo
!     
!     check whether node is on axis
!     
      if(nodei.eq.noded) then
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
     &     (tolloc.le.dsqrt((rp-rcs0(node))**2+(zp-zcs0(node))**2)))
     &     then
!     
!     the nodal positions on the dependent and independent
!     sides of the mpc's do no agree: interpolation is
!     necessary. 
!     
         write(*,*) '*WARNING in generatecycmpcs: no cyclic'
         write(*,*) '         symmetric partner found for'
         write(*,*) '         dependent node ',noded,'.'
         write(*,*) '         allowed tolerance:',tolloc
         write(*,*) '         best partner node number:',nodei
         write(*,*) '         actual distance in a radial plane: ',
     &           dsqrt((rp-rcs0(node))**2+(zp-zcs0(node))**2)
         write(*,*) '         Remedy: the node is connected to an'
         write(*,*) '         independent element side.'
         write(*,*)
!     
         interpolation=.true.
!     
         if(.not.triangulation) then
c            write(*,*) 'beforetria ',rcs0(101),zcs0(101)
            call triangulate(ics,rcs0,zcs0,ncsnodes,
     &           rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,
     &           straight,ne,ipkon,kon,lakon,lcs,netri,ifacetet,
     &           inodface)
            triangulation=.true.
c            write(*,*) 'aftertria ',rcs0(101),zcs0(101)
!
            open(70,file='triangulation.frd',status='unknown')
            c='C'
            m1=' -1'
            m2=' -2'
            m3=' -3'
            p0='    0'
            p1='    1'
            p2='    2'
            p3='    3'
            p7='    7'
            p9999=' 9999'
            one=1
            write(70,'(a5,a1)') p1,c
            write(70,'(a5,a1,67x,i1)') p2,c,one
            do i=1,nk
               write(70,'(a3,i10,1p,3e12.5)') m1,i,(co(j,i),j=1,3)
            enddo
            write(70,'(a3)') m3
            write(70,'(a5,a1,67x,i1)') p3,c,one
            do i=1,netri
               write(70,'(a3,i10,2a5)')m1,i,p7,p0
               write(70,'(a3,3i10)') m2,(kontri(j,i),j=1,3)
            enddo
            write(70,'(a3)') m3
            write(70,'(a5)') p9999
            close(70)
!     
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
         nodei=nk+1
         co(1,nodei)=csab(1)+zp*xn+rp*(x2*dcos(phi)+x3*dsin(phi))
         co(2,nodei)=csab(2)+zp*yn+rp*(y2*dcos(phi)+y3*dsin(phi))
         co(3,nodei)=csab(3)+zp*zn+rp*(z2*dcos(phi)+z3*dsin(phi))
!     
c            write(*,*) 'beforeli ',rcs0(101),zcs0(101)
         call linkdissimilar(co,nk,ics,csab,ncsnodes,
     &        rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,straight,
     &        lcs,nodef,ratio,nterms,rp,zp,netri,
     &        nodesonaxis,nodei,ifacetet,inodface,noded,tolloc,xn,yn,
     &        zn)
c            write(*,*) 'afterli ',rcs0(101),zcs0(101)
!     
         xap=co(1,nodei)-csab(1)
         yap=co(2,nodei)-csab(2)
         zap=co(3,nodei)-csab(3)
!     
         zp=xap*xn+yap*yn+zap*zn
         rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
         x2=(xap-zp*xn)/rp
         y2=(yap-zp*yn)/rp
         z2=(zap-zp*zn)/rp
         x3=yn*z2-y2*zn
         y3=x2*zn-xn*z2
         z3=xn*y2-x2*yn
!     
         co(1,noded)=csab(1)+zp*xn+rp*(x2*dcos(phi)-x3*dsin(phi))
         co(2,noded)=csab(2)+zp*yn+rp*(y2*dcos(phi)-y3*dsin(phi))
         co(3,noded)=csab(3)+zp*zn+rp*(z2*dcos(phi)-z3*dsin(phi))
      else
         if(ics(node).lt.0) return
!
!        moving the dependent node such that is corresponds exactly
!        to the independent node
!     
         xap=co(1,nodei)-csab(1)
         yap=co(2,nodei)-csab(2)
         zap=co(3,nodei)-csab(3)
!     
         zp=xap*xn+yap*yn+zap*zn
         rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
         x2=(xap-zp*xn)/rp
         y2=(yap-zp*yn)/rp
         z2=(zap-zp*zn)/rp
         x3=yn*z2-y2*zn
         y3=x2*zn-xn*z2
         z3=xn*y2-x2*yn
!     
         co(1,noded)=csab(1)+zp*xn+rp*(x2*dcos(phi)-x3*dsin(phi))
         co(2,noded)=csab(2)+zp*yn+rp*(y2*dcos(phi)-y3*dsin(phi))
         co(3,noded)=csab(3)+zp*zn+rp*(z2*dcos(phi)-z3*dsin(phi))
      endif
!     
!     generating the mechanical MPC's; the generated MPC's are for
!     nodal diameter 0. For other nodal diameters the MPC's are
!     changed implicitly in mastructcs and mafillsmcs
!     
      call transformatrix(csab,co(1,noded),al)
      call transformatrix(csab,co(1,nodei),ar)
!     
!     checking for latin hypercube positions in matrix al none of
!     which are zero
!     
      do inum=1,6
         if((dabs(al(lathyp(1,inum),1)).gt.1.d-3).and.
     &      (dabs(al(lathyp(2,inum),2)).gt.1.d-3).and.
     &      (dabs(al(lathyp(3,inum),3)).gt.1.d-3)) exit
      enddo
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
         number=lathyp(ndir,inum)
         idof=8*(noded-1)+number
         call nident(ikmpc,idof,nmpc-1,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
               write(*,*) '*WARNING in generatecycmpcs: cyclic MPC in no
     &de'
               write(*,*) '         ',noded,' and direction ',ndir
               write(*,*) '         cannot be created: the'
               write(*,*) '         DOF in this node is already used'
               cycle
            endif
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
            if(dabs(al(number,ndir)).lt.1.d-5) cycle
            nodempc(1,mpcfree)=noded
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
            if(dabs(ar(number,ndir)).lt.1.d-5) cycle
            if(.not.interpolation) then
               nodempc(1,mpcfree)=nodei
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
      idof=8*(noded-1)
      call nident(ikmpc,idof,nmpc-1,id)
      if(id.gt.0) then
         if(ikmpc(id).eq.idof) then
            write(*,*) '*ERROR in generatecycmpcs: temperature'
            write(*,*) '       in node',noded,'is already used'
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
      nodempc(1,mpcfree)=noded
      nodempc(2,mpcfree)=0
      coefmpc(mpcfree)=1.d0
      mpcfree=nodempc(3,mpcfree)
      if(mpcfree.eq.0) then
         write(*,*)'*ERROR in generatecycmpcs: increase nmpc_'
         stop
      endif
      if(.not.interpolation) then
         nodempc(1,mpcfree)=nodei
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

