!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
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
      subroutine cycsymmods(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,nr,nz,
     &  rcs0,zcs0,ncs_,cs,labmpc,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ntie,mcs,lprev,ithermal,rcscg,rcs0cg,zcscg,
     &  zcs0cg,nrcg,nzcg,jcs,kontri,straight,ne,ipkon,kon,
     &  lakon,lcs,ifacetet,inodface,ipoinpc,maxsectors,
     &  trab,ntrans,ntrans_,jobnamec,vold,cfd,mi)
!
!     reading the input deck: *CYCLIC SYMMETRY MODEL
!
!     several cyclic symmetry parts can be defined for one and the
!     same model; for each part there must be a *CYCLIC SYMMETRY MODEL
!     card
!
!     cs(1,mcs): # segments in 360 degrees
!     cs(2,mcs): minimum node diameter
!     cs(3,mcs): maximum node diameter
!     cs(4,mcs): # nodes on the independent side
!     cs(5,mcs): # sectors to be plotted
!     cs(6,mcs) up to cs(12,mcs): cyclic symmetry axis
!     cs(13,mcs): number of the element set
!     cs(14,mcs): sum of previous independent nodes
!     cs(15,mcs): cos(angle); angle = 2*pi/cs(1,mcs)
!     cs(16,mcs): sin(angle)
!     cs(17,mcs): number of tie constraint
!
      implicit none
!
      logical triangulation,calcangle,nodesonaxis,check
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*80 tie
      character*81 set(*),depset,indepset,tieset(3,*),elset
      character*132 textpart(16),jobnamec(*)
!
      integer istartset(*),iendset(*),ialset(*),ipompc(*),
     &  nodempc(3,*),
     &  nset,istep,istat,n,key,i,j,k,nk,nmpc,nmpc_,mpcfree,ics(*),
     &  nr(*),nz(*),jdep,jindep,l,noded,ikmpc(*),ilmpc(*),lcs(*),
     &  kflag,node,ncsnodes,ncs_,iline,ipol,inl,ipoinp(2,*),nneigh,
     &  inp(3,*),itie,iset,ipos,mcs,lprev,ntie,ithermal,ncounter,
     &  nrcg(*),nzcg(*),jcs(*),kontri(3,*),ne,ipkon(*),kon(*),nodei,
     &  ifacetet(*),inodface(*),ipoinpc(0:*),maxsectors,
     &  noden(2),ntrans,ntrans_,cfd,mi(*)
!
      real*8 tolloc,co(3,*),coefmpc(*),rcs(*),zcs(*),rcs0(*),zcs0(*),
     &  csab(7),xn,yn,zn,dd,xap,yap,zap,tietol(2,*),cs(17,*),xsectors,
     &  gsectors,x3,y3,z3,phi,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),
     &  straight(9,*),x1,y1,z1,x2,y2,z2,zp,rp,dist,trab(7,*),
     &  vold(0:mi(2),*),calculated_angle,user_angle
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '       *CYCLIC SYMMETRY MODEL should'
         write(*,*) '       be placed before all step definitions'
         stop
      endif
!
      check=.true.
      gsectors=1
      elset='
     &                      '
      tie='
     &                   '
      do i=2,n
         if(textpart(i)(1:2).eq.'N=') then
            read(textpart(i)(3:22),'(f20.0)',iostat=istat) xsectors
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:8).eq.'CHECK=NO') then
            check=.false.
         elseif(textpart(i)(1:7).eq.'NGRAPH=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) gsectors
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:4).eq.'TIE=') then
            read(textpart(i)(5:84),'(a80)',iostat=istat) tie
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            read(textpart(i)(7:86),'(a80)',iostat=istat) elset
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &                 '*WARNING reading *CYCLIC SYMMETRY MODEL:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
      if(xsectors.le.0) then
         write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '       the required parameter N'
         write(*,*) '       is lacking on the *CYCLIC SYMMETRY MODEL'
         write(*,*) '       keyword card or has a value <=0'
         stop
      endif
      if(gsectors.lt.1) then
         write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '         cannot plot less than'
         write(*,*) '         one sector: one sector will be plotted'
         gsectors=1
      endif
      if(gsectors.gt.xsectors) then
         write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '         cannot plot more than'
         write(*,*) '         ',xsectors,'sectors;',
     &           xsectors,' sectors will'
         write(*,*) '       be plotted'
         gsectors=xsectors
      endif
!
      maxsectors=max(maxsectors,int(xsectors+0.5d0))
!
      mcs=mcs+1
      cs(2,mcs)=-0.5
      cs(3,mcs)=-0.5
      cs(14,mcs)=lprev+0.5
!
!     determining the tie constraint
!
      itie=0
      do i=1,ntie
         if((tieset(1,i)(1:80).eq.tie).and.
     &      (tieset(1,i)(81:81).ne.'C').and.
     &      (tieset(1,i)(81:81).ne.'T')) then
            itie=i
            exit
         endif
      enddo
      if(itie.eq.0) then
         if(ntie.eq.1) then
            itie=1
         else
            write(*,*)
     &                 '*ERROR reading *CYCLIC SYMMETRY MODEL:'
            write(*,*) '       tie constraint is nonexistent'
            call inputerror(inpc,ipoinpc,iline)
         endif
      endif
!
      cs(1,mcs)=xsectors
      cs(5,mcs)=gsectors+0.5
      cs(17,mcs)=itie+0.5
      depset=tieset(2,itie)
      indepset=tieset(3,itie)
      tolloc=tietol(1,itie)
!
!     determining the element set
!
      iset=0
      if(elset.eq.'                     ') then
         write(*,*) '*INFO reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '      no element set given'
         call inputinfo(inpc,ipoinpc,iline)
      else
         do i=1,nset
            if(set(i).eq.elset) then
               iset=i
               exit
            endif
         enddo
         if(iset.eq.0) then
            write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
            write(*,*) '       element set does not'
            write(*,*) '       exist; '
            call inputerror(inpc,ipoinpc,iline)
         endif
      endif
      cs(13,mcs)=iset+0.5
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '      definition of the cyclic'
         write(*,*) '      symmetry model is not complete'
         stop
      endif
!
      ntrans=ntrans+1
      if(ntrans.gt.ntrans_) then
         write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '       increase ntrans_'
         stop
      endif
!
      do i=1,6
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) csab(i)
         trab(i,ntrans)=csab(i)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
      enddo
!
!     cyclic coordinate system
!
      csab(7)=-1.d0
!
!     marker for cyclic symmetry axis
!
      trab(7,ntrans)=2
!
!     check whether depset and indepset exist
!
      do i=1,nset
         if(set(i).eq.depset) exit
      enddo
      if(i.gt.nset) then
         write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '       surface ',depset
         write(*,*) '       has not yet been defined.' 
         stop
      endif
      jdep=i
!
      do i=1,nset
         if(set(i).eq.indepset) exit
      enddo
      if(i.gt.nset) then
         write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '       surface ',indepset
         write(*,*) '       has not yet been defined.' 
         stop
      endif
      jindep=i
!
!     unit vector along the rotation axis (xn,yn,zn)
!
      xn=csab(4)-csab(1)
      yn=csab(5)-csab(2)
      zn=csab(6)-csab(3)
      dd=dsqrt(xn*xn+yn*yn+zn*zn)
      xn=xn/dd
      yn=yn/dd
      zn=zn/dd
!
!     defining the indepset as a 2-D data field (axes: r=radial
!     coordinate, z=axial coordinate): needed to allocate a node
!     of the depset to a node of the indepset for the cyclic
!     symmetry equations
!
      l=0
      do j=istartset(jindep),iendset(jindep)
         if(ialset(j).gt.0) then
            l=l+1
            if(lprev+l.gt.ncs_) then
               write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
               write(*,*) '       increase ncs_'
               stop
            endif
            node =ialset(j)
!
            xap=co(1,node)-csab(1)
            yap=co(2,node)-csab(2)
            zap=co(3,node)-csab(3)
!
            ics(l)=node
            zcs(l)=xap*xn+yap*yn+zap*zn
            rcs(l)=dsqrt((xap-zcs(l)*xn)**2+
     &                   (yap-zcs(l)*yn)**2+
     &                   (zap-zcs(l)*zn)**2)
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               l=l+1
               if(l.gt.ncs_) then
                  write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
                  write(*,*) '       increase ncs_'
                  stop
               endif
               node=k
!
               xap=co(1,node)-csab(1)
               yap=co(2,node)-csab(2)
               zap=co(3,node)-csab(3)
!
               ics(l)=node
               zcs(l)=xap*xn+yap*yn+zap*zn
               rcs(l)=dsqrt((xap-zcs(l)*xn)**2+
     &                      (yap-zcs(l)*yn)**2+
     &                      (zap-zcs(l)*zn)**2)
            enddo
         endif
      enddo
!
      ncsnodes=l
!
!     initialization of near2d
!
      do i=1,ncsnodes
         nr(i)=i
         nz(i)=i
         rcs0(i)=rcs(i)
         zcs0(i)=zcs(i)
      enddo
      kflag=2
      call dsort(rcs,nr,ncsnodes,kflag)
      call dsort(zcs,nz,ncsnodes,kflag)
!
!     check whether a tolerance was defined. If not, a tolerance
!     is calculated as 0.5 % of the mean of the distance of every
!     independent node to its nearest neighbour
!
      if(tolloc.lt.0.d0) then
         nneigh=2
         dist=0.d0
         do i=1,ncsnodes
            nodei=ics(i)
!
            xap=co(1,nodei)-csab(1)
            yap=co(2,nodei)-csab(2)
            zap=co(3,nodei)-csab(3)
!     
            zp=xap*xn+yap*yn+zap*zn
            rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
            call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,noden,
     &            nneigh)
!
            dist=dist+dsqrt((co(1,nodei)-co(1,noden(2)))**2+
     &                     (co(2,nodei)-co(2,noden(2)))**2+
     &                     (co(3,nodei)-co(3,noden(2)))**2)
         enddo
         tolloc=1.d-10*dist/ncsnodes
         write(*,*) '*INFO reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '      no tolerance was defined'
         write(*,*) '      in the *TIE option; a tolerance of ',
     &       tolloc
         write(*,*) '      will be used'
         write(*,*)
      endif
!
!     calculating the angle between dependent and independent
!     side and check for nodes on the axis
!
!     this angle may be different from 2*pi/xsectors: in that way
!     the user can simulate fractional nodal diameters
!
!     (x2,y2,z2): unit vector on the dependent side and orthogonal
!                 to the rotation axis
!     (x3,y3,z3): unit vector on the independent side and orthogonal
!                 to the rotation axis
!     (x1,y1,z1)=(x2,y2,z2)x(x3,y3,z3)
!                points in the same direction of xn if the independent
!                side is on the clockwise side of the dependent side if
!                looking in the direction of xn
!
      calcangle=.false.
      nodesonaxis=.false.
!
      nneigh=1
      do i=istartset(jdep),iendset(jdep)
         if(ialset(i).gt.0) then
            if(i.gt.istartset(jdep)) then
               if(ialset(i).eq.ialset(i-1)) cycle
            endif
            noded=ialset(i)
!
            xap=co(1,noded)-csab(1)
            yap=co(2,noded)-csab(2)
            zap=co(3,noded)-csab(3)
!     
            zp=xap*xn+yap*yn+zap*zn
            rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
            if((.not.calcangle).and.(rp.gt.1.d-10)) then
               x2=(xap-zp*xn)/rp
               y2=(yap-zp*yn)/rp
               z2=(zap-zp*zn)/rp
            endif
!
            call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,node,
     &            nneigh)
!
            nodei=ics(node)
            if(nodei.lt.0) cycle
            if(nodei.eq.noded) then
               ics(node)=-nodei
               nodesonaxis=.true.
               cycle
            endif
!
            xap=co(1,nodei)-csab(1)
            yap=co(2,nodei)-csab(2)
            zap=co(3,nodei)-csab(3)
!     
            zp=xap*xn+yap*yn+zap*zn
            rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
            if((.not.calcangle).and.(rp.gt.1.d-10)) then
               x3=(xap-zp*xn)/rp
               y3=(yap-zp*yn)/rp
               z3=(zap-zp*zn)/rp
!
               x1=y2*z3-y3*z2
               y1=x3*z2-x2*z3
               z1=x2*y3-x3*y2
!
               phi=(x1*xn+y1*yn+z1*zn)/dabs(x1*xn+y1*yn+z1*zn)*
     &              dacos(x2*x3+y2*y3+z2*z3)
               if(check) then
                  calculated_angle=dacos(x2*x3+y2*y3+z2*z3)
                  user_angle=6.28318531d0/cs(1,mcs)
                  if(dabs(calculated_angle-user_angle)/calculated_angle
     &                 .gt.0.01d0) then
                     write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL'
                     write(*,*) '       number of segments does not'
                     write(*,*) '       agree with the geometry'
                     write(*,*) '       angle based on N:',
     &                  user_angle*57.29577951d0
                     write(*,*) '       angle based on the geometry:',
     &                       calculated_angle*57.29577951d0
                     stop
                  endif
               endif
               calcangle=.true.
            endif
!
         else
            k=ialset(i-2)
            do
               k=k-ialset(i)
               if(k.ge.ialset(i-1)) exit
               noded=k
!
               xap=co(1,noded)-csab(1)
               yap=co(2,noded)-csab(2)
               zap=co(3,noded)-csab(3)
!     
               zp=xap*xn+yap*yn+zap*zn
               rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
               if((.not.calcangle).and.(rp.gt.1.d-10)) then
                  x2=(xap-zp*xn)/rp
                  y2=(yap-zp*yn)/rp
                  z2=(zap-zp*zn)/rp
               endif
!     
               call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,node,
     &              nneigh)
!     
               nodei=ics(node)
               if(nodei.lt.0) cycle
               if(nodei.eq.noded) then
                  ics(node)=-nodei
                  nodesonaxis=.true.
                  cycle
               endif
!
               xap=co(1,nodei)-csab(1)
               yap=co(2,nodei)-csab(2)
               zap=co(3,nodei)-csab(3)
!     
               zp=xap*xn+yap*yn+zap*zn
               rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
               if((.not.calcangle).and.(rp.gt.1.d-10)) then
                  x3=(xap-zp*xn)/rp
                  y3=(yap-zp*yn)/rp
                  z3=(zap-zp*zn)/rp
!     
                  x1=y2*z3-y3*z2
                  y1=x3*z2-x2*z3
                  z1=x2*y3-x3*y2
!     
                  phi=(x1*xn+y1*yn+z1*zn)/dabs(x1*xn+y1*yn+z1*zn)*
     &              dacos(x2*x3+y2*y3+z2*z3)
                  if(check) then
                     calculated_angle=dacos(x2*x3+y2*y3+z2*z3)
                     user_angle=6.28318531d0/cs(1,mcs)
                     if(dabs(calculated_angle-user_angle)
     &                    /calculated_angle.gt.0.01d0) then
                        write(*,*) 
     &                     '*ERROR reading *CYCLIC SYMMETRY MODEL'
                        write(*,*) '       number of segments does not'
                        write(*,*) '       agree with the geometry'
                        write(*,*) '       angle based on N:',
     &                    user_angle*57.29577951d0
                        write(*,*) '       angle based on the geometry:'
     &                       ,calculated_angle*57.29577951d0
                        stop
                     endif
                  endif
                  calcangle=.true.
               endif
!     
            enddo
         endif
!
      enddo
!
!     allocating a node of the depset to each node of the indepset 
!
      ncounter=0
      triangulation=.false.
!
      do i=istartset(jdep),iendset(jdep)
         if(ialset(i).gt.0) then
            if(i.gt.istartset(jdep)) then
               if(ialset(i).eq.ialset(i-1)) cycle
            endif
            noded=ialset(i)
!
            call generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &         coefmpc,nmpc,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,
     &         nr,nz,rcs0,zcs0,labmpc,
     &         mcs,triangulation,csab,xn,yn,zn,phi,noded,
     &         ncsnodes,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,
     &         nzcg,jcs,lcs,kontri,straight,ne,ipkon,kon,lakon,
     &         ifacetet,inodface,ncounter,jobnamec,vold,cfd,mi,
     &         indepset)
!
         else
            k=ialset(i-2)
            do
               k=k-ialset(i)
               if(k.ge.ialset(i-1)) exit
               noded=k
!
               call generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &           coefmpc,nmpc,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,
     &           nr,nz,rcs0,zcs0,labmpc,
     &           mcs,triangulation,csab,xn,yn,zn,phi,noded,
     &           ncsnodes,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,
     &           nzcg,jcs,lcs,kontri,straight,ne,ipkon,kon,lakon,
     &           ifacetet,inodface,ncounter,jobnamec,vold,cfd,mi,
     &           indepset)
            enddo
         endif
!
      enddo
!
      if(ncounter.ne.0) then
         write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
         write(*,*) '        for at least one dependent'
         write(*,*) '        node in a cyclic symmetry definition no '
         write(*,*) '        independent counterpart was found'
c     next line was commented on 19/04/2012
c         stop
      endif
!
!     sorting ics
!     ics contains the master (independent) nodes
!
      kflag=1
      call isortii(ics,nr,ncsnodes,kflag)
      cs(4,mcs)=ncsnodes+0.5
      lprev=lprev+ncsnodes
!
!     check orientation of (xn,yn,zn) (important for copying of base
!     sector in arpackcs)
!
      if(phi.lt.0.d0) then
         csab(4)=2.d0*csab(1)-csab(4)
         csab(5)=2.d0*csab(2)-csab(5)
         csab(6)=2.d0*csab(3)-csab(6)
      endif
!
      do i=1,7
         cs(5+i,mcs)=csab(i)
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

