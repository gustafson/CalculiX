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
      subroutine umat_walker(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab)
!
!     calculates stiffness and stresses for a user defined material
!     law
!
!     icmd=3: calcutates stress at mechanical strain
!     else: calculates stress at mechanical strain and the stiffness
!           matrix
!
!     INPUT:
!
!     amat               material name
!     iel                element number
!     iint               integration point number
!
!     kode               material type (-100-#of constants entered
!                        under *USER MATERIAL): can be used for materials
!                        with varying number of constants
!
!     elconloc(21)       user defined constants defined by the keyword
!                        card *USER MATERIAL (max. 21, actual # =
!                        -kode-100), interpolated for the
!                        actual temperature t1l
!
!     emec(6)            Lagrange mechanical strain tensor (component order:
!                        11,22,33,12,13,23) at the end of the increment
!                        (thermal strains are subtracted)
!     emec0(6)           Lagrange mechanical strain tensor at the start of the
!                        increment (thermal strains are subtracted)
!     beta(6)            residual stress tensor (the stress entered under
!                        the keyword *INITIAL CONDITIONS,TYPE=STRESS)
!
!     xokl(3,3)          deformation gradient at the start of the increment
!     voj                Jacobian at the start of the increment
!     xkl(3,3)           deformation gradient at the end of the increment
!     vj                 Jacobian at the end of the increment
!
!     ithermal           0: no thermal effects are taken into account
!                        1: thermal effects are taken into account (triggered
!                        by the keyword *INITIAL CONDITIONS,TYPE=TEMPERATURE)
!     t1l                temperature at the end of the increment
!     dtime              time length of the increment
!     time               step time at the end of the current increment
!     ttime              total time at the start of the current increment
!
!     icmd               not equal to 3: calculate stress and stiffness
!                        3: calculate only stress
!     ielas              0: no elastic iteration: irreversible effects
!                        are allowed
!                        1: elastic iteration, i.e. no irreversible
!                           deformation allowed
!
!     mint_              max. # of integration points per element in the
!                        model
!     nstate_            max. # of state variables in the model
!
!     xstateini(nstate_,mint_,# of elements)
!                        state variables at the start of the increment
!     xstate(nstate_,mint_,# of elements)
!                        state variables at the end of the increment
!
!     stre(6)            Piola-Kirchhoff stress of the second kind
!                        at the start of the increment
!
!     iorien             number of the local coordinate axis system
!                        in the integration point at stake (takes the value
!                        0 if no local system applies)
!     pgauss(3)          global coordinates of the integration point
!     orab(7,*)          description of all local coordinate systems.
!                        If a local coordinate system applies the global 
!                        tensors can be obtained by premultiplying the local
!                        tensors with skl(3,3). skl is  determined by calling
!                        the subroutine transformatrix: 
!                        call transformatrix(orab(1,iorien),pgauss,skl)
!
!
!     OUTPUT:
!
!     xstate(nstate_,mint_,# of elements)
!                        updated state variables at the end of the increment
!     stre(6)            Piola-Kirchhoff stress of the second kind at the
!                        end of the increment
!     stiff(21):         consistent tangent stiffness matrix in the material
!                        frame of reference at the end of the increment. In
!                        other words: the derivative of the PK2 stress with
!                        respect to the Lagrangian strain tensor. The matrix
!                        is supposed to be symmetric, only the upper half is
!                        to be given in the same order as for a fully
!                        anisotropic elastic material (*ELASTIC,TYPE=ANISO).
!                        Notice that the matrix is an integral part of the 
!                        fourth order material tensor, i.e. the Voigt notation
!                        is not used.
!
      implicit none
!
      logical plastic,outofrange
!
      character*80 amat
!
      integer ithermal,icmd,kode,ielas,iel,iint,nstate_,mint_,iorien,
     &  i,j,n,nrhs,lda,ldb,ipiv(7),info
!
      real*8 elconloc(21),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,xstate(nstate_,mint_,*),xstateini(nstate_,mint_,*),
     &  epq0,ep0(6),bs0(6),bl0(6),b0(6),um,un,aa,cc,f,hh0,h0,h1,xm,
     &  xn,qq,tm,delta,tt,dd0,ugas,um2,dc,s0(6),xi(6),dxi,xk,theta,hd,
     &  ta,dg,ddg,ep(6),s(6),ds(6),b(6),de(6),h,ddc,cddc,sddc,dgxi,
     &  r(7),rg,xide,xibs,xibl,xmm,residual,xll,xls,dlldd,dlsdd,
     &  bs1(6),bl1(6),xls1,xll1,dij(7,7),eta(6,6),a(7,7),constant,
     &  constant1,constant2,dbdd(6),dbdg(6),fv(7),bb(7,2),um2i,
     &  dedbdd,xidbdd,xidbsdd,xidbldd,bsdbdd,bldbdd,gg,
     &  tan(6,6),dedbdg,xidbdg,xidbsdg,xidbldg,bsdbdg,bldbdg,
     &  sbd,csbd,ssbd,bl(6),bs(6),hv(7),cg,dbsdd(6),dbldd(6),
     &  dbsdg(6),dbldg(6),ak,ehydro,constant3,af(7),ar(7),edev(6),
     &  ehydro0,ainv(7,7),acp(7,7),d,d0,dbds,umxls1,umxll1,al,
     &  xlsi,xlli,dummy
!
      data dij /1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &          0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &          0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     &          0.d0,0.d0,0.d0,.5d0,0.d0,0.d0,0.d0,
     &          0.d0,0.d0,0.d0,0.d0,.5d0,0.d0,0.d0,
     &          0.d0,0.d0,0.d0,0.d0,0.d0,.5d0,0.d0,
     &          0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0/
!
      if(ithermal.eq.0) then
         write(*,*) '*ERROR in umat_walker: no temperature defined'
         stop
      endif
!
!     material constants
!
      um=elconloc(1)
      un=elconloc(2)
      aa=elconloc(3)
      cc=elconloc(4)
      f=elconloc(5)
      hh0=1.d0/elconloc(6)
      h0=elconloc(7)
      h1=elconloc(8)
      xm=elconloc(9)
      xn=elconloc(10)
      qq=elconloc(11)
      tm=elconloc(12)
!
      delta=0.035d0
      tt=tm/2.d0
      dd0=cc/100.d0
      ugas=8.314d0
      um2=2.d0*um
      um2i=1.d0/um2
      dc=1.d0/(delta*cc)
!
!     internal variables at the start of the increment
!
!     equivalent plastic strain
!
      epq0=xstateini(1,iint,iel)
!
!     plastic strain
!
      do i=1,6
         ep0(i)=xstateini(1+i,iint,iel)
      enddo
!
!     short range back stress
!
      do i=1,6
         bs0(i)=xstateini(7+i,iint,iel)
      enddo
!
!     long range back stress
!
      do i=1,6
         bl0(i)=xstateini(13+i,iint,iel)
      enddo
!
!     drag
!
      d0=xstateini(20,iint,iel)
      if(d0.le.0.d0) d0=dd0+1.d-10
!
!     total back stress
!     
      do i=1,6
         b0(i)=bs0(i)+bl0(i)
      enddo
!
!     hydrostatic strain
!
      ehydro0=(emec0(1)+emec0(2)+emec0(3))/3.d0
      ehydro=(emec(1)+emec(2)+emec(3))/3.d0
      ak=3.d0*(2.d0*um*(1.d0+un))/(3.d0*(1.d0-2.d0*un))
!
!     deviatoric strain
!
      do i=1,3
         edev(i)=emec(i)-ehydro
      enddo
      do i=4,6
         edev(i)=emec(i)
      enddo
!
!     relative stress vector
!
      do i=1,6
         s0(i)=um2*(edev(i)-ep0(i))
         xi(i)=s0(i)-b0(i)
      enddo
      dxi=dsqrt((xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3))/2.d0+
     &           xi(4)*xi(4)+xi(5)*xi(5)+xi(6)*xi(6))
      if(dxi.lt.1.d-10) then
         do i=1,3
            stre(i)=s0(i)+ak*ehydro
         enddo
         do i=4,6
            stre(i)=s0(i)
         enddo
!
         if(icmd.ne.3) then
            al=2.d0*un*um/(1.d0-2.d0*un)
            stiff(1)=al+2.d0*um
            stiff(2)=al
            stiff(3)=al+2.d0*um
            stiff(4)=al
            stiff(5)=al
            stiff(6)=al+2.d0*um
            stiff(7)=0.d0
            stiff(8)=0.d0
            stiff(9)=0.d0
            stiff(10)=um
            stiff(11)=0.d0
            stiff(12)=0.d0
            stiff(13)=0.d0
            stiff(14)=0.d0
            stiff(15)=um
            stiff(16)=0.d0
            stiff(17)=0.d0
            stiff(18)=0.d0
            stiff(19)=0.d0
            stiff(20)=0.d0
            stiff(21)=um
         endif
!
         return
      endif
!
!     plastic limit
!      
      xk=0.25d0*(cc-dd0)*(cc+dd0)*dc
!
!     plastic or viscoplastic?
!
      if(dxi.lt.xk) then
         plastic=.false.
      else
         plastic=.true.
      endif
!
!     determining the temperature-dependent constants
!
      if(t1l.ge.tm) then
         write(*,*) '*ERROR in umat_walker: temperature ',t1l
         write(*,*) '       exceeds melting point'
         stop
      elseif(t1l.ge.tt) then
         theta=dexp(-qq/(ugas*t1l))
         hd=h1
      else
         theta=dexp(-qq*(dlog(tt/t1l)+1.d0)/(ugas*tt))
         hd=h0-(h0-h1)*t1l/tt
      endif
      ta=theta*aa
!
!     initialization of the fields
!
      dg=0.d0
      ddg=0.d0
      do i=1,6
         ep(i)=ep0(i)
         s(i)=s0(i)
         de(i)=emec(i)-emec0(i)
!
         bs(i)=um2*de(i)
         bl(i)=bl0(i)+bs(i)*hh0
         bs(i)=bs0(i)+bs(i)
         b(i)=bs(i)+bl(i)
      enddo
!
!     total deviatoric strain increment
!
      constant=ehydro-ehydro0
      do i=1,3
         de(i)=de(i)-constant
      enddo
!
      d=d0
!
!     determining the ls and ll parameters
!
      ddc=(d-dd0)*dc
      xll=(cc-d)*ddc
      xls=f*xll
      xll=xll-xls
!
!     determining auxiliary variables
!
      do i=1,6
         bs1(i)=um2*de(i)-ds(i)
         bl1(i)=bl0(i)+bs1(i)*hh0
         bs1(i)=bs0(i)+bs1(i)
      enddo
      xls1=1.d0/xls
      xll1=1.d0/xll
      umxls1=um*xls1*xls1
      umxll1=um*xll1*xll1*hh0
!
!     starting the loop to determine the consistency parameter
!
      do
!
!        actual relative stress tensor
!
         do i=1,6
            xi(i)=s(i)-b(i)
         enddo
         dxi=dsqrt((xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3))/2.d0+
     &        xi(4)*xi(4)+xi(5)*xi(5)+xi(6)*xi(6))
!
!        actual value of h
!
         cddc=dcosh(ddc)
         sddc=dsinh(ddc)
         if(d.le.dd0) then
            h=hd
         else
            h=hd*(ddc/sddc)**xm
         endif
!
!        determining the residuals
!
         dgxi=dg/(2.d0*dxi)
         do i=1,6
            r(i)=ep0(i)-ep(i)+xi(i)*dgxi
         enddo
         r(7)=d0-d+h*(dg-ta*sddc**xn)
         if(.not.plastic) then
            sbd=dxi/d
            if((sbd.ge.(23.026d0-dlog(ta))/xn+0.693d0).and.
     &         (sbd.ge.10.d0)) then
               outofrange=.true.
               rg=dg-1.d10
            else
               outofrange=.false.
               ssbd=dsinh(sbd)
               rg=dg-ta*ssbd**xn
            endif
         else
            xlsi=1.d0/xls
            xlli=1.d0/xll
            xide=xi(1)*de(1)+xi(2)*de(2)+xi(3)*de(3)+
     &         2.d0*(xi(4)*de(4)+xi(5)*de(5)+xi(6)*de(6))
            xibs=xi(1)*bs(1)+xi(2)*bs(2)+xi(3)*bs(3)+
     &         2.d0*(xi(4)*bs(4)+xi(5)*bs(5)+xi(6)*bs(6))
            xibl=xi(1)*bl(1)+xi(2)*bl(2)+xi(3)*bl(3)+
     &         2.d0*(xi(4)*bl(4)+xi(5)*bl(5)+xi(6)*bl(6))
            xmm=1.d0/(4.d0*dxi-xlsi*xibs+(2.d0*dxi-xlli*xibl)*hh0)
            rg=dg-2.d0*xide*xmm
         endif
!
!        check convergence
!
         residual=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)+r(4)*r(4)+
     &     r(5)*r(5)+r(6)*r(6)+r(7)*r(7)+rg*rg
         dummy=1.e-10-residual
         if((residual.le.(1.e-10)).or.(dabs(ddg).lt.(1.d-3*dabs(dg)))) 
     &         exit
!
!        determining the derivatives of the ls and ll parameters
!        w.r.t. d
!
         dlldd=(cc+dd0-2.d0*d)*dc
         dlsdd=f*dlldd
         dlldd=dlldd-dlsdd
!
!        calculation of the eta tensor
!
         constant=1.d0/(2.d0*dxi*dxi)
         do i=1,6
            do j=1,6
               eta(i,j)=dij(i,j)-xi(i)*xi(j)*constant
            enddo
         enddo
!
!        the upper left 6x6 tensor in a
!
         dbds=1.d0+xls*xls1+xll*xll1*hh0
         constant2=dgxi*dbds
         do i=1,6
            do j=1,6
               a(i,j)=dij(i,j)*um2i+eta(i,j)*constant2
            enddo
         enddo
!
!        derivative of vector B with respect to scalar D
!
         constant1=umxls1*dg*dlsdd
         constant2=umxll1*dg*dlldd
         do i=1,6
            dbsdd(i)=bs1(i)*constant1
            dbldd(i)=bl1(i)*constant2
            dbdd(i)=dbsdd(i)+dbldd(i)
         enddo
         xidbdd=xi(1)*dbdd(1)+xi(2)*dbdd(2)+xi(3)*dbdd(3)+
     &        2.d0*(xi(4)*dbdd(4)+xi(5)*dbdd(5)+xi(6)*dbdd(6))
!
         do i=1,6
            a(i,7)=-(eta(i,1)*dbdd(1)+eta(i,2)*dbdd(2)+
     &              eta(i,3)*dbdd(3)+2.d0*(eta(i,4)*dbdd(4)+
     &              eta(i,5)*dbdd(5)+eta(i,6)*dbdd(6)))*dgxi
            a(7,i)=0.d0
         enddo
!
!        diagonal element in the lower right corner of a
!
         a(7,7)=-1.d0-h*ta*xn*dc*sddc**(xn-1)*cddc
         if(d.gt.dd0) then
            a(7,7)=a(7,7)+hd*xm*(ddc/sddc)**(xm-1)*
     &             (sddc*dc-ddc*cddc*dc)/(sddc*sddc)*
     &             (dg-ta*sddc**xn)
         endif
!
!        copying a
!
         do i=1,7
            do j=1,7
               ainv(i,j)=a(j,i)
               acp(i,j)=a(i,j)
            enddo
         enddo
!
!        -1 * derivative of vector B w.r.t. the scalar dg
!
         constant1=umxls1*xls
         constant2=umxll1*xll
         do i=1,6
            dbsdg(i)=bs1(i)*constant1
            dbldg(i)=bl1(i)*constant2
            dbdg(i)=dbsdg(i)+dbldg(i)
         enddo
         xidbdg=xi(1)*dbdg(1)+xi(2)*dbdg(2)+xi(3)*dbdg(3)+
     &        2.d0*(xi(4)*dbdg(4)+xi(5)*dbdg(5)+xi(6)*dbdg(6))
!
!        vector f
!
         constant=1.d0/(2.d0*dxi)
         do i=1,6
            fv(i)=(eta(i,1)*dbdg(1)+eta(i,2)*dbdg(2)+
     &             eta(i,3)*dbdg(3)+2.d0*(eta(i,4)*dbdg(4)+
     &             eta(i,5)*dbdg(5)+eta(i,6)*dbdg(6)))*dgxi
     &           +xi(i)*constant
         enddo
         fv(7)=h
!
!        solving for A:R and A:F
!
         do i=1,7
            bb(i,1)=r(i)
            bb(i,2)=fv(i)
         enddo
         n=7
         nrhs=2
         lda=7
         ldb=7
         call dgesv(n,nrhs,a,lda,ipiv,bb,ldb,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in umat_walker:'
            write(*,*) '       singular system of equations'
            stop
         endif
         do i=1,7
            ar(i)=bb(i,1)
            af(i)=bb(i,2)
         enddo
!
!        determination of vector field h and the constant cg
!
         if(.not.plastic) then
!
!           auxiliary constants
!
            if(outofrange) then
               constant1=1.d10*xn
            else
               csbd=dcosh(sbd)
               constant1=ta*xn*ssbd**(xn-1)*csbd
            endif
            constant2=1.d0/(2.d0*d*dxi)
            constant=-constant1*dbds*constant2
!
            do i=1,6
               hv(i)=xi(i)*constant
            enddo
            hv(7)=constant1*(xidbdd*constant2+sbd/d)
            cg=1.d0+constant1*xidbdg*constant2
         else
!
!           auxiliary constants
!
            dedbdd=de(1)*dbdd(1)+de(2)*dbdd(2)+de(3)*dbdd(3)+
     &         2.d0*(de(4)*dbdd(4)+de(5)*dbdd(5)+de(6)*dbdd(6))
            xidbsdd=xi(1)*dbsdd(1)+xi(2)*dbsdd(2)+xi(3)*dbsdd(3)+
     &         2.d0*(xi(4)*dbsdd(4)+xi(5)*dbsdd(5)+xi(6)*dbsdd(6))
            xidbldd=xi(1)*dbldd(1)+xi(2)*dbldd(2)+xi(3)*dbldd(3)+
     &         2.d0*(xi(4)*dbldd(4)+xi(5)*dbldd(5)+xi(6)*dbldd(6))
            bsdbdd=bs(1)*dbdd(1)+bs(2)*dbdd(2)+bs(3)*dbdd(3)+
     &         2.d0*(bs(4)*dbdd(4)+bs(5)*dbdd(5)+bs(6)*dbdd(6))
            bldbdd=bl(1)*dbdd(1)+bl(2)*dbdd(2)+bl(3)*dbdd(3)+
     &         2.d0*(bl(4)*dbdd(4)+bl(5)*dbdd(5)+bl(6)*dbdd(6))
!
            dedbdg=de(1)*dbdg(1)+de(2)*dbdg(2)+de(3)*dbdg(3)+
     &         2.d0*(de(4)*dbdg(4)+de(5)*dbdg(5)+de(6)*dbdg(6))
            xidbsdg=xi(1)*dbsdg(1)+xi(2)*dbsdg(2)+xi(3)*dbsdg(3)+
     &         2.d0*(xi(4)*dbsdg(4)+xi(5)*dbsdg(5)+xi(6)*dbsdg(6))
            xidbldg=xi(1)*dbldg(1)+xi(2)*dbldg(2)+xi(3)*dbldg(3)+
     &         2.d0*(xi(4)*dbldg(4)+xi(5)*dbldg(5)+xi(6)*dbldg(6))
            bsdbdg=bs(1)*dbdg(1)+bs(2)*dbdg(2)+bs(3)*dbdg(3)+
     &         2.d0*(bs(4)*dbdg(4)+bs(5)*dbdg(5)+bs(6)*dbdg(6))
            bldbdg=bl(1)*dbdg(1)+bl(2)*dbdg(2)+bl(3)*dbdg(3)+
     &         2.d0*(bl(4)*dbdg(4)+bl(5)*dbdg(5)+bl(6)*dbdg(6))
!
            constant=dbds/dxi
            constant1=-2.d0*dbds*xmm
            constant2=2.d0*xide*xmm**2
            constant3=constant2*((2.d0+hh0)*constant
     &          +xls1+xll1*hh0)
!            
            do i=1,6
               hv(i)=constant1*de(i)+constant3*xi(i)+
     &               constant2*dbds*(xlsi*bs(i)+xlli*bl(i)*hh0)
            enddo
            hv(7)=2.d0*xmm*dedbdd+constant2*(-(2.d0+hh0)*xidbdd/dxi
     &           -dlsdd*xibs*xlsi*xlsi-dlldd*xibl*xlli*xlli*hh0
     &           -xlsi*xidbsdd-xlli*xidbldd*hh0
     &           +xlsi*bsdbdd+xlli*bldbdd*hh0)
!
            cg=1.d0+2.d0*dedbdg*xmm+constant2*
     &         (-(2.d0+hh0)*xidbdg/dxi
     &         -xlsi*xidbsdg-xlli*xidbldg*hh0
     &         +xlsi*bsdbdg+xlli*bldbdg*hh0)
!
         endif
!
!        calculating ddg
!
         gg=(hv(1)*af(1)+hv(2)*af(2)+hv(3)*af(3)
     &        +hv(4)*af(4)+hv(5)*af(5)+hv(6)*af(6)+hv(7)*af(7))-cg
         ddg=(rg-(hv(1)*ar(1)+hv(2)*ar(2)+hv(3)*ar(3)
     &        +hv(4)*ar(4)+hv(5)*ar(5)+hv(6)*ar(6)+hv(7)*ar(7)))/gg
!
         dg=dg+ddg
!
!        updated value of the stress
!         
         do i=1,3
            s(i)=s(i)-ar(i)-ddg*af(i)
         enddo
         do i=4,6
            s(i)=s(i)-0.5d0*(ar(i)+ddg*af(i))
         enddo
!
!        updated value of the drag
!
         d=d-ar(7)-ddg*af(7)
!
!        updated value of the plastic strain
!
         do i=1,6
            ep(i)=edev(i)-s(i)/um2
         enddo
!
!        determining the ls and ll parameters
!
         ddc=(d-dd0)*dc
         xll=(cc-d)*ddc
         xls=f*xll
         xll=xll-xls
!
!        determining auxiliary variables
!
         do i=1,6
            bs1(i)=um2*de(i)-ds(i)
            bl1(i)=bl0(i)+bs1(i)*hh0
            bs1(i)=bs0(i)+bs1(i)
         enddo
         xls1=1.d0/(xls+um*dg)
         xll1=1.d0/(xll+um*dg*hh0)
         umxls1=um*xls1*xls1
         umxll1=um*xll1*xll1*hh0
!
!        calculate the back stress
!
         do i=1,6
            bs(i)=xls*bs1(i)*xls1
            bl(i)=xll*bl1(i)*xll1
            b(i)=bs(i)+bl(i)
         enddo
!
      enddo
!
!     convergence: calculate the stress
!
      do i=1,3
         stre(i)=s(i)+ak*ehydro
      enddo
      do i=4,6
         stre(i)=s(i)
      enddo
!
      if(icmd.ne.3) then
!
!     tangent matrix
!
         nrhs=1
         call dgesv(n,nrhs,ainv,lda,ipiv,hv,ldb,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in umat_walker:'
            write(*,*) '       singular system of equations'
            stop
         endif
!
         do i=1,6
            do j=1,6
               tan(i,j)=acp(i,1)*(dij(1,j)-fv(1)*hv(j)/gg)+
     &                  acp(i,2)*(dij(2,j)-fv(2)*hv(j)/gg)+
     &                  acp(i,3)*(dij(3,j)-fv(3)*hv(j)/gg)+
     &                  acp(i,4)*(dij(4,j)-fv(4)*hv(j)/gg)+
     &                  acp(i,5)*(dij(5,j)-fv(5)*hv(j)/gg)+
     &                  acp(i,6)*(dij(6,j)-fv(6)*hv(j)/gg)+
     &                  acp(i,7)*(dij(7,j)-fv(7)*hv(j)/gg)
               if(j.gt.3) tan(i,j)=0.5d0*tan(i,j)
            enddo
         enddo
!
!
!        symmatrizing the stiffness matrix
!
         stiff(1)=tan(1,1)
         stiff(2)=(tan(1,2)+tan(2,1))/2.d0
         stiff(3)=tan(2,2)
         stiff(4)=(tan(1,3)+tan(3,1))/2.d0
         stiff(5)=(tan(2,3)+tan(3,2))/2.d0
         stiff(6)=tan(3,3)
         stiff(7)=(tan(1,4)+tan(4,1))/2.d0
         stiff(8)=(tan(2,4)+tan(4,2))/2.d0
         stiff(9)=(tan(3,4)+tan(4,3))/2.d0
         stiff(10)=tan(4,4)
         stiff(11)=(tan(1,5)+tan(5,1))/2.d0
         stiff(12)=(tan(2,5)+tan(5,2))/2.d0
         stiff(13)=(tan(3,5)+tan(5,3))/2.d0
         stiff(14)=(tan(4,5)+tan(5,4))/2.d0
         stiff(15)=tan(5,5)
         stiff(16)=(tan(1,6)+tan(6,1))/2.d0
         stiff(17)=(tan(2,6)+tan(6,2))/2.d0
         stiff(18)=(tan(3,6)+tan(6,3))/2.d0
         stiff(19)=(tan(4,6)+tan(6,4))/2.d0
         stiff(20)=(tan(5,6)+tan(6,5))/2.d0
         stiff(21)=tan(6,6)
!
      endif
!
!     internal variables at the end of the increment
!
!     equivalent plastic strain
!
      xstate(1,iint,iel)=epq0+dg
!
!     plastic strain
!
      do i=1,6
         xstate(1+i,iint,iel)=ep(i)
      enddo
!
!     short range back stress
!
      do i=1,6
         xstate(7+i,iint,iel)=bs(i)
      enddo
!
!     long range back stress
!
      do i=1,6
         xstate(13+i,iint,iel)=bl(i)
      enddo
!
!     drag
!
      xstate(20,iint,iel)=d
!
      return
      end


