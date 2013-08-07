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
      subroutine umat_gurson(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab)
!
!     calculates stiffness and stresses for a Gurson-type material
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
!                        >0: thermal effects are taken into account (triggered
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
!     mi(1)              max. # of integration points per element in the
!                        model
!     nstate_            max. # of state variables in the model
!
!     xstateini(nstate_,mi(1),# of elements)
!                        state variables at the start of the increment
!     xstate(nstate_,mi(1),# of elements)
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
!     xstate(nstate_,mi(1),# of elements)
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
      character*80 amat
!
      integer ithermal,icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &  i,j,n,nrhs,lda,ldb,ipiv(7),info
!
      real*8 elconloc(21),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  um,un,aa,f,um2,dg,ddg,ep(6),s(6),de(6),h,r(7),rg,residual,
     &  dij(6,6),a(4,4),constant,constant1,constant2,fv(4),bb(4,2),
     &  ra(4),fa(4),gg,tan(6,6),hv(4),cg,ak,ehydro,af(7),ar(7),edev(6),
     &  ainv(6,6),acp(6,6),d,al,dij4(4,4),xd(6),s23,s32,c1,c2,c3,fc,ff,
     &  fn,en,sn,r0,um3,q0,f0,el(6),p,ds,xn(6),svm,svm2,dp,dsvm,dei,den,
     &  aki,arg,cco,sco,c1c2,c1c2f,c1c2fp,q2,q3,q4,depr,devm,dhs,dhss,
     &  dhsq,dhp,dhpp,dhpq,dhpf,dhq,dhf,b11,b12,b21,b22,q,ep0(6),fnsn,
     &  eplm
!
      data dij /1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &          0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     &          0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,
     &          0.d0,0.d0,0.d0,.5d0,0.d0,0.d0,
     &          0.d0,0.d0,0.d0,0.d0,.5d0,0.d0,
     &          0.d0,0.d0,0.d0,0.d0,0.d0,.5d0/
      data dij4 /1.d0,0.d0,0.d0,0.d0,
     &           0.d0,1.d0,0.d0,0.d0,
     &           0.d0,0.d0,1.d0,0.d0,
     &           0.d0,0.d0,0.d0,1.d0/
      data xd /1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
!
      s23=dsqrt(2.d0/3.d0)
      s32=dsqrt(1.5d0)
!
!     material constants
!
      um=elconloc(1)
      un=elconloc(2)
      c1=elconloc(3)
      c2=elconloc(4)
      c3=elconloc(5)
      fc=elconloc(6)
      ff=elconloc(7)
      fn=elconloc(8)
      en=elconloc(9)
      sn=elconloc(10)
      r0=elconloc(11)
!
      um2=2.d0*um
      um3=1.d0/(3.d0*um)
!
!     internal variables at the start of the increment
!
!     yield stress of the fully dense material
!
      q0=xstateini(1,iint,iel)
!
!     plastic strain
!
      do i=1,6
         ep0(i)=xstateini(1+i,iint,iel)
      enddo
!
!     void volume fraction
!
      f0=xstateini(8,iint,iel)
!
!     elastic strain in the assumption that no plasticity
!     occurs in the present increment
!
      do i=1,6
         el(i)=emec(i)-ep0(i)
      enddo
!
!     hydrostatic strain
!
      ehydro=(el(1)+el(2)+el(3))/3.d0
!
!     deviatoric strain
!
      do i=1,3
         edev(i)=el(i)-ehydro
      enddo
      do i=4,6
         edev(i)=el(i)
      enddo
!
!     deviatoric trial stress
!
      do i=1,6
         s(i)=um2*edev(i)
      enddo
!
!     trial pressure
!
      ak=3.d0*(2.d0*um*(1.d0+un))/(3.d0*(1.d0-2.d0*un))
      p=ak*ehydro
!
!     radial vector
!
      ds=dsqrt(s(1)*s(1)+s(2)*s(2)+s(3)*s(3)+
     &     2.d0*(s(4)*s(4)+s(5)*s(5)+s(6)*s(6)))
      do i=1,6
         xn(i)=s(i)/ds
      enddo
!
!     von Mises stress
!
      svm=s32*ds
!
!     yield criterion
!
      h=(svm*svm)/(q0*q0)+2.d0*c1*f0*dcosh(3.d0*p*c2/(2.d0*q0))
     & -(1.d0+c3*f0*f0)
!
      if(h.le.0.d0) then
         do i=1,3
            stre(i)=s(i)+ak*ehydro
         enddo
         do i=4,6
            stre(i)=s(i)
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
!     plasticity; initialization of the fields
!
      dg=0.d0
      ddg=0.d0
      dp=0.d0
      dsvm=0.d0
      q=q0
      f=f0
!
!     total strain increment
!
      do i=1,6
         de(i)=emec(i)-emec0(i)
      enddo
      dei=de(1)+de(2)+de(3)
      den=de(1)*xn(1)+de(2)*xn(2)+de(3)*xn(3)+
     &    2.d0*(de(4)*xn(4)+de(5)*xn(5)+de(6)*xn(6))
!
!     auxiliary variables
!
      aki=1.d0/ak
      fnsn=fn/(sn*dsqrt(8.d0*datan(1.d0)))
!
!     starting the loop to determine the consistency parameter
!
      do
!
!     inverse of the tangent hardening modulus (to complete!)
!
         d=1.
         eplm=1.
!
!        void nucleation constant
!
         aa=fnsn*dexp(-((eplm-en)/sn)**2/2.d0)
!
!        auxiliary variables
!
         arg=3.d0*c2*p/(2.d0*q)
         cco=dcosh(arg)
         sco=dsinh(arg)
         c1c2=c1*c2
         c1c2f=c1c2*f
         q2=q*q
!
         depr=dei+dp*aki
         devm=dsvm*um3-s23*den
         svm2=svm*svm
!
!        determining the residuals
!
         r(1)=-depr+dg*3.d0*c1c2f*sco/q
         r(2)=-devm+dg*2.d0*svm/q2
         r(3)=(1-f)*d*q*(q0-q)+devm*svm-depr*p
         r(4)=f0-f+(1.d0-f)*depr-aa*d*(q-q0)
!
         rg=(svm2)/q2+2.d0*c1*f*cco+1.d0+c3*f*f
!
!        check convergence
!
         residual=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)+r(4)*r(4)+rg*rg
         if((residual.le.1.d-10).or.(dabs(ddg).lt.1.d-3*dabs(dg))) exit
!
!        auxiliary variables
!
         c1c2fp=c1c2f*p
         q3=q2*q
         q4=q3*q
!
!        derivatives of the yield function
!         
         dhs=2.d0*svm/q2
         dhss=2.d0/q2
         dhsq=-4.d0*svm/q
!
         dhp=3.d0*c1c2f*sco/q
         dhpp=9.d0*c1c2f*c2*cco/(2.d0*q2)
         dhpq=-3.d0*c1c2f*sco/q2-9.d0*c1c2fp*c2*cco/(2.d0*q3)
         dhpf=3.d0*c1c2*sco/q
!
         dhq=-2.d0*svm2/q3-3.d0*c1c2fp*sco/q2
!
         dhf=2.d0*c1*cco+2.d0*c3*f
!
         a(1,1)=dg*dhpp-1.d0*aki
         a(1,2)=0.d0
         a(1,3)=dg*dhpq
         a(4,4)=dg*dhpf
         a(2,1)=0.d0
         a(2,2)=dg*dhss+um3
         a(2,3)=dg*dhsq
         a(2,4)=0.d0
         a(3,1)=-2.d0*p*aki
         a(3,2)=-2.d0*svm*um3
         a(3,3)=-(1.d0-f)*d*(2.d0*q-q0)
         a(3,4)=d*q*(q-q0)
         a(4,1)=-f*aki
         a(4,2)=0.d0
         a(4,3)=-aa*d
         a(4,4)=-depr-1.d0
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
!        vector f
!
         fv(1)=dhp
         fv(2)=dhs
         fv(3)=0.d0
         fv(4)=0.d0
!
!        solving for A:R and A:F
!
         do i=1,4
            bb(i,1)=r(i)
            bb(i,2)=fv(i)
         enddo
         n=4
         nrhs=2
         lda=4
         ldb=4
         call dgesv(n,nrhs,a,lda,ipiv,bb,ldb,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in umat_gurson:'
            write(*,*) '       singular system of equations'
            stop
         endif
         do i=1,4
            ra(i)=bb(i,1)
            fa(i)=bb(i,2)
         enddo
!
!        determination of vector field h and the constant cg
!
         hv(1)=dhp
         hv(2)=dhs
         hv(3)=dhq
         hv(4)=dhf
         cg=0.d0
!
!        calculating ddg
!
         gg=(hv(1)*af(1)+hv(2)*af(2)+hv(3)*af(3)+hv(4)*af(4))-cg
         ddg=(rg-(hv(1)*ar(1)+hv(2)*ar(2)+hv(3)*ar(3)+hv(4)*ar(4)))/gg
!
         dg=dg+ddg
!
!        update p,svm,q and f
!
         dp=-ar(1)-ddg*af(1)
         p=p+dp
         dsvm=-ar(2)-ddg*af(2)
         svm=svm+dsvm
         q=q-ar(3)-ddg*af(3)
         f=f-ar(4)-ddg*af(4)
!
      enddo
!
!     convergence: calculate the plastic strain
!
      devm=s23*dei-dsvm*um3
      depr=dei+dp*aki
      constant1=s32*devm
      constant2=depr/3.d0
      do i=1,6
         ep(i)=ep0(i)+constant1*xn(i)
      enddo
      do i=1,3
         ep(i)=ep(i)+constant2
      enddo
!
      if(icmd.ne.3) then
!
!     tangent matrix
!
         nrhs=1
         call dgesv(n,nrhs,ainv,lda,ipiv,hv,ldb,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in umat_gurson:'
            write(*,*) '       singular system of equations'
            stop
         endif
!
         do i=1,2
            do j=1,2
               tan(i,j)=acp(i,1)*(dij4(1,j)-fv(1)*hv(j)/gg)+
     &                  acp(i,2)*(dij4(2,j)-fv(2)*hv(j)/gg)+
     &                  acp(i,3)*(dij4(3,j)-fv(3)*hv(j)/gg)+
     &                  acp(i,4)*(dij4(4,j)-fv(4)*hv(j)/gg)
            enddo
         enddo
!
         constant=s23*svm+um2/ds
         b11=-tan(1,1)-constant/3.d0
         b12=-s23*tan(1,2)
         b21=s23*tan(2,1)
         b22=s23*s23*tan(2,2)-constant
!
         do i=1,6
            do j=1,6
               tan(i,j)=dij(i,j)*constant+b11*xd(i)*xd(j)+
     &                  b12*xn(i)*xd(j)+b21*xd(i)*xn(j)+
     &                  b22*xn(i)*xn(j)
            enddo
         enddo
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
!     yield stress of the fully dense material
!
      xstate(1,iint,iel)=q
!
!     plastic strain
!
      do i=1,6
         xstate(1+i,iint,iel)=ep(i)
      enddo
!
!     void volume fraction
!
      xstate(8,iint,iel)=f
!
      return
      end


