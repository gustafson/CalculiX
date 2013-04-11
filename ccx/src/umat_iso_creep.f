!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine umat_iso_creep(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,
     &        mi,nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,
     &        orab)
!
!     calculates stiffness and stresses for an elastically isotropic
!     material with isotropic creep
!
!     icmd=3: calculates stress at mechanical strain
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
!                                        at mechanical strain
!                        3: calculate only stress at mechanical strain
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
!
      character*20 amat
!
      integer ithermal,icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien
!
      integer i
!
      real*8 elconloc(21),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,arg
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      real*8 c1,c2,c3,ep0(6),eei(6),e,un,al,am1,um2,ep(6),dg,
     &  ddg,stri(6),p,eeq0,eeq,um,dstri,c4,c5,f,df
!
      INTEGER LEXIMP,LEND,NSTATV,KSPT,KSTEP,KINC,LAYER
      REAL*8 DECRA(5),DESWA,STATEV,SERD,EC0,ESW0,DTEMP,PREDEF,DPRED,
     &  DUMMY,COORDS
!
      data c1 /0.8164965809277260d0/
      data c2 /0.6666666666666666d0/
      data leximp /1/
!
!     state variables
!
      eeq0=xstateini(1,iint,iel)
      do i=1,6
         ep0(i)=xstateini(i+1,iint,iel)
      enddo
!
!     elastic strains
!
      do i=1,6
         eei(i)=emec(i)-ep0(i)
      enddo
!
!     elastic constants
!
      e=elconloc(1)
      un=elconloc(2)
!
      um2=e/(1.d0+un)
      al=un*um2/(1.d0-2.d0*un)
      am1=al+um2
      um=um2/2.d0
!
      if(ielas.eq.1) then
!
         stre(1)=am1*eei(1)+al*(eei(2)+eei(3))
         stre(2)=am1*eei(2)+al*(eei(1)+eei(3))
         stre(3)=am1*eei(3)+al*(eei(1)+eei(2))
         stre(4)=um2*eei(4)
         stre(5)=um2*eei(5)
         stre(6)=um2*eei(6)
!
         if(icmd.ne.3) then
            stiff(1)=am1
            stiff(2)=al
            stiff(3)=am1
            stiff(4)=al
            stiff(5)=al
            stiff(6)=am1
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
         return
      endif
!
!     creep
!
      stri(1)=am1*eei(1)+al*(eei(2)+eei(3))
      stri(2)=am1*eei(2)+al*(eei(1)+eei(3))
      stri(3)=am1*eei(3)+al*(eei(1)+eei(2))
      stri(4)=um2*eei(4)
      stri(5)=um2*eei(5)
      stri(6)=um2*eei(6)
!
      p=-(stri(1)+stri(2)+stri(3))/3.d0
      do i=1,3
         stri(i)=stri(i)+p
      enddo
!
      dstri=dsqrt(stri(1)*stri(1)+stri(2)*stri(2)+stri(3)*stri(3)+
     &      2.d0*(stri(4)*stri(4)+stri(5)*stri(5)+stri(6)*stri(6)))
!
!     unit trial vector
!
      do i=1,6
         stri(i)=stri(i)/dstri
      enddo
!
      dg=0.d0
      eeq=eeq0+c1*dg
!
!     determination of the consistency parameter
!
      do
         arg=(dstri-um2*dg)/c1
         call CREEP( DECRA, DESWA, STATEV, SERD, EC0, ESW0, p, arg,
     &        t1l, DTEMP, PREDEF, DPRED, DUMMY, dtime, amat,
     &        leximp, LEND, COORDS, NSTATV, iel, iint, LAYER,
     &        KSPT, KSTEP, KINC )
         f=decra(1)
         df=decra(5)
         ddg=(c1*f-c2*dg)/(um2*df+c2)
         dg=dg+ddg
         eeq=eeq0+c1*dg
         if((ddg.lt.dg*1.d-4).or.(ddg.lt.1.d-10)) exit
      enddo
!
      do i=1,6
         ep(i)=dg*stri(i)
         eei(i)=eei(i)-ep(i)
         ep(i)=ep0(i)+ep(i)
      enddo
!
!     stress values
!
      stre(1)=am1*eei(1)+al*(eei(2)+eei(3))
      stre(2)=am1*eei(2)+al*(eei(1)+eei(3))
      stre(3)=am1*eei(3)+al*(eei(1)+eei(2))
      stre(4)=um2*eei(4)
      stre(5)=um2*eei(5)
      stre(6)=um2*eei(6)
!
!     stiffness matrix
!
      if(icmd.ne.3) then
!
         c3=um2*um2
         c4=c3*dg/dstri
         c3=c4-c3*df/(um2*df+c2)
         c5=c4/3.d0
!
         stiff(1)=am1+c3*stri(1)*stri(1)+c5-c4
         stiff(2)=al+c3*stri(1)*stri(2)+c5
         stiff(3)=am1+c3*stri(2)*stri(2)+c5-c4
         stiff(4)=al+c3*stri(1)*stri(3)+c5
         stiff(5)=al+c3*stri(2)*stri(3)+c5
         stiff(6)=am1+c3*stri(3)*stri(3)+c5-c4
         stiff(7)=0.d0+c3*stri(1)*stri(4)
         stiff(8)=0.d0+c3*stri(2)*stri(4)
         stiff(9)=0.d0+c3*stri(3)*stri(4)
         stiff(10)=um+c3*stri(4)*stri(4)-c4/2.d0
         stiff(11)=0.d0+c3*stri(1)*stri(5)
         stiff(12)=0.d0+c3*stri(2)*stri(5)
         stiff(13)=0.d0+c3*stri(3)*stri(5)
         stiff(14)=0.d0+c3*stri(4)*stri(5)
         stiff(15)=um+c3*stri(5)*stri(5)-c4/2.d0
         stiff(16)=0.d0+c3*stri(1)*stri(6)
         stiff(17)=0.d0+c3*stri(2)*stri(6)
         stiff(18)=0.d0+c3*stri(3)*stri(6)
         stiff(19)=0.d0+c3*stri(4)*stri(6)
         stiff(20)=0.d0+c3*stri(5)*stri(6)
         stiff(21)=um+c3*stri(6)*stri(6)-c4/2.d0
      endif
!
!     state variables
!
      xstate(1,iint,iel)=eeq
      do i=1,6
         xstate(i+1,iint,iel)=ep(i)
      enddo
!
      return
      end
