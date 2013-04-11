!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine resultstherm(co,kon,ipkon,lakon,ne,v,
     &  elcon,nelcon,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &  ntmat_,t0,iperturb,fn,
     &  iout,qa,vold,ipompc,nodempc,coefmpc,nmpc,
     &  dtime,time,ttime,plicon,nplicon,xstateini,xstiff,xstate,npmat_,
     &  matname,mi,ncmat_,nstate_,cocon,ncocon,
     &  qfx,ikmpc,ilmpc,istep,iinc,springarea,
     &  calcul_fn,calcul_qa,nal,nea,neb)
!
!     calculates the heat flux and the material tangent at the integration
!     points and the internal concentrated flux at the nodes
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(20),iperm(20),ikmpc(*),ilmpc(*),mi(*),
     &  nelcon(2,*),nrhcon(*),ielmat(mi(3),*),ielorien(mi(3),*),
     &  ntmat_,ipkon(*),ipompc(*),nodempc(3,*),
     &  ncocon(2,*),iflag,nshcon,istep,iinc,mt,ne,mattyp,
     &  i,j,k,m1,jj,i1,m3,indexe,nope,norien,iperturb(*),iout,
     &  nal,nmpc,kode,imat,mint3d,iorien,istiff,ncmat_,nstate_,
     &  nplicon(0:ntmat_,*),npmat_,calcul_fn,calcul_qa,nea,neb
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,20),
     &  xl(3,20),vl(0:mi(2),20),elcon(0:ncmat_,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),qfx(3,mi(1),*),orab(7,*),
     &  rho,fn(0:mi(2),*),tnl(9),timeend(2),q(0:mi(2),20),
     &  vkl(0:3,3),t0(*),vold(0:mi(2),*),coefmpc(*),
     &  springarea(2,*),elconloc(21),cocon(0:6,ntmat_,*),
     &  shcon,sph,c1,xi,et,ze,xsj,qa(3),t0l,t1l,dtime,
     &  weight,pgauss(3),coconloc(6),qflux(3),time,ttime,
     &  t1lold,plicon(0:2*npmat_,ntmat_,*),xstiff(27,mi(1),*),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      include "gauss.f"
!
c      data iflag /3/
c      data iperm /5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/
      iflag=3
      iperm=(/5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/)
!
      mt=mi(2)+1
!
!     calculation of temperatures and thermal flux
!
      nal=0
!
      do i=nea,neb
!
         if(ipkon(i).lt.0) cycle
         imat=ielmat(1,i)
         amat=matname(imat)
         if(norien.gt.0) then
            iorien=ielorien(1,i)
         else
            iorien=0
         endif
!
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(1:1).eq.'E') then
c            read(lakon(i)(8:8),'(i1)') nope
            nope=ichar(lakon(i)(8:8))-48
!
!           local contact spring number
!
            if(lakon(i)(7:7).eq.'C') konl(nope+1)=kon(indexe+nope+1)
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            if(lakon(i)(6:7).eq.'RA') then
               mint3d=4
            else
               mint3d=8
            endif
         elseif(lakon(i)(4:4).eq.'2') then
            mint3d=27
         elseif(lakon(i)(4:5).eq.'10') then
            mint3d=4
         elseif(lakon(i)(4:4).eq.'4') then
            mint3d=1
         elseif(lakon(i)(4:5).eq.'15') then
            mint3d=9
         elseif(lakon(i)(4:4).eq.'6') then
            mint3d=2
         elseif(lakon(i)(1:1).eq.'E') then
            mint3d=0
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
            enddo
            vl(0,j)=v(0,konl(j))
         enddo
c         write(*,*) ((xl(k,j),k=1,3),j=1,nope)
c         write(*,*) ((vl(k,j),k=1,3),j=1,nope)
!
!        q contains the nodal forces per element; initialisation of q
!
         if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &          then
            do m1=1,nope
               q(0,m1)=fn(0,konl(m1))
            enddo
         endif
!
!        calculating the concentrated flux for the contact elements
!
         if(mint3d.eq.0) then
!
            lakonl=lakon(i)
!
!           spring elements (including contact springs)
!     
            if(lakonl(2:2).eq.'S') then
!
!              velocity may be needed for contact springs
!
               kode=nelcon(1,imat)
               if(kode.eq.-51) then
                  timeend(1)=time
                  timeend(2)=ttime+dtime
                  call springforc_th(xl,vl,imat,elcon,nelcon,
     &              tnl,ncmat_,ntmat_,nope,kode,elconloc,
     &              plicon,nplicon,npmat_,mi,springarea(1,konl(nope+1)),
     &              timeend,matname,konl(nope),i,istep,iinc)
               endif
!
               do j=1,nope
                     fn(0,konl(j))=fn(0,konl(j))+tnl(j)
               enddo
            endif
         endif
!
         do jj=1,mint3d
            if(lakon(i)(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,jj)
               et=gauss3d2(2,jj)
               ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            elseif(lakon(i)(4:4).eq.'2') then
               xi=gauss3d3(1,jj)
               et=gauss3d3(2,jj)
               ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif(lakon(i)(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif(lakon(i)(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakon(i)(4:5).eq.'15') then
               xi=gauss3d8(1,jj)
               et=gauss3d8(2,jj)
               ze=gauss3d8(3,jj)
               weight=weight3d8(jj)
            elseif(lakon(i)(4:4).eq.'6') then
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif
!
            if(nope.eq.20) then
               if(lakon(i)(7:7).eq.'A') then
                  call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
               elseif((lakon(i)(7:7).eq.'E').or.
     &                (lakon(i)(7:7).eq.'S')) then
                  call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               endif
            elseif(nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
            c1=xsj*weight
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m3=1,3
               vkl(0,m3)=0.d0
            enddo
!
            do m1=1,nope
               do m3=1,3
                  vkl(0,m3)=vkl(0,m3)+shp(m3,m1)*vl(0,m1)
               enddo
            enddo
!
            kode=ncocon(1,imat)
!
!              calculating the temperature difference in
!              the integration point
!
            t1lold=0.d0
            t1l=0.d0
            if(lakon(i)(4:5).eq.'8 ') then
               do i1=1,nope
                  t1lold=t1lold+vold(0,konl(i1))/8.d0
                  t1l=t1l+v(0,konl(i1))/8.d0
               enddo
            elseif(lakon(i)(4:6).eq.'20 ') then
               call lintemp_th(t0,vold,konl,nope,jj,t0l,t1lold,mi)
               call lintemp_th(t0,v,konl,nope,jj,t0l,t1l,mi)
            else
               do i1=1,nope
                  t1lold=t1lold+shp(4,i1)*vold(0,konl(i1))
                  t1l=t1l+shp(4,i1)*v(0,konl(i1))
               enddo
            endif
!
!           calculating the coordinates of the integration point
!           for material orientation purposes (for cylindrical
!           coordinate systems)
!
            if((iorien.gt.0).or.(kode.le.-100)) then
               do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                     pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data; for linear elastic materials
!                 this includes the calculation of the stiffness
!                 matrix
!
            istiff=0
!
            call materialdata_th(cocon,ncocon,imat,iorien,pgauss,orab,
     &           ntmat_,coconloc,mattyp,t1l,rhcon,nrhcon,rho,shcon,
     &           nshcon,sph,xstiff,jj,i,istiff,mi(1))
!
            call thermmodel(amat,i,jj,kode,coconloc,vkl,dtime,
     &           time,ttime,mi(1),nstate_,xstateini,xstate,qflux,xstiff,
     &           iorien,pgauss,orab,t1l,t1lold,vold,co,lakon(i),konl,
     &           ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc)
! 
            qfx(1,jj,i)=qflux(1)
            qfx(2,jj,i)=qflux(2)
            qfx(3,jj,i)=qflux(3)
            if(lakon(i)(6:7).eq.'RA') then
               qfx(1,jj+4,i)=qflux(1)
               qfx(2,jj+4,i)=qflux(2)
               qfx(3,jj+4,i)=qflux(3)
            endif
!
!           calculation of the nodal flux
!
            if(calcul_fn.eq.1)then
!
!                    calculating fn using skl
!
               if(lakon(i)(6:7).eq.'RA') then
                  do m1=1,nope
                     fn(0,konl(m1))=fn(0,konl(m1))
     &                  -c1*(qflux(1)*(shp(1,m1)+shp(1,iperm(m1)))
     &                      +qflux(2)*(shp(2,m1)+shp(2,iperm(m1)))
     &                      +qflux(3)*(shp(3,m1)+shp(3,iperm(m1))))
                  enddo
               else
                  do m1=1,nope
                     do m3=1,3
                        fn(0,konl(m1))=fn(0,konl(m1))-
     &                       c1*qflux(m3)*shp(m3,m1)
                     enddo
                  enddo
               endif
            endif
         enddo
!
!        q contains the contributions to the nodal force in the nodes
!        belonging to the element at stake from other elements (elements
!        already treated). These contributions have to be
!        subtracted to get the contributions attributable to the element
!        at stake only
!
         if(calcul_qa.eq.1) then
            do m1=1,nope
               qa(2)=qa(2)+dabs(fn(0,konl(m1))-q(0,m1))
            enddo
            nal=nal+nope
         endif
      enddo
!
      return
      end
