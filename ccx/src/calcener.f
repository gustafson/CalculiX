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
      subroutine calcener(co,kon,ipkon,lakon,ne,
     &  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,
     &  eme,iperturb,fn,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mint_,ielas,icmd,
     &  ncmat_,nstate_)
!
!     calculates and prints the displacements, temperatures and forces 
!     at the nodes and the stress and  strain  at the reduced integration 
!     points and at the nodes
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(20),nelcon(2,*),nrhcon(*),nalcon(2,*),
     &  ielmat(*),ielorien(*),ntmat_,ipkon(*),iflag
!
      integer ne,mattyp,ithermal,i,j,k,m1,m2,jj,
     &  i1,m3,kk,indexe,nope,norien,iperturb,icmd,ihyper,nmethod,
     &  kode,imat,mint3d,mint_,iorien,ielas,istiff,ncmat_,nstate_
!
      integer nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),shp(4,20),xl(3,20),vl(0:3,20),stre(6),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),elas(21),rho,fn(0:3,*),
     &  fnl(3,9),beta(6),vkl(0:3,3),t0(*),t1(*),
     &  eme(6,mint_,*),ckl(3,3),vold(0:4,*),eloc(9),veold(0:3,*),
     &  elconloc(21),eth(6),xkl(3,3),voldl(0:3,20),xikl(3,3),emec(6),
     &  emec0(6),vel(1:3,20),veoldl(3,9)
!
      real*8 xi,et,ze,tt,xsj,vj,t0l,t1l,dtime,weight,pgauss(3),vij,time,
     &  ttime,potenergy
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mint_,*),xstate(nstate_,mint_,*),plconloc(82),
     &  xstateini(nstate_,mint_,*)
!
      include "gauss.f"
!
      data iflag /3/
!
      potenergy=0.d0
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         imat=ielmat(i)
         amat=matname(imat)
         if(norien.gt.0) then
            iorien=ielorien(i)
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
            read(lakon(i)(8:8),'(i1)') nope
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            mint3d=8
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
               voldl(k,j)=vold(k,konl(j))
               veoldl(k,j)=veold(k,konl(j))
            enddo
         enddo
!
!        check for hyperelastic material
!
         ihyper=0
!
!        calculating the forces for the contact elements
!
         if(mint3d.eq.0) then
            if(lakon(i)(7:7).eq.'A') then
               kode=nelcon(1,imat)
               lakonl=lakon(i)
               t0l=0.d0
               t1l=0.d0
               if(ithermal.eq.1) then
                  t0l=(t0(konl(1))+t0(konl(2)))/2.d0
                  t1l=(t1(konl(1))+t1(konl(2)))/2.d0
               elseif(ithermal.ge.2) then
                  t0l=(t0(konl(1))+t0(konl(2)))/2.d0
                  t1l=(vold(0,konl(1))+vold(0,konl(2)))/2.d0
               endif
            endif
            if(lakon(i)(2:2).eq.'S') then
!
!              velocity may be needed for contact springs
!
               if(lakon(i)(7:7).eq.'C') then
                  do j=1,nope
                     do k=1,3
                        veoldl(k,j)=veold(k,konl(j))
                     enddo
                  enddo
               endif
               call springforc(xl,konl,vl,imat,elcon,nelcon,elas,
     &              fnl,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &              plicon,nplicon,npmat_,veoldl)
               do j=1,nope
                  do k=1,3
                     fn(k,konl(j))=fn(k,konl(j))+fnl(k,j)
                  enddo
               enddo
            elseif((nmethod.eq.4).or.
     &             ((nmethod.eq.1).and.(iperturb.ge.2))) then
               do j=1,nope
                  konl(j)=kon(indexe+j)
                  do k=1,3
                     vel(k,j)=veold(k,konl(j))
                  enddo
               enddo
               call dashforc(xl,konl,vl,imat,elcon,nelcon,
     &              elas,fn,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,
     &              elconloc,plicon,nplicon,npmat_,vel,time,nmethod)
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
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m2=1,3
               do m3=1,3
                  vkl(m2,m3)=0.d0
               enddo
            enddo
!
            do m1=1,nope
               do m2=1,3
                  do m3=1,3
                     vkl(m2,m3)=vkl(m2,m3)+shp(m3,m1)*voldl(m2,m1)
                  enddo
               enddo
            enddo
!
!           calculating the strain
!
            eloc(1)=vkl(1,1)
            eloc(2)=vkl(2,2)
            eloc(3)=vkl(3,3)
            eloc(4)=(vkl(1,2)+vkl(2,1))/2.d0
            eloc(5)=(vkl(1,3)+vkl(3,1))/2.d0
            eloc(6)=(vkl(2,3)+vkl(3,2))/2.d0
!
!                 prestress values
!
            do kk=1,6
               beta(kk)=0.d0
            enddo
!
            if(ithermal.ge.1) then
!
!              calculating the temperature difference in
!              the integration point
!
               t0l=0.d0
               t1l=0.d0
               if(ithermal.eq.1) then
                  if(lakon(i)(4:5).eq.'8 ') then
                     do i1=1,nope
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+t1(konl(i1))/8.d0
                     enddo
                  elseif(lakon(i)(4:6).eq.'20 ') then
                     call lintemp(t0,t1,konl,nope,jj,t0l,t1l)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*t1(konl(i1))
                     enddo
                  endif
               elseif(ithermal.ge.2) then
                  if(lakon(i)(4:5).eq.'8 ') then
                     do i1=1,nope
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+vold(0,konl(i1))/8.d0
                     enddo
                  elseif(lakon(i)(4:6).eq.'20 ') then
                     call lintemp_th(t0,vold,konl,nope,jj,t0l,t1l)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                     enddo
                  endif
               endif
               tt=t1l-t0l
            endif
!
!                 calculating the coordinates of the integration point
!                 for material orientation purposes (for cylindrical
!                 coordinate systems)
!
            if(iorien.gt.0) then
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
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &           nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &           elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &           istiff,elconloc,eth,kode,plicon,nplicon,
     &           plkcon,nplkcon,npmat_,plconloc,mint_,dtime,i,jj,
     &           xstiff,ncmat_)
!
!           determining the mechanical strain
!
            if(ithermal.ne.0) then
               do m1=1,6
                  emec(m1)=eloc(m1)-eth(m1)
                  emec0(m1)=eme(m1,jj,i)
               enddo
            else
               do m1=1,6
                  emec(m1)=eloc(m1)
                  emec0(m1)=eme(m1,jj,i)
               enddo
            endif
!
!           calculating the local stiffness and stress
!
            call mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &           icmd,beta,stre,xkl,ckl,vj,xikl,vij,
     &           plconloc,xstate,xstateini,ielas,
     &           amat,t1l,dtime,time,ttime,i,jj,nstate_,mint_,
     &           iorien,pgauss,orab,eloc,mattyp)
! 
!           updating the energy
!
            potenergy=potenergy+
     &         ((((eloc(1)-eth(1))*stre(1)+
     &            (eloc(2)-eth(2))*stre(2)+
     &            (eloc(3)-eth(3))*stre(3)))/2.d0+
     &            (eloc(4)-eth(4))*stre(4)+
     &            (eloc(5)-eth(5))*stre(5)+
     &            (eloc(6)-eth(6))*stre(6))*weight
!
         enddo
!
      enddo
!
      return
      end
