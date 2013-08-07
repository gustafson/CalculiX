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
      subroutine remeshcontactel(tieset,ntie,set,nset,istartset,
     &  iendset,ialset,ipkon,kon,nkon,lakon,nodface,ipoface,
     &  nk,ipompc,nodempc,ikmpc,ilmpc,nmpc,nmpc_,labmpc,coefmpc,
     &  mpcfree,nalset,co,ithermal,nk0,ne,ielmat,ielorien,mi,t0,
     &  vold,veold,iponoel,inoel,xstate,nstate_,prestr,iprestr)
!
!     remeshing quadratic elements adjacent to contact surfaces
!     using appropriate linear elements (C3D8I,C3D4 and C3D6)
!
      implicit none
!
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 tieset(3,*),surfset,set(*)
!
      integer ipoface(*),nodface(9,*),nodes(8),nk,jface,iaux,
     &  ne,ipkon(*),kon(*),indexe,ifaceq(9,6),index1,ialset(*),
     &  ifacew(8,5),kflag,i,j,k,l,m,n,istartset(*),iendset(*),
     &  ntie,ipos,ij,ifour,nset,konl(27),ipompc(*),nodempc(3,*),
     &  ikmpc(*),ilmpc(*),nmpc,nmpc_,mpcfree,nalset,ielface10(4,4),
     &  ielface15(4,5),ielface20(4,6),idof,ithermal(2),jmin,jmax,
     &  kon10(4,8),kon15(6,8),kon20(8,8),nkon,is,ie,nk0,mpcfreeold,
     &  mi(*),ielmat(mi(3),*),ielorien(mi(3),*),iponoel(*),inoel(2,*),
     &  ifacet(7,4),ik,node,indexe1,konl1(10),nstate_,iprestr,ll,nope,
     &  iflag,mm,nf
!
      real*8 coefmpc(*),co(3,*),xstate(nstate_,mi(1),*),shp(4,20),
     &  prestr(6,mi(1),*),field(6,20),fieldst(nstate_,20),a8(8,8),
     &  a4(4,4),a27(20,27),a9(6,9),xi,et,ze,g10(3,8),g15(3,8),g20(3,8),
     &  xl(3,20),xsj,vold(0:mi(2),*),veold(0:mi(2),*),t0(*)
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data ielface10 /1,2,3,4,
     &                1,2,7,8,
     &                2,3,5,8,
     &                1,3,6,8/
      data ielface15 /1,2,3,4,
     &                5,6,7,8,
     &                1,2,5,6,
     &                2,3,6,7,
     &                1,3,5,7/
      data ielface20 /1,2,3,4,
     &                5,6,7,8,
     &                1,2,5,6,
     &                2,4,6,8,
     &                3,4,7,8,
     &                1,3,5,7/
      data kon10 /1,5,7,8,
     &            5,2,6,9,
     &            7,6,3,10,
     &            5,6,7,10,
     &            5,6,10,9,
     &            7,5,10,8,
     &            5,9,10,8,
     &            8,9,10,4/
      data kon15 /1,7,9,13,16,18,
     &            7,2,8,16,14,17,
     &            9,8,3,18,17,15,
     &            7,8,9,16,17,18,
     &            13,16,18,4,10,12,
     &            16,14,17,10,5,11,
     &            18,17,15,12,11,6,
     &            16,17,18,10,11,12/
      data kon20 /1,9,21,12,17,23,27,26,
     &            9,2,10,21,23,18,24,27,
     &            12,21,11,4,26,27,25,20,
     &            21,10,3,11,27,24,19,25,
     &            17,23,27,26,5,13,22,16,
     &            23,18,24,27,13,6,14,22,
     &            26,27,25,20,16,22,15,8,
     &            27,24,19,25,22,14,7,15/
!
      data a4 /  1.92705, -0.30902, -0.30902, -0.30902,
     &          -0.30902,  1.92705, -0.30902, -0.30902,
     &          -0.30902, -0.30902,  1.92705, -0.30902,
     &          -0.30902, -0.30902, -0.30902,  1.92705/
      data a9 / 1.63138,-0.32628,-0.32628,-0.52027, 0.10405, 0.10405,
     &         -0.32628, 1.63138,-0.32628, 0.10405,-0.52027, 0.10405,
     &         -0.32628,-0.32628, 1.63138, 0.10405, 0.10405,-0.52027,
     &          0.55556,-0.11111,-0.11111, 0.55556,-0.11111,-0.11111,
     &         -0.11111, 0.55556,-0.11111,-0.11111, 0.55556,-0.11111,
     &         -0.11111,-0.11111, 0.55556,-0.11111,-0.11111, 0.55556,
     &         -0.52027, 0.10405, 0.10405, 1.63138,-0.32628,-0.32628,
     &          0.10405,-0.52027, 0.10405,-0.32628, 1.63138,-0.32628,
     &          0.10405, 0.10405,-0.52027,-0.32628,-0.32628, 1.63138/
      data a8 /2.549,-.683,.183,-.683,-.683,.183,
     &        -.04904,.183,-.683,2.549,-.683,.183,
     &        .183,-.683,.183,-.04904,-.683,.183,
     &        -.683,2.549,.183,-.04904,.183,-.683,
     &        .183,-.683,2.549,-.683,-.04904,.183,
     &        -.683,.183,-.683,.183,-.04904,.183,
     &        2.549,-.683,.183,-.683,.183,-.683,
     &        .183,-.04904,-.683,2.549,-.683,.183,
     &        .183,-.04904,.183,-.683,-.683,.183,
     &        -.683,2.549,-.04904,.183,-.683,.183,
     &        .183,-.683,2.549,-.683/      
      data a27 /
     &  2.37499,-0.12559,-0.16145,-0.12559,-0.12559,-0.16145, 0.11575,
     & -0.16145, 0.32628, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     & -0.10405, 0.11111, 0.32628, 0.11111,-0.10405, 0.11111,-0.31246,
     & -0.31246, 0.31481, 0.31481, 0.31481, 0.31481,-0.16902,-0.16902,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.19444, 0.15961,-0.00661,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.12559, 2.37499,
     & -0.12559,-0.16145,-0.16145,-0.12559,-0.16145, 0.11575, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111,-0.10405,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.31481,-0.16902,-0.16902, 0.31481,-0.27072,-0.19444,
     & -0.27072, 1.28439, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.22898, 0.22898, 0.22898, 0.22898, 0.05556, 0.05556, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556,-0.22222,-0.22222,
     & -0.22222,-0.22222, 0.31481,-0.31246,-0.31246, 0.31481,-0.16902,
     &  0.31481, 0.31481,-0.16902,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.19444, 0.15961,-0.00661, 0.15961,-0.27072,-0.27072,
     &  0.15961,-0.12559,-0.16145,-0.12559, 2.37499,-0.16145, 0.11575,
     & -0.16145,-0.12559, 0.11111, 0.11111, 0.32628, 0.32628,-0.10405,
     & -0.10405, 0.11111, 0.11111, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.16902,-0.16902, 0.31481,
     &  0.31481,-0.19444,-0.27072, 1.28439,-0.27072,-0.00661, 0.15961,
     & -0.19444, 0.15961, 0.15961, 0.15961,-0.27072,-0.27072,-0.16145,
     & -0.12559, 2.37499,-0.12559, 0.11575,-0.16145,-0.12559,-0.16145,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.11111,
     & -0.10405,-0.10405, 0.11111, 0.32628, 0.11111,-0.31246, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     &  1.28439,-0.19444,-0.00661,-0.19444,-0.48824,-0.48824, 0.22898,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898, 0.05556,-0.22222,
     &  0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.31481,-0.31246, 0.31481,-0.16902,
     &  0.31481,-0.31246, 0.31481,-0.16902,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.19444, 1.28439,
     & -0.19444,-0.00661,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898, 0.22898,-0.48824,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,
     & -0.29630,-0.29630,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     & -0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.31246,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     & -0.27072,-0.19444,-0.00661,-0.19444, 1.28439, 0.22898, 0.22898,
     & -0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222,
     &  0.05556, 0.05556, 0.05556, 0.05556,-0.16902, 0.31481,-0.31246,
     &  0.31481,-0.16902, 0.31481,-0.31246, 0.31481, 0.15961,-0.27072,
     & -0.27072, 0.15961, 0.15961,-0.27072,-0.27072, 0.15961,-0.00661,
     & -0.19444, 1.28439,-0.19444,-0.12559,-0.16145, 0.11575,-0.16145,
     &  2.37499,-0.12559,-0.16145,-0.12559, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.31481, 0.31481,-0.16902,-0.16902,-0.31246,
     & -0.31246, 0.31481, 0.31481,-0.19444, 0.15961,-0.00661, 0.15961,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.16145,-0.12559,-0.16145, 0.11575,-0.12559, 2.37499,
     & -0.12559,-0.16145, 0.11111, 0.11111,-0.10405,-0.10405, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     &  0.31481,-0.16902,-0.16902, 0.31481,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,-0.19444,
     & -0.27072, 1.28439,-0.27072, 0.15961, 0.15961,-0.27072, 0.22898,
     &  0.22898, 0.22898, 0.22898,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.22222,-0.22222,-0.22222,-0.22222,-0.16902, 0.31481,
     &  0.31481,-0.16902, 0.31481,-0.31246,-0.31246, 0.31481, 0.15961,
     & -0.19444, 0.15961,-0.00661,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.27072,-0.27072, 0.15961,-0.16145, 0.11575,-0.16145,
     & -0.12559,-0.12559,-0.16145,-0.12559, 2.37499,-0.10405,-0.10405,
     &  0.11111, 0.11111, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.32628,-0.16902,-0.16902, 0.31481, 0.31481,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.00661, 0.15961,-0.19444,
     &  0.15961,-0.19444,-0.27072, 1.28439,-0.27072, 0.15961, 0.15961,
     & -0.27072,-0.27072, 0.11575,-0.16145,-0.12559,-0.16145,-0.16145,
     & -0.12559, 2.37499,-0.12559,-0.10405, 0.11111, 0.11111,-0.10405,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.11111/
!
      data g10 /.125d0,.125d0,.125d0,
     &          .625d0,.125d0,.125d0,
     &          .125d0,.625d0,.125d0,
     &          .250d0,.375d0,.125d0,
     &          .375d0,.250d0,.250d0,
     &          .125d0,.250d0,.250d0,
     &          .250d0,.125d0,.375d0,
     &          .125d0,.125d0,.625d0/
!
      data g15 /.1666666666d0,.1666666666d0,-.5d0,
     &          .6666666666d0,.1666666666d0,-.5d0,
     &          .1666666666d0,.6666666666d0,-.5d0,
     &          .3333333333d0,.3333333333d0,-.5d0,
     &          .1666666666d0,.1666666666d0,0.5d0,
     &          .6666666666d0,.1666666666d0,0.5d0,
     &          .1666666666d0,.6666666666d0,0.5d0,
     &          .3333333333d0,.3333333333d0,0.5d0/
!
      data g20 /-.5d0,-.5d0,-.5d0,
     &          0.5d0,-.5d0,-.5d0,
     &          -.5d0,0.5d0,-.5d0,
     &          0.5d0,0.5d0,-.5d0,
     &          -.5d0,-.5d0,0.5d0,
     &          0.5d0,-.5d0,0.5d0,
     &          -.5d0,0.5d0,0.5d0,
     &          0.5d0,0.5d0,0.5d0/
!
      include "gauss.f"
!
      ifour=4
      kflag=1
!
!     degrees of freedom for the MPC's
!     the thermal equations generated for ithermal(2)=1 are
!     used in tempload.f
!
c      if(ithermal(2).le.1) then
      if(ithermal(2).eq.0) then
         jmin=1
         jmax=3
      elseif(ithermal(2).eq.2) then
         jmin=0
         jmax=0
      else
         jmin=0
         jmax=3
      endif
!
!     generating additional nodes in the middle of the 
!     8-node faces
!
      do i=1,nk0
         index1=ipoface(i)
         do
            if(index1.eq.0) exit
            nk=nk+1
            nodes(1)=i
            do j=2,8
               nodes(j)=nodface(j-1,index1)
            enddo
!
!           coordinates of the new node
!
            do j=1,3
               co(j,nk)=0.d0
               do k=1,4
                  co(j,nk)=co(j,nk)-co(j,nodes(k))
               enddo
               do k=5,8
                  co(j,nk)=co(j,nk)+2.d0*co(j,nodes(k))
               enddo
               co(j,nk)=co(j,nk)/4.d0
            enddo
!
!           initial conditions
!
            if(ithermal(1).gt.0) then
               t0(nk)=0.d0
               do k=1,4
                  t0(nk)=t0(nk)-t0(nodes(k))
               enddo
               do k=5,8
                  t0(nk)=t0(nk)+2.d0*t0(nodes(k))
               enddo
               t0(nk)=t0(nk)/4.d0
            endif
!
            do j=0,mi(2)
               vold(j,nk)=0.d0
               veold(j,nk)=0.d0
               do k=1,4
                  vold(j,nk)=vold(j,nk)-vold(j,nodes(k))
                  veold(j,nk)=veold(j,nk)-veold(j,nodes(k))
               enddo
               do k=5,8
                  vold(j,nk)=vold(j,nk)+2.d0*vold(j,nodes(k))
                  veold(j,nk)=veold(j,nk)+2.d0*veold(j,nodes(k))
               enddo
               vold(j,nk)=vold(j,nk)/4.d0
               veold(j,nk)=veold(j,nk)/4.d0
            enddo
!
!           generating MPC's if the node belongs to an interface,
!           i.e. is external to the volume of elements adjacent to 
!           the contact surfaces and does not belong to the contact
!           surfaces itself
!
            if(nodface(8,index1).gt.0) then
               do j=jmin,jmax
                  idof=8*(nk-1)+j
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) '*ERROR in remeshcontactel:'
                     write(*,*) '       increase nmpc_'
                     stop
                  endif
                  ipompc(nmpc)=mpcfree
                  labmpc(nmpc)='                    '
                  ikmpc(nmpc)=idof
                  ilmpc(nmpc)=nmpc
                  nodempc(1,mpcfree)=nk
                  nodempc(2,mpcfree)=j
                  coefmpc(mpcfree)=4.d0
                  mpcfree=nodempc(3,mpcfree)
                  do k=1,4
                     nodempc(1,mpcfree)=nodes(k)
                     nodempc(2,mpcfree)=j
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                  enddo
                  do k=5,8
                     nodempc(1,mpcfree)=nodes(k)
                     nodempc(2,mpcfree)=j
                     coefmpc(mpcfree)=-2.d0
                     if(k.eq.8) mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                  enddo
                  nodempc(3,mpcfreeold)=0
               enddo
            endif
!
!           storing the additional node in field nodface
!
            nodface(8,index1)=nk
            index1=nodface(9,index1)
         enddo
      enddo
!
!     remeshing the elements
!
      do ll=1,ntie
!     
!     check for contact conditions
!     
         if((tieset(1,ll)(81:81).eq.'C').or.
     &        (tieset(1,ll)(81:81).eq.'-')) then
!     
!     contact constraint
!     
            do m=2,3
               surfset=tieset(m,ll)
!     
!     check whether facial surface
!     
               ipos=index(surfset,' ')-1
c               if(surfset(ipos:ipos).eq.'S') cycle
!     
               do n=1,nset
                  if(set(n).eq.surfset) exit
               enddo
!
               if(n.le.nset) then
                  if(surfset(ipos:ipos).eq.'S') cycle
               else
                  do n=1,nset
                     if((set(n)(1:ipos-1).eq.surfset(1:ipos-1)).and.
     &                  (set(n)(ipos:ipos).eq.'T')) exit
                  enddo
               endif
!
!              storing the actual starting and ending values
!
               is=istartset(n)
               ie=iendset(n)
!     
               do ij=is,ie
!     
                  i=int(ialset(ij)/10.d0)
                  jface=ialset(ij)-10*i
                  indexe=ipkon(i)
                  if(indexe.lt.0) cycle
!     
!     quadratic hexahedral element (beam and shell inclusive)
!     
                  if((lakon(i)(4:4).eq.'2').and.
     &               ((lakon(i)(7:7).eq.' ').or.
     &                (lakon(i)(7:7).eq.'L').or.
     &                (lakon(i)(7:7).eq.'B'))) then
!
!                    storing the nodes of the C3D20 element
!
                     do j=1,20
                        konl(j)=kon(indexe+j)
                     enddo
!
!                    storing the nodes in the middle of the faces
!
                     do j=1,6
                        do k=1,4
                           nodes(k)=kon(indexe+ifaceq(k,j))
                        enddo
                        call isortii(nodes,iaux,ifour,kflag)
                        index1=ipoface(nodes(1))
                        do
                           if(index1.eq.0) then
                              write(*,*) '*ERROR in remeshcontactel'
                              write(*,*) '       face',ialset(ij),
     &                          '       is not catalogued'
                              stop
                           endif
                           if((nodface(1,index1).eq.nodes(2)).and.
     &                          (nodface(2,index1).eq.nodes(3)).and.
     &                          (nodface(3,index1).eq.nodes(4))) then
                              konl(20+j)=nodface(8,index1)
                              exit
                           endif
                           index1=nodface(9,index1)
                        enddo
                     enddo
!
!                    generating a new node in the middle of the
!                    element
!
                     nk=nk+1
                     konl(27)=nk
!
                     do j=1,3
                        co(j,nk)=0.d0
                        do k=1,8
                           co(j,nk)=co(j,nk)-co(j,konl(k))
                        enddo
                        do k=9,20
                           co(j,nk)=co(j,nk)+co(j,konl(k))
                        enddo
                        co(j,nk)=co(j,nk)/4.d0
                     enddo
!
!     initial conditions
!     
                     if(ithermal(1).gt.0) then
                        t0(nk)=0.d0
                        do k=1,4
                           t0(nk)=t0(nk)-t0(nodes(k))
                        enddo
                        do k=5,8
                           t0(nk)=t0(nk)+2.d0*t0(nodes(k))
                        enddo
                        t0(nk)=t0(nk)/4.d0
                     endif
!     
                     do j=0,mi(2)
                        vold(j,nk)=0.d0
                        veold(j,nk)=0.d0
                        do k=1,4
                           vold(j,nk)=vold(j,nk)-vold(j,nodes(k))
                           veold(j,nk)=veold(j,nk)-veold(j,nodes(k))
                        enddo
                        do k=5,8
                           vold(j,nk)=vold(j,nk)+2.d0*vold(j,nodes(k))
                           veold(j,nk)=
     &                         veold(j,nk)+2.d0*veold(j,nodes(k))
                        enddo
                        vold(j,nk)=vold(j,nk)/4.d0
                        veold(j,nk)=veold(j,nk)/4.d0
                     enddo
!
!                    deactivating the remeshed element
!
                     ipkon(i)=-ipkon(i)-2
!
!                    first node of the original element is 
!                    replaced by the first subelement
!
                     kon(indexe+1)=ne+1
!
!                    extrapolating the initial integration point
!                    values to the nodes
!
                     if((iprestr.gt.0).or.(nstate_.gt.0)) then
                        iflag=1
                        if(lakon(i)(6:6).eq.'R') then
                           if(iprestr.gt.0) then
                              do j=1,8
                                 do k=1,6
                                    field(k,j)=0.d0
                                    do l=1,8
                                       field(k,j)=field(k,j)+
     &                                            a8(j,l)*prestr(k,l,i)
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(nstate_.gt.0) then
                              do j=1,8
                                 do k=1,nstate_
                                    fieldst(k,j)=0.d0
                                    do l=1,8
                                       fieldst(k,j)=fieldst(k,j)+
     &                                            a8(j,l)*xstate(k,l,i)
                                    enddo
                                 enddo
                              enddo
                           endif
                        else
                           if(iprestr.gt.0) then
                              do j=1,20
                                 do k=1,6
                                    field(k,j)=0.d0
                                    do l=1,27
                                       field(k,j)=field(k,j)+
     &                                            a27(j,l)*prestr(k,l,i)
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(nstate_.gt.0) then
                              do j=1,20
                                 do k=1,nstate_
                                    fieldst(k,j)=0.d0
                                    do l=1,27
                                       fieldst(k,j)=fieldst(k,j)+
     &                                            a27(j,l)*xstate(k,l,i)
                                    enddo
                                 enddo
                              enddo
                           endif
                        endif
                     endif
!
!                    generating the new elements
!
                     do j=1,8
                        ne=ne+1
                        ipkon(ne)=nkon
                        lakon(ne)='C3D8IC  '
                        do k=1,mi(3)
                           ielmat(k,ne)=ielmat(k,i)
                           ielorien(k,ne)=ielorien(k,i)
                        enddo
                        do k=1,8
                           kon(nkon+k)=konl(kon20(k,j))
                        enddo
!
!                       internal nodes for the C3D8I element
!
                        do k=9,11
                           nk=nk+1
                           kon(nkon+k)=nk
                           do l=1,3
                              co(l,nk)=0.d0
                           enddo
                           do l=0,mi(2)
                              vold(l,nk)=0.d0
                              veold(l,nk)=0.d0
                              if(ithermal(1).gt.0) then
                                 t0(nk)=0.d0
                              endif
                           enddo
                        enddo
                        nkon=nkon+11
!
!                       interpolation of initial nodal values to the
!                       new integration points
!
                        if((iprestr.gt.0).or.(nstate_.gt.0)) then
                           do k=1,8
                              xi=g20(1,k)+gauss3d2(1,k)/2.d0
                              et=g20(1,k)+gauss3d2(2,k)/2.d0
                              ze=g20(1,k)+gauss3d2(3,k)/2.d0
                              if(lakon(i)(6:6).eq.'R') then
                                 call shape8h(xi,et,ze,xl,xsj,shp,iflag)
                                 nope=8
                              else
                                 call shape20h(xi,et,ze,xl,xsj,shp,
     &                                         iflag)
                                 nope=20
                              endif
                              if(iprestr.gt.0) then
                                 do l=1,6
                                    do mm=1,nope
                                       prestr(l,k,ne)=prestr(l,k,ne)+
     &                                    shp(4,mm)*field(l,mm)
                                    enddo
                                 enddo
                              endif
                              if(nstate_.gt.0) then
                                 do l=1,nstate_
                                    do mm=1,nope
                                       xstate(l,k,ne)=xstate(l,k,ne)+
     &                                    shp(4,mm)*fieldst(l,mm)
                                    enddo
                                 enddo
                              endif
                           enddo
                        endif
!     
                     enddo
!     
!     quadratic tetrahedral element
!     
                  elseif((lakon(i)(4:5).eq.'10')) then
!     
!     storing the nodes of the C3D10 element
!     
                     do j=1,10
                        konl(j)=kon(indexe+j)
                     enddo
!     
!                    deactivating the remeshed element
!     
                     ipkon(i)=-ipkon(i)-2
!     
!                    first node of the original element is 
!                    replaced by the first subelement
!     
                     kon(indexe+1)=ne+1
!
!                    extrapolating the initial integration point
!                    values to the nodes
!
                     if((iprestr.gt.0).or.(nstate_.gt.0)) then
                        iflag=1
                        if(iprestr.gt.0) then
                           do j=1,4
                              do k=1,6
                                 field(k,j)=0.d0
                                 do l=1,4
                                    field(k,j)=field(k,j)+
     &                                   a4(j,l)*prestr(k,l,i)
                                 enddo
                              enddo
                           enddo
                        endif
                        if(nstate_.gt.0) then
                           do j=1,4
                              do k=1,nstate_
                                 fieldst(k,j)=0.d0
                                 do l=1,4
                                    fieldst(k,j)=fieldst(k,j)+
     &                                   a4(j,l)*xstate(k,l,i)
                                 enddo
                              enddo
                           enddo
                        endif
                     endif
!     
!                    generating the new elements
!     
                     do j=1,8
                        ne=ne+1
                        ipkon(ne)=nkon
                        lakon(ne)='C3D4 C  '
                        do k=1,mi(3)
                           ielmat(k,ne)=ielmat(k,i)
                           ielorien(k,ne)=ielorien(k,i)
                        enddo
                        do k=1,4
                           kon(nkon+k)=konl(kon10(k,j))
                        enddo
                        nkon=nkon+4
!
!                       interpolation of initial nodal values to the
!                       new integration points
!
                        if((iprestr.gt.0).or.(nstate_.gt.0)) then
                           do k=1,1
                              xi=g10(1,k)+gauss3d4(1,k)/2.d0
                              et=g10(1,k)+gauss3d4(2,k)/2.d0
                              ze=g10(1,k)+gauss3d4(3,k)/2.d0
                              call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
                              nope=4
                              if(iprestr.gt.0) then
                                 do l=1,6
                                    do mm=1,nope
                                       prestr(l,k,ne)=prestr(l,k,ne)+
     &                                    shp(4,mm)*field(l,mm)
                                    enddo
                                 enddo
                              endif
                              if(nstate_.gt.0) then
                                 do l=1,nstate_
                                    do mm=1,nope
                                       xstate(l,k,ne)=xstate(l,k,ne)+
     &                                    shp(4,mm)*fieldst(l,mm)
                                    enddo
                                 enddo
                              endif
                           enddo
                        endif
                     enddo
!     
!                    remeshing all C3D10 elements containing an end node
!                    of the contact surface but which do not contain a
!                    face of the contact surface
!
                     do ik=1,3
                        node=konl(ifacet(ik,jface))
                        index1=iponoel(node)
                        do
                           if(index1.eq.0) exit
                           i=inoel(1,index1)
!
!                          next lines are for C3D10 elements which
!                          have a face in common with the contact
!                          surface
!
                           if(i.eq.0) then
                              index1=inoel(2,index1)
                              cycle
                           endif
                           indexe1=ipkon(i)
                           if(indexe1.ge.0) then
!
!                              storing the nodes of the C3D10 element
!     
                              do j=1,10
                                 konl1(j)=kon(indexe1+j)
                              enddo
!
!                             deactivating the remeshed element
!     
                              ipkon(i)=-ipkon(i)-2
!     
!                             first node of the original element is 
!                             replaced by the first subelement
!     
                              kon(indexe1+1)=ne+1
!
!                    extrapolating the initial integration point
!                    values to the nodes
!
                              if((iprestr.gt.0).or.(nstate_.gt.0)) then
                                 iflag=1
                                 if(iprestr.gt.0) then
                                    do j=1,4
                                       do k=1,6
                                          field(k,j)=0.d0
                                          do l=1,4
                                             field(k,j)=field(k,j)+
     &                                            a4(j,l)*prestr(k,l,i)
                                          enddo
                                       enddo
                                    enddo
                                 endif
                                 if(nstate_.gt.0) then
                                    do j=1,4
                                       do k=1,nstate_
                                          fieldst(k,j)=0.d0
                                          do l=1,4
                                             fieldst(k,j)=fieldst(k,j)+
     &                                            a4(j,l)*xstate(k,l,i)
                                          enddo
                                       enddo
                                    enddo
                                 endif
                              endif
!     
!                             generating the new elements
!     
                              do j=1,8
                                 ne=ne+1
                                 ipkon(ne)=nkon
                                 lakon(ne)='C3D4 C  '
                                 do k=1,mi(3)
                                    ielmat(k,ne)=ielmat(k,i)
                                    ielorien(k,ne)=ielorien(k,i)
                                 enddo
                                 do k=1,4
                                    kon(nkon+k)=konl1(kon10(k,j))
                                 enddo
                                 nkon=nkon+4
!
!                                interpolation of initial nodal values to the
!                                new integration points
!
                                 if((iprestr.gt.0).or.(nstate_.gt.0)) 
     &                                  then
                                    do k=1,1
                                       xi=g10(1,k)+gauss3d4(1,k)/2.d0
                                       et=g10(1,k)+gauss3d4(2,k)/2.d0
                                       ze=g10(1,k)+gauss3d4(3,k)/2.d0
                                       call shape4tet(xi,et,ze,xl,xsj,
     &                                      shp,iflag)
                                       nope=4
                                       if(iprestr.gt.0) then
                                          do l=1,6
                                             do mm=1,nope
                                                prestr(l,k,ne)=
     &                                            prestr(l,k,ne)+
     &                                            shp(4,mm)*field(l,mm)
                                             enddo
                                          enddo
                                       endif
                                       if(nstate_.gt.0) then
                                          do l=1,nstate_
                                             do mm=1,nope
                                                xstate(l,k,ne)=
     &                                           xstate(l,k,ne)+
     &                                           shp(4,mm)*fieldst(l,mm)
                                             enddo
                                          enddo
                                       endif
                                    enddo
                                 endif
                              enddo
                           endif
                           index1=inoel(2,index1)
                        enddo
                     enddo
!     
!     quadratic wedge element (beam and shell inclusive)
!     
                  elseif((lakon(i)(4:5).eq.'15').and.
     &               ((lakon(i)(7:7).eq.' ').or.
     &                (lakon(i)(7:7).eq.'L').or.
     &                (lakon(i)(7:7).eq.'B'))) then
!
!                    storing the nodes of the C3D15 element
!
                     do j=1,15
                        konl(j)=kon(indexe+j)
                     enddo
!
!                    storing the nodes in the middle of the faces
!
                     do j=3,5
                        do k=1,4
                           nodes(k)=kon(indexe+ifacew(k,j))
                        enddo
                        call isortii(nodes,iaux,ifour,kflag)
                        index1=ipoface(nodes(1))
                        do
                           if(index1.eq.0) then
                              write(*,*) '*ERROR in remeshcontactel'
                              write(*,*) '       face',ialset(ij),
     &                           '       is not catalogued'
                              stop
                           endif
                           if((nodface(1,index1).eq.nodes(2)).and.
     &                          (nodface(2,index1).eq.nodes(3)).and.
     &                          (nodface(3,index1).eq.nodes(4))) then
                              konl(13+j)=nodface(8,index1)
                              exit
                           endif
                           index1=nodface(9,index1)
                        enddo
                     enddo
!
!                    deactivating the remeshed element
!
                     ipkon(i)=-ipkon(i)-2
!
!                    first node of the original element is 
!                    replaced by the first subelement
!
                     kon(indexe+1)=ne+1
!
!                    extrapolating the initial integration point
!                    values to the nodes
!
                     if((iprestr.gt.0).or.(nstate_.gt.0)) then
                        iflag=1
                        if(iprestr.gt.0) then
                           do j=1,6
                              do k=1,6
                                 field(k,j)=0.d0
                                 do l=1,9
                                    field(k,j)=field(k,j)+
     &                                   a9(j,l)*prestr(k,l,i)
                                 enddo
                              enddo
                           enddo
                        endif
                        if(nstate_.gt.0) then
                           do j=1,6
                              do k=1,nstate_
                                 fieldst(k,j)=0.d0
                                 do l=1,9
                                    fieldst(k,j)=fieldst(k,j)+
     &                                   a9(j,l)*xstate(k,l,i)
                                 enddo
                              enddo
                           enddo
                        endif
                     endif
!
!                    generating the new elements
!
                     do j=1,8
                        ne=ne+1
                        ipkon(ne)=nkon
                        lakon(ne)='C3D6 C  '
                        do k=1,mi(3)
                           ielmat(k,ne)=ielmat(k,i)
                           ielorien(k,ne)=ielorien(k,i)
                        enddo
                        do k=1,6
                           kon(nkon+k)=konl(kon15(k,j))
                        enddo
                        nkon=nkon+6
!
!                       interpolation of initial nodal values to the
!                       new integration points
!
                        if((iprestr.gt.0).or.(nstate_.gt.0)) then
                           do k=1,2
                              xi=g15(1,k)+gauss3d7(1,k)/2.d0
                              et=g15(1,k)+gauss3d7(2,k)/2.d0
                              ze=g15(1,k)+gauss3d7(3,k)/2.d0
                              call shape6w(xi,et,ze,xl,xsj,shp,iflag)
                              nope=6
                              if(iprestr.gt.0) then
                                 do l=1,6
                                    do mm=1,nope
                                       prestr(l,k,ne)=prestr(l,k,ne)+
     &                                    shp(4,mm)*field(l,mm)
                                    enddo
                                 enddo
                              endif
                              if(nstate_.gt.0) then
                                 do l=1,nstate_
                                    do mm=1,nope
                                       xstate(l,k,ne)=xstate(l,k,ne)+
     &                                    shp(4,mm)*fieldst(l,mm)
                                    enddo
                                 enddo
                              endif
                           enddo
                        endif
                     enddo
                  endif
               enddo
!
            enddo
         endif
      enddo
!
!     updating the contact faces
!
      do ll=1,ntie
!     
!     check for contact conditions
!     
         if((tieset(1,ll)(81:81).eq.'C').or.
     &        (tieset(1,ll)(81:81).eq.'-')) then
!     
!     contact constraint
!     
            do m=2,3
               surfset=tieset(m,ll)
!     
!     check whether facial surface
!     
               ipos=index(surfset,' ')-1
!     
               do n=1,nset
                  if(set(n).eq.surfset) exit
               enddo
!
               if(n.le.nset) then
                  if(surfset(ipos:ipos).eq.'S') cycle
               else
                  do n=1,nset
                     if((set(n)(1:ipos-1).eq.surfset(1:ipos-1)).and.
     &                  (set(n)(ipos:ipos).eq.'T')) exit
                  enddo
               endif
!
!              storing the actual starting and ending values
!
               is=istartset(n)
               ie=iendset(n)
!
               istartset(n)=nalset+1
!     
               do ij=is,ie
!     
                  i=int(ialset(ij)/10.d0)
                  jface=ialset(ij)-10*i
                  indexe=ipkon(i)
!
!                 nf is the number of the element preceding the
!                 elements used for the remeshing of element i
!
                  if(indexe.lt.-1) then
                     nf=kon(-indexe-1)-1
                  elseif(indexe.eq.-1) then
                     cycle
                  endif
!     
!     quadratic hexahedral element (beam and shell inclusive)
!     
                  if((lakon(i)(4:4).eq.'2').and.
     &               ((lakon(i)(7:7).eq.' ').or.
     &                (lakon(i)(7:7).eq.'L').or.
     &                (lakon(i)(7:7).eq.'B'))) then
!
!                    catalogueing the new contact faces
!
                     do j=1,4
                        nalset=nalset+1
                        ialset(nalset)=10*(nf+ielface20(j,jface))+jface
                     enddo
!     
!     quadratic tetrahedral element
!     
                  elseif((lakon(i)(4:5).eq.'10')) then
!     
!                    catalogueing the new contact faces
!     
                     do j=1,4
                        nalset=nalset+1
                        ialset(nalset)=10*(nf+ielface10(j,jface))+jface
                     enddo
!     
!     quadratic wedge element (beam and shell inclusive)
!     
                  elseif((lakon(i)(4:5).eq.'15').and.
     &               ((lakon(i)(7:7).eq.' ').or.
     &                (lakon(i)(7:7).eq.'L').or.
     &                (lakon(i)(7:7).eq.'B'))) then
!
!                    catalogueing the new contact faces
!
                     do j=1,4
                        nalset=nalset+1
                        ialset(nalset)=10*(nf+ielface15(j,jface))+jface
                     enddo
!     
                  else
!
!                    linear elements: storing the contact face
!
                     nalset=nalset+1
                     ialset(nalset)=ialset(ij)
                  endif
               enddo
!
               iendset(n)=nalset
            enddo
         endif
      enddo
!     
      return
      end
