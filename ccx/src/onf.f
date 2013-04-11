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
      subroutine onf(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
     &  kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,xstaten,
     &  nstate_,istep,iinc,iperturb,ener,mi)
!
!     stores the results in frd format
!
      implicit none
!
      character*5 m1
      character*8 lakon(*)
      character*80 matname(*)
      character*87 filab(*)
!
      integer kon(*),inum(*),nk,ne,nmethod,kode,i,j,ipkon(*),indexe,
     &  one,ielmat(*),nstate_,istep,iinc,ianatyp,iperturb,
     &  konl(20),jj,mint3d,k,nope,mi(2),n,kflag,iy(3),iflag,
     &  nktrue,netrue
!
      real*8 co(3,*),v(3,*),stn(6,*),een(6,*),t1(*),fn(3,*),time,
     &  epn(*),enern(*),xstaten(nstate_,*),zero,stnprin(3),str(6),
     &  str2(6),v1,v2,v3,ener(mi(1),*),xi,et,ze,weight,volume,
     &  energy,totenergy,xsj,xl(3,20),shp(4,20),tt,cm,cn,
     &  bb,cc,pi
!
      include "gauss.f"
!
      data iflag /2/
!
      one=1
      zero=0.d0
      m1='   -1'
      pi=4.d0*datan(1.d0)
      n=3
      kflag=-1
!
c         open(12,file='beam_520.onf',status='unknown')
c         open(13,file='beam.onf',status='unknown')
      if(nmethod.eq.1) then
         if(iperturb.gt.1) then
            ianatyp=5
         else
            ianatyp=1
         endif
      elseif(nmethod.eq.2) then
         ianatyp=2
      elseif(nmethod.eq.3) then
         ianatyp=6
      elseif(nmethod.eq.4) then
         if(iperturb.gt.1) then
            ianatyp=2
         else
            ianatyp=4
         endif
      endif
!
!     calculating the true number of nodes and elements
!
      nktrue=0
      do i=1,nk
         if(inum(i).ne.0) nktrue=nktrue+1
      enddo
      netrue=0
      do i=1,ne
         if(ipkon(i).ge.0) netrue=netrue+1
      enddo
!
!     storing the frequency and/or the buckling eigenvalue
!
      if((nmethod.eq.2).or.(nmethod.eq.3)) then
         write(11,'(a5)') m1
         write(11,'(a3)') '500'
         write(11,'(i1)') one
         write(11,'(i5,",",i5,",",i5)') ianatyp,istep,iinc
         write(11,'(i10)') one
         write(11,99) time,zero,zero,zero
      endif
!
!     storing the displacements of the nodes
!
      if(filab(1)(1:4).eq.'U   ') then
!
         write(11,'(a5)') m1
         write(11,'(a3)') '510'
         write(11,'(i1)') one
         write(11,'(i5,",",i5,",",i5)') ianatyp,istep,iinc
         write(11,'(i10)') nktrue
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(11,100) i,(v(j,i),j=1,3),zero,zero,zero
         enddo
!
         write(11,'(a5)') m1
      endif
!
!     storing the stresses in the nodes
!
      if(filab(3)(1:4).eq.'S   ') then
!
!        calculating the nodal principal stress
!
         write(11,'(a5)') m1
         write(11,'(a3)') '520'
         write(11,'(i1)') one
         write(11,'(i5,",",i5,",",i5)') ianatyp,istep,iinc
         write(11,'(i10)') nktrue
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            do j=1,6
               str(j)=stn(j,i)
            enddo
            str2(1)=str(1)*str(1)+str(4)*str(4)+str(5)*str(5)
            str2(2)=str(4)*str(4)+str(2)*str(2)+str(6)*str(6)
            str2(3)=str(5)*str(5)+str(6)*str(6)+str(3)*str(3)
c            str2(4)=str(1)*str(4)+str(4)*str(2)+str(5)*str(6)
c            str2(5)=str(1)*str(5)+str(4)*str(6)+str(5)*str(3)
c            str2(6)=str(4)*str(5)+str(2)*str(6)+str(6)*str(3)
            v1=str(1)+str(2)+str(3)
            v2=(v1*v1-str2(1)-str2(2)-str2(3))/2.d0
            v3=str(1)*(str(2)*str(3)-str(6)*str(6))
     &        -str(4)*(str(4)*str(3)-str(5)*str(6))
     &        +str(5)*(str(4)*str(6)-str(5)*str(2))
            bb=v2-v1*v1/3.d0
            cc=-2.d0*v1**3/27.d0+v1*v2/3.d0-v3
            if(dabs(bb).le.1.d-10) then
               stnprin(1)=0.d0
               stnprin(2)=0.d0
               stnprin(3)=0.d0
            else
               cm=2.d0*dsqrt(-bb/3.d0)
               cn=3.d0*cc/(cm*bb)
               if(dabs(cn).gt.1.d0) then
                  if(cn.gt.1.d0) then
                     cn=1.d0
                  else
                     cn=-1.d0
                  endif
               endif
               tt=datan2(dsqrt(1.d0-cn*cn),cn)/3.d0
               stnprin(1)=cm*dcos(tt)
               stnprin(2)=cm*dcos(tt+2.d0*pi/3.d0)
               stnprin(3)=cm*dcos(tt+4.d0*pi/3.d0)
            endif
            do j=1,3
               stnprin(j)=stnprin(j)+v1/3.d0
            enddo
            call dsort(stnprin,iy,n,kflag)
            write(11,101) i,one,(stnprin(j),j=1,3)
         enddo
         write(11,'(a5)') m1
!
!        calculating the elemental principal stress
!
         write(11,'(a5)') m1
         write(11,'(a3)') '530'
         write(11,'(i1)') one
         write(11,'(i5,",",i5,",",i5)') ianatyp,istep,iinc
         write(11,'(i10)') netrue
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
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
            else
               nope=6
            endif
!
!           calculating the eigenvalues of the mean stress
!           tensor. The mean tensor is taken over the nodes
!           belonging to the element
!
            do j=1,6
               str(j)=0.d0
               do k=1,nope
                  str(j)=str(j)+stn(j,kon(indexe+k))
               enddo
               str(j)=str(j)/nope
            enddo
            str2(1)=str(1)*str(1)+str(4)*str(4)+str(5)*str(5)
            str2(2)=str(4)*str(4)+str(2)*str(2)+str(6)*str(6)
            str2(3)=str(5)*str(5)+str(6)*str(6)+str(3)*str(3)
c            str2(4)=str(1)*str(4)+str(4)*str(2)+str(5)*str(6)
c            str2(5)=str(1)*str(5)+str(4)*str(6)+str(5)*str(3)
c            str2(6)=str(4)*str(5)+str(2)*str(6)+str(6)*str(3)
            v1=str(1)+str(2)+str(3)
            v2=(v1*v1-str2(1)-str2(2)-str2(3))/2.d0
            v3=str(1)*(str(2)*str(3)-str(6)*str(6))
     &        -str(4)*(str(4)*str(3)-str(5)*str(6))
     &        +str(5)*(str(4)*str(6)-str(5)*str(2))
            bb=v2-v1*v1/3.d0
            cc=-2.d0*v1**3/27.d0+v1*v2/3.d0-v3
            if(dabs(bb).le.1.d-10) then
               stnprin(1)=0.d0
               stnprin(2)=0.d0
               stnprin(3)=0.d0
            else
               cm=2.d0*dsqrt(-bb/3.d0)
               cn=3.d0*cc/(cm*bb)
               if(dabs(cn).gt.1.d0) then
                  if(cn.gt.1.d0) then
                     cn=1.d0
                  else
                     cn=-1.d0
                  endif
               endif
               tt=datan2(dsqrt(1.d0-cn*cn),cn)/3.d0
               stnprin(1)=cm*dcos(tt)
               stnprin(2)=cm*dcos(tt+2.d0*pi/3.d0)
               stnprin(3)=cm*dcos(tt+4.d0*pi/3.d0)
            endif
            call dsort(stnprin,iy,n,kflag)
            do j=1,3
               stnprin(j)=stnprin(j)+v1/3.d0
            enddo
            write(11,101) i,one,(stnprin(j),j=1,3)
         enddo
!
         write(11,'(a5)') m1
      endif
!
      if(filab(7)(1:4).eq.'ENER') then
!
!        calculating the energy
!
         write(11,'(a5)') m1
         write(11,'(a3)') '540'
         write(11,'(i1)') one
         write(11,'(i5,",",i5,",",i5)') ianatyp,istep,iinc
         write(11,'(i10)') netrue
!
!        calculating the total energy
!
         totenergy=0.d0
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            indexe=ipkon(i)
!
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
            else
               nope=6
            endif
!
            do j=1,nope
               konl(j)=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl(j))
               enddo
            enddo
!
            if(lakon(i)(4:5).eq.'8R') then
               mint3d=1
            elseif((lakon(i)(4:4).eq.'8').or.
     &              (lakon(i)(4:6).eq.'20R')) then
               mint3d=8
            elseif(lakon(i)(4:4).eq.'2') then
               mint3d=27
            elseif(lakon(i)(4:5).eq.'10') then
               mint3d=4
            elseif(lakon(i)(4:4).eq.'4') then
               mint3d=1
            elseif(lakon(i)(4:5).eq.'15') then
               mint3d=9
            else
               mint3d=2
            endif
!
            do jj=1,mint3d
               if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,jj)
                  et=gauss3d1(2,jj)
                  ze=gauss3d1(3,jj)
                  weight=weight3d1(jj)
               elseif((lakon(i)(4:4).eq.'8').or.
     &                 (lakon(i)(4:6).eq.'20R'))
     &                 then
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
               else
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
     &                    (lakon(i)(7:7).eq.'S')) then
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
               totenergy=totenergy+weight*ener(jj,i)*xsj
            enddo
         enddo
!
!        calculating the element energy...
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            indexe=ipkon(i)
!
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
            else
               nope=6
            endif
!
            do j=1,nope
               konl(j)=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl(j))
               enddo
            enddo
!
            energy=0.d0
            volume=0.d0
!
            if(lakon(i)(4:5).eq.'8R') then
               mint3d=1
            elseif((lakon(i)(4:4).eq.'8').or.
     &              (lakon(i)(4:6).eq.'20R')) then
               mint3d=8
            elseif(lakon(i)(4:4).eq.'2') then
               mint3d=27
            elseif(lakon(i)(4:5).eq.'10') then
               mint3d=4
            elseif(lakon(i)(4:4).eq.'4') then
               mint3d=1
            elseif(lakon(i)(4:5).eq.'15') then
               mint3d=9
            else
               mint3d=2
            endif
!
            do jj=1,mint3d
               if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,jj)
                  et=gauss3d1(2,jj)
                  ze=gauss3d1(3,jj)
                  weight=weight3d1(jj)
               elseif((lakon(i)(4:4).eq.'8').or.
     &                 (lakon(i)(4:6).eq.'20R'))
     &                 then
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
               else
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
     &                    (lakon(i)(7:7).eq.'S')) then
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
               energy=energy+weight*xsj*ener(jj,i)
               volume=volume+weight*xsj
            enddo
            write(11,102) i,energy,energy/totenergy,energy/volume
         enddo
         write(11,'(a5)') m1
      endif
  99  format(e12.5,3(",",e12.5))
 100  format(i10,6(",",e12.5))
 101  format(i10,",",i5,3(",",e12.5))
 102  format(i10,3(",",e12.5))
!
      return
      end
