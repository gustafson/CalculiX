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
      subroutine allocation(nload,nforc,nboun,nk,ne,nmpc,
     &  nset,nalset,nmat,ntmat,npmat,norien,nam,nprint,
     &  mi,ntrans,set,meminset,rmeminset,ncs,
     &  namtot,ncmat,memmpc,ne1d,ne2d,nflow,jobnamec,irstrt,
     &  ithermal,nener,nstate,irestartstep,inpc,ipoinp,inp,
     &  ntie,nbody,nprop,ipoinpc,nevdamp,npt,nslavs,nkon)
!
!     calculates a conservative estimate of the size of the 
!     fields to be allocated
!
!     the underscores were dropped since they caused problems in the
!     DDD debugger. 
!
!     meminset=total # of terms in sets
!     rmeminset=total # of reduced terms (due to use of generate) in
!               sets
!
      implicit none
!
      logical igen,lin,frequency,isochoric,cyclicsymmetry,composite
!
      character*1 selabel,sulabel,inpc(*)
      character*5 llab
      character*8 label
      character*20 mpclabel
      character*81 set(*),noset,elset,slavset,mastset,noelset,
     &  surface
      character*132 jobnamec(*),textpart(16)
!
      integer nload,nforc,nboun,nk,ne,nmpc,nset,nalset,
     &  nmat,ntmat,npmat,norien,nam,nprint,kode,iline,
     &  istat,n,key,meminset(*),i,js,inoset,mi(*),ii,ipol,inl,
     &  ibounstart,ibounend,ibound,ntrans,ntmatl,npmatl,ityp,l,
     &  ielset,nope,nteller,nterm,ialset(16),ncs,rmeminset(*),
     &  islavset,imastset,namtot,ncmat,nconstants,memmpc,j,ipos,
     &  maxrmeminset,ne1d,ne2d,necper,necpsr,necaxr,nesr,
     &  neb32,nn,nflow,nradiate,irestartread,irestartstep,icntrl,
     &  irstrt,ithermal(2),nener,nstate,ipoinp(2,*),inp(3,*),
     &  ntie,nbody,nprop,ipoinpc(0:*),nevdamp,npt,
     &  iposs,iposm,nslavs,nlayer,nkon,nopeexp
!
      real*8 temperature,tempact,xfreq,tpinc,tpmin,tpmax
!
      integer nentries
      parameter(nentries=14)
!
!     in the presence of mechanical steps the highest number
!     of DOF is at least 3
!
      if(ithermal(2).ne.2) mi(2)=3
!
!     initialisation of ipoinp
!
      do i=1,nentries
         if(ipoinp(1,i).ne.0) then
            ipol=i
            inl=ipoinp(1,i)
            iline=inp(1,inl)-1
            exit
         endif
      enddo
!
      istat=0
!
      nset=0
      maxrmeminset=0
      necper=0
      necpsr=0
      necaxr=0
      nesr=0
      neb32=0
      nradiate=0
      nkon=0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      do
         if(istat.lt.0) then
            exit
         endif
!
         if(textpart(1)(1:10).eq.'*AMPLITUDE') then
            nam=nam+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               namtot=namtot+4
            enddo
         elseif(textpart(1)(1:12).eq.'*BEAMSECTION') then
            mi(3)=max(mi(3),2)
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:9).eq.'*BOUNDARY') then
            nam=nam+1
            namtot=namtot+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
!
               read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
               if(textpart(3)(1:1).eq.' ') then
                  ibounend=ibounstart
               else
                  read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               endif
               ibound=ibounend-ibounstart+1
               ibound=max(1,ibound)
               ibound=min(3,ibound)
!
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.eq.0) then
                  nboun=nboun+ibound
                  if(ntrans.gt.0) then
                     nmpc=nmpc+ibound
                     memmpc=memmpc+4*ibound
                     nk=nk+1
                  endif
               else
                  read(textpart(1)(1:80),'(a80)',iostat=istat) noset
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  do i=1,nset
                     if(set(i).eq.noset) then
                        nboun=nboun+ibound*meminset(i)
                        if(ntrans.gt.0)then
                           nmpc=nmpc+ibound*meminset(i)
                           memmpc=memmpc+4*ibound*meminset(i)
                           nk=nk+meminset(i)
                        endif
                        exit
                     endif
                  enddo
               endif
            enddo
         elseif(textpart(1)(1:6).eq.'*CFLUX') then
            nam=nam+1
            namtot=namtot+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
!
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.eq.0) then
                  nforc=nforc+1
               else
                  read(textpart(1)(1:80),'(a80)',iostat=istat) noset
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  do i=1,nset
                     if(set(i).eq.noset) then
                        nforc=nforc+meminset(i)
                        exit
                     endif
                  enddo
               endif
            enddo
         elseif(textpart(1)(1:6).eq.'*CLOAD') then
            nam=nam+1
            namtot=namtot+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
!
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.eq.0) then
                  if(ntrans.eq.0) then
                     nforc=nforc+1
                  else
                     nforc=nforc+3
                  endif
               else
                  read(textpart(1)(1:80),'(a80)',iostat=istat) noset
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  do i=1,nset
                     if(set(i).eq.noset) then
                        if(ntrans.eq.0) then
                           nforc=nforc+meminset(i)
                        else
                           nforc=nforc+3*meminset(i)
                        endif
                        exit
                     endif
                  enddo
               endif
            enddo
         elseif((textpart(1)(1:13).eq.'*CONDUCTIVITY').or.
     &          (textpart(1)(1:8).eq.'*DENSITY').or.
     &          (textpart(1)(1:10).eq.'*EXPANSION').or.
     &          (textpart(1)(1:15).eq.'*FLUIDCONSTANTS').or.
     &          (textpart(1)(1:13).eq.'*SPECIFICHEAT')) then
            ntmatl=0
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               ntmatl=ntmatl+1
               ntmat=max(ntmatl,ntmat)
            enddo
         elseif(textpart(1)(1:15).eq.'*CONTACTDAMPING') then
            ncmat=max(5,ncmat)
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:12).eq.'*CONTACTPAIR') then
            ntie=ntie+1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:13).eq.'*CONTACTPRINT') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               nprint=nprint+n
            enddo
         elseif(textpart(1)(1:6).eq.'*CREEP') then
            ncmat=max(9,ncmat)
            npmat=max(2,npmat)
            if(ncmat.ge.9) ncmat=max(19,ncmat)
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:16).eq.'*CYCLICHARDENING') then
            ntmatl=0
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &                 temperature
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if(ntmatl.eq.0) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               elseif(temperature.ne.tempact) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               endif
               npmatl=npmatl+1
               npmat=max(npmatl,npmat)
            enddo
         elseif(textpart(1)(1:20).eq.'*CYCLICSYMMETRYMODEL') then
!
!           possible MPC's: static temperature, displacements(velocities)
!           and static pressure
!
            nk=nk+1
            nmpc=nmpc+5*ncs
            memmpc=memmpc+125*ncs
            ntrans=ntrans+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
            enddo
         elseif(textpart(1)(1:8).eq.'*DASHPOT') then
            nmat=nmat+1
            frequency=.false.
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &           xfreq
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            if(xfreq.gt.0.d0) frequency=.true.
            iline=iline-1
            if(.not.frequency) then
               ntmatl=0
               ncmat=max(2,ncmat)
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
               enddo
            else
               ntmatl=0
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &                 temperature
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  if(ntmatl.eq.0) then
                     npmatl=0
                     ntmatl=ntmatl+1
                     ntmat=max(ntmatl,ntmat)
                     tempact=temperature
                  elseif(temperature.ne.tempact) then
                     npmatl=0
                     ntmatl=ntmatl+1
                     ntmat=max(ntmatl,ntmat)
                     tempact=temperature
                  endif
                  npmatl=npmatl+1
                  npmat=max(npmatl,npmat)
               enddo
               if(ncmat.ge.9) ncmat=max(19,ncmat)
            endif
         elseif(textpart(1)(1:22).eq.'*DEFORMATIONPLASTICITY') then
            ncmat=max(5,ncmat)
            ntmatl=0
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               ntmatl=ntmatl+1
               ntmat=max(ntmatl,ntmat)
            enddo
         elseif(textpart(1)(1:7).eq.'*DEPVAR') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.lt.0) exit
               nstate=max(l,nstate)
            enddo
         elseif(textpart(1)(1:20).eq.'*DISTRIBUTEDCOUPLING') then
            do i=2,n
               if(textpart(i)(1:8).eq.'SURFACE=') then
                  surface=textpart(i)(9:88)
                  ipos=index(surface,' ')
                  surface(ipos:ipos)='T'
                  exit
               endif
            enddo
            do i=1,nset
               if(set(i).eq.surface) then
!
!                 worst case: 8 nodes per element face
!
                  nk=nk+8*meminset(i)
!
!                 1 distributed coupling MPC
!                  
                  nmpc=nmpc+1
!
!                 3 terms * # of nodes +1 parallel to coupling
!                 direction
!
                  memmpc=memmpc+24*meminset(i)+1
                  exit
!
               endif
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif((textpart(1)(1:6).eq.'*DLOAD').or.
     &           (textpart(1)(1:6).eq.'*DFLUX').or.
     &           (textpart(1)(1:5).eq.'*FILM')) then
            if(textpart(1)(1:5).ne.'*FILM') then
               nam=nam+1
               namtot=namtot+1
            else
               nam=nam+2
               namtot=namtot+2
            endif
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(2)(1:5),'(a5)',iostat=istat) llab
               if((llab.eq.'GRAV ').or.(llab.eq.'CENTR').or.
     &            (llab.eq.'NEWTO')) then
                  nbody=nbody+1
                  cycle
               endif
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.eq.0) then
                  nload=nload+1
               else
                  read(textpart(1)(1:80),'(a80)',iostat=istat) elset
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  do i=1,nset
                     if(set(i).eq.elset) then
                        nload=nload+meminset(i)
                        exit
                     endif
                  enddo
               endif
            enddo
         elseif((textpart(1)(1:8).eq.'*DYNAMIC').or.
     &      (textpart(1)(1:32).eq.'*COUPLEDTEMPERATURE-DISPLACEMENT'))
     &          then
            if((mi(1).eq.1).or.(mi(1).eq.8).or.(mi(1).eq.27)) then
               mi(1)=27
            elseif(mi(1).eq.4) then
               mi(1)=15
            else
               mi(1)=18
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:8).eq.'*ELPRINT') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               nprint=nprint+n
            enddo
         elseif(textpart(1)(1:8).eq.'*ELASTIC') then
            ntmatl=0
            ityp=2
            ncmat=max(2,ncmat)
            do i=2,n
               if(textpart(i)(1:5).eq.'TYPE=') then
                  if(textpart(i)(6:8).eq.'ISO') then
                     ityp=2
                     ncmat=max(2,ncmat)
                  elseif((textpart(i)(6:10).eq.'ORTHO').or.
     &                   (textpart(i)(6:10).eq.'ENGIN')) then
                     ityp=9
                     ncmat=max(9,ncmat)
                  elseif(textpart(i)(6:10).eq.'ANISO') then
                     ityp=21
                     ncmat=max(21,ncmat)
                  endif
                  exit
               endif
            enddo
            if(ityp.eq.2) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
               enddo
            elseif(ityp.eq.9) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  iline=iline+1
               enddo
            elseif(ityp.eq.21) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  iline=iline+2
               enddo
            endif
         elseif((textpart(1)(1:8).eq.'*ELEMENT').and.
     &          (textpart(1)(1:14).ne.'*ELEMENTOUTPUT')) then
            ielset=0
            isochoric=.false.
!
            loop1: do i=2,n
               if(textpart(i)(1:6).eq.'ELSET=') then
                  elset=textpart(i)(7:86)
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  ielset=1
                  do js=1,nset
                     if(set(js).eq.elset) exit
                  enddo
                  if(js.gt.nset) then
                     nset=nset+1
                     set(nset)=elset
                  endif
               elseif(textpart(i)(1:5).eq.'TYPE=') then
                  read(textpart(i)(6:13),'(a8)') label
                  if(label.eq.'        ') then
                     write(*,*) 
     &                   '*ERROR in allocation: element type is lacking'
                     write(*,*) '       '
                     call inputerror(inpc,ipoinpc,iline)
                     stop
                  endif
                  if(label(1:2).eq.'DC') then
                     label(1:7)=label(2:8)
                     label(8:8)=' '
                  endif
                  if((label.eq.'C3D20   ').or.
     &               (label.eq.'F3D20   ')) then
                     mi(1)=max(mi(1),27)
                     nope=20
                     nopeexp=20
                  elseif((label(1:8).eq.'C3D20R  ').or.
     &                   (label(1:8).eq.'F3D20R  ')) then
                     mi(1)=max(mi(1),8)
                     nope=20
                     nopeexp=20
                  elseif(label(1:8).eq.'C3D20RI ') then
                     mi(1)=max(mi(1),8)
                     nope=20
                     nopeexp=20
                     isochoric=.true.
                  elseif((label.eq.'C3D8R   ').or.
     &                   (label.eq.'F3D8R   ')) then
                     mi(1)=max(mi(1),1)
                     nope=8
                     nopeexp=8
                  elseif((label.eq.'C3D10   ').or.
     &                   (label.eq.'F3D10   ')) then
                     mi(1)=max(mi(1),4)
                     nope=10
                     nopeexp=10
                  elseif((label.eq.'C3D4    ').or.
     &                   (label.eq.'F3D4    ')) then
                     mi(1)=max(mi(1),1)
                     nope=4
                     nopeexp=4
                  elseif((label.eq.'C3D15   ').or.
     &                   (label.eq.'F3D15   ')) then
                     mi(1)=max(mi(1),9)
                     nope=15
                     nopeexp=15
                  elseif((label.eq.'C3D6    ').or.
     &                   (label.eq.'F3D6    ')) then
                     mi(1)=max(mi(1),2)
                     nope=6
                     nopeexp=6
                  elseif(label.eq.'F3D6R   ') then
                     mi(1)=max(mi(1),1)
                     nope=6
                     nopeexp=6
                  elseif((label.eq.'CPE8R   ').or.
     &                    (label.eq.'CPS8R   ').or.
     &                    (label.eq.'CAX8R   ').or.
     &                    (label.eq.'S8R     ')) then
                     mi(1)=max(mi(1),8)
                     nope=8
                     nopeexp=28
                  elseif((label.eq.'C3D8    ').or.
     &                    (label.eq.'F3D8    ')) then
                     mi(1)=max(mi(1),8)
                     nope=8
                     nopeexp=8
c    Bernhardi start
                  elseif(label.eq.'C3D8I   ') then
                     mi(1)=max(mi(1),8)
                     nope=8
                     nopeexp=11
c    Bernhardi end
                  elseif((label.eq.'CPE8    ').or.
     &                    (label.eq.'CPS8    ').or.
     &                    (label.eq.'CAX8    ').or.
     &                    (label.eq.'S8      ')) then
                     mi(1)=max(mi(1),27)
                     nope=8
                     nopeexp=28
                  elseif((label.eq.'CPE6    ').or.
     &                    (label.eq.'CPS6    ').or.
     &                    (label.eq.'CAX6    ').or.
     &                    (label.eq.'S6      ')) then
                     mi(1)=max(mi(1),9)
                     nope=6
                     nopeexp=21
                  elseif(label.eq.'SB4     ') then
                     mi(2)=max(mi(2),7)
                     nope=4
                  elseif(label.eq.'B32     ') then
                     mi(1)=max(mi(1),27)
                     nope=3
                     nopeexp=23
                  elseif(label.eq.'B32R    ') then
                     mi(1)=max(mi(1),8)
                     nope=3
                     nopeexp=23
                  elseif(label(1:8).eq.'DASHPOTA') then
                     label='EDSHPTA2'
                     nope=2
                     nopeexp=2
                  elseif(label(1:1).eq.'D') then
                     nope=3
                     nopeexp=3
                     mi(2)=max(3,mi(2))
                  elseif(label(1:7).eq.'SPRINGA') then
                     mi(1)=max(mi(1),1)
                     label='ESPRNGA2'
                     nope=2
                     nopeexp=2
                  elseif(label.eq.'GAPUNI  ') then
                     nope=2
                     nopeexp=2
                  endif
                  if(label(1:1).eq.'F') mi(2)=max(mi(2),4)
               endif
            enddo loop1
!
            loop2:do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(1)(1:10),'(i10)',iostat=istat) i
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
c    Bernhardi start
c              space for incompatible mode nodes
               if(label(1:5).eq.'C3D8I') then
                  nk=nk+3
               endif
c    Bernhardi end
               if(label(1:2).ne.'C3') then
                  if(label(1:3).eq.'CPE') then
                     necper=necper+1
                  elseif(label(1:2).eq.'CP') then
                     necpsr=necpsr+1
                  elseif(label(1:1).eq.'C') then
                     necaxr=necaxr+1
                  elseif(label(1:1).eq.'S') then
                     nesr=nesr+1
                  elseif(label(1:1).eq.'B') then
                     neb32=neb32+1
                  elseif(label(1:1).eq.'D') then
                     nflow=nflow+1
                  endif
               endif
               nteller=n-1
               if(nteller.lt.nope) then
                  do
                     call getnewline(inpc,textpart,istat,n,key,iline,
     &                    ipol,inl,ipoinp,inp,ipoinpc)
                     if((istat.lt.0).or.(key.eq.1)) exit loop2
                     if(nteller+n.gt.nope) n=nope-nteller
                     nteller=nteller+n
                     if(nteller.eq.nope) exit
                  enddo
               endif
               ne=max(ne,i)
               nkon=nkon+nopeexp
               if(ielset.eq.1) then
                  meminset(js)=meminset(js)+1
                  rmeminset(js)=rmeminset(js)+1
               endif
!
!              up to 8 new mpc's with 22 terms in each mpc
!              (21 = 7 nodes x 3 dofs + inhomogeneous term)
!
               if(isochoric) then
                  nmpc=nmpc+8
                  nk=nk+8
                  nboun=nboun+8
                  memmpc=memmpc+176
               endif
            enddo loop2
         elseif((textpart(1)(1:5).eq.'*NSET').or.
     &           (textpart(1)(1:6).eq.'*ELSET')) then
            if(textpart(1)(1:5).eq.'*NSET')
     &           then
               noelset=textpart(2)(6:85)
               noelset(81:81)=' '
               ipos=index(noelset,' ')
               noelset(ipos:ipos)='N'
               kode=0
            else
               noelset=textpart(2)(7:86)
               noelset(81:81)=' '
               ipos=index(noelset,' ')
               noelset(ipos:ipos)='E'
               kode=1
            endif
!
!              check whether new set name or old one
!
            do js=1,nset
               if(set(js).eq.noelset) exit
            enddo
            if(js.gt.nset) then
               nset=nset+1
               set(nset)=noelset
               nn=nset
            else
               nn=js
            endif
!
            if((n.gt.2).and.(textpart(3)(1:8).eq.'GENERATE')) then
               igen=.true.
            else
               igen=.false.
            endif
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               if(igen) then
                  if(textpart(2)(1:1).eq.' ')
     &                 textpart(2)=textpart(1)
                  if(textpart(3)(1:1).eq.' ')
     &                 textpart(3)='1        '
                  do i=1,3
                     read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &                         ialset(i)
                     if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  enddo
                  meminset(nn)=meminset(nn)+
     &                 (ialset(2)-ialset(1))/ialset(3)+1
                  rmeminset(nn)=rmeminset(nn)+3
               else
                  do i=1,n
                     read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &                         ialset(i)
                     if(istat.gt.0) then
                        noelset=textpart(i)(1:80)
                        noelset(81:81)=' '
                        ipos=index(noelset,' ')
                        if(kode.eq.0) then
                           noelset(ipos:ipos)='N'
                        else
                           noelset(ipos:ipos)='E'
                        endif
                        do j=1,nset
                           if(noelset.eq.set(j)) then
                              meminset(nn)=meminset(nn)+
     &                             meminset(j)
                              rmeminset(nn)=rmeminset(nn)+
     &                             rmeminset(j)
                              exit
                           endif
                        enddo
                     else
                        meminset(nn)=meminset(nn)+1
                        rmeminset(nn)=rmeminset(nn)+1
                     endif
                  enddo
               endif
            enddo
         elseif(textpart(1)(1:9).eq.'*EQUATION') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit 
               read(textpart(1)(1:10),'(i10)',iostat=istat) nterm
               if(ntrans.eq.0) then
                  nmpc=nmpc+1
                  memmpc=memmpc+nterm
               else
                  nmpc=nmpc+3
                  memmpc=memmpc+3*nterm
               endif
               ii=0
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ii=ii+n/3
                  if(ii.eq.nterm) exit
               enddo
            enddo
         elseif(textpart(1)(1:13).eq.'*FACEPRINT') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               nprint=nprint+n
            enddo
         elseif(textpart(1)(1:5).eq.'*GAP ') then
            elset='                    '
            do i=2,n
               if(textpart(i)(1:6).eq.'ELSET=') then
                  elset=textpart(i)(7:86)
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  exit
               endif
            enddo
            if(elset.ne.'                     ') then
               do i=1,nset
                  if(set(i).eq.elset) then
                     nk=nk+2+2*meminset(i)
                     nmpc=nmpc+meminset(i)
                     nboun=nboun+2*meminset(i)
                     memmpc=memmpc+8*meminset(i)
                  endif
               enddo
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:13).eq.'*FLUIDSECTION') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               nprop=nprop+8
            enddo
         elseif(textpart(1)(1:9).eq.'*FRICTION') then
            ncmat=max(7,ncmat)
            ntmat=max(1,ntmat)
c            nstate=max(nstate,9)
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:15).eq.'*GAPCONDUCTANCE') then
            nmat=nmat+1
            ntmatl=0
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &              inl,ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &              temperature
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if(ntmatl.eq.0) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               elseif(temperature.ne.tempact) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               endif
               npmatl=npmatl+1
               npmat=max(npmatl,npmat)
            enddo
         elseif(textpart(1)(1:13).eq.'*HYPERELASTIC') then
            ntmatl=0
            ityp=-7
            do i=2,n
               if(textpart(i)(1:12).eq.'ARRUDA-BOYCE') then
                  ityp=-1
                  ncmat=max(3,ncmat)
               elseif(textpart(i)(1:13).eq.'MOONEY-RIVLIN') then
                  ityp=-2
                  ncmat=max(3,ncmat)
               elseif(textpart(i)(1:8).eq.'NEOHOOKE') then
                  ityp=-3
                  ncmat=max(2,ncmat)
               elseif(textpart(i)(1:5).eq.'OGDEN') then
                  ityp=-4
                  ncmat=max(3,ncmat)
               elseif(textpart(i)(1:10).eq.'POLYNOMIAL') then
                  ityp=-7
                  ncmat=max(3,ncmat)
               elseif(textpart(i)(1:17).eq.'REDUCEDPOLYNOMIAL')
     &                 then
                  ityp=-10
                  ncmat=max(2,ncmat)
               elseif(textpart(i)(1:11).eq.'VANDERWAALS') then
                  ityp=-13
                  ncmat=max(5,ncmat)
               elseif(textpart(i)(1:4).eq.'YEOH') then
                  ityp=-14
                  ncmat=max(6,ncmat)
               elseif(textpart(i)(1:2).eq.'N=') then
                  if(textpart(i)(3:3).eq.'1') then
                  elseif(textpart(i)(3:3).eq.'2') then
                     if(ityp.eq.-4) then
                        ityp=-5
                        ncmat=max(6,ncmat)
                     elseif(ityp.eq.-7) then
                        ityp=-8
                        ncmat=max(7,ncmat)
                     elseif(ityp.eq.-10) then
                        ityp=-11
                        ncmat=max(4,ncmat)
                     endif
                  elseif(textpart(i)(3:3).eq.'3') then
                     if(ityp.eq.-4) then
                        ityp=-6
                        ncmat=max(9,ncmat)
                     elseif(ityp.eq.-7) then
                        ityp=-9
                        ncmat=max(12,ncmat)
                     elseif(ityp.eq.-10) then
                        ityp=-12
                        ncmat=max(6,ncmat)
                     endif
                  endif
               endif
            enddo
            if((ityp.ne.-6).and.(ityp.ne.-9)) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmat,ntmatl)
               enddo
            else
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmat,ntmatl)
                  iline=iline+1
               enddo
            endif
         elseif(textpart(1)(1:10).eq.'*HYPERFOAM') then
            ntmatl=0
            ityp=-15
            ncmat=max(3,ncmat)
            do i=2,n
               if(textpart(i)(1:2).eq.'N=') then
                  if(textpart(i)(3:3).eq.'1') then
                  elseif(textpart(i)(3:3).eq.'2') then
                     ityp=-16
                     ncmat=max(6,ncmat)
                  elseif(textpart(i)(3:3).eq.'3') then
                     ityp=-17
                     ncmat=max(9,ncmat)
                  endif
               endif
            enddo
            if(ityp.ne.-17) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmat,ntmatl)
               enddo
            else
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmat,ntmatl)
                  iline=iline+1
               enddo
            endif
         elseif(textpart(1)(1:9).eq.'*MATERIAL') then
            nmat=nmat+1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:13).eq.'*MODALDAMPING') then
            if(textpart(2)(1:8).ne.'RAYLEIGH') then
               nevdamp=0
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  read(textpart(1)(1:10),'(i10)',iostat=istat) i
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  nevdamp = max(nevdamp,i)
                  read(textpart(2)(1:10),'(i10)',iostat=istat) i
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  nevdamp = max(nevdamp,i)
               enddo
            else
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
            endif
         elseif(textpart(1)(1:4).eq.'*MPC') then
            mpclabel='                    '
            if((mpclabel(1:8).ne.'STRAIGHT').and.
     &           (mpclabel(1:4).ne.'PLANE')) then
               nk=nk+1
               nmpc=nmpc+1
               nboun=nboun+1
               memmpc=memmpc+1
            endif
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               do i=1,n
                  read(textpart(i)(1:10),'(i10)',iostat=istat) ialset(i)
                  if(mpclabel.eq.'                    ')
     &                 mpclabel=textpart(i)(1:20)
                  if(istat.gt.0) then
                     noelset=textpart(i)(1:80)
                     noelset(81:81)=' '
                     ipos=index(noelset,' ')
                     noelset(ipos:ipos)='N'
                     do j=1,nset
                        if(noelset.eq.set(j)) then
                           if(mpclabel(1:8).eq.'STRAIGHT') then
                              nk=nk+2*meminset(j)
                              nmpc=nmpc+2*meminset(j)
                              nboun=nboun+2*meminset(j)
                              memmpc=memmpc+14*meminset(j)
                           elseif(mpclabel(1:5).eq.'PLANE') then
                              nk=nk+meminset(j)
                              nmpc=nmpc+meminset(j)
                              nboun=nboun+meminset(j)
                              memmpc=memmpc+13*meminset(j)
                           else
                              memmpc=memmpc+meminset(j)
                           endif
                           exit
                        endif
                     enddo
                  else
                     if(mpclabel(1:8).eq.'STRAIGHT') then
                        nk=nk+2
                        nmpc=nmpc+2
                        nboun=nboun+2
                        memmpc=memmpc+14
                     elseif(mpclabel(1:5).eq.'PLANE') then
                        nk=nk+1
                        nmpc=nmpc+1
                        nboun=nboun+1
                        memmpc=memmpc+13
                     else
                        memmpc=memmpc+1
                     endif
                  endif
               enddo
            enddo
         elseif((textpart(1)(1:5).eq.'*NODE').and.
     &           (textpart(1)(1:10).ne.'*NODEPRINT').and.
     &           (textpart(1)(1:9).ne.'*NODEFILE').and.
     &           (textpart(1)(1:11).ne.'*NODEOUTPUT')) then
            inoset=0
            loop3: do i=2,n
               if(textpart(i)(1:5).eq.'NSET=') then
                  noset=textpart(i)(6:85)
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  inoset=1
                  do js=1,nset
                     if(set(js).eq.noset) exit
                  enddo
                  if(js.gt.nset) then
                     nset=nset+1
                     set(nset)=noset
                  endif
               endif
            enddo loop3
!
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(1)(1:10),'(i10)',iostat=istat) i
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               nk=max(nk,i)
               if(inoset.eq.1) then
                  meminset(js)=meminset(js)+1
                  rmeminset(js)=rmeminset(js)+1
               endif
            enddo
         elseif(textpart(1)(1:10).eq.'*NODEPRINT') then
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               nprint=nprint+n
            enddo
         elseif(textpart(1)(1:12).eq.'*ORIENTATION') then
            norien=norien+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
            enddo
         elseif(textpart(1)(1:8).eq.'*PLASTIC') then
            ntmatl=0
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &                temperature
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if(ntmatl.eq.0) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               elseif(temperature.ne.tempact) then
                  npmatl=0
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
                  tempact=temperature
               endif
               npmatl=npmatl+1
               npmat=max(npmatl,npmat)
            enddo
            if(ncmat.ge.9) ncmat=max(19,ncmat)
         elseif(textpart(1)(1:19).eq.'*PRE-TENSIONSECTION') then
            do i=2,n
               if(textpart(i)(1:8).eq.'SURFACE=') then
                  surface=textpart(i)(9:88)
                  ipos=index(surface,' ')
                  surface(ipos:ipos)='T'
                  exit
               endif
            enddo
            do i=1,nset
               if(set(i).eq.surface) then
!
!                 worst case: 8 nodes per element face
!
                  nk=nk+8*meminset(i)
                  npt=npt+8*meminset(i)
!
!                 2 MPC's per node perpendicular to tension direction
!                 + 1 MPC in tension direction
!                  
                  nmpc=nmpc+16*meminset(i)+1
!
!                 6 terms per MPC perpendicular to tension direction
!                 + 6 terms * # of nodes +1 parallel to tension
!                 direction
!
                  memmpc=memmpc+96*meminset(i)+48*meminset(i)+1
                  exit
!
               endif
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:8).eq.'*RADIATE') then
            nam=nam+2
            namtot=namtot+2
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(2)(1:5),'(a5)',iostat=istat) llab
               if((llab.eq.'GRAV ').or.(llab.eq.'CENTR')) exit
               read(textpart(1)(1:10),'(i10)',iostat=istat) l
               if(istat.eq.0) then
                  nload=nload+1
                  nradiate=nradiate+1
               else
                  read(textpart(1)(1:80),'(a80)',iostat=istat) elset
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  do i=1,nset
                     if(set(i).eq.elset) then
                        nload=nload+meminset(i)
                        nradiate=nradiate+meminset(i)
                        exit
                     endif
                  enddo
               endif
            enddo
         elseif(textpart(1)(1:8).eq.'*RESTART') then
            irestartread=0
            irestartstep=0
            do i=1,n
               if(textpart(i)(1:4).eq.'READ') then
                  irestartread=1
c                  if(irestartstep.eq.0) irestartstep=1
               endif
               if(textpart(i)(1:5).eq.'STEP=') then
                  read(textpart(i)(6:15),'(i10)',iostat=istat) 
     &                  irestartstep
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               endif
            enddo
            if(irestartread.eq.1) then
               icntrl=1
               call restartshort(nset,nload,nbody,nforc,nboun,nk,ne,
     &              nmpc,nalset,nmat,ntmat,npmat,norien,nam,nprint,
     &              mi,ntrans,ncs,namtot,ncmat,memmpc,
     &              ne1d,ne2d,nflow,set,meminset,rmeminset,jobnamec,
     &              irestartstep,icntrl,ithermal,nener,nstate,ntie,
     &              nslavs,nkon)
               irstrt=-1
            else
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:10).eq.'*RIGIDBODY') then
            noset='
     &                           '
            elset='
     &                           '
            do i=2,n
               if(textpart(i)(1:5).eq.'NSET=')
     &              then
                  noset=textpart(i)(6:85)
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  exit
               elseif(textpart(i)(1:6).eq.'ELSET=')
     &                 then
                  elset=textpart(i)(7:86)
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  exit
               endif
            enddo
            if(noset(1:1).ne.' ') then
               do i=1,nset
                  if(set(i).eq.noset) then
                     nk=nk+2+meminset(i)
                     nmpc=nmpc+3*meminset(i)
                     memmpc=memmpc+18*meminset(i)
                     nboun=nboun+3*meminset(i)
                  endif
               enddo
            elseif(elset(1:1).ne.' ') then
               do i=1,nset
                  if(set(i).eq.elset) then
                     nk=nk+2+20*meminset(i)
                     nmpc=nmpc+60*meminset(i)
                     memmpc=memmpc+360*meminset(i)
                     nboun=nboun+60*meminset(i)
                  endif
               enddo
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:13).eq.'*SHELLSECTION') then
            composite=.false.
            do i=2,n
               if(textpart(i)(1:9).eq.'COMPOSITE') then
                  composite=.true.
                  nlayer=0
               elseif(textpart(i)(1:6).eq.'ELSET=') then
                  elset=textpart(i)(7:86)
                  elset(81:81)=' '
                  ipos=index(elset,' ')
                  elset(ipos:ipos)='E'
                  do js=1,nset
                     if(set(js).eq.elset) exit
                  enddo
               endif
            enddo
            if(composite) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) then
                     mi(1)=max(mi(1),8*nlayer)
                     mi(3)=max(mi(3),nlayer)
                     if(js.lt.nset) then
                        nk=nk+20*nlayer*meminset(js)
                        nkon=nkon+20*nlayer*meminset(js)
                     endif
                     exit
                  endif
                  nlayer=nlayer+1
               enddo
            else
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &              inl,ipoinp,inp,ipoinpc)
            endif
         elseif(textpart(1)(1:7).eq.'*SPRING') then
            nmat=nmat+1
            lin=.true.
            do i=2,n
               if(textpart(i)(1:9).eq.'NONLINEAR') then
                  lin=.false.
                  exit
               endif
            enddo
            if(lin) then
               ntmatl=0
               ncmat=max(2,ncmat)
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  ntmatl=ntmatl+1
                  ntmat=max(ntmatl,ntmat)
               enddo
            else
               ntmatl=0
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) exit
                  read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &                 temperature
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  if(ntmatl.eq.0) then
                     npmatl=0
                     ntmatl=ntmatl+1
                     ntmat=max(ntmatl,ntmat)
                     tempact=temperature
                  elseif(temperature.ne.tempact) then
                     npmatl=0
                     ntmatl=ntmatl+1
                     ntmat=max(ntmatl,ntmat)
                     tempact=temperature
                  endif
                  npmatl=npmatl+1
                  npmat=max(npmatl,npmat)
               enddo
               if(ncmat.ge.9) ncmat=max(19,ncmat)
            endif
         elseif(textpart(1)(1:9).eq.'*SURFACE ') then
            nset=nset+1
            sulabel='T'
            do i=2,n
               if(textpart(i)(1:5).eq.'NAME=')
     &              then
                  set(nset)=textpart(i)(6:85)
                  set(nset)(81:81)=' '
               elseif(textpart(i)(1:9).eq.'TYPE=NODE') then
                  sulabel='S'
               endif
            enddo
            ipos=index(set(nset),' ')
            set(nset)(ipos:ipos)=sulabel
            if(sulabel.eq.'S') then
               selabel='N'
            else
               selabel='E'
            endif
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               read(textpart(1)(1:10),'(i10)',iostat=istat) ialset(1)
               if(istat.gt.0) then
                  noset=textpart(1)(1:80)
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)=selabel
                  do i=1,nset-1
                     if(set(i).eq.noset) then
                        meminset(nset)=meminset(nset)+meminset(i)
!
!                       surfaces are stored in expanded form 
!                       (no equivalent to generate)
!
                        rmeminset(nset)=rmeminset(nset)+meminset(i)
                     endif
                  enddo
               else
                  meminset(nset)=meminset(nset)+1
                  rmeminset(nset)=rmeminset(nset)+1
               endif
            enddo
         elseif(textpart(1)(1:16).eq.'*SURFACEBEHAVIOR') then
            ncmat=max(2,ncmat)
            ntmat=max(1,ntmat)
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:19).eq.'*SURFACEINTERACTION') then
            nmat=nmat+1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:11).eq.'*TEMPERATURE') then
            nam=nam+1
            namtot=namtot+1
         elseif(textpart(1)(1:4).eq.'*TIE') then
            ntie=ntie+1
            cyclicsymmetry=.false.
            do i=1,n
               if((textpart(i)(1:14).eq.'CYCLICSYMMETRY').or.
     &            (textpart(i)(1:10).eq.'MULTISTAGE')) then
                  cyclicsymmetry=.true.
               endif
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if(.not.cyclicsymmetry) cycle
            if((istat.lt.0).or.(key.eq.1)) cycle
            slavset=textpart(1)(1:80)
            slavset(81:81)=' '
            iposs=index(slavset,' ')
            slavset(iposs:iposs)='S'
            mastset=textpart(2)(1:80)
            mastset(81:81)=' '
            iposm=index(mastset,' ')
            mastset(iposm:iposm)='S'
            islavset=0
            imastset=0
            do i=1,nset
               if(set(i).eq.slavset) then
                  islavset=i
               elseif(set(i).eq.mastset) then
                  imastset=i
               endif
            enddo
            if((islavset.ne.0).and.(imastset.ne.0)) then
               ncs=ncs+max(meminset(islavset),meminset(imastset))
            else
               write(*,*) '*ERROR in allocation: either the slave'
               write(*,*) '       surface or the master surface in a'
               write(*,*) '       cyclic symmetry *TIE option or both'
               write(*,*) '       do not exist or are no nodal surfaces'
               write(*,*) '       slave set:',slavset(iposs-1:iposs-1)
               write(*,*) '       master set:',mastset(iposm-1:iposm-1)
               stop
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         elseif(textpart(1)(1:11).eq.'*TIMEPOINTS') then
            igen=.false.
            nam=nam+1
            do i=2,n
               if(textpart(i)(1:8).eq.'GENERATE') then
                  igen=.true.
                  exit
               endif
            enddo
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               if(igen)then
                  if(n.lt.3)then
                     write(*,*)'*ERROR in allocation: *TIMEPOINTS'
                     call inputerror(inpc,ipoinpc,iline)
                     stop
                  else
                     read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &                 tpmin
                     if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                     read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &                 tpmax
                     if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                     read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &                 tpinc
                     if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
                     if((tpinc.le.0).or.(tpmin.ge.tpmax)) then
                        write(*,*) '*ERROR in allocation: *TIMEPOINTS'
                        call inputerror(inpc,ipoinpc,iline)
                        stop
                     else
                        namtot=namtot+2+INT((tpmax-tpmin)/tpinc)
                     endif

                  endif
               else
                  namtot=namtot+8
               endif
            enddo
         elseif(textpart(1)(1:10).eq.'*TRANSFORM') then
            ntrans=ntrans+1
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
            enddo
         elseif(textpart(1)(1:13).eq.'*USERMATERIAL') then
            ntmatl=0
            do i=2,n
               if(textpart(i)(1:10).eq.'CONSTANTS=') then
                  read(textpart(i)(11:20),'(i10)',iostat=istat) 
     %                  nconstants
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
                  ncmat=max(nconstants,ncmat)
                  exit
               endif
            enddo
            do
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               ntmatl=ntmatl+1
               ntmat=max(ntmatl,ntmat)
               do i=2,(nconstants-1)/8+1
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
               enddo
            enddo
         else
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         endif
      enddo
!
      do i=1,nset
         nalset=nalset+rmeminset(i)
         maxrmeminset=max(maxrmeminset,rmeminset(i))
      enddo
!
!        extra space needed for rearrangement in elements.f and
!        noelsets.f
!
      nalset=nalset+maxrmeminset
!
      nmpc=nmpc+1
      memmpc=memmpc+1
!
      if(irstrt.eq.0) then
         ne1d=neb32
         ne2d=necper+necpsr+necaxr+nesr
      endif
!
!     introducing a fake tie for axisymmetric elements
!     (needed for cavity radiation)
!
      if(necaxr.gt.0) ntie=max(1,ntie)
!
!        providing space for the expansion of shell and beam elements
!        to genuine volume elements
!
      nk=nk+3*8*ne2d+8*3*ne1d
      if(ne1d.gt.0) then
         nboun=nboun*9
         nforc=nforc*9
      elseif(ne2d.gt.0) then
         nboun=nboun*4
         nforc=nforc*4
      endif
!
!     providing for rigid nodes (knots)
!
!     number of knots: 8*ne2d+3*ne1d
!     number of expanded nodes: 3*8*ne2d+8*3*ne1d
!
!     number of extra nodes (1 rotational node and
!     1 expansion node per knot
!     and one inhomogeneous term node per expanded node)
!
      nk=nk+(2+3)*8*ne2d+(2+8)*3*ne1d
!
!     number of equations (3 per expanded node)
!
      nmpc=nmpc+3*(3*8*ne2d+8*3*ne1d)
!
!     number of terms: 7 per equation
!
      memmpc=memmpc+7*3*(3*8*ne2d+8*3*ne1d)
!
!     number of SPC's: 1 per DOF per expanded node
!
      nboun=nboun+3*(3*8*ne2d+8*3*ne1d)
!
!     temperature DOF in knots
!
      nmpc=nmpc+(3*8*ne2d+8*3*ne1d)
      memmpc=memmpc+2*(3*8*ne2d+8*3*ne1d)
!
!        extra MPCs to avoid undefinid rotation of rigid body nodes
!        lying on a line
!
      nmpc=nmpc+3*8*ne2d+8*3*ne1d
      memmpc=memmpc+3*(3*8*ne2d+8*3*ne1d)
!
!        expanding the MPCs: 2-node MPC link (2D elements) or
!        5-node MPC link (1D elements) between nodes defined by
!        the user and generated mid-nodes
!
c      nmpc=nmpc+3*ne1d+8*ne2d
c      memmpc=memmpc+15*ne1d+24*ne2d
!
!        extra nodes for the radiation boundary conditions
!
      nk=nk+nradiate
!
      write(*,*)
      write(*,*) ' The numbers below are estimated upper bounds'
      write(*,*)
      write(*,*) ' number of:'
      write(*,*)
      write(*,*) '  nodes: ',nk
      write(*,*) '  elements: ',ne
      write(*,*) '  one-dimensional elements: ',ne1d
      write(*,*) '  two-dimensional elements: ',ne2d
      write(*,*) '  integration points per element: ',mi(1)
      write(*,*)
      write(*,*) '  distributed facial loads: ',nload
      write(*,*) '  distributed volumetric loads: ',nbody
      write(*,*) '  concentrated loads: ',nforc
      write(*,*) '  single point constraints: ',nboun
      write(*,*) '  multiple point constraints: ',nmpc
      write(*,*) '  terms in all multiple point constraints: ',memmpc
      write(*,*) '  tie constraints: ',ntie
      write(*,*) '  dependent nodes tied by cyclic constraints: ',ncs
      write(*,*) '  dependent nodes in pre-tension constraints: ',npt
      write(*,*)
      write(*,*) '  sets: ',nset
      write(*,*) '  terms in all sets: ',nalset
      write(*,*)
      write(*,*) '  materials: ',nmat
      write(*,*) '  constants per material and temperature: ',ncmat
      write(*,*) '  temperature points per material: ',ntmat
      write(*,*) '  plastic data points per material: ',npmat
      write(*,*)
      write(*,*) '  orientations: ',norien
      write(*,*) '  amplitudes: ',nam
      write(*,*) '  data points in all amplitudes: ',namtot
      write(*,*) '  print requests: ',nprint
      write(*,*) '  transformations: ',ntrans
      write(*,*) '  property cards: ',nprop
      write(*,*)
!
      return
      end











