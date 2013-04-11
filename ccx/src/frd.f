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
      subroutine frd(co,nk,kon,ipkon,lakon,ne0,v,stn,inum,nmethod,
     &  kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,xstaten,
     &  nstate_,istep,iinc,ithermal,qfn,mode,noddiam,trab,inotr,
     &  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
     &  mi,stx,vr,vi,stnr,stni,vmax,stnmax,ngraph,veold,ener,ne,
     &  cs,set,nset,istartset,iendset,ialset,eenmax)
!
!     stores the results in frd format
!
!     iselect selects which nodes are to be stored:
!          iselect=-1 means only those nodes for which inum negative
!                     ist, i.e. network nodes
!          iselect=+1 means only those nodes for which inum positive
!                     ist, i.e. structural nodes
!          iselect=0  means both of the above
!
      implicit none
!
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*5 p0,p1,p2,p3,p4,p5,p6,p8,p10,p11,p12
      character*8 lakon(*),date,newclock,fmat
      character*10 clock
      character*12 description
      character*20 newdate
      character*80 matname(*)
      character*81 set(*)
      character*87 filab(*)
      character*132 text
!
      integer kon(*),inum(*),nk,ne0,nmethod,kode,i,j,ipkon(*),indexe,
     &  one,ielmat(*),nstate_,l,ithermal,mode,mi(2),norien,
     &  noddiam,null,icounter,inotr(2,*),ntrans,ipneigh(*),neigh(2,*),
     &  ielorien(*),iinc,istep,nkcoords,ngraph,k,nodes,nope,ne,
     &  nout,nset,istartset(*),iendset(*),ialset(*),iset,m,
     &  noutloc,ncomp,nksegment,iselect,noutplus,noutmin,ncomma
!
      real*8 co(3,*),v(0:mi(2),*),stn(6,*),een(6,*),t1(*),fn(0:mi(2),*),
     &  time,epn(*),enern(*),xstaten(nstate_,*),pi,qfn(3,*),oner,
     &  trab(7,*),stx(6,mi(1),*),orab(7,*),vr(0:mi(2),*),
     &  vi(0:mi(2),*),stnr(6,*),stni(6,*),vmax(0:3,*),stnmax(0:6,*),
     &  veold(0:mi(2),*),ener(mi(1),*),cs(17,*),eenmax(0:6,*)
!
      data icounter /0/
      save icounter,nkcoords,nout,noutmin,noutplus
!
      pi=4.d0*datan(1.d0)
!
      c='C'
!
      m1=' -1'
      m2=' -2'
      m3=' -3'
      m4=' -4'
      m5=' -5'
!
      p0='    0'
      p1='    1'
      p2='    2'
      p3='    3'
      p4='    4'
      p5='    5'
      p6='    6'
      p8='    8'
      p10='   10'
      p11='   11'
      p12='   12'
!
      if((time.le.0.d0).or.(nmethod.eq.2)) then
         fmat(1:8)='(e12.5) '
      elseif((dlog10(time).ge.0.d0).and.(dlog10(time).lt.10.d0)) then
         fmat(1:5)='(f12.'
         ncomma=10-int(dlog10(time)+1.d0)
         write(fmat(6:6),'(i1)') ncomma
         fmat(7:8)=') '
      else
         fmat(1:8)='(e12.5) '
      endif
!
      null=0
      one=1
      oner=1.d0
!
      if(kode.eq.1) then
!
        write(7,'(a5,a1)') p1,c
        call date_and_time(date,clock)
        newdate(1:20)='                    '
        newdate(1:2)=date(7:8)
        newdate(3:3)='.'
        if(date(5:6).eq.'01') then
           newdate(4:11)='january.'
           newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'02') then
           newdate(4:12)='february.'
           newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'03') then
           newdate(4:9)='march.'
           newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'04') then
           newdate(4:9)='april.'
           newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'05') then
           newdate(4:7)='may.'
           newdate(8:11)=date(1:4)
        elseif(date(5:6).eq.'06') then
           newdate(4:8)='june.'
           newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'07') then
           newdate(4:8)='july.'
           newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'08') then
           newdate(4:10)='august.'
           newdate(11:14)=date(1:4)
        elseif(date(5:6).eq.'09') then
           newdate(4:13)='september.'
           newdate(14:17)=date(1:4)
        elseif(date(5:6).eq.'10') then
           newdate(4:11)='october.'
           newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'11') then
           newdate(4:12)='november.'
           newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'12') then
           newdate(4:12)='december.'
           newdate(13:16)=date(1:4)
        endif
        newclock(1:2)=clock(1:2)
        newclock(3:3)=':'
        newclock(4:5)=clock(3:4)
        newclock(6:6)=':'
        newclock(7:8)=clock(5:6)
        write(7,'(a5,''UUSER'')') p1
        write(7,'(a5,''UDATE'',14x,a20)') p1,newdate
        write(7,'(a5,''UTIME'',14x,a8)') p1,newclock
        write(7,'(a5,''UHOST'')') p1
        write(7,'(a5,''UPGM               CalculiX'')') p1
        write(7,'(a5,''UDIR'')') p1
        write(7,'(a5,''UDBN'')') p1
!
!       storing the coordinates of the nodes
!
        write(7,'(a5,a1,67x,i1)') p2,c,one
!
        if(nmethod.ne.0) then
          nout=0
          noutplus=0
          noutmin=0
          do i=1,nk
             if(inum(i).eq.0) cycle
             write(7,101) m1,i,(co(j,i),j=1,3)
             nout=nout+1
             if(inum(i).gt.0) noutplus=noutplus+1
             if(inum(i).lt.0) noutmin=noutmin+1
          enddo
        else
          do i=1,nk
             write(7,101) m1,i,(co(j,i),j=1,3)
          enddo
          nout=nk
        endif
!
!       nkcoords is the number of nodes at the time when
!       the nodal coordinates are stored in the frd file.
!
        nkcoords=nk
!
        write(7,'(a3)') m3
!
!       storing the element topology
!
        write(7,'(a5,a1,67x,i1)') p3,c,one
!
        do i=1,ne0
!
           if(ipkon(i).lt.0) cycle
           indexe=ipkon(i)
           if(lakon(i)(4:4).eq.'2') then
              if((lakon(i)(7:7).eq.' ').or.(filab(1)(5:5).eq.'E').or.
     &           (lakon(i)(7:7).eq.'I')) then
              write(7,'(a3,i10,3a5)') m1,i,p4,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=1,10)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=11,12),
     &             (kon(indexe+j),j=17,19),kon(indexe+20),
     &             (kon(indexe+j),j=13,16)
              elseif(lakon(i)(7:7).eq.'B') then
              write(7,'(a3,i10,3a5)')m1,i,p12,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,3i10)') m2,kon(indexe+21),kon(indexe+23),
     &               kon(indexe+22)
              else
              write(7,'(a3,i10,3a5)')m1,i,p10,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,8i10)') m2,(kon(indexe+20+j),j=1,8)
              endif
           elseif(lakon(i)(4:4).eq.'8') then
              write(7,'(a3,i10,3a5)') m1,i,p1,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,8i10)') m2,(kon(indexe+j),j=1,8)
           elseif(lakon(i)(4:5).eq.'10') then
              write(7,'(a3,i10,3a5)') m1,i,p6,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=1,10)
           elseif(lakon(i)(4:4).eq.'4') then
              write(7,'(a3,i10,3a5)') m1,i,p3,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,4i10)') m2,(kon(indexe+j),j=1,4)
           elseif(lakon(i)(4:5).eq.'15') then
              if((lakon(i)(7:7).eq.' ').or.(filab(1)(5:5).eq.'E')) then
              write(7,'(a3,i10,3a5)') m1,i,p5,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=1,9),
     &          kon(indexe+13)
              write(7,'(a3,5i10)') m2,(kon(indexe+j),j=14,15),
     &          (kon(indexe+j),j=10,12)
              else
              write(7,'(a3,i10,3a5)') m1,i,p8,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,6i10)') m2,(kon(indexe+15+j),j=1,6)
              endif
           elseif(lakon(i)(4:4).eq.'6') then
              write(7,'(a3,i10,3a5)') m1,i,p2,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,6i10)') m2,(kon(indexe+j),j=1,6)
           elseif(lakon(i)(1:1).eq.'D') then
              if((kon(indexe+1).eq.0).or.(kon(indexe+3).eq.0)) cycle
              write(7,'(a3,i10,3a5)')m1,i,p12,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,3i10)') m2,kon(indexe+1),kon(indexe+3),
     &                 kon(indexe+2)
           elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).eq.'A'))then
              write(7,'(a3,i10,3a5)')m1,i,p11,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,2i10)') m2,(kon(indexe+j),j=1,2)
           endif
!
        enddo
!
        write(7,'(a3)') m3
!
        if(nmethod.eq.0) return
      endif
!
!     for cyclic symmetry frequency calculations only results
!     for even numbers (= odd modes, numbering starts at 0)are stored
!
      if((nmethod.eq.2).and.(((mode/2)*2.ne.mode).and.
     &                        (noddiam.ge.0))) return
!
!     storing the displacements of the nodes
!
      if(filab(1)(1:4).eq.'U   ') then
!
         iselect=1
         call frdset(filab(1),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  DISP        4    1'
         write(7,'(a132)') text
         text=' -5  D1          1    2    1    0'
         write(7,'(a132)') text
         text=' -5  D2          1    2    2    0'
         write(7,'(a132)') text
         text=' -5  D3          1    2    3    0'
         write(7,'(a132)') text
         text=' -5  ALL         1    2    0    0    1ALL'
         write(7,'(a132)') text
!
         call frdvector(v,iset,ntrans,filab(1),nkcoords,inum,m1,inotr,
     &        trab,co,istartset,iendset,ialset,mi,ngraph)
!
         write(7,'(a3)') m3
      endif
!
!     storing the imaginary part of displacements of the nodes
!     for the odd modes of cyclic symmetry calculations
!
      if(noddiam.ge.0) then
         if(filab(1)(1:4).eq.'U   ') then
!
            call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
            text=' -4  DISP        4    1'
            write(7,'(a132)') text
            text=' -5  D1          1    2    1    0'
            write(7,'(a132)') text
            text=' -5  D2          1    2    2    0'
            write(7,'(a132)') text
            text=' -5  D3          1    2    3    0'
            write(7,'(a132)') text
            text=' -5  ALL         1    2    0    0    1ALL'
            write(7,'(a132)') text
!     
c         call frdvector(v((mi(2)+1)*nk,1),iset,ntrans,filab,nkcoords,
c     &        inum,m1,inotr,trab,co,istartset,iendset,ialset,mi,ngraph)
         call frdvector(v(0,nk+1),iset,ntrans,filab(1),nkcoords,
     &        inum,m1,inotr,trab,co,istartset,iendset,ialset,mi,ngraph)
!     
            write(7,'(a3)') m3
         endif
      endif
!
!     storing the velocities of the nodes
!
      if(filab(21)(1:4).eq.'V   ') then
!
         iselect=1
         call frdset(filab(21),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  VELO        4    1'
         write(7,'(a132)') text
         text=' -5  V1          1    2    1    0'
         write(7,'(a132)') text
         text=' -5  V2          1    2    2    0'
         write(7,'(a132)') text
         text=' -5  V3          1    2    3    0'
         write(7,'(a132)') text
         text=' -5  ALL         1    2    0    0    1ALL'
         write(7,'(a132)') text
!
         call frdvector(veold,iset,ntrans,filab(21),nkcoords,inum,m1,
     &        inotr,trab,co,istartset,iendset,ialset,mi,ngraph)
!     
         write(7,'(a3)') m3
      endif
!
!     storing the temperatures in the nodes
!
      if(filab(2)(1:4).eq.'NT  ') then
!
         iselect=0
         call frdset(filab(2),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  NDTEMP      1    1'
         write(7,'(a132)') text
         text=' -5  T           1    1    0    0'
         write(7,'(a132)') text
!
         if(ithermal.le.1) then
            call frdscalar(t1,iset,nkcoords,inum,m1,
     &           istartset,iendset,ialset,ngraph,iselect)
         else
            ncomp=0
            call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &           istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
         endif
!     
         write(7,'(a3)') m3
      endif
!
!     storing the stresses in the nodes
!
      if(filab(3)(1:4).eq.'S   ') then
!     
         iselect=1
         call frdset(filab(3),set,iset,istartset,iendset,ialset,
     &        inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!     
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &        noutloc,description,kode,nmethod,fmat)
!     
         text=' -4  STRESS      6    1'
         write(7,'(a132)') text
         text=' -5  SXX         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  SYY         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  SZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  SXY         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  SYZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  SZX         1    4    3    1'
         write(7,'(a132)') text
!     
         call frdtensor(stn,iset,nkcoords,inum,m1,istartset,iendset,
     &          ialset,ngraph)
!     
         write(7,'(a3)') m3
      endif
!
!     storing the imaginary part of the stresses in the nodes
!     for the odd modes of cyclic symmetry calculations
!
      if(noddiam.ge.0) then
         if(filab(3)(1:4).eq.'S   ') then
!
            call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
            text=' -4  STRESS      6    1'
            write(7,'(a132)') text
            text=' -5  SXX         1    4    1    1'
            write(7,'(a132)') text
            text=' -5  SYY         1    4    2    2'
            write(7,'(a132)') text
            text=' -5  SZZ         1    4    3    3'
            write(7,'(a132)') text
            text=' -5  SXY         1    4    1    2'
            write(7,'(a132)') text
            text=' -5  SYZ         1    4    2    3'
            write(7,'(a132)') text
            text=' -5  SZX         1    4    3    1'
            write(7,'(a132)') text
!     
            call frdtensor(stn(1,nk+1),iset,nkcoords,inum,m1,istartset,
     &           iendset,ialset,ngraph)
!     
            write(7,'(a3)') m3
         endif
      endif
!     
!     storing the strains in the nodes
!
      if(filab(4)(1:4).eq.'E   ') then
!
         iselect=1
         call frdset(filab(4),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  STRAIN      6    1'
         write(7,'(a132)') text
         text=' -5  EXX         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  EYY         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  EZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  EXY         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  EYZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  EZX         1    4    3    1'
         write(7,'(a132)') text
!
         call frdtensor(een,iset,nkcoords,inum,m1,istartset,iendset,
     &          ialset,ngraph)
!
         write(7,'(a3)') m3
      endif
!
!     storing the forces in the nodes
!
      if(filab(5)(1:4).eq.'RF  ') then
!
         iselect=1
         call frdset(filab(5),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  FORC        4    1'
         write(7,'(a132)') text
         text=' -5  F1          1    2    1    0'
         write(7,'(a132)') text
         text=' -5  F2          1    2    2    0'
         write(7,'(a132)') text
         text=' -5  F3          1    2    3    0'
         write(7,'(a132)') text
         text=' -5  ALL         1    2    0    0    1ALL'
         write(7,'(a132)') text
!
         call frdvector(fn,iset,ntrans,filab(5),nkcoords,inum,m1,inotr,
     &        trab,co,istartset,iendset,ialset,mi,ngraph)
!     
         write(7,'(a3)') m3
      endif
!
!     storing the equivalent plastic strains in the nodes
!
      if(filab(6)(1:4).eq.'PEEQ') then
!
         iselect=1
         call frdset(filab(6),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  PE          1    1'
         write(7,'(a132)') text
         text=' -5  PE          1    1    0    0'
         write(7,'(a132)') text
!
         call frdscalar(epn,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the energy in the nodes
!
      if(filab(7)(1:4).eq.'ENER') then
!
         iselect=1
         call frdset(filab(7),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  ENER        1    1'
         write(7,'(a132)') text
         text=' -5  ENER        1    1    0    0'
         write(7,'(a132)') text
!
         call frdscalar(enern,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the contact informations at the nodes
!     with CDIS,CSTR
! 
      if(filab(26)(1:4).eq.'CONT') then
!
         do i=ne,1,-1
            if((lakon(i)(2:2).ne.'S').or.
     &           (lakon(i)(7:7).ne.'C')) exit
         enddo
         noutloc=ne-i
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  CONTACT     6    1'
         write(7,'(a132)') text
         text=' -5  COPEN       1    4    1    1'
         write(7,'(a132)') text
         text=' -5  CSLIP1      1    4    2    2'
         write(7,'(a132)') text
         text=' -5  CSLIP2      1    4    3    3'
         write(7,'(a132)') text
         text=' -5  CPRESS      1    4    1    2'
         write(7,'(a132)') text
         text=' -5  CSHEAR1     1    4    2    3'
         write(7,'(a132)') text
         text=' -5  CSHEAR2     1    4    3    1'
         write(7,'(a132)') text
!
         do i=ne,1,-1
            if((lakon(i)(2:2).ne.'S').or.
     &           (lakon(i)(7:7).ne.'C')) exit
            read(lakon(i)(8:8),'(i1)') nope
            nodes=kon(ipkon(i)+nope)
            write(7,101) m1,nodes,(stx(j,1,i),j=1,6)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the contact energy in the nodes
!
      if(filab(27)(1:4).eq.'CELS') then
!
         do i=ne,1,-1
            if((lakon(i)(2:2).ne.'S').or.
     &           (lakon(i)(7:7).ne.'C')) exit
         enddo
         noutloc=ne-i
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  CELS        1    1'
         write(7,'(a132)') text
         text=' -5  CELS        1    1    0    0'
         write(7,'(a132)') text
!
         do i=ne,1,-1
            if((lakon(i)(2:2).ne.'S').or.
     &           (lakon(i)(7:7).ne.'C')) exit
            read(lakon(i)(8:8),'(i1)') nope
            nodes=kon(ipkon(i)+nope)
            write(7,101) m1,nodes,ener(1,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the internal state variables in the nodes
!
      if(filab(8)(1:4).eq.'SDV ') then
!
         iselect=1
         call frdset(filab(8),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  SDV         6    1'
         if(nstate_.le.9) then
            write(text(18:18),'(i1)') nstate_
         else
            write(text(17:18),'(i2)') nstate_
         endif
         write(7,'(a132)') text
         do j=1,nstate_
            text=' -5  SDV         1    1    0    0'
            if(j.le.9) then
               write(text(9:9),'(i1)') j
            else
               write(text(9:10),'(i2)') j
            endif
            write(7,'(a132)') text
         enddo
!     
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               do k=1,int((nstate_+5)/6)
                  if(k.eq.1) then
                     write(7,101) m1,i,(xstaten(j,i),j=1,min(6,nstate_))
                  else
                     write(7,102) m2,(xstaten(j,i),j=(k-1)*6+1,
     &                    min(k*6,nstate_))
                  endif
               enddo
            enddo
         else
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  i=ialset(k)
                  if(inum(i).le.0) cycle
                  do l=1,int((nstate_+5)/6)
                     if(l.eq.1) then
                        write(7,101) m1,i,
     &                    (xstaten(j,i),j=1,min(6,nstate_))
                     else
                        write(7,102) m2,(xstaten(j,i),j=(l-1)*6+1,
     &                       min(l*6,nstate_))
                     endif
                  enddo
               else
                  i=ialset(k-2)
                  do
                     i=i-ialset(k)
                     if(i.ge.ialset(k-1)) exit
                     if(inum(i).le.0) cycle
                     do l=1,int((nstate_+5)/6)
                        if(l.eq.1) then
                           write(7,101) m1,i,
     &                       (xstaten(j,i),j=1,min(6,nstate_))
                        else
                           write(7,102) m2,(xstaten(j,i),j=(l-1)*6+1,
     &                          min(l*6,nstate_))
                        endif
                     enddo
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
!     storing the heat flux in the nodes
!
      if((filab(9)(1:4).eq.'HFL ').and.(ithermal.gt.1)) then
!
         iselect=1
         call frdset(filab(9),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
c         text=
c     & '  100CL       .00000E+00                                 3    1'
c         write(text(25:36),'(i12)') nout
c         text(37:48)=description
c         text(75:75)='1'
c         write(text(8:12),'(i5)') 100+kode
c         write(text(13:24),fmat) time
c         write(text(59:63),'(i5)') kode
c         write(7,'(a132)') text
         text=' -4  FLUX        4    1'
         write(7,'(a132)') text
         text=' -5  F1          1    2    1    0'
         write(7,'(a132)') text
         text=' -5  F2          1    2    2    0'
         write(7,'(a132)') text
         text=' -5  F3          1    2    3    0'
         write(7,'(a132)') text
         text=' -5  ALL         1    2    0    0    1ALL'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               write(7,101) m1,i,(qfn(j,i),j=1,3)
            enddo
         else
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  i=ialset(k)
                  if(inum(i).le.0) cycle
                  write(7,101) m1,i,(qfn(j,i),j=1,3)
               else
                  i=ialset(k-2)
                  do
                     i=i-ialset(k)
                     if(i.ge.ialset(k-1)) exit
                     if(inum(i).le.0) cycle
                     write(7,101) m1,i,(qfn(j,i),j=1,3)
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
!     storing the heat generation in the nodes
!
      if((filab(10)(1:4).eq.'RFL ').and.(ithermal.gt.1)) then
!
         iselect=1
         call frdset(filab(10),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
c         text=
c     & '  100CL       .00000E+00                                 3    1'
c         write(text(25:36),'(i12)') nout
c         text(37:48)=description
c         text(75:75)='1'
c         write(text(8:12),'(i5)') 100+kode
c         write(text(13:24),fmat) time
c         write(text(59:63),'(i5)') kode
c         write(7,'(a132)') text
         text=' -4  RFL         1    1'
         write(7,'(a132)') text
         text=' -5  RFL         1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=0
         call frdvectorcomp(fn,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the Zienkiewicz-Zhu improved stresses in the nodes
!
      if(filab(13)(1:3).eq.'ZZS') then
! 
         call estimator(co,nk,kon,ipkon,lakon,ne0,stn,
     &            ipneigh,neigh,stx,mi(1))
!
         iselect=1
         call frdset(filab(13),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  ZZSTR       6    1'
         write(7,'(a132)') text
         text=' -5  SXX         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  SYY         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  SZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  SXY         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  SYZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  SZX         1    4    3    1'
         write(7,'(a132)') text
!
         call frdtensor(stn,iset,nkcoords,inum,m1,istartset,iendset,
     &                 ialset,ngraph)
!
         write(7,'(a3)') m3
      endif
!
!     storing the imaginary part of the Zienkiewicz-Zhu stresses in the nodes
!     for the odd modes of cyclic symmetry calculations
!
      if(noddiam.ge.0) then
         if(filab(13)(1:3).eq.'ZZS') then
!     
            call estimator(co,nk,kon,ipkon,lakon,ne0,stn,
     &           ipneigh,neigh,stx(1,1,ne+1),mi(1))
!     
            iselect=1
            call frdset(filab(13),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!     
            call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &           noutloc,description,kode,nmethod,fmat)
!     
            text=' -4  ZZSTR       6    1'
            write(7,'(a132)') text
            text=' -5  SXX         1    4    1    1'
            write(7,'(a132)') text
            text=' -5  SYY         1    4    2    2'
            write(7,'(a132)') text
            text=' -5  SZZ         1    4    3    3'
            write(7,'(a132)') text
            text=' -5  SXY         1    4    1    2'
            write(7,'(a132)') text
            text=' -5  SYZ         1    4    2    3'
            write(7,'(a132)') text
            text=' -5  SZX         1    4    3    1'
            write(7,'(a132)') text
!     
            call frdtensor(stn,iset,nkcoords,inum,m1,istartset,iendset,
     &           ialset,ngraph)
!     
            write(7,'(a3)') m3
         endif
      endif
!
!     storing the total temperature in the fluid nodes
!
      if(filab(14)(1:4).eq.'TT  ') then
!
         iselect=-1
         call frdset(filab(14),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  TOTEMP      1    1'
         write(7,'(a132)') text
         text=' -5  TT          1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=0
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the mass flow in the fluid nodes
!
      if(filab(15)(1:4).eq.'MF  ') then
!
         iselect=-1
         call frdset(filab(15),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  MAFLOW      1    1'
         write(7,'(a132)') text
         text=' -5  MF          1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=1
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the total pressure in the gas network nodes
!
      if(filab(16)(1:4).eq.'PT  ') then
!
         iselect=-1
         call frdset(filab(16),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  TOPRES      1    1'
         write(7,'(a132)') text
         text=' -5  PT          1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=2
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the static pressure in the liquid network nodes
!
      if(filab(22)(1:4).eq.'PT  ') then
!
         iselect=-1
         call frdset(filab(22),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  STPRES      1    1'
         write(7,'(a132)') text
         text=' -5  PS          1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=2
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the liquid depth in the channel nodes
!
      if(filab(28)(1:4).eq.'DEPT') then
!
         iselect=-1
         call frdset(filab(28),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  DEPTH       1    1'
         write(7,'(a132)') text
         text=' -5  DEPTH       1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=2
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the liquid depth in the channel nodes
!
      if(filab(29)(1:4).eq.'HCRI') then
!
         iselect=-1
         call frdset(filab(29),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  HCRIT       1    1'
         write(7,'(a132)') text
         text=' -5  HCRIT       1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=3
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
!     storing the static temperature in the fluid nodes
!
      if(filab(17)(1:4).eq.'TS  ') then
!
         iselect=-1
         call frdset(filab(17),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  STTEMP      1    1'
         write(7,'(a132)') text
         text=' -5  TS          1    1    0    0'
         write(7,'(a132)') text
!
         ncomp=3
         call frdvectorcomp(v,iset,nkcoords,inum,m1,
     &        istartset,iendset,ialset,ncomp,mi,ngraph,iselect)
!
         write(7,'(a3)') m3
      endif
!
c      if((nmethod.ne.2).and.(nmethod.lt.4)) return
!
!     the remaining lines only apply to frequency calculations
!     with cyclic symmetry and steady state calculations
!
      if((nmethod.ne.2).and.(nmethod.ne.5)) return
      if((nmethod.eq.5).and.(mode.eq.-1)) return
!
!     storing the displacements of the nodes (magnitude, phase)
!
      if(filab(11)(1:4).eq.'PU  ') then
!
         iselect=1
         call frdset(filab(11),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  PDISP       6    1'
         write(7,'(a132)') text
         text=' -5  MAG1        1   12    1    0'
         write(7,'(a132)') text
         text=' -5  MAG2        1   12    2    0'
         write(7,'(a132)') text
         text=' -5  MAG3        1   12    3    0'
         write(7,'(a132)') text
         text=' -5  PHA1        1   12    4    0'
         write(7,'(a132)') text
         text=' -5  PHA2        1   12    5    0'
         write(7,'(a132)') text
         text=' -5  PHA3        1   12    6    0'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).eq.0) cycle
               write(7,101) m1,i,(vr(j,i),j=1,3),(vi(j,i),j=1,3)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).eq.0) cycle
                     write(7,101) m1,i,(vr(j,i),j=1,3),(vi(j,i),j=1,3)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).eq.0) cycle
                        write(7,101) m1,i,(vr(j,i),j=1,3),
     &                      (vi(j,i),j=1,3)
                     enddo
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
!     storing the temperatures of the nodes (magnitude,phase)
!
      if(filab(12)(1:4).eq.'PNT ') then
!
         iselect=1
         call frdset(filab(12),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  PNDTEMP     2    1'
         write(7,'(a132)') text
         text=' -5  MAG1        1    1    1    0'
         write(7,'(a132)') text
         text=' -5  PHA1        1    1    2    0'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).eq.0) cycle
               write(7,101) m1,i,vr(0,i),vi(0,i)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).eq.0) cycle
                     write(7,101) m1,i,vr(0,i),vi(0,i)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).eq.0) cycle
                        write(7,101) m1,i,vr(0,i),vi(0,i)
                     enddo
                  enddo
               endif
            enddo
         endif
!
         write(7,'(a3)') m3
      endif
!
!     storing the stresses in the nodes (magnitude,phase)
!
      if(filab(18)(1:4).eq.'PHS ') then
!
         iselect=1
         call frdset(filab(18),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  PSTRESS    12    1'
         write(7,'(a132)') text
         text=' -5  MAGXX       1    4    1    1'
         write(7,'(a132)') text
         text=' -5  MAGYY       1    4    2    2'
         write(7,'(a132)') text
         text=' -5  MAGZZ       1    4    3    3'
         write(7,'(a132)') text
         text=' -5  MAGXY       1    4    1    2'
         write(7,'(a132)') text
         text=' -5  MAGYZ       1    4    2    3'
         write(7,'(a132)') text
         text=' -5  MAGZX       1    4    3    1'
         write(7,'(a132)') text
         text=' -5  PHAXX       1    4    1    1'
         write(7,'(a132)') text
         text=' -5  PHAYY       1    4    2    2'
         write(7,'(a132)') text
         text=' -5  PHAZZ       1    4    3    3'
         write(7,'(a132)') text
         text=' -5  PHAXY       1    4    1    2'
         write(7,'(a132)') text
         text=' -5  PHAYZ       1    4    2    3'
         write(7,'(a132)') text
         text=' -5  PHAZX       1    4    3    1'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               write(7,101) m1,i,(stnr(j,i),j=1,4),
     &              stnr(6,i),stnr(5,i)
               write(7,101) m2,i,(stni(j,i),j=1,4),
     &              stni(6,i),stni(5,i)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).le.0) cycle
                     write(7,101) m1,i,(stnr(j,i),j=1,4),
     &                    stnr(6,i),stnr(5,i)
                     write(7,101) m2,i,(stni(j,i),j=1,4),
     &                    stni(6,i),stni(5,i)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).le.0) cycle
                        write(7,101) m1,i,(stnr(j,i),j=1,4),
     &                       stnr(6,i),stnr(5,i)
                        write(7,101) m2,i,(stni(j,i),j=1,4),
     &                       stni(6,i),stni(5,i)
                     enddo
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
      if(nmethod.ne.2) return
!
!     storing the maximum displacements of the nodes 
!     in the basis sector
!         (magnitude, components)
!
      if(filab(19)(1:4).eq.'MAXU') then
!
         iselect=1
         call frdset(filab(19),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &           ngraph)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  MDISP       4    1'
         write(7,'(a132)') text
         text=' -5  DX          1    4    1    0'
         write(7,'(a132)') text
         text=' -5  DY          1    4    2    0'
         write(7,'(a132)') text
         text=' -5  DZ          1    4    3    0'
         write(7,'(a132)') text
         text=' -5  ANG         1    4    4    0'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).eq.0) cycle
               write(7,101) m1,i,(vmax(j,i),j=1,3),vmax(0,i)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).eq.0) cycle
                     write(7,101) m1,i,(vmax(j,i),j=1,3),vmax(0,i)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).eq.0) cycle
                        write(7,101) m1,i,(vmax(j,i),j=1,3),vmax(0,i)
                     enddo
                  enddo
               endif
            enddo
         endif
!
         write(7,'(a3)') m3
      endif
!
!     storing the worst principal stress at the nodes
!     in the basis sector (components, magnitude)
! 
!     the worst principal stress is the maximum of the
!     absolute value of all principal stresses, times
!     its original sign
!
      if(filab(20)(1:4).eq.'MAXS') then
!
         iselect=1
         call frdset(filab(20),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  MSTRESS     7    1'
         write(7,'(a132)') text
         text=' -5  SXX         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  SYY         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  SZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  SXY         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  SYZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  SZX         1    4    3    1'
         write(7,'(a132)') text
         text=' -5  MAG         1    4    0    0'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               write(7,101) m1,i,(stnmax(j,i),j=1,4),
     &              stnmax(6,i),stnmax(5,i)
               write(7,101) m2,i,stnmax(0,i)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).le.0) cycle
                     write(7,101) m1,i,(stnmax(j,i),j=1,4),
     &                    stnmax(6,i),stnmax(5,i)
                     write(7,101) m2,i,stnmax(0,i)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).le.0) cycle
                        write(7,101) m1,i,(stnmax(j,i),j=1,4),
     &                       stnmax(6,i),stnmax(5,i)
                        write(7,101) m2,i,stnmax(0,i)
                     enddo
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
!     storing the worst principal strain at the nodes
!     in the basis sector (components, magnitude)
! 
!     the worst principal strain is the maximum of the
!     absolute value of all principal strains, times
!     its original sign
!
      if(filab(30)(1:4).eq.'MAXE') then
!
         iselect=1
         call frdset(filab(30),set,iset,istartset,iendset,ialset,
     &           inum,noutloc,nout,nset,noutmin,noutplus,iselect)
!
         call frdheader(icounter,oner,time,pi,noddiam,cs,null,mode,
     &                  noutloc,description,kode,nmethod,fmat)
!
         text=' -4  MSTRAIN     7    1'
         write(7,'(a132)') text
         text=' -5  EXX         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  EYY         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  EZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  EXY         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  EYZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  EZX         1    4    3    1'
         write(7,'(a132)') text
         text=' -5  MAG         1    4    0    0'
         write(7,'(a132)') text
!
         if(iset.eq.0) then
            do i=1,nkcoords
               if(inum(i).le.0) cycle
               write(7,101) m1,i,(eenmax(j,i),j=1,4),
     &              eenmax(6,i),eenmax(5,i)
               write(7,101) m2,i,eenmax(0,i)
            enddo
         else
            nksegment=nkcoords/ngraph
            do k=istartset(iset),iendset(iset)
               if(ialset(k).gt.0) then
                  do l=0,ngraph-1
                     i=ialset(k)+l*nksegment
                     if(inum(i).le.0) cycle
                     write(7,101) m1,i,(eenmax(j,i),j=1,4),
     &                    eenmax(6,i),eenmax(5,i)
                     write(7,101) m2,i,eenmax(0,i)
                  enddo
               else
                  l=ialset(k-2)
                  do
                     l=l-ialset(k)
                     if(l.ge.ialset(k-1)) exit
                     do m=0,ngraph-1
                        i=l+m*nksegment
                        if(inum(i).le.0) cycle
                        write(7,101) m1,i,(eenmax(j,i),j=1,4),
     &                       eenmax(6,i),eenmax(5,i)
                        write(7,101) m2,i,eenmax(0,i)
                     enddo
                  enddo
               endif
            enddo
         endif
!     
         write(7,'(a3)') m3
      endif
!
 101  format(a3,i10,1p,6e12.5)
 102  format(a3,10x,1p,6e12.5)
!
      return
      end


