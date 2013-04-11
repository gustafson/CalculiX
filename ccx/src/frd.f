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
      subroutine frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
     &  kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,xstaten,
     &  nstate_,istep,iinc,ithermal,qfn,mode,noddiam,trab,inotr,
     &  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
     &  mint_,sti)
!
!     stores the results in frd format
!
      implicit none
!
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*5 p0,p1,p2,p3,p4,p5,p6,p8,p10,p11,p12
      character*6 filab(*)
      character*8 lakon(*),date,newclock,fmat
      character*10 clock
      character*12 description
      character*20 newdate
      character*80 matname(*)
      character*132 text
!
      integer kon(*),inum(*),nk,ne,nmethod,kode,i,j,ipkon(*),indexe,
     &  one,ielmat(*),lb,nterms,nstate_,l,ithermal,mode,mint_,norien,
     &  noddiam,null,icounter,inotr(2,*),ntrans,ipneigh(*),neigh(2,*),
     &  ielorien(*),iinc,istep
!
      real*8 co(3,*),v(0:3,*),stn(6,*),een(6,*),t1(*),fn(0:3,*),time,
     &  epn(*),enern(*),xstaten(nstate_,*),pi,qfn(3,*),oner,trab(7,*),
     &  a(3,3),sti(6,mint_,*),orab(7,*)
!
      data icounter /0/
      save icounter
!
      icounter=icounter+1
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
      if(time.le.0.d0) then
         fmat(1:8)='(e12.5) '
      elseif((dlog10(time).ge.0.d0).and.(dlog10(time).lt.11.d0)) then
         fmat(1:5)='(f12.'
         write(fmat(6:7),'(i2)') 11-int(dlog10(time)+1.d0)
         fmat(8:8)=')'
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
          do i=1,nk
             if(inum(i).eq.0) cycle
             write(7,100) m1,i,(co(j,i),j=1,3)
          enddo
        else
          do i=1,nk
             write(7,100) m1,i,(co(j,i),j=1,3)
          enddo
        endif
!
        write(7,'(a3)') m3
!
!       storing the element topology
!
        write(7,'(a5,a1,67x,i1)') p3,c,one
!
        do i=1,ne
!
           if(ipkon(i).lt.0) cycle
           indexe=ipkon(i)
           if(lakon(i)(4:4).eq.'2') then
              if((lakon(i)(7:7).eq.' ').or.(filab(1)(5:5).eq.'E')) then
              write(7,'(a3,i10,3a5)') m1,i,p4,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=1,10)
              write(7,'(a3,10i10)') m2,(kon(indexe+j),j=11,12),
     &             (kon(indexe+j),j=17,19),kon(indexe+20),
     &             (kon(indexe+j),j=13,16)
              elseif(lakon(i)(7:7).eq.'B') then
              write(7,'(a3,i10,3a5)')m1,i,p12,p0,matname(ielmat(i))(1:5)
              write(7,'(a3,3i10)') m2,(kon(indexe+20+j),j=1,3)
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
              write(7,'(a3,3i10)') m2,(kon(indexe+j),j=1,3)
           elseif(lakon(i)(1:1).eq.'E') then
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
!     storing the displacements of the nodes
!
      if(filab(1)(1:4).eq.'U   ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         if((ntrans.eq.0).or.(filab(1)(6:6).eq.'G')) then
            do i=1,nk
               if(inum(i).eq.0) cycle
               write(7,100) m1,i,(v(j,i),j=1,3)
            enddo
         else
            do i=1,nk
               if(inum(i).eq.0) cycle
               if(inotr(1,i).eq.0) then
                  write(7,100) m1,i,(v(j,i),j=1,3)
               else
                  call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                  write(7,100) m1,i,
     &               v(1,i)*a(1,1)+v(2,i)*a(2,1)+v(3,i)*a(3,1),
     &               v(1,i)*a(1,2)+v(2,i)*a(2,2)+v(3,i)*a(3,2),
     &               v(1,i)*a(1,3)+v(2,i)*a(2,3)+v(3,i)*a(3,3)
               endif
            enddo
         endif
!
         write(7,'(a3)') m3
      endif
!
!     storing the temperatures in the nodes
!
      if(filab(2)(1:4).eq.'NT  ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  NDTEMP      1    1'
         write(7,'(a132)') text
         text=' -5  T           1    1    0    0'
         write(7,'(a132)') text
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            if(ithermal.le.1) then
               write(7,100) m1,i,t1(i)
            else
               write(7,100) m1,i,v(0,i)
            endif
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the stresses in the nodes
!
      if(filab(3)(1:4).eq.'S   ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode      
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,101) m1,i,(stn(j,i),j=1,4),
     &           stn(6,i),stn(5,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the strains in the nodes
!
      if(filab(4)(1:4).eq.'E   ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode      
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,101) m1,i,(een(j,i),j=1,4),
     &           een(6,i),een(5,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the forces in the nodes
!
      if(filab(5)(1:4).eq.'RF  ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         if((ntrans.eq.0).or.(filab(5)(6:6).eq.'G')) then
            do i=1,nk
               if(inum(i).eq.0) cycle
               write(7,100) m1,i,(fn(j,i),j=1,3)
            enddo
         else
            do i=1,nk
               if(inum(i).eq.0) cycle
               if(inotr(1,i).eq.0) then
                  write(7,100) m1,i,(fn(j,i),j=1,3)
               else
                  call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                  write(7,100) m1,i,
     &                 fn(1,i)*a(1,1)+fn(2,i)*a(2,1)+fn(3,i)*a(3,1),
     &                 fn(1,i)*a(1,2)+fn(2,i)*a(2,2)+fn(3,i)*a(3,2),
     &                 fn(1,i)*a(1,3)+fn(2,i)*a(2,3)+fn(3,i)*a(3,3)
               endif
            enddo
         endif
!
         write(7,'(a3)') m3
      endif
!
!     storing the equivalent plastic strains in the nodes
!
      if(filab(6)(1:4).eq.'PE  ') then
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  PE          1    1'
         write(7,'(a132)') text
         text=' -5  PE          1    1    0    0'
         write(7,'(a132)') text
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,epn(i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the energy in the nodes
!
      if(filab(7)(1:4).eq.'ENER') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  ENER        1    1'
         write(7,'(a132)') text
         text=' -5  ENER        1    1    0    0'
         write(7,'(a132)') text
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,enern(i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the internal state variables in the nodes
!
      if(filab(8)(1:4).eq.'SDV ') then
         do l=1,(nstate_+5)/6
            lb=(l-1)*6
            text=
     & '  100CL       .00000E+00                                 3    1'
            text(37:48)=description
            text(75:75)='1'
            write(text(8:10),'(i3)') 100+kode      
            write(text(13:24),fmat) time
            write(text(59:63),'(i5)') kode
            write(7,'(a132)') text
            if(l.eq.(nstate_+5)/6) then
               nterms=nstate_-lb
            else
               nterms=6
            endif
            text=' -4  SDV         6    1'
            write(text(18:18),'(i1)') nterms
            if(lb+1.le.9) then
               write(text(9:9),'(i1)') lb+1
            else
               write(text(9:10),'(i2)') lb+1
            endif
            write(7,'(a132)') text
            do j=1,nterms
               text=' -5  SDV         1    1    0    0'
               if(lb+j.le.9) then
                  write(text(9:9),'(i1)') lb+j
               else
                  write(text(9:10),'(i2)') lb+j
               endif
               write(7,'(a132)') text
            enddo
!
            if(l.eq.(nstate_+5)/6) then
               do i=1,nk
                  if(inum(i).eq.0) cycle
                  write(7,101) m1,i,(xstaten(lb+j,i),j=1,nstate_-lb)
               enddo
            else
               do i=1,nk
                  if(inum(i).eq.0) cycle
                  write(7,101) m1,i,(xstaten(lb+j,i),j=1,6)
               enddo
            endif
!
            write(7,'(a3)') m3
         enddo
      endif
!
!     storing the heat flux in the nodes
!
      if((filab(9)(1:4).eq.'HFL ').and.(ithermal.gt.1)) then
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,(qfn(j,i),j=1,3)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the heat generation in the nodes
!
      if((filab(10)(1:4).eq.'RFL ').and.(ithermal.gt.1)) then
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  RFL         1    1'
         write(7,'(a132)') text
         text=' -5  RFL         1    1    0    0'
         write(7,'(a132)') text
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,fn(0,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the residual forces in the nodes (only in case
!     of divergence)
!
      if(filab(13)(1:5).eq.'RFRES') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  RFRES       4    1'
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
         if((ntrans.eq.0).or.(filab(5)(6:6).eq.'G')) then
            do i=1,nk
               if(inum(i).eq.0) cycle
               write(7,100) m1,i,(fn(j,i),j=1,3)
            enddo
         else
            do i=1,nk
               if(inum(i).eq.0) cycle
               if(inotr(1,i).eq.0) then
                  write(7,100) m1,i,(fn(j,i),j=1,3)
               else
                  call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                  write(7,100) m1,i,
     &                 fn(1,i)*a(1,1)+fn(2,i)*a(2,1)+fn(3,i)*a(3,1),
     &                 fn(1,i)*a(1,2)+fn(2,i)*a(2,2)+fn(3,i)*a(3,2),
     &                 fn(1,i)*a(1,3)+fn(2,i)*a(2,3)+fn(3,i)*a(3,3)
               endif
            enddo
         endif
!
         write(7,'(a3)') m3
      endif
!
!     storing the residual heat generation in the nodes (only in the
!     case of divergence)
!
      if((filab(14)(1:6).eq.'RFLRES').and.(ithermal.gt.1)) then
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
         text=' -4  RFLRES      1    1'
         write(7,'(a132)') text
         text=' -5  RFLRES      1    1    0    0'
         write(7,'(a132)') text
!
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,fn(0,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the stress errors in the nodes
!
      if(filab(15)(1:3).eq.'ZZS') then
! 
         call estimator(co,nk,kon,ipkon,lakon,ne,stn,
     &            ipneigh,neigh,sti,mint_)
!
         text='    1PSTEP'
         write(text(25:36),'(i12)') icounter
         write(7,'(a132)') text
         if(nmethod.eq.2) then
            text='    1PGM'
            write(text(25:36),'(e12.6)') oner
            write(7,'(a132)') text
            text='    1PGK'
            write(text(25:36),'(e12.6)') (time*2.d0*pi)**2
            write(7,'(a132)') text
            text='    1PHID'
            write(text(25:36),'(i12)') noddiam
            write(7,'(a132)') text
            text='    1PSUBC'
            write(text(25:36),'(i12)') null
            write(7,'(a132)') text
            text='    1PMODE'
            write(text(25:36),'(i12)') mode+1
            write(7,'(a132)') text
         endif
!
         text=
     & '  100CL       .00000E+00                                 3    1'
         text(37:48)=description
         if(nmethod.eq.2) text(64:68)='MODAL'
         text(75:75)='1'
         write(text(8:10),'(i3)') 100+kode      
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         write(7,'(a132)') text
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
         do i=1,nk
            if(inum(i).eq.0) cycle
            write(7,101) m1,i,(stn(j,i),j=1,4),
     &           stn(6,i),stn(5,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
 100  format(a3,i10,1p,3e12.5)
 101  format(a3,i10,1p,6e12.5)
!
      return
      end


