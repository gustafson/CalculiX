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
      subroutine frdphase(kode,time,nk,inum,vr,vi,stnr,stni,filab,
     &   mode,noddiam,nmethod,vmax,stnmax,nkcoords)
!
!     stores the results in frd format
!
      implicit none
!
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*6 filab(*)
      character*5 p0,p1,p2,p3,p4,p5,p6
      character*8 fmat
      character*132 text
!
      integer inum(*),nk,kode,i,j,one,null,mode,noddiam,nmethod,
     &   nkcoords
!
      real*8 vr(0:4,*),time,oner,vi(0:4,*),stnr(6,*),stni(6,*),pi,
     &  vmax(0:3,*),stnmax(0:6,*)
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
!
      pi=4.d0*datan(1.d0)
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
!     storing the displacements of the nodes (magnitude, phase)
!
      if(filab(11)(1:4).eq.'PU  ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') kode
         write(7,'(a132)') text
!
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
         text(75:75)='1'
         write(text(8:12),'(i5)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         if(nmethod.eq.2) text(64:68)='MODAL'
         write(7,'(a132)') text
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
         do i=1,nkcoords
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,(vr(j,i),j=1,3),(vi(j,i),j=1,3)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the temperatures of the nodes (magnitude,phase)
!
      if(filab(12)(1:4).eq.'PNT ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') kode
         write(7,'(a132)') text
!
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
         text(75:75)='1'
         write(text(8:12),'(i5)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         if(nmethod.eq.2) text(64:68)='MODAL'
         write(7,'(a132)') text
!
         text=' -4  PNDTEMP     7    1'
         write(7,'(a132)') text
         text=' -5  MAG1        1    1    1    0'
         write(7,'(a132)') text
         text=' -5  PHA1        1    1    2    0'
         write(7,'(a132)') text
         text=' -5  ALL         1   12    0    0    1ALL'
         write(7,'(a132)') text
!
         do i=1,nkcoords
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,vr(0,i),vr(0,i),vr(0,i),
     &            vi(0,i),vi(0,i),vi(0,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
      if(nmethod.ne.2) return
!
!     storing the stresses in the nodes (magnitude,phase)
!
      if(filab(18)(1:4).eq.'PHS ') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') kode
         write(7,'(a132)') text
!
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
         text(75:75)='1'
         write(text(8:12),'(i5)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         if(nmethod.eq.2) text(64:68)='MODAL'
         write(7,'(a132)') text
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
         do i=1,nkcoords
            if(inum(i).le.0) cycle
            write(7,100) m1,i,(stnr(j,i),j=1,4),
     &           stnr(6,i),stnr(5,i)
            write(7,100) m2,i,(stni(j,i),j=1,4),
     &           stni(6,i),stni(5,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the maximum displacements of the nodes 
!     in the basis sector
!         (magnitude, components)
!
      if(filab(19)(1:4).eq.'MAXU') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') kode
         write(7,'(a132)') text
!
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
         text(75:75)='1'
         write(text(8:12),'(i5)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         if(nmethod.eq.2) text(64:68)='MODAL'
         write(7,'(a132)') text
!
         text=' -4  MDISP       4    1'
         write(7,'(a132)') text
         text=' -5  D1          1    4    2    0'
         write(7,'(a132)') text
         text=' -5  D2          1    4    3    0'
         write(7,'(a132)') text
         text=' -5  D3          1    4    4    0'
         write(7,'(a132)') text
         text=' -5  MAG         1    4    1    0'
         write(7,'(a132)') text
!
         do i=1,nkcoords
            if(inum(i).eq.0) cycle
            write(7,100) m1,i,(vmax(j,i),j=1,3),vmax(0,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
!     storing the worst principal stress at the nodes
!     in the basis sector (components, magnitude)
! 
!     the worst principal stress is the maximum of the
!     absolute value of all principal stresses
!
      if(filab(20)(1:4).eq.'MAXS') then
         text='    1PSTEP'
         write(text(25:36),'(i12)') kode
         write(7,'(a132)') text
!
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
         text(75:75)='1'
         write(text(8:12),'(i5)') 100+kode
         write(text(13:24),fmat) time
         write(text(59:63),'(i5)') kode
         if(nmethod.eq.2) text(64:68)='MODAL'
         write(7,'(a132)') text
!
         text=' -4  MSTRESS     7    1'
         write(7,'(a132)') text
         text=' -5  SRR         1    4    1    1'
         write(7,'(a132)') text
         text=' -5  STT         1    4    2    2'
         write(7,'(a132)') text
         text=' -5  SZZ         1    4    3    3'
         write(7,'(a132)') text
         text=' -5  SRT         1    4    1    2'
         write(7,'(a132)') text
         text=' -5  STZ         1    4    2    3'
         write(7,'(a132)') text
         text=' -5  SZR         1    4    3    1'
         write(7,'(a132)') text
         text=' -5  MAG         1    4    0    0'
         write(7,'(a132)') text
!
         do i=1,nkcoords
            if(inum(i).le.0) cycle
            write(7,100) m1,i,(stnmax(j,i),j=1,4),
     &           stnmax(6,i),stnmax(5,i)
            write(7,100) m2,i,stnmax(0,i)
         enddo
!
         write(7,'(a3)') m3
      endif
!
 100  format(a3,i10,1p,6e12.5)
!
      return
      end


