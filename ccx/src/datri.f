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
      subroutine datri(al,au,ad,jp,neq,flg)
      implicit none
c....triangular decomposition of a matrix stored in profile form
      logical flg
      integer jp(*),ior,iow,j,neq,jd,jr,jh,is,ie,i,jrh,idh,ifig,id,ih
      real*8 al(*),au(*),ad(*),zero,one,tol,dd,dimx,dimn,dfig,daval,dot
      common /iofile/ ior,iow
c....n.b. tol should be set to approximate half-word precision
      data zero,one/0.0d0,1.0d0/, tol/0.5d-07/
!
      dd=0.d0
!
c....set initial values for conditioning check
      dimx = zero
      dimn = zero
      do 50 j = 1,neq
      dimn = max(dimn,abs(ad(j)))
   50 continue
      dfig = zero
c....loop trough the columns to perform the triangular decomposition
      jd = 1
      do 200 j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1
c....if diagonal is zero compute a norm for singularity test
          if(ad(j).eq.zero) call datest(au(jr),jh,daval)
          do 100 i = is,ie
            jr =jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0)then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr) = au(jr) - dot(au(jrh),al(idh),ih)
              if(flg) al(jr) = al(jr) - dot(al(jrh),au(idh),ih)
            endif
  100     continue
        endif
c....reduce the diagonal
        if(jh.ge.0) then
          dd = ad(j)
          jr = jd - jh
          jrh = j - jh - 1
          call dredu(al(jr),au(jr),ad(jrh),jh+1,flg  ,ad(j))
c....check for possible errors and print warnings
          if(abs(ad(j)).lt.tol*abs(dd))  write(iow,2000) j
          if(dd.lt.zero.and.ad(j).gt.zero) write(iow,2001) j
          if(dd.gt.zero.and.ad(j).lt.zero) write(iow,2001) j
          if(ad(j).eq.zero)                write(iow,2002) j
          if(dd.eq.zero.and.jh.gt.0) then
            if(abs(ad(j)).lt.tol*daval)   write(iow,2003)j
          endif
          if(ior.lt.0) then
            if(abs(ad(j)).lt.tol*abs(dd))  write(*,2000)j
            if(dd.lt.zero.and.ad(j).gt.zero) write(*,2001)j
            if(dd.gt.zero.and.ad(j).lt.zero) write(*,2001)j
            if(ad(j).eq.zero)                write(*,2002)j
            if(dd.eq.zero.and.jh.gt.0) then
              if(abs(ad(j)).lt.tol*daval)   write(*,2003)j
            endif
          endif
        endif
c....store reciprocal of diagonal, compute condition checks
        if(ad(j).ne.zero) then
          dimx  = max(dimx,abs(ad(j)))   
          dimn  = min(dimn,abs(ad(j)))   
          dfig  = max(dfig,abs(dd/ad(j))) 
          ad(j) = one/ad(j)
        endif
  200 continue  
c....print conditioning information
      dd=zero
      if(dimn.ne.zero) dd = dimx/dimn
      ifig = dlog10(dfig) + 0.6
!      write(iow,2004) dimx,dimn,dd,ifig
!      if(ior.lt.0) write(*,2004) dimx,dimn,dd,ifig
      return
c....formats
 2000 format(' ***DATRI WARNING 1*** Loss of at least 7 digits in',
     1 ' reducing diagonal of equation',i5) 
 2001 format(' ***DATRI WARNING 2*** Sign of diagonal changed when',
     1 ' reducing equation',i5) 
 2002 format(' ***DATRI WARNING 3*** Reduced diagonal is zero for',
     1 ' equation',i5) 
 2003 format(' ***DATRI WARNING 4*** Rank failure for zero unreduced',
     1 ' diagonal in equation',i5) 
 2004 format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     1 '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)
      end
