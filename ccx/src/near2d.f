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
      subroutine near2d(xo,yo,x,y,nx,ny,xp,yp,n,node,k)
C
c     determines the k closest nodes out of n with coordinates in
c     xo,yo to the point with coordinates (xp,yp);
c
c     order xo and yo before the first call to near2d and
c     store the results with the corresponding permutation array
c     in x,y and nx,ny, respectively
c
      implicit none
C
      integer node(*),nx(*),ny(*),ni(24),n,m,idx,idy,it,j,iz,kflag,i,
     &  k,ks,nboundary
c
      real*8 x(*),y(*),xo(*),yo(*),r(24),xr,yr,aan,aas,aaw,aao,
     &  aano,aaso,aasw,aanw,aamin,xx,yy,rr,xp,yp
c
      if(k.gt.n) then
c         write(*,*) '*ERROR in near2d: neighbors requested'
c         write(*,*) '       supersedes number of nodes'
c         write(*,*) k,n
c         stop
         k=n
      endif
c
      if(k.gt.20) then
         write(*,*) '*ERROR in near2d: no more than 20 neighbors'
         write(*,*) '       can be identified'
         stop
      endif
c
      CALL IDENT(X,XP,N,IDX)
      CALL IDENT(Y,YP,N,IDY)
C
      kflag=2
C
      do i=1,k
         XR=XO(i)-XP
         YR=YO(i)-YP
         R(i)=dSQRT(XR*XR+YR*YR)
         NI(i)=i
      enddo
      call dsort(r,ni,k,kflag)
C
      i=1
      ks=0
c
      loop: do
c
         IZ=0
C
C     WEST
C
         M=IDX-i+1
         IF(M.LE.0) THEN
            AAW=R(k)
            GO TO 4
         END IF
C
         XX=XO(NX(M))
         YY=YO(NX(M))
         AAW=XX-XP
         YR=YY-YP
         RR=SQRT(AAW*AAW+YR*YR)
C
         IF(RR.GE.R(k)) GO TO 4
         IT=IZ+k
         DO 8 J=1,IT
            IF(NX(M).EQ.NI(J)) GO TO 4
 8       CONTINUE
         IZ=IZ+1
         R(IZ+k)=RR
         NI(IZ+k)=NX(M)
C
C     EAST
C
 4       M=IDX+i
         IF(M.GT.N) THEN
            AAO=R(k)
            GO TO 6
         END IF
C
         XX=XO(NX(M))
         YY=YO(NX(M))
         AAO=XX-XP
         YR=YY-YP
         RR=SQRT(AAO*AAO+YR*YR)
C
         IF(RR.GE.R(k)) GO TO 6
         IT=IZ+k
         DO 10 J=1,IT
            IF(NX(M).EQ.NI(J)) GO TO 6
 10      CONTINUE
         IZ=IZ+1
         R(IZ+k)=RR
         NI(IZ+k)=NX(M)
C
C     SOUTH
C
 6       M=IDY-i+1
         IF(M.LE.0) THEN
            AAS=R(k)
            GO TO 5
         END IF
C
         XX=XO(NY(M))
         YY=YO(NY(M))
         XR=XX-XP
         AAS=YY-YP
         RR=SQRT(XR*XR+AAS*AAS)
C
         IF(RR.GE.R(k)) GO TO 5
         IT=IZ+k
         DO 9 J=1,IT
            IF(NY(M).EQ.NI(J)) GO TO 5
 9       CONTINUE
         IZ=IZ+1
         R(IZ+k)=RR
         NI(IZ+k)=NY(M)
C
C     NORTH
C
 5       M=IDY+i
         IF(M.GT.N) THEN
            AAN=R(k)
            GO TO 7
         END IF
C
         XX=XO(NY(M))
         YY=YO(NY(M))
         XR=XX-XP
         AAN=YY-YP
         RR=SQRT(XR*XR+AAN*AAN)
C
         IF(RR.GE.R(k)) GO TO 7
         IT=IZ+k
         DO 11 J=1,IT
            IF(NY(M).EQ.NI(J)) GO TO 7
 11      CONTINUE
         IZ=IZ+1
         R(IZ+k)=RR
         NI(IZ+k)=NY(M)
C
 7       AANO=SQRT(AAN*AAN+AAO*AAO)
         AASO=SQRT(AAS*AAS+AAO*AAO)
         AASW=SQRT(AAS*AAS+AAW*AAW)
         AANW=SQRT(AAN*AAN+AAW*AAW)
         AAMIN=MIN(AANO,AASO,AASW,AANW)
C
c         write(*,*) iz+k,(ni(j),j=1,iz+k)
         IF(IZ.NE.0) THEN
            IZ=IZ+k
            CALL dsort(R,NI,IZ,kflag)
         ENDIF
C
         nboundary=ks
         do j=nboundary+1,k
            IF(R(j).LE.AAMIN) THEN
               ks=j
               if(ks.eq.k) exit loop
            ELSE
               i=i+1
               cycle loop
            ENDIF
         enddo
C
      enddo loop
C
      do i=1,k
          node(i)=ni(i)
c          write(*,*) 'near2d ',i,node(i),r(i),
c     &         xo(ni(i)),yo(ni(i)),xp,yp
      enddo
c
      RETURN
      END
