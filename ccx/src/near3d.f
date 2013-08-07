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
      subroutine near3d(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,n,node,k)
C
c     determines the k closest nodes out of n with coordinates in
c     xo,yo,zo to the point with coordinates (xp,yp,zp);
c
c     order xo, yo and zo before the first call to near2d and
c     store the results with the corresponding permutation array
c     in x,y,z and nx,ny,nz, respectively
c
      implicit none
C
      integer node(*),nx(*),ny(*),nz(*),i,j,k,m,n,kflag,ni(26),ks,
     &  it,iz,idx,idy,idz,nboundary
c     
      real*8 xo(*),yo(*),zo(*),x(*),y(*),z(*),xp,yp,zp,r(26),xx,yy,
     &  zz,rr,aaw,aao,aas,aan,aau,aad,aanod,aanou,aasod,aasou,aaswd,
     &  aaswu,aanwd,aanwu,aamin,xr,yr,zr
c
      if(k.gt.n) then
         k=n
      endif
c
      if(k.gt.20) then
         write(*,*) '*ERROR in near3d: no more than 20 neighbors'
         write(*,*) '       can be identified'
         stop
      endif
c
      CALL IDENT(X,XP,N,IDX)
      CALL IDENT(Y,YP,N,IDY)
      CALL IDENT(Z,ZP,N,IDZ)
C
      kflag=2
C
C     DETERMINATION OF THE k NEAREST POINTS
C
      DO I=1,K
        XR=XO(I)-XP
        YR=YO(I)-YP
        zr=zo(i)-zp
        R(I)=DSQRT(XR*XR+YR*YR+zr*zr)
        NI(I)=I
      ENDDO
C
      CALL dsort(R,NI,K,kflag)
c
      i=1
      ks=0
C
      loop: do
c
      IZ=0
C
C     WEST
C
      M=IDX-I+1
      IF(M.LE.0) THEN
        AAW=R(K)
        GO TO 4
      END IF
C
      XX=XO(NX(M))
      YY=YO(NX(M))
      zz=zo(nx(m))
      AAW=XX-XP
      YR=YY-YP
      zr=zz-zp
      RR=DSQRT(AAW*AAW+YR*YR+zr*zr)
C
      IF(RR.GE.R(K)) GO TO 4
      IT=IZ+K
      DO 8 J=1,IT
        IF(NX(M).EQ.NI(J)) GO TO 4
8     CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=NX(M)
C
C     EAST
C
4     M=IDX+I
      IF(M.GT.N) THEN
        AAO=R(K)
        GO TO 6
      END IF
C
      XX=XO(NX(M))
      YY=YO(NX(M))
      zz=zo(nx(m))
      AAO=XX-XP
      YR=YY-YP
      zr=zz-zp
      RR=DSQRT(AAO*AAO+YR*YR+zr*zr)
C
      IF(RR.GE.R(K)) GO TO 6
      IT=IZ+K
      DO 10 J=1,IT
        IF(NX(M).EQ.NI(J)) GO TO 6
10    CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=NX(M)
C
C     SOUTH
C
6     M=IDY-I+1
      IF(M.LE.0) THEN
        AAS=R(K)
        GO TO 5
      END IF
C
      XX=XO(NY(M))
      YY=YO(NY(M))
      zz=zo(ny(m))
      XR=XX-XP
      AAS=YY-YP
      zr=zz-zp
      RR=DSQRT(XR*XR+AAS*AAS+zr*zr)
C
      IF(RR.GE.R(K)) GO TO 5
      IT=IZ+K
      DO 9 J=1,IT
        IF(NY(M).EQ.NI(J)) GO TO 5
9     CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=NY(M)
C
C     NORTH
C
5     M=IDY+I
      IF(M.GT.N) THEN
        AAN=R(K)
        GO TO 7
      END IF
C
      XX=XO(NY(M))
      YY=YO(NY(M))
      zz=zo(ny(m))
      XR=XX-XP
      AAN=YY-YP
      zr=zz-zp
      RR=DSQRT(XR*XR+AAN*AAN+zr*zr)
C
      IF(RR.GE.R(K)) GO TO 7
      IT=IZ+K
      DO 11 J=1,IT
      IF(NY(M).EQ.NI(J)) GO TO 7
11    CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=NY(M)
C
C     up
C
7     M=IDz-I+1
      IF(M.LE.0) THEN
        AAu=R(K)
        GO TO 20
      END IF
C
      XX=XO(Nz(M))
      YY=YO(Nz(M))
      zz=zo(nz(m))
      XR=XX-XP
      yr=YY-YP
      aau=zz-zp
      RR=DSQRT(XR*XR+yr*yr+aau*aau)
C
      IF(RR.GE.R(K)) GO TO 20
      IT=IZ+K
      DO 21 J=1,IT
        IF(Nz(M).EQ.NI(J)) GO TO 20
21    CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=Nz(M)
C
C     down
C
20    M=IDz+I
      IF(M.GT.N) THEN
        AAd=R(K)
        GO TO 22
      END IF
C
      XX=XO(Nz(M))
      YY=YO(Nz(M))
      zz=zo(nz(m))
      XR=XX-XP
      yr=YY-YP
      aad=zz-zp
      RR=DSQRT(XR*XR+yr*yr+aad*aad)
C
      IF(RR.GE.R(K)) GO TO 22
      IT=IZ+K
      DO 23 J=1,IT
      IF(Nz(M).EQ.NI(J)) GO TO 22
23    CONTINUE
      IZ=IZ+1
      R(IZ+K)=RR
      NI(IZ+K)=Nz(M)
C
22    AANOd=DSQRT(AAN*AAN+AAO*AAO+aad*aad)
      AANOu=DSQRT(AAN*AAN+AAO*AAO+aau*aau)
      AASOd=DSQRT(AAS*AAS+AAO*AAO+aad*aad)
      AASOu=DSQRT(AAS*AAS+AAO*AAO+aau*aau)
      AASWd=DSQRT(AAS*AAS+AAW*AAW+aad*aad)
      AASWu=DSQRT(AAS*AAS+AAW*AAW+aau*aau)
      AANWd=DSQRT(AAN*AAN+AAW*AAW+aad*aad)
      AANWu=DSQRT(AAN*AAN+AAW*AAW+aau*aau)
      AAMIN=MIN(AANOd,aanou,AASOd,aasou,AASWd,aaswu,AANWd,aanwu)
C
      IF(IZ.NE.0) THEN
        IZ=IZ+K
        CALL dsort(R,NI,IZ,kflag)
      ENDIF
C
      NBOUNDARY=KS
      DO J=NBOUNDARY+1,K
        IF(R(J).LE.AAMIN) THEN
          KS=J
          IF(KS.EQ.K) exit loop
        ELSE
          i=i+1
          cycle loop
        ENDIF
      ENDDO
C
      enddo loop
C
      do i=1,k
        node(i)=ni(i)
      enddo
c
      RETURN
      END
