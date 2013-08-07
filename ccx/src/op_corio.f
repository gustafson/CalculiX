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
C
C-----MATRIX-VECTOR MULTIPLY FOR REAL SPARSE ANTISYMMETRIC MATRICES-----
C     i.e. the transpose is the negative matrix: A^T=-A
C
      SUBROUTINE OP_corio(n,p,W,U,ad,asd,icol,irow,nzl)
!
      implicit none
!
C-----------------------------------------------------------------------
      INTEGER  IROW(*),ICOL(*),n,nzl,l,lfirst,i,j,llast
      real*8   U(*),W(*),Asd(*),AD(*),p(*)
C-----------------------------------------------------------------------
C    SPARSE MATRIX-VECTOR MULTIPLY FOR LANCZS  U = A*W
c    the vector p is not needed but is kept for compatibility reasons
c    with the calling program
C-----------------------------------------------------------------------
C
C     COMPUTE THE DIAGONAL TERMS
      DO 10 I = 1,N
         U(I) = AD(I)*W(I)
 10   CONTINUE
C
C     COMPUTE BY COLUMN
      LLAST = 0
      DO 30 J = 1,NZL
C
         IF (ICOL(J).EQ.0) GO TO 30
         LFIRST = LLAST + 1
         LLAST = LLAST + ICOL(J)
C
         DO 20 L = LFIRST,LLAST
            I = IROW(L)
C
            U(I) = U(I) + Asd(L)*W(J)
            U(J) = U(J) - Asd(L)*W(I)
C
 20      CONTINUE
C
 30   CONTINUE
C
      RETURN
      END




