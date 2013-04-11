*DECK ISORT
      SUBROUTINE ISORTIDDC2 (IX1,ix2, DY1,DY2,CY, N, KFLAG)
!
!     modified to sort in addition a double (dy) and char*20 (cy) array!
!
C***BEGIN PROLOGUE  ISORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2A
C***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   ISORT sorts array IX1 and optionally makes the same interchanges in
C   array IY.  The array IX1 may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IX1 - integer array of values to be sorted
C      IY - integer array to be (optionally) carried along
C      N  - number of values in integer array IX1 to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX1 in increasing order and carry IY along.
C            =  1  means sort IX1 in increasing order (ignoring IY)
C            = -1  means sort IX1 in decreasing order (ignoring IY)
C            = -2  means sort IX1 in decreasing order and carry IY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   810801  Modified by David K. Kahaner.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to IX1,IY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  ISORT
C     .. Scalar Arguments ..
      implicit none
c
      INTEGER KFLAG, N,iside,istat
C     .. Array Arguments ..
      INTEGER IX1(2,*),ix2(2,*)
      real*8 DY1(2,*),DY2(2,*)
      character*20 CY(*)
C     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, J, K, KK, L, M, NN, T, TT,tx21,tx12,tx22,
     &  ttx21,ttx12,ttx22
      real*8 TTY11,TTY12,TY11,TY12,TTY21,TTY22,TY21,TY22
      character*20 UUY,UY
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
!      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  ISORT
!
      do i=1,n
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix1(1,i)=10*ix1(1,i)+iside
      enddo
!
      NN = N
      IF (NN .LT. 1) THEN
!         CALL XERMSG ('SLATEC', 'ISORT',
!     +      'The number of values to be sorted is not positive.', 1, 1)
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
!         CALL XERMSG ('SLATEC', 'ISORT',
!     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
!     +      1)
         RETURN
      ENDIF
C
C     Alter array IX1 to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IX1(1,I) = -IX1(1,I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort IX1 only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX1(1,IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX1(1,I) .GT. T) THEN
         IX1(1,IJ) = IX1(1,I)
         IX1(1,I) = T
         T = IX1(1,IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (IX1(1,J) .LT. T) THEN
         IX1(1,IJ) = IX1(1,J)
         IX1(1,J) = T
         T = IX1(1,IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX1(1,I) .GT. T) THEN
            IX1(1,IJ) = IX1(1,I)
            IX1(1,I) = T
            T = IX1(1,IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (IX1(1,L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (IX1(1,K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX1(1,L)
         IX1(1,L) = IX1(1,K)
         IX1(1,K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = IX1(1,I+1)
      IF (IX1(1,I) .LE. T) GO TO 80
      K = I
C
   90 IX1(1,K+1) = IX1(1,K)
      K = K-1
      IF (T .LT. IX1(1,K)) GO TO 90
      IX1(1,K+1) = T
      GO TO 80
C
C     Sort IX1 and carry IY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX1(1,IJ)
      TY11 = DY1(1,IJ)
      TY21 = DY1(2,IJ)
      TY12 = DY2(1,IJ)
      TY22 = DY2(2,IJ)
      TX21 = IX1(2,IJ)
      tx12=ix2(1,ij)
      tx22=ix2(2,ij)
      uy = cy(ij)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX1(1,I) .GT. T) THEN
         IX1(1,IJ) = IX1(1,I)
         IX1(1,I) = T
         T = IX1(1,IJ)
         DY1(1,IJ) = DY1(1,I)
         DY1(2,IJ) = DY1(2,I)
         DY2(1,IJ) = DY2(1,I)
         DY2(2,IJ) = DY2(2,I)
         IX1(2,IJ) = IX1(2,I)
         ix2(1,ij)=ix2(1,i)
         ix2(2,ij)=ix2(2,i)
         cy(ij) = cy(i)
         DY1(1,I) = TY11
         DY1(2,I) = TY21
         DY2(1,I) = TY12
         DY2(2,I) = TY22
         IX1(2,I) = TX21
         ix2(1,i)=tx12
         ix2(2,i)=tx22
         cy(i) = uy
         TY11 = DY1(1,IJ)
         TY21 = DY1(2,IJ)
         TY12 = DY2(1,IJ)
         TY22 = DY2(2,IJ)
         TX21 = IX1(2,IJ)
         tx12=ix2(1,ij)
         tx22=ix2(2,ij)
         uy = cy(ij)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IX1(1,J) .LT. T) THEN
         IX1(1,IJ) = IX1(1,J)
         IX1(1,J) = T
         T = IX1(1,IJ)
         DY1(1,IJ) = DY1(1,J)
         DY1(2,IJ) = DY1(2,J)
         DY2(1,IJ) = DY2(1,J)
         DY2(2,IJ) = DY2(2,J)
         IX1(2,IJ) = IX1(2,J)
         ix2(1,ij)=ix2(1,j)
         ix2(2,ij)=ix2(2,j)
         cy(ij) = cy(j)
         DY1(1,J) = TY11
         DY1(2,J) = TY21
         DY2(1,J) = TY12
         DY2(2,J) = TY22
         IX1(2,J) = TX21
         ix2(1,j)=tx12
         ix2(2,j)=tx22
         cy(j) = uy
         TY11 = DY1(1,IJ)
         TY21 = DY1(2,IJ)
         TY12 = DY2(1,IJ)
         TY22 = DY2(2,IJ)
         TX21 = IX1(2,IJ)
         tx12=ix2(1,ij)
         tx22=ix2(2,ij)
         uy = cy(ij)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX1(1,I) .GT. T) THEN
            IX1(1,IJ) = IX1(1,I)
            IX1(1,I) = T
            T = IX1(1,IJ)
            DY1(1,IJ) = DY1(1,I)
            DY1(2,IJ) = DY1(2,I)
            DY2(1,IJ) = DY2(1,I)
            DY2(2,IJ) = DY2(2,I)
            IX1(2,IJ) = IX1(2,I)
            ix2(1,ij)=ix2(1,i)
            ix2(2,ij)=ix2(2,i)
            cy(ij) = cy(i)
            DY1(1,I) = TY11
            DY1(2,I) = TY21
            DY2(1,I) = TY12
            DY2(2,I) = TY22
            IX1(2,I) = TX21
            ix2(1,i)=tx12
            ix2(2,i)=tx22
            cy(i) = uy
            TY11 = DY1(1,IJ)
            TY21 = DY1(2,IJ)
            TY12 = DY2(1,IJ)
            TY22 = DY2(2,IJ)
            TX21 = IX1(2,IJ)
            tx12=ix2(1,ij)
            tx22=ix2(2,ij)
            uy = cy(ij)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IX1(1,L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IX1(1,K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX1(1,L)
         IX1(1,L) = IX1(1,K)
         IX1(1,K) = TT
         TTY11 = DY1(1,L)
         TTY21 = DY1(2,L)
         TTY12 = DY2(1,L)
         TTY22 = DY2(2,L)
         TTX21 = IX1(2,L)
         ttx12=ix2(1,l)
         ttx22=ix2(2,l)
         uuy = cy(l)
         DY1(1,L) = DY1(1,K)
         DY1(2,L) = DY1(2,K)
         DY2(1,L) = DY2(1,K)
         DY2(2,L) = DY2(2,K)
         IX1(2,L) = IX1(2,K)
         ix2(1,l)=ix2(1,k)
         ix2(2,l)=ix2(2,k)
         cy(l) = cy(k)
         DY1(1,K) = TTY11
         DY1(2,K) = TTY21
         DY2(1,K) = TTY12
         DY2(2,K) = TTY22
         IX1(2,K) = TTX21
         ix2(1,k)=ttx12
         ix2(2,k)=ttx22
         cy(k) = uuy
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = IX1(1,I+1)
      TY11 = DY1(1,I+1)
      TY21 = DY1(2,I+1)
      TY12 = DY2(1,I+1)
      TY22 = DY2(2,I+1)
      TX21 = IX1(2,I+1)
      tx12=ix2(1,i+1)
      tx22=ix2(2,i+1)
      uy = cy(i+1)
      IF (IX1(1,I) .LE. T) GO TO 170
      K = I
C
  180 IX1(1,K+1) = IX1(1,K)
      DY1(1,K+1) = DY1(1,K)
      DY1(2,K+1) = DY1(2,K)
      DY2(1,K+1) = DY2(1,K)
      DY2(2,K+1) = DY2(2,K)
      IX1(2,K+1) = IX1(2,K)
      ix2(1,k+1)=ix2(1,k)
      ix2(2,k+1)=ix2(2,k)
      cy(k+1) = cy(k)
      K = K-1
      IF (T .LT. IX1(1,K)) GO TO 180
      IX1(1,K+1) = T
      DY1(1,K+1) = TY11
      DY1(2,K+1) = TY21
      DY2(1,K+1) = TY12
      DY2(2,K+1) = TY22
      IX1(2,K+1) = TX21
      ix2(1,k+1)=tx12
      ix2(2,k+1)=tx22
      cy(k+1) = uy
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IX1(1,I) = -IX1(1,I)
  200    CONTINUE
      ENDIF
!
      do i=1,nn
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix1(1,i)=(ix1(1,i)-iside)/10
      enddo
!
      RETURN
      END
