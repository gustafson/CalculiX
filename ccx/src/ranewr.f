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
      REAL FUNCTION RANEWR ()
C
C ERZEUGUNG GLEICHVERTEILTER ZUFALLSZAHLEN ZWISCHEN 0 UND 1.
C PORTABLER ZUFALLSZAHLENGENERATOR IN STANDARD F77
C
C AUTOR: H. PFOERTNER
C
C AENDERUNGSSTAND :
C 26.08.95 EXTERNAL RAEWIN ENGEFUEGT, UM UEBER BLOCKDATA-LINK
C          STARTBELEGUNG AUCH OHNE INIRAN-AUFRUF ZU ERZWINGEN
C 07.12.92 BASISVERSION
C
C LITERATUR: WICHMANN AND HILL: APPL. STATIST. (JRSSC),
C                               (31) 188-190, (1982)
C
C GEDAECHTNIS:
C MUSS VOR DEM ERSTEN AUFRUF VON RANEWR DURCH EINEN AUFRUF VON
C INIRAN VORBELEGT WERDEN.
      INTEGER IX, IY, IZ
      COMMON /XXXRAN/ IX, IY, IZ
C
      EXTERNAL RAEWIN
C
C MODULO-OPERATIONEN
      IX = 171 * MOD ( IX, 177) -  2 * ( IX / 177 )
      IY = 172 * MOD ( IY, 176) - 35 * ( IY / 176 )
      IZ = 170 * MOD ( IZ, 178) - 63 * ( IZ / 178 )
C
C AUF POSITIVEN BEREICH BRINGEN
      IF ( IX .LT. 0 ) IX = IX + 30269
      IF ( IY .LT. 0 ) IY = IY + 30307
      IF ( IZ .LT. 0 ) IZ = IZ + 30323
C
C ZAHL ZWISCHEN 0 UND 1 ERZEUGEN
      RANEWR = MOD ( REAL(IX) / 30269.0
     &             + REAL(IY) / 30307.0
     &             + REAL(IZ) / 30323.0,  1.0 )
C
      RETURN
C ENDE DER FUNCTION RANEWR
      END
C *******************************************************************
      SUBROUTINE INIRAN(i1,i2,i3)
C
C STARTBELEGUNG FUER DEN ZUFALLSZAHLENGENERATOR RANEWR
C
C AUTOR: H. PFOERTNER
C
C AENDERUNGSSTAND :
C 07.12.92 BASISVERSION
C
C LITERATUR: WICHMANN AND HILL: APPL. STATIST. (JRSSC),
C                               (31) 188-190, (1982)
C
C GEDAECHTNIS:
      INTEGER IX, IY, IZ
      COMMON /XXXRAN/ IX, IY, IZ
C
C VORBELEGUNG
      IX = i1
      IY = i2
      IZ = i3
C
      RETURN
C ENDE DES UP. INIRAN
      END
C *******************************************************************
      BLOCKDATA RAEWIN
C
C ERZWINGUNG EINER STARTBELEGUNG (Z.B. BEI VERGESSENEM INIRAN-AUFRUF)
C
C HUGO PFOERTNER / OBERHACHING
C 26.08.95 BASISVERSION
C
      INTEGER IX, IY, IZ
      COMMON /XXXRAN/ IX, IY, IZ
      DATA IX, IY, IZ / 1974, 235, 337 /
      END
