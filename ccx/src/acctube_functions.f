!
!     user subroutine for acc_tube.
!
      real*8 function calc_next_hole_abs_position(x,tube_length,
     &   hole_dist_option,hole_number_distribution,
     &   hole_interval_factor)
!     
      implicit none
!      
      integer
     &hole_dist_option,
     &section,
     &section_length,
     &find_section
!      
      real*8 
     &x,
     &tube_length,
     &hole_number_distribution(50,2),
     &hole_interval_factor,
     &interval,
     &next_hole_position,
     &gradient,
     &calc_section_rel_length
!
      calc_next_hole_abs_position=0.d0
!     
      return
!
      end
!
!=======================================================================
!
      subroutine MACHPI (MACH, PI,kappa, rgas)
!
!-----------------------------------------------------------------------
!                                                                      |
!     Dieses Unterprogramm berechnet die Mach-Zahl fuer das            |
!     eingegebene Druckverhaeltnis PI.                                 |
!                                                                      |
!     Eingabe-Groessen:                                                |
!       PI     = Druckverhaeltnis PS/PT                                |
!                                                                      |
!     Ausgabe-Groessen:                                                |
!       MACH   = Mach-Zahl                                             |
!                                                                      |
!-----------------------------------------------------------------------
!
      IMPLICIT CHARACTER*1 (A-Z)
      real*8    PI, MACH, MA2, kappa, rgas, kappam,kappax,pikrit
!
!-----------------------------------------------------------------------
!
      kappax = (kappa-1)/kappa
      KAPPAM = 2. / (KAPPA - 1.)
      PIKRIT = (2./(KAPPA+1.)) ** (KAPPA/(KAPPA-1.))
!
      IF (PI.GE.1.) THEN
!       Druckverhaeltnis groesser gleich 1
        MACH = 0.
      ELSEIF (PI.GT.PIKRIT) THEN
!       Druckverhaeltnis unterkritisch
        MA2  = KAPPAM * (PI**(-KAPPAX) - 1.)
        IF (MA2.GT.0) THEN
          MACH = SQRT (MA2)
        ELSE
          MACH = 0.
        ENDIF
      ELSEIF (PI.GT.0.) THEN
!       Druckverhaeltnis ueberkritisch
        MACH = 1.
      ELSE
!       Druckverhaeltnis ungueltig
        MACH = 1.E20
      ENDIF
!
      RETURN
      END
!
!=======================================================================
!
      subroutine WPI (W, PI, Q, SQTT,kappa,RGAS)
!
!-----------------------------------------------------------------------
!                                                                      |
!     Dieses Unterprogramm berechnet die Stroemungs-Geschwindigkeit    |
!     fuer das eingegebene Druckverhaeltnis PI.                        |
!                                                                      |
!     Eingabe-Groessen:                                                |
!       PI     = Druckverhaeltnis PS/PT                                |
!       Q      = reduzierter Durchsatz                                 |
!       SQTT   = SQRT (Totaltemperatur)                                |
!                                                                      |
!     Ausgabe-Groessen:                                                |
!       W      = Stroemungs-Geschwindigkeit                            |
!                                                                      |
!-----------------------------------------------------------------------
!
      IMPLICIT CHARACTER*1 (A-Z)
!       INCLUDE 'comkapfk.inc'
      real*8    W, PI, Q, SQTT,kappaq,kappa,RGAS,pikrit,kappah,wkritf
!
!-----------------------------------------------------------------------
!
      kappaq = 1/kappa
      PIKRIT = (2./(KAPPA+1.)) ** (KAPPA/(KAPPA-1.))
!
      KAPPAH = 2. * KAPPA / (KAPPA + 1.)
      WKRITF = SQRT( KAPPAH * RGAS )
!
      IF (PI.GE.1.) THEN
!       Druckverhaeltnis groesser gleich 1
        W    = 0.
      ELSEIF (PI.GT.PIKRIT) THEN
!       Druckverhaeltnis unterkritisch
        IF (Q.GT.0.) THEN
          W    = Q * RGAS * SQTT * PI**(-KAPPAQ)
        ELSE
          W    = 0.
        ENDIF
      ELSEIF (PI.GT.0.) THEN
!       Druckverhaeltnis ueberkritisch
        W    = WKRITF * SQTT
      ELSE
!       Druckverhaeltnis ungueltig
        W    = 1.E20
      ENDIF
!
      RETURN
      END
      
      
      
