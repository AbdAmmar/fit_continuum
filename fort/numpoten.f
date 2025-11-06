!C
!C   PROGRAM NUMPOTEN (demo for subroutine package RADIAL)
!C
!C     This program solves the radial wave equation for a numerical
!C  potential V(R) such that R*V(R) is finite for all R.
!C
!C     The potential is read from an input file [2 columns with values
!C  of R and R*V(R), in free format]. The radial grid R(I) _must_ include
!C  the origin (R=0), and extend up to radial distances for which the
!C  function R*V(R) reaches its (constant) asymptotic value. The input
!C  grid points must be in non-decreasing order, i.e. R(I+1).GE.R(I).
!C
      INCLUDE 'radial.f'  ! File included to simplify compilation.
!C
!C  NB: RADIAL must be compiled first to make global parameters available
!C      to other subprograms.
!C
!C  *********************************************************************
!C                       MAIN PROGRAM
!C  *********************************************************************
!C
!C     Solves the radial wave equations for the input potential.
!C
      USE CONSTANTS
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER FILEN*20,VFNAME*20
      PARAMETER (PI=3.1415926535897932D0)
      DIMENSION DR0(NDIM)  ! Output from SGRID.
!C  ****  Potential.
      DIMENSION R0(NDIM),RV0(NDIM)
!C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
!C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
!C
!C  ****  Input potential.
!C
      WRITE(6,'(2X,2A)') 'Enter the (path-)filename of the input',
     1  ' potential file ...'
      WRITE(6,'(4X,2A)') '(2 columns with values of R and R*V(R),',
     1  ' free format)'
      READ(5,'(A20)') VFNAME
      WRITE(6,'(/A,A)') ' # Input potential file: ',VFNAME
!C
      OPEN(3,FILE=VFNAME)
      DO I=1,NDIM
 1      CONTINUE
        READ(3,*,ERR=1,END=10) R0(I),RV0(I)
        NV=I
      ENDDO
 10   CONTINUE
      CLOSE(3)
!C
      WRITE(6,'(/A,I5)') ' # Potential grid. Number of radii =',NV
      CALL SPLERR(R0,RV0,0.0D0,0.0D0,ERR,NV,1)
      WRITE(6,'(A,1P,E9.1)') ' # Spline interpolation error =',ERR
!C
!C  ****  Spline interpolation of the potential.
      CALL VINT(R0,RV0,NV)
      RANGE=VRANGE()
      WRITE(6,'(A,1P,E9.1)') ' # Range of the potential =',RANGE
!C
!C
!C -\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-\/-
!C  ****  Print tables of the potential function and its interpolating
!C        spline.
!C
      OPEN(8,FILE='potential.dat')
      WRITE(8,'(A)') '# Potential function.'
      WRITE(8,'(''#'',2X,''r'',13X,''r*V(r)'')')
      DO I=1,NV
        WRITE(8,'(1X,1P,2E14.6,E10.1)') R0(I),RV0(I)
      ENDDO
      CLOSE(8)
!C
      OPEN(8,FILE='pot-spline.dat')
      WRITE(8,'(A)') '# Interpolated potential function.'
      WRITE(8,'(''#'',2X,''r'',13X,''r*V(r)'')')
      NINT=4
      DO I=1,NV
        IF(I.LT.NV) THEN
          DR=(R0(I+1)-R0(I))/DBLE(NINT)
        ELSE
          DR=(R0(NV)-R0(NV-1))/DBLE(NINT)
        ENDIF
        IF(DR.GT.1.0D-12*ABS(R0(I))) THEN
          DO J=1,NINT
            RINT=R0(I)+(J-1)*DR
            IF(J.EQ.1.AND.I.GT.1) THEN
              IF(R0(I)-R0(I-1).LT.1.0D-12*ABS(R0(I)))
     1          RINT=R0(I)*(1.0D0+1.0D-12)
            ENDIF
            RVS=RVSPL(RINT)
            WRITE(8,'(1X,1P,2E14.6)') RINT,RVS
          ENDDO
        ELSE
          RINT=R0(I)*(1.0D0-1.0D-12)
          RVS=RVSPL(RINT)
          WRITE(8,'(1X,1P,2E14.6,''  #'')') RINT,RVS
        ENDIF
      ENDDO
      CLOSE(8)
!C -/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-
!C
!C
      OPEN(8,FILE='resn.dat')
      WRITE(8,'(/A,A)') ' # Input potential file: ',VFNAME
      WRITE(8,'(A,I5)') ' # Potential grid. Number of radii =',NV
      WRITE(8,'(A,1P,E9.1)') ' # Spline interpolation error =',ERR
      WRITE(8,'(A,1P,E9.1)') ' # Range of the potential =',RANGE
!C
!C  ****  High-energy limit of the Dirac inner phase shift.
!C
      CALL DELINF(HEDEL)
      WRITE(6,'(/,'' # Delta_infty = '',1P,E22.15)') HEDEL
      Z=RV0(NV)
!C
!C  ****  Pre-defined user grid.
!C
      WRITE(6,1003)
 1003 FORMAT(/2X,'The user grid can be read from a file (single',
     1  ' column, increasing radii)'/2X,'or determined automati',
     2  'cally. If you wish to use a prepared grid, enter'/2X,
     1  'the file name, otherwise type ''n'' ...')
      READ(5,'(A12)') FILEN
      IF(FILEN.NE.'N'.AND.FILEN.NE.'n') THEN
        IGRID=1
        OPEN(7,FILE=FILEN)
        DO I=1,NDIM
          READ(7,*,END=20) RAD(I)
          NGP=I
        ENDDO
        CLOSE(7)
      ELSE
        IGRID=0
      ENDIF
!C
!C  ************  Calculation of radial functions.
!C
 20   CONTINUE
      WRITE(6,*) '  '
      WRITE(6,*) ' Select one option ...'
      WRITE(6,*) '   1: Schrodinger bound state,   2: Sch',
     1  'rodinger free state,'
      WRITE(6,*) '   3: Dirac bound state,         4: Dir',
     2  'ac free state.'
      IF(Z.LT.-0.5D0)
     1  WRITE(6,*) '   5: Quantum defect.'
      READ(5,*,ERR=20) IOPT
      IF(IOPT.LT.0.OR.IOPT.GT.5) STOP
      IF(IOPT.EQ.5.AND.Z.GT.-0.5D0) THEN
        WRITE(6,'(/2X,'' Quantum defects are not defined for attrac'',
     1    ''tive potentials,'',/2X,''Press any key to continue...'')')
        READ(5,*)
        GO TO 20
      ENDIF
      WRITE(6,*) '  '
      WRITE(8,*) '  '
!C
!C  ****  Schrodinger equation. Bound state.
!C
      IF(IOPT.EQ.1) THEN
        WRITE(6,*) ' Enter N, L and EPS ...'
        READ(5,*,ERR=20) N,L,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(N.LT.1) THEN
          WRITE(6,'(A,I5)') '  N =',N
          WRITE(6,'(A)') '  N must be >0.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
        IF(L.LT.0.OR.L.GE.N) THEN
          WRITE(6,'(A,2I5)') '  L,N =',L,N
          WRITE(6,'(A)') '  L must be >0 and <N.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
!C
        IF(IGRID.EQ.0) THEN
          NGP=4000
          RN=3500.0D0
          CALL SGRID(RAD,DR0,RN,1.0D-6,1.0D0,NGP,NDIM,IERS)
          IF(IERS.NE.0) STOP 'Error in the grid definition (SB).'
        ENDIF
        E=-Z**2/(2.0D0*N*N)
        CALL SBOUND(E,EPS,N,L)
        IF(IER.NE.0) GO TO 20
!C
        WRITE(6,1101) VFNAME,N,L,EPS,E
        WRITE(8,1101) VFNAME,N,L,EPS,E
 1101   FORMAT(1X,1P,'# **** Schrodinger Eq. Numerical potential.',
     1    /' #',6X,'Input potential file: ',A20,
     2    /' #',6X,'Bound state: N=',I4,', L=',I4,'  (EPS=',E8.1,')'
     3    /' #',6X,'Binding energy=',E22.15)
        FILEN='schrodinger.dat'
!C
!C  ****  Schrodinger equation. Free state.
!C
      ELSE IF(IOPT.EQ.2) THEN
        WRITE(6,*) ' Enter E, L and EPS ...'
        READ(5,*,ERR=20) E,L,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(E.LT.0.0D0) THEN
          WRITE(6,'(A,1P,E14.6)') '  E =',E
          WRITE(6,'(A)') '  The energy must be positive.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
        IF(L.LT.0) THEN
          WRITE(6,'(A,I5)') '  L =',L
          WRITE(6,'(A)') '  L must be a non-negative integer.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
!C
        IF(IGRID.EQ.0) THEN
          NGP=2000
          WAVEL=2.0D0*PI/SQRT(E+E)
          DRN=WAVEL/40.0D0
          RN=DRN*DBLE(NGP-300)
          CALL SGRID(RAD,DR0,RN,1.0D-6,DRN,NGP,NDIM,IERS)
          IF(IERS.NE.0) STOP 'Error in the grid definition (SF).'
        ENDIF
        CALL SFREE(E,EPS,PHASE,L,1)
        IF(IER.NE.0) THEN
          WRITE(6,'(A,I3)') 'Error in SFREE. IER =',IER
          GO TO 20
        ENDIF
!C
        WRITE(6,1201) VFNAME,E,L,EPS,PHASE,DELTA,ETA
        WRITE(8,1201) VFNAME,E,L,EPS,PHASE,DELTA,ETA
 1201   FORMAT(1X,1P,'# **** Schrodinger Eq. Numerical potential.',
     1    /' #',6X,'Input potential file: ',A20,
     2    /' #',6X,'Free state: E=',E13.6,', L=',I4,'  (EPS=',E8.1,')'
     3    /' #',6X,'  Inner phase shift=',E22.15,
     4    /' #',6X,'Coulomb phase shift=',E22.15,'  (ETA=',E13.6,')')
        WRITE(6,1202) ILAST,RAD(ILAST)
 1202   FORMAT(' #',6X,'Matching radius: RAD(',I5,')=',1P,E22.15)
        FILEN='schrodinger.dat'
!C
!C  ****  Dirac equation. Bound state.
!C
      ELSE IF(IOPT.EQ.3) THEN
        WRITE(6,*) ' Enter N, K and EPS ...'
        READ(5,*,ERR=20) N,K,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(N.LT.1) THEN
          WRITE(6,'(A,I5)') '  N =',N
          WRITE(6,'(A)') '  N must be >0.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
        IF(K.EQ.0.OR.K.GE.N.OR.K.LT.-N) THEN
          WRITE(6,'(A,2I5)') '  K,N =',K,N
          WRITE(6,'(A)') '  K must be between -N and N-1 and .NE.0'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
!C
        IF(IGRID.EQ.0) THEN
          NGP=4000
          RN=3500.0D0
          CALL SGRID(RAD,DR0,RN,1.0D-6,1.0D0,NGP,NDIM,IERS)
          IF(IERS.NE.0) STOP 'Error in the grid definition (DB).'
        ENDIF
        E=-Z**2/(2.0D0*N*N)
        CALL DBOUND(E,EPS,N,K)
        IF(IER.NE.0) GO TO 20
!C
        WRITE(6,1301) VFNAME,N,K,EPS,E
        WRITE(8,1301) VFNAME,N,K,EPS,E
 1301   FORMAT(1X,1P,'# **** Dirac equation. Numerical potential.',
     1    /' #',6X,'Input potential file: ',A20,
     2    /' #',6X,'Bound state: N=',I4,', K=',I4,'  (EPS=',E8.1,')'
     3    /' #',6X,'Binding energy=',E22.15)
        FILEN='dirac.dat'
!C
!C  ****  Dirac equation. Free state.
!C
      ELSE IF(IOPT.EQ.4) THEN
        WRITE(6,*) ' Enter E, K and EPS ...'
        READ(5,*,ERR=20) E,K,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(E.LT.0.0D0) THEN
          WRITE(6,'(A,1P,E14.6)') '  E =',E
          WRITE(6,'(A)') '  The energy must be positive.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
        IF(K.EQ.0) THEN
          WRITE(6,'(A,I5)') '  K =',K
          WRITE(6,'(A)') '  K must be an integer different from 0.'
          WRITE(6,'(/2X,''Press any key to continue...'')')
          READ(5,*)
          GO TO 20
        ENDIF
!C
        IF(K.LT.0) THEN
          L=-K-1
        ELSE
          L=K
        ENDIF
        IF(IGRID.EQ.0) THEN
          NGP=2000
          WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
          DRN=WAVEL/40.0D0
          RN=DRN*DBLE(NGP-300)
          CALL SGRID(RAD,DR0,RN,1.0D-6,DRN,NGP,NDIM,IERS)
          IF(IERS.NE.0) STOP 'Error in the grid definition (DF).'
        ENDIF
        CALL DFREE(E,EPS,PHASE,K,1)
        IF(IER.NE.0) THEN
          WRITE(6,'(A,I3)') 'Error in DFREE. IER =',IER
          GO TO 20
        ENDIF
!C
        WRITE(6,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
        WRITE(8,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
 1401   FORMAT(1X,1P,'# **** Dirac equation. Numerical potential.',
     1    /' #',6X,'Input potential file: ',A20,
     2    /' #',6X,'Free state: E=',E13.6,', K=',I4,'  (EPS=',E8.1,')'
     3    /' #',6X,'  Inner phase shift=',E22.15,
     4    /' #',6X,'Coulomb phase shift=',E22.15,'  (ETA=',E13.6,')')
        WRITE(6,1202) ILAST,RAD(ILAST)
        FILEN='dirac.dat'
!C
!C  ************  Quantum-defect function (Dirac).
!C
      ELSE IF(IOPT.EQ.5) THEN
        IF(Z.GT.-0.5D0) THEN
          WRITE(6,*) ' Quantum defects are defined only for attractive'
          WRITE(6,*) ' Coulomb fields (Z < -0.5).'
          GO TO 20
        ENDIF
        WRITE(6,*) ' Enter K, EPS ...'
        READ(5,*,ERR=20) K,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(K.EQ.0) THEN
          WRITE(6,'(A,I5)') ' K =',K
          WRITE(6,'(A)') ' K must be an integer different from 0.'
          GO TO 20
        ENDIF
        CALL QNTDEF(K,QD0,QDA,QDB,EPS,ERRM)
        WRITE(6,'('' #'',4X,''QD0 = '',1P,E13.6)') QD0
        WRITE(6,'('' #'',4X,''  A = '',1P,E13.6)') QDA
        WRITE(6,'('' #'',4X,''  B = '',1P,E13.6)') QDB
        WRITE(6,'(/'' #    Largest error (%) ='',1P,E13.6)') ERRM
        WRITE(6,*) ' '
        GO TO 20
      ELSE
        GO TO 10
      ENDIF
!C
!C  ****  Radial wave functions printed in output file.
!C
      READ(5,'(A)') FILEN  ! Uncomment to change the output filename.
      OPEN(10,FILE=FILEN)
 1501 FORMAT(1X,'# Radial wave functions calculated by RADIAL.')
      WRITE(10,1501)
      IF(IOPT.EQ.1) THEN
        WRITE(10,1101) VFNAME,N,L,EPS,E
      ELSE IF(IOPT.EQ.2) THEN
        WRITE(10,1201) VFNAME,E,L,EPS,PHASE,DELTA,ETA
        WRITE(10,1202) ILAST,RAD(ILAST)
      ELSE IF(IOPT.EQ.3) THEN
        WRITE(10,1301) VFNAME,N,K,EPS,E
      ELSE IF(IOPT.EQ.4) THEN
        WRITE(10,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
        WRITE(10,1202) ILAST,RAD(ILAST)
      ENDIF
!C
      NTAB=NGP
      DO I=NGP,1,-1
        IF(ABS(P(I)).GT.1.0D-35) THEN
          NTAB=I
          GO TO 30
        ENDIF
      ENDDO
 30   CONTINUE
!C
      WRITE(10,1502)
 1502 FORMAT(1X,'#',7X,'R',14X,'P(R)',12X,'Q(R)')
      DO I=1,NTAB
!C  ----  Do not print values less than 1.0D-99  ------------------------
        IF(ABS(P(I)).LT.1.0D-98) P(I)=0.0D0
        IF(ABS(Q(I)).LT.1.0D-98) Q(I)=0.0D0
!C  ---------------------------------------------------------------------
        WRITE(10,'(1X,1P,3E16.8)') RAD(I),P(I),Q(I)
      ENDDO
      CLOSE(10)
!C
      GO TO 20
      END
