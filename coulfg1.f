      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     1                  MODE1,KFN,IFAIL,M1)
CCCCCC
C
C    REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD
C
C    A. R. BARNETT           MANCHESTER  MARCH   1981
C
C    ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395
C                   + 'RCWFF'      IN    CPC 11 (1976) 141-142
C    FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
C    THIS VERSION WRITTEN UP       IN    CPC XX (1981) YYY-ZZZ
C
C    COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0),
C     AND REAL LAMDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMDA VALUES
C     THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
C     EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
C     THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS
C
C    FOR A RANGE OF LAMDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,
C    STARTING ARRAY ELEMENT IS M1 = MAX0(INT(XLMIN+ACCUR),0) + 1
C
C    IF 'MODE' = 1  GET F,G,F',G'     FOR INTEGER-SPACED LAMDA VALUES
C              = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN
C              = 3      F               CALL TO AT LEAST LENGTH (1)
C    IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED
C              = 1 SPHERICAL   BESSEL      "      "     "
C              = 2 CYLINDRICAL BESSEL      "      "     "
C    THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT
C
C    PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'
C       IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))
CCCCCC
C
      implicit real*8(a-h,o-z) 
      COMPLEX*16   C1,AA,BB,DD ,DL,PQ,TWOI
c      DIMENSION    FC(99),GC( 99),FCP(99),GCP(99)
      DIMENSION    FC(int(XLMAX-XLMIN)),GC(int(XLMAX-XLMIN))
      DIMENSION    FCP(int(XLMAX-XLMIN)),GCP(int(XLMAX-XLMIN))
      LOGICAL      ETANE0
C
      DATA ZERO,ONE,TWO,TEN4,ABORTL /0.E0,1.0E0,2.0E0,10000.E0,20000.E0/
      DATA HALF,TM30 / 0.5E0,1.0E-30 /
      DATA RT2DPI /0.79788 45608 02865 /
C ***      THIS CONSTANT IS  SQRT(TWO/PI)
      DATA C1,TWOI / ( 1.0E0,0.0E0 ), ( 0.0E0,2.0E0 ) /
      CXAMOD(AA) = ABS(DBLE(AA)) + ABS(IMAG(AA))
C
                        ACCUR = 1.0D-12
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION AS ABOVE
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      NFP   = 0
      NPQ   = 0
      ETA   = ETA1
      IF(KFN .NE. 0) ETA = ZERO
      ETANE0= ETA .NE. ZERO
      ACC   = ACCUR
      ACC4  = ACC*TEN4
      ACCH  = SQRT(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE.SQRT(ACCUR)
C
      IF(ABS(XX) .LE. ACCH)                    GO TO 100
      X     =   XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
C     IF(ABS(MOD(DELL,ONE)) .GT. ACCH) WRITE(6,2040)XLMAX,XLMIN,DELL
      LXTRA = INT(DELL)
      XLL   = XLM + FLOAT(LXTRA)
C ***         LXTRA IS NUMBER OF ADDITIONAL LAMDA VALUES TO BE COMPUTED
C ***         XLL  IS MAX LAMDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(INT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,X,ETA)/F(XL,X,ETA)
C
      XI  = ONE/X
      FCL = ONE
      RL  = ONE
      PK  = XLL + ONE
      PX  = PK + ABORTL
    2 EK  = ETA/PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 = PK + ONE
C ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
             IF(ABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK
      GO TO 2
    3 D   = ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
           IF(FCL .NE. ONE ) FCL = -ONE
           IF(D   .LT. ZERO) FCL = -FCL
      F   = F  + DF
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMDA + 1
C
    4 PK    = PK1
        PK1 = PK1 + ONE
        EK  = ETA / PK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   = TK - D*(ONE + EK*EK)
              IF(ABS(D) .GT. ACCH)             GO TO 5
              WRITE (6,1000) D,DF,ACCH,PK,EK,X,ETA
              RL = RL + ONE
              IF(  RL .GT. TWO)                 GO TO 110
    5 D     = ONE/D
              IF (D .LT. ZERO) FCL = -FCL
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
              IF(PK .GT. PX)                    GO TO 110
      IF(ABS(DF) .GE. ABS(F)*ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMDA = XLM.  ARRAY GC,IF PRESENT,STORES RL
C
C                TM30 = ONE
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = SQRT(ONE + EL*EL)
         SL    = EL + XL*XI
         L     = L1 - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
    6 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMDA = XLMIN = XLM
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C
    7 TA = TWO*ABORTL
      WI = TWO*ETA
      PK = ZERO
      PQ = CMPLX(ZERO,ONE - ETA*XI)
      AA = CMPLX(-E2MM1,ETA)
      BB = CMPLX(TWO*(X - ETA),TWO)
      DD = C1/BB
      DL = AA*DD*CMPLX(ZERO,XI)
    8 PQ    = PQ + DL
         PK = PK + TWO
         AA = AA + CMPLX(PK,WI)
         BB = BB + TWOI
         DD = C1/(AA*DD + BB)
         DL = DL*(BB*DD - C1)
            IF(PK .GT. TA)                      GO TO 120
      IF(CXAMOD(DL) .GE. CXAMOD(PQ)*ACC)        GO TO 8
      P     = DBLE(PQ + DL)
      Q     = IMAG(PQ + DL)
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/MIN(ABS(Q),ONE)
                      IF(ABS(P) .GT. ABS(Q)) PACCQ = PACCQ*ABS(P)
C
C *** SOLVE FOR FCM = F AT LAMDA = XLM,THEN FIND NORM FACTOR W=W/FCM
C
C???            IF(Q .LE. ACC4*ABS(P))             GO TO 130
      GAM = (F - P)/Q
      W   = ONE/SQRT((F - P)*GAM + Q)
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
                      ALPHA = ZERO
      IF(KFN  .EQ. 1) ALPHA = XI
      IF(KFN  .EQ. 2) ALPHA = XI*HALF
      IF(KFN  .EQ. 1) W     = W*XI
      IF(KFN  .EQ. 2) W     = W*SQRT(XI)*RT2DPI
    9 FCM = SIGN(W,FCL)
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 10
           GCL     = FCM*GAM
                      IF(KFN .NE. 0 ) GCL = -GCL
           GC(M1)  =  GCL
           GPL     = GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 10
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   10 IF(LXTRA.EQ. 0 .OR. IFAIL .NE. 0) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = W/ABS(FCL)
         LMAX = L1 - 1
      DO 11 L = M1,LMAX
                    IF(MODE .EQ. 3)             GO TO 11
                    XL = XL + ONE
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = GC (L+1)
                    SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                    IF(MODE .EQ. 2)             GO TO 11
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   11 FC(L+1)     = W*FC(L+1)
      RETURN
   20 W    = ZERO
      GO TO 9
 1000 FORMAT(/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,X,ETA = ',1P7D9.2/)
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1PD12.3,' TRY SMALL-X ASYMPTOTIC SOLUTIONS'/
     *      ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     *1P3D15.6//)
      RETURN
  110 IFAIL =  1
      WRITE (6,2010) ABORTL,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',1P5D12.3//)
      GO TO 20
  120 IFAIL =  2
      WRITE (6,2020) ABORTL,PQ,DL,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DL,ACCUR =  ',1P4D17.7,D12.3//)
      PRINT *,XX,ETA1,XLMIN,XLMAX, MODE1,KFN,IFAIL,M1
C     call flush(6)
C     call abortp('COULFG')
      GO TO 20
  130 IFAIL =  3
      WRITE (6,2030) PQ,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE.ABS(P)*ACC*10**4 , P,Q,ACC = ',1P3D12.3,4X,
     *' DELL,LXTRA,M1 = ',D12.3,2I5 /)
      GO TO 20
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3D20.10/)
      END
