
DOUBLE PRECISION FUNCTION DUMACH ()
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
DOUBLE PRECISION U, COMP
!***FIRST EXECUTABLE STATEMENT  DUMACH
U = 1.0D0
10   U = U*0.5D0
CALL DUMSUM(1.0D0, U, COMP)
IF (COMP .NE. 1.0D0) GOTO 10
DUMACH = U*2.0D0
RETURN
!----------------------- End of Function DUMACH ------------------------
END
SUBROUTINE DUMSUM(A,B,C)
!     Routine to force normal storing of A + B, for DUMACH.
DOUBLE PRECISION A, B, C
C = A + B
RETURN
END

SUBROUTINE DCFODE (METH, ELCO, TESCO)
!***BEGIN PROLOGUE  DCFODE
!***SUBSIDIARY
!***PURPOSE  Set ODE integrator coefficients.
!***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DCFODE is called by the integrator routine to set coefficients
!  needed there.  The coefficients for the current method, as
!  given by the value of METH, are set for all orders and saved.
!  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
!  (A smaller value of the maximum order is also allowed.)
!  DCFODE is called once at the beginning of the problem,
!  and is not called again unless and until METH is changed.
!
!  The ELCO array contains the basic method coefficients.
!  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
!  order nq are stored in ELCO(i,nq).  They are given by a genetrating
!  polynomial, i.e.,
!      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
!  For the implicit Adams methods, l(x) is given by
!      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
!  For the BDF methods, l(x) is given by
!      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
!  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
!  The TESCO array contains test constants used for the
!  local error test and the selection of step size and/or order.
!  At order nq, TESCO(k,nq) is used for the selection of step
!  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
!  nq + 1 if k = 3.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DCFODE
!**End
INTEGER METH
INTEGER I, IB, NQ, NQM1, NQP1
DOUBLE PRECISION ELCO, TESCO
DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,RQFAC, RQ1FAC, TSIGN, XPIN
DIMENSION ELCO(13,12), TESCO(3,12)
DIMENSION PC(12)
!
!***FIRST EXECUTABLE STATEMENT  DCFODE
GOTO (100, 200), METH
!
100  ELCO(1,1) = 1.0D0
ELCO(2,1) = 1.0D0
TESCO(1,1) = 0.0D0
TESCO(2,1) = 2.0D0
TESCO(1,2) = 1.0D0
TESCO(3,12) = 0.0D0
PC(1) = 1.0D0
RQFAC = 1.0D0
DO 140 NQ = 2,12
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
RQ1FAC = RQFAC
RQFAC = RQFAC/NQ
NQM1 = NQ - 1
FNQM1 = NQM1
NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq-1). ----------------------------------
PC(NQ) = 0.0D0
DO 110 IB = 1,NQM1
I = NQP1 - IB
110      PC(I) = PC(I-1) + FNQM1*PC(I)
PC(1) = FNQM1*PC(1)
! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
PINT = PC(1)
XPIN = PC(1)/2.0D0
TSIGN = 1.0D0
DO 120 I = 2,NQ
TSIGN = -TSIGN
PINT = PINT + TSIGN*PC(I)/I
120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
! Store coefficients in ELCO and TESCO. --------------------------------
ELCO(1,NQ) = PINT*RQ1FAC
ELCO(2,NQ) = 1.0D0
DO 130 I = 2,NQ
130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
AGAMQ = RQFAC*XPIN
RAGQ = 1.0D0/AGAMQ
TESCO(2,NQ) = RAGQ
IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
TESCO(3,NQM1) = RAGQ
140    CONTINUE
RETURN
!
200  PC(1) = 1.0D0
RQ1FAC = 1.0D0
DO 230 NQ = 1,5
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
FNQ = NQ
NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq). ------------------------------------
PC(NQP1) = 0.0D0
DO 210 IB = 1,NQ
I = NQ + 2 - IB
210      PC(I) = PC(I-1) + FNQ*PC(I)
PC(1) = FNQ*PC(1)
! Store coefficients in ELCO and TESCO. --------------------------------
DO 220 I = 1,NQP1
220      ELCO(I,NQ) = PC(I)/PC(2)
ELCO(2,NQ) = 1.0D0
TESCO(1,NQ) = RQ1FAC
TESCO(2,NQ) = NQP1/ELCO(1,NQ)
TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
RQ1FAC = RQ1FAC/FNQ
230    CONTINUE
RETURN
!----------------------- END OF SUBROUTINE DCFODE ----------------------
END

SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
!***BEGIN PROLOGUE  DINTDY
!***SUBSIDIARY
!***PURPOSE  Interpolate solution derivatives.
!***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DINTDY computes interpolated values of the K-th derivative of the
!  dependent variable vector y, and stores it in DKY.  This routine
!  is called within the package with K = 0 and T = TOUT, but may
!  also be called by the user for any K up to the current order.
!  (See detailed instructions in the usage documentation.)
!
!  The computed values in DKY are gotten by interpolation using the
!  Nordsieck history array YH.  This array corresponds uniquely to a
!  vector-valued polynomial of degree NQCUR or less, and DKY is set
!  to the K-th derivative of this polynomial at T.
!  The formula for DKY is:
!               q
!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
!              j=K
!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
!  communicated by COMMON.  The above sum is done in reverse order.
!  IFLAG is returned negative if either K or T is out of bounds.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  XERRWD
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   050427  Corrected roundoff decrement in TP. (ACH)
!***END PROLOGUE  DINTDY
!**End
INTEGER K, NYH, IFLAG
DOUBLE PRECISION T, YH, DKY
DIMENSION YH(NYH,*), DKY(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF
INTEGER LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6),&
ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,&
MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
DOUBLE PRECISION C, R, S, TP
CHARACTER*80 MSG
!
!***FIRST EXECUTABLE STATEMENT  DINTDY
IFLAG = 0
IF (K .LT. 0 .OR. K .GT. NQ) GOTO 80
TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
IF ((T-TP)*(T-TN) .GT. 0.0D0) GOTO 90
!
S = (T - TN)/H
IC = 1
IF (K .EQ. 0) GOTO 15
JJ1 = L - K
DO 10 JJ = JJ1,NQ
10     IC = IC*JJ
15   C = IC
DO 20 I = 1,N
20     DKY(I) = C*YH(I,L)
IF (K .EQ. NQ) GOTO 55
JB2 = NQ - K
DO 50 JB = 1,JB2
J = NQ - JB
JP1 = J + 1
IC = 1
IF (K .EQ. 0) GOTO 35
JJ1 = JP1 - K
DO 30 JJ = JJ1,J
30       IC = IC*JJ
35     C = IC
DO 40 I = 1,N
40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
50     CONTINUE
IF (K .EQ. 0) RETURN
55   R = H**(-K)
DO 60 I = 1,N
60     DKY(I) = R*DKY(I)
RETURN
!
80   MSG = 'DINTDY-  K (=I1) illegal      '
CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
IFLAG = -1
RETURN
90   MSG = 'DINTDY-  T (=R1) illegal      '
CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
IFLAG = -2
RETURN
!----------------------- END OF SUBROUTINE DINTDY ----------------------
END

SUBROUTINE DPREPJ (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, F, JAC)
!***BEGIN PROLOGUE  DPREPJ
!***SUBSIDIARY
!***PURPOSE  Compute and process Newton iteration matrix.
!***TYPE      DOUBLE PRECISION (SPREPJ-S, DPREPJ-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DPREPJ is called by DSTODE to compute and process the matrix
!  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
!  Here J is computed by the user-supplied routine JAC if
!  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
!  If MITER = 3, a diagonal approximation to J is used.
!  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
!  subjected to LU decomposition in preparation for later solution
!  of linear systems with P as coefficient matrix.  This is done
!  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!
!  In addition to variables described in DSTODE and DLSODE prologues,
!  communication with DPREPJ uses the following:
!  Y     = array containing predicted values on entry.
!  FTEM  = work array of length N (ACOR in DSTODE).
!  SAVF  = array containing f evaluated at predicted y.
!  WM    = real work space for matrices.  On output it contains the
!          inverse diagonal matrix if MITER = 3 and the LU decomposition
!          of P if MITER is 1, 2 , 4, or 5.
!          Storage of matrix elements starts at WM(3).
!          WM also contains the following matrix-related data:
!          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
!          WM(2) = H*EL0, saved for later use if MITER = 3.
!  IWM   = integer work space containing pivot information, starting at
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!  EL0   = EL(1) (input).
!  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!          P matrix found to be singular.
!  JCUR  = output flag = 1 to indicate that the Jacobian matrix
!          (or approximation) is now current.
!  This routine also uses the COMMON variables EL0, H, TN, UROUND,
!  MITER, N, NFE, and NJE.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DGBFA, DGEFA, DVNORM
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890504  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DPREPJ
!**End
EXTERNAL F, JAC
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR
INTEGER LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6),&
ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,&
MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP, MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
DOUBLE PRECISION CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ, DVNORM
!
!***FIRST EXECUTABLE STATEMENT  DPREPJ
NJE = NJE + 1
IERPJ = 0
JCUR = 1
HL0 = H*EL0
GOTO (100, 200, 300, 400, 500), MITER
! If MITER = 1, call JAC and multiply by scalar. -----------------------
100  LENP = N*N
DO 110 I = 1,LENP
110    WM(I+2) = 0.0D0
CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
CON = -HL0
DO 120 I = 1,LENP
120    WM(I+2) = WM(I+2)*CON
GOTO 240
! If MITER = 2, make N calls to F to approximate J. --------------------
200  FAC = DVNORM (N, SAVF, EWT)
R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
SRUR = WM(1)
J1 = 2
DO 230 J = 1,N
YJ = Y(J)
R = MAX(SRUR*ABS(YJ),R0/EWT(J))
Y(J) = Y(J) + R
FAC = -HL0/R
CALL F (NEQ, TN, Y, FTEM)
DO 220 I = 1,N
220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
Y(J) = YJ
J1 = J1 + N
230    CONTINUE
NFE = NFE + N
! Add identity matrix. -------------------------------------------------
240  J = 3
NP1 = N + 1
DO 250 I = 1,N
WM(J) = WM(J) + 1.0D0
250    J = J + NP1
! Do LU decomposition on P. --------------------------------------------
CALL DGEFA (WM(3), N, N, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
300  WM(2) = HL0
R = EL0*0.1D0
DO 310 I = 1,N
310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
CALL F (NEQ, TN, Y, WM(3))
NFE = NFE + 1
DO 320 I = 1,N
R0 = H*SAVF(I) - YH(I,2)
DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
WM(I+2) = 1.0D0
IF (ABS(R0) .LT. UROUND/EWT(I)) GOTO 320
IF (ABS(DI) .EQ. 0.0D0) GOTO 330
WM(I+2) = 0.1D0*R0/DI
320    CONTINUE
RETURN
330  IERPJ = 1
RETURN
! If MITER = 4, call JAC and multiply by scalar. -----------------------
400  ML = IWM(1)
MU = IWM(2)
ML3 = ML + 3
MBAND = ML + MU + 1
MEBAND = MBAND + ML
LENP = MEBAND*N
DO 410 I = 1,LENP
410    WM(I+2) = 0.0D0
CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
CON = -HL0
DO 420 I = 1,LENP
420    WM(I+2) = WM(I+2)*CON
GOTO 570
! If MITER = 5, make MBAND calls to F to approximate J. ----------------
500  ML = IWM(1)
MU = IWM(2)
MBAND = ML + MU + 1
MBA = MIN(MBAND,N)
MEBAND = MBAND + ML
MEB1 = MEBAND - 1
SRUR = WM(1)
FAC = DVNORM (N, SAVF, EWT)
R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
DO 560 J = 1,MBA
DO 530 I = J,N,MBAND
YI = Y(I)
R = MAX(SRUR*ABS(YI),R0/EWT(I))
530      Y(I) = Y(I) + R
CALL F (NEQ, TN, Y, FTEM)
DO 550 JJ = J,N,MBAND
Y(JJ) = YH(JJ,1)
YJJ = Y(JJ)
R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
FAC = -HL0/R
I1 = MAX(JJ-MU,1)
I2 = MIN(JJ+ML,N)
II = JJ*MEB1 - ML + 2
DO 540 I = I1,I2
540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
550      CONTINUE
560    CONTINUE
NFE = NFE + MBA
! Add identity matrix. -------------------------------------------------
570  II = MBAND + 2
DO 580 I = 1,N
WM(II) = WM(II) + 1.0D0
580    II = II + MEBAND
! Do LU decomposition of P. --------------------------------------------
CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
!----------------------- END OF SUBROUTINE DPREPJ ----------------------
END

SUBROUTINE DSOLSY (WM, IWM, X, TEM)
!***BEGIN PROLOGUE  DSOLSY
!***SUBSIDIARY
!***PURPOSE  ODEPACK linear system solver.
!***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This routine manages the solution of the linear system arising from
!  a chord iteration.  It is called if MITER .ne. 0.
!  If MITER is 1 or 2, it calls DGESL to accomplish this.
!  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
!  matrix, and then computes the solution.
!  If MITER is 4 or 5, it calls DGBSL.
!  Communication with DSOLSY uses the following variables:
!  WM    = real work space containing the inverse diagonal matrix if
!          MITER = 3 and the LU decomposition of the matrix otherwise.
!          Storage of matrix elements starts at WM(3).
!          WM also contains the following matrix-related data:
!          WM(1) = SQRT(UROUND) (not used here),
!          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
!  IWM   = integer work space containing pivot information, starting at
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!  X     = the right-hand side vector on input, and the solution vector
!          on output, of length N.
!  TEM   = vector of work space of length N, not used in this version.
!  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
!          IERSL = 1 if a singular matrix arose with MITER = 3.
!  This routine also uses the COMMON variables EL0, H, MITER, and N.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DGBSL, DGESL
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSOLSY
!**End
INTEGER IWM
DOUBLE PRECISION WM, X, TEM
DIMENSION WM(*), IWM(*), X(*), TEM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, MEBAND, ML, MU
DOUBLE PRECISION DI, HL0, PHL0, R
!
!***FIRST EXECUTABLE STATEMENT  DSOLSY
IERSL = 0
GOTO (100, 100, 300, 400, 400), MITER
100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
RETURN
!
300  PHL0 = WM(2)
HL0 = H*EL0
WM(2) = HL0
IF (HL0 .EQ. PHL0) GOTO 330
R = HL0/PHL0
DO 320 I = 1,N
DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
IF (ABS(DI) .EQ. 0.0D0) GOTO 390
320    WM(I+2) = 1.0D0/DI
330  DO 340 I = 1,N
340    X(I) = WM(I+2)*X(I)
RETURN
390  IERSL = 1
RETURN
!
400  ML = IWM(1)
MU = IWM(2)
MEBAND = 2*ML + MU + 1
CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
RETURN
!----------------------- END OF SUBROUTINE DSOLSY ----------------------
END

SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
!***BEGIN PROLOGUE  DSRCOM
!***SUBSIDIARY
!***PURPOSE  Save/restore ODEPACK COMMON blocks.
!***TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This routine saves or restores (depending on JOB) the contents of
!  the COMMON block DLS001, which is used internally
!  by one or more ODEPACK solvers.
!
!  RSAV = real array of length 218 or more.
!  ISAV = integer array of length 37 or more.
!  JOB  = flag indicating to save or restore the COMMON blocks:
!         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
!         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
!         A call with JOB = 2 presumes a prior call with JOB = 1.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   921116  Deleted treatment of block /EH0001/.  (ACH)
!   930801  Reduced Common block length by 2.  (ACH)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced Common block length by 209+12. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   031112  Added SAVE statement for data-loaded constants.
!***END PROLOGUE  DSRCOM
!**End
INTEGER ISAV, JOB
INTEGER ILS
INTEGER I, LENILS, LENRLS
DOUBLE PRECISION RSAV,   RLS
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS
COMMON /DLS001/ RLS(218), ILS(37)
DATA LENRLS/218/, LENILS/37/
!
!***FIRST EXECUTABLE STATEMENT  DSRCOM
IF (JOB .EQ. 2) GOTO 100
!
DO 10 I = 1,LENRLS
10     RSAV(I) = RLS(I)
DO 20 I = 1,LENILS
20     ISAV(I) = ILS(I)
RETURN
!
100  CONTINUE
DO 110 I = 1,LENRLS
110     RLS(I) = RSAV(I)
DO 120 I = 1,LENILS
120     ILS(I) = ISAV(I)
RETURN
!----------------------- END OF SUBROUTINE DSRCOM ----------------------
END

SUBROUTINE DSTODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM, F, JAC, PJAC, SLVS)
!***BEGIN PROLOGUE  DSTODE
!***SUBSIDIARY
!***PURPOSE  Performs one step of an ODEPACK integration.
!***TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DSTODE performs one step of the integration of an initial value
!  problem for a system of ordinary differential equations.
!  Note:  DSTODE is independent of the value of the iteration method
!  indicator MITER, when this is .ne. 0, and hence is independent
!  of the type of chord method used, or the Jacobian structure.
!  Communication with DSTODE is done with the following variables:
!
!  NEQ    = integer array containing problem size in NEQ(1), and
!           passed as the NEQ argument in all calls to F and JAC.
!  Y      = an array of length .ge. N used as the Y argument in
!           all calls to F and JAC.
!  YH     = an NYH by LMAX array containing the dependent variables
!           and their approximate scaled derivatives, where
!           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!           j-th derivative of y(i), scaled by h**j/factorial(j)
!           (j = 0,1,...,NQ).  on entry for the first step, the first
!           two columns of YH must be set from the initial values.
!  NYH    = a constant integer .ge. N, the first dimension of YH.
!  YH1    = a one-dimensional array occupying the same space as YH.
!  EWT    = an array of length N containing multiplicative weights
!           for local error measurements.  Local errors in Y(i) are
!           compared to 1.0/EWT(i) in various error tests.
!  SAVF   = an array of working storage, of length N.
!           Also used for input of YH(*,MAXORD+2) when JSTART = -1
!           and MAXORD .lt. the current order NQ.
!  ACOR   = a work array of length N, used for the accumulated
!           corrections.  On a successful return, ACOR(i) contains
!           the estimated one-step local error in Y(i).
!  WM,IWM = real and integer work arrays associated with matrix
!           operations in chord iteration (MITER .ne. 0).
!  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
!           and P = I - h*el0*JAC, if a chord method is being used.
!  SLVS   = name of routine to solve linear system in chord iteration.
!  CCMAX  = maximum relative change in h*el0 before PJAC is called.
!  H      = the step size to be attempted on the next step.
!           H is altered by the error control algorithm during the
!           problem.  H can be either positive or negative, but its
!           sign must remain constant throughout the problem.
!  HMIN   = the minimum absolute value of the step size h to be used.
!  HMXI   = inverse of the maximum absolute value of h to be used.
!           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
!           HMIN and HMXI may be changed at any time, but will not
!           take effect until the next change of h is considered.
!  TN     = the independent variable. TN is updated on each step taken.
!  JSTART = an integer used for input only, with the following
!           values and meanings:
!                0  perform the first step.
!            .gt.0  take a new step continuing from the last.
!               -1  take the next step with a new value of H, MAXORD,
!                     N, METH, MITER, and/or matrix parameters.
!               -2  take the next step with a new value of H,
!                     but with other inputs unchanged.
!           On return, JSTART is set to 1 to facilitate continuation.
!  KFLAG  = a completion code with the following meanings:
!                0  the step was succesful.
!               -1  the requested error could not be achieved.
!               -2  corrector convergence could not be achieved.
!               -3  fatal error in PJAC or SLVS.
!           A return with KFLAG = -1 or -2 means either
!           abs(H) = HMIN or 10 consecutive failures occurred.
!           On a return with KFLAG negative, the values of TN and
!           the YH array are as of the beginning of the last
!           step, and H is the last step size attempted.
!  MAXORD = the maximum order of integration method to be allowed.
!  MAXCOR = the maximum number of corrector iterations allowed.
!  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
!  MXNCF  = maximum number of convergence failures allowed.
!  METH/MITER = the method flags.  See description in driver.
!  N      = the number of first-order differential equations.
!  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
!  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DCFODE, DVNORM
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSTODE
!**End
EXTERNAL F, JAC, PJAC, SLVS
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L
INTEGER LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM

COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU

!
!***FIRST EXECUTABLE STATEMENT  DSTODE
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set to 2
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
MEO = METH
NSLP = 0
IPUP = MITER
IRET = 3
GOTO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MEO) GOTO 110
CALL DCFODE (METH, ELCO, TESCO)
MEO = METH
IF (NQ .GT. MAXORD) GOTO 120
IALTH = L
IRET = 1
GOTO 150
110  IF (NQ .LE. MAXORD) GOTO 160
120  NQ = MAXORD
L = LMAX
DO 125 I = 1,L
125    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
EXDN = 1.0D0/L
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
RH = MIN(RHDN,1.0D0)
IREDO = 3
IF (H .EQ. HOLD) GOTO 170
RH = MIN(RH,ABS(H/HOLD))
H = HOLD
GOTO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140  CALL DCFODE (METH, ELCO, TESCO)
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal Triangle matrix.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called, if a Jacobian is involved.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the R.M.S. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220  M = 0
DO 230 I = 1,N
230    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
IF (IPUP .LE. 0) GOTO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*el(1)*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
IPUP = 0
RC = 1.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .NE. 0) GOTO 430
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
270  IF (MITER .NE. 0) GOTO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
DO 290 I = 1,N
SAVF(I) = H*SAVF(I) - YH(I,2)
290    Y(I) = SAVF(I) - ACOR(I)
DEL = DVNORM (N, Y, EWT)
DO 300 I = 1,N
Y(I) = YH(I,1) + EL(1)*SAVF(I)
300    ACOR(I) = SAVF(I)
GOTO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
350  DO 360 I = 1,N
360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
CALL SLVS (WM, IWM, Y, SAVF)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
DEL = DVNORM (N, Y, EWT)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + Y(I)
380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M.gt.0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
IF (DCON .LE. 1.0D0) GOTO 450
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
DELP = DEL
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
GOTO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GOTO 430
ICF = 1
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
RMAX = 2.0D0
TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 680
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 670
IF (NCF .EQ. MXNCF) GOTO 670
RH = 0.25D0
IPUP = MITER
IREDO = 1
GOTO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450  JCUR = 0
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
HU = H
NQU = NQ
DO 470 J = 1,L
DO 470 I = 1,N
470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -3) GOTO 640
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 560
DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GOTO 610
IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, l, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640  IF (KFLAG .EQ. -10) GOTO 660
RH = 0.1D0
RH = MAX(HMIN/ABS(H),RH)
H = H*RH
DO 645 I = 1,N
645    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
DO 650 I = 1,N
650    YH(I,2) = H*SAVF(I)
IPUP = MITER
IALTH = 5
IF (NQ .EQ. 1) GOTO 200
NQ = 1
L = 2
IRET = 3
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -3
GOTO 720
690  RMAX = 10.0D0
700  R = 1.0D0/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- END OF SUBROUTINE DSTODE ----------------------
END

SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
!***BEGIN PROLOGUE  DEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DEWSET
!**End
INTEGER N, ITOL
INTEGER I
DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
!
!***FIRST EXECUTABLE STATEMENT  DEWSET
GOTO (10, 20, 30, 40), ITOL
10   CONTINUE
DO 15 I = 1,N
15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
RETURN
20   CONTINUE
DO 25 I = 1,N
25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
RETURN
30   CONTINUE
DO 35 I = 1,N
35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
RETURN
40   CONTINUE
DO 45 I = 1,N
45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
RETURN
!----------------------- END OF SUBROUTINE DEWSET ----------------------
END

DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
!***BEGIN PROLOGUE  DVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the array V, with weights
!  contained in the array W of length N:
!    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DVNORM
!**End
INTEGER N,   I
DOUBLE PRECISION V, W,   SUM
DIMENSION V(N), W(N)
!
!***FIRST EXECUTABLE STATEMENT  DVNORM
SUM = 0.0D0
DO 10 I = 1,N
10     SUM = SUM + (V(I)*W(I))**2
DVNORM = SQRT(SUM/N)
RETURN
!----------------------- END OF FUNCTION DVNORM ------------------------
END

SUBROUTINE DIPREP (NEQ, Y, RWORK, IA, JA, IPFLAG, F, JAC)
EXTERNAL F, JAC
INTEGER NEQ, IA, JA, IPFLAG
DOUBLE PRECISION Y, RWORK
DIMENSION NEQ(*), Y(*), RWORK(*), IA(*), JA(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMAX, LEWTN, LYHD, LYHN
!-----------------------------------------------------------------------
! This routine serves as an interface between the driver and
! Subroutine DPREP.  It is called only if MITER is 1 or 2.
! Tasks performed here are:
!  * call DPREP,
!  * reset the required WM segment length LENWK,
!  * move YH back to its final location (following WM in RWORK),
!  * reset pointers for YH, SAVF, EWT, and ACOR, and
!  * move EWT to its new position if ISTATE = 1.
! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
! no trouble, and IPFLAG is the value of the DPREP error flag IPPER
! if there was trouble in Subroutine DPREP.
!-----------------------------------------------------------------------
IPFLAG = 0
! Call DPREP to do matrix preprocessing operations. --------------------
CALL DPREP (NEQ, Y, RWORK(LYH), RWORK(LSAVF), RWORK(LEWT), RWORK(LACOR), IA, JA, RWORK(LWM), RWORK(LWM), IPFLAG, F, JAC)
LENWK = MAX(LREQ,LWMIN)
IF (IPFLAG .LT. 0) RETURN
! If DPREP was successful, move YH to end of required space for WM. ----
LYHN = LWM + LENWK
IF (LYHN .GT. LYH) RETURN
LYHD = LYH - LYHN
IF (LYHD .EQ. 0) GOTO 20
IMAX = LYHN - 1 + LENYHM
DO 10 I = LYHN,IMAX
10     RWORK(I) = RWORK(I+LYHD)
LYH = LYHN
! Reset pointers for SAVF, EWT, and ACOR. ------------------------------
20   LSAVF = LYH + LENYH
LEWTN = LSAVF + N
LACOR = LEWTN + N
IF (ISTATC .EQ. 3) GOTO 40
! If ISTATE = 1, move EWT (left) to its new position. ------------------
IF (LEWTN .GT. LEWT) RETURN
DO 30 I = 1,N
30     RWORK(I+LEWTN-1) = RWORK(I+LEWT-1)
40   LEWT = LEWTN
RETURN
!----------------------- End of Subroutine DIPREP ----------------------
END

SUBROUTINE DPREP (NEQ, Y, YH, SAVF, EWT, FTEM, IA, JA, WK, IWK, IPPER, F, JAC)
EXTERNAL F,JAC
INTEGER NEQ, IA, JA, IWK, IPPER
DOUBLE PRECISION Y, YH, SAVF, EWT, FTEM, WK
DIMENSION NEQ(*), Y(*), YH(*), SAVF(*), EWT(*), FTEM(*), IA(*), JA(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN,&
IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, &
NLU, NNZ, NSP, NZL, NZU
INTEGER I, IBR, IER, IPIL, IPIU, IPTT1, IPTT2, J, JFOUND, K, KNEW, KMAX, KMIN, LDIF, LENIGP, LIWK, MAXG, NP1, NZSUT
DOUBLE PRECISION DQ, DYJ, ERWT, FAC, YJ
!-----------------------------------------------------------------------
! This routine performs preprocessing related to the sparse linear
! systems that must be solved if MITER = 1 or 2.
! The operations that are performed here are:
!  * compute sparseness structure of Jacobian according to MOSS,
!  * compute grouping of column indices (MITER = 2),
!  * compute a new ordering of rows and columns of the matrix,
!  * reorder JA corresponding to the new ordering,
!  * perform a symbolic LU factorization of the matrix, and
!  * set pointers for segments of the IWK/WK array.
! In addition to variables described previously, DPREP uses the
! following for communication:
! YH     = the history array.  Only the first column, containing the
!          current Y vector, is used.  Used only if MOSS .ne. 0.
! SAVF   = a work array of length NEQ, used only if MOSS .ne. 0.
! EWT    = array of length NEQ containing (inverted) error weights.
!          Used only if MOSS = 2 or if ISTATE = MOSS = 1.
! FTEM   = a work array of length NEQ, identical to ACOR in the driver,
!          used only if MOSS = 2.
! WK     = a real work array of length LENWK, identical to WM in
!          the driver.
! IWK    = integer work array, assumed to occupy the same space as WK.
! LENWK  = the length of the work arrays WK and IWK.
! ISTATC = a copy of the driver input argument ISTATE (= 1 on the
!          first call, = 3 on a continuation call).
! IYS    = flag value from ODRV or CDRV.
! IPPER  = output error flag with the following values and meanings:
!          0  no error.
!         -1  insufficient storage for internal structure pointers.
!         -2  insufficient storage for JGROUP.
!         -3  insufficient storage for ODRV.
!         -4  other error flag from ODRV (should never occur).
!         -5  insufficient storage for CDRV.
!         -6  other error flag from CDRV.
!-----------------------------------------------------------------------
IBIAN = LRAT*2
IPIAN = IBIAN + 1
NP1 = N + 1
IPJAN = IPIAN + NP1
IBJAN = IPJAN - 1
LIWK = LENWK*LRAT
IF (IPJAN+N-1 .GT. LIWK) GOTO 210
IF (MOSS .EQ. 0) GOTO 30
!
IF (ISTATC .EQ. 3) GOTO 20
! ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
DO 10 I = 1,N
ERWT = 1.0D0/EWT(I)
FAC = 1.0D0 + 1.0D0/(I + 1.0D0)
Y(I) = Y(I) + FAC*SIGN(ERWT,Y(I))
10     CONTINUE
GOTO (70, 100), MOSS
!
20   CONTINUE
! ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
DO 25 I = 1,N
25     Y(I) = YH(I)
GOTO (70, 100), MOSS
!
! MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
30   KNEW = IPJAN
KMIN = IA(1)
IWK(IPIAN) = 1
DO 60 J = 1,N
JFOUND = 0
KMAX = IA(J+1) - 1
IF (KMIN .GT. KMAX) GOTO 45
DO 40 K = KMIN,KMAX
I = JA(K)
IF (I .EQ. J) JFOUND = 1
IF (KNEW .GT. LIWK) GOTO 210
IWK(KNEW) = I
KNEW = KNEW + 1
40       CONTINUE
IF (JFOUND .EQ. 1) GOTO 50
45     IF (KNEW .GT. LIWK) GOTO 210
IWK(KNEW) = J
KNEW = KNEW + 1
50     IWK(IPIAN+J) = KNEW + 1 - IPJAN
KMIN = KMAX + 1
60     CONTINUE
GOTO 140
!
! MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
70   CONTINUE
! A dummy call to F allows user to create temporaries for use in JAC. --
CALL F (NEQ, TN, Y, SAVF)
K = IPJAN
IWK(IPIAN) = 1
DO 90 J = 1,N
IF (K .GT. LIWK) GOTO 210
IWK(K) = J
K = K + 1
DO 75 I = 1,N
75       SAVF(I) = 0.0D0
CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), SAVF)
DO 80 I = 1,N
IF (ABS(SAVF(I)) .LE. SETH) GOTO 80
IF (I .EQ. J) GOTO 80
IF (K .GT. LIWK) GOTO 210
IWK(K) = I
K = K + 1
80       CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
90     CONTINUE
GOTO 140
!
! MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
100  K = IPJAN
IWK(IPIAN) = 1
CALL F (NEQ, TN, Y, SAVF)
DO 120 J = 1,N
IF (K .GT. LIWK) GOTO 210
IWK(K) = J
K = K + 1
YJ = Y(J)
ERWT = 1.0D0/EWT(J)
DYJ = SIGN(ERWT,YJ)
Y(J) = YJ + DYJ
CALL F (NEQ, TN, Y, FTEM)
Y(J) = YJ
DO 110 I = 1,N
DQ = (FTEM(I) - SAVF(I))/DYJ
IF (ABS(DQ) .LE. SETH) GOTO 110
IF (I .EQ. J) GOTO 110
IF (K .GT. LIWK) GOTO 210
IWK(K) = I
K = K + 1
110      CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
120    CONTINUE
!
140  CONTINUE
IF (MOSS .EQ. 0 .OR. ISTATC .NE. 1) GOTO 150
! If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
DO 145 I = 1,N
145    Y(I) = YH(I)
150  NNZ = IWK(IPIAN+N) - 1
LENIGP = 0
IPIGP = IPJAN + NNZ
IF (MITER .NE. 2) GOTO 160
!
! Compute grouping of column indices (MITER = 2). ----------------------
MAXG = NP1
IPJGP = IPJAN + NNZ
IBJGP = IPJGP - 1
IPIGP = IPJGP + N
IPTT1 = IPIGP + NP1
IPTT2 = IPTT1 + N
LREQ = IPTT2 + N - 1
IF (LREQ .GT. LIWK) GOTO 220
CALL JGROUP (N, IWK(IPIAN), IWK(IPJAN), MAXG, NGP, IWK(IPIGP), IWK(IPJGP), IWK(IPTT1), IWK(IPTT2), IER)
IF (IER .NE. 0) GOTO 220
LENIGP = NGP + 1
!
! Compute new ordering of rows/columns of Jacobian. --------------------
160  IPR = IPIGP + LENIGP
IPC = IPR
IPIC = IPC + N
IPISP = IPIC + N
IPRSP = (IPISP - 2)/LRAT + 2
IESP = LENWK + 1 - IPRSP
IF (IESP .LT. 0) GOTO 230
IBR = IPR - 1
DO 170 I = 1,N
170    IWK(IBR+I) = I
NSP = LIWK + 1 - IPISP
CALL ODRV (N, IWK(IPIAN), IWK(IPJAN), WK, IWK(IPR), IWK(IPIC), NSP, IWK(IPISP), 1, IYS)
IF (IYS .EQ. 11*N+1) GOTO 240
IF (IYS .NE. 0) GOTO 230
!
! Reorder JAN and do symbolic LU factorization of matrix. --------------
IPA = LENWK + 1 - NNZ
NSP = IPA - IPRSP
LREQ = MAX(12*N/LRAT, 6*N/LRAT+2*N+NNZ) + 3
LREQ = LREQ + IPRSP - 1 + NNZ
IF (LREQ .GT. LENWK) GOTO 250
IBA = IPA - 1
DO 180 I = 1,NNZ
180    WK(IBA+I) = 0.0D0
IPISP = LRAT*(IPRSP - 1) + 1
CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)
LREQ = LENWK - IESP
IF (IYS .EQ. 10*N+1) GOTO 250
IF (IYS .NE. 0) GOTO 260
IPIL = IPISP
IPIU = IPIL + 2*N + 1
NZU = IWK(IPIL+N) - IWK(IPIL)
NZL = IWK(IPIU+N) - IWK(IPIU)
IF (LRAT .GT. 1) GOTO 190
CALL ADJLR (N, IWK(IPISP), LDIF)
LREQ = LREQ + LDIF
190  CONTINUE
IF (LRAT .EQ. 2 .AND. NNZ .EQ. N) LREQ = LREQ + 1
NSP = NSP + LREQ - LENWK
IPA = LREQ + 1 - NNZ
IBA = IPA - 1
IPPER = 0
RETURN
!
210  IPPER = -1
LREQ = 2 + (2*N + 1)/LRAT
LREQ = MAX(LENWK+1,LREQ)
RETURN
!
220  IPPER = -2
LREQ = (LREQ - 1)/LRAT + 1
RETURN
!
230  IPPER = -3
CALL CNTNZU (N, IWK(IPIAN), IWK(IPJAN), NZSUT)
LREQ = LENWK - IESP + (3*N + 4*NZSUT - 1)/LRAT + 1
RETURN
!
240  IPPER = -4
RETURN
!
250  IPPER = -5
RETURN
!
260  IPPER = -6
LREQ = LENWK
RETURN
!----------------------- End of Subroutine DPREP -----------------------
END

SUBROUTINE JGROUP (N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
INTEGER N, IA, JA, MAXG, NGRP, IGP, JGP, INCL, JDONE, IER
DIMENSION IA(*), JA(*), IGP(*), JGP(*), INCL(*), JDONE(*)
!-----------------------------------------------------------------------
! This subroutine constructs groupings of the column indices of
! the Jacobian matrix, used in the numerical evaluation of the
! Jacobian by finite differences.
!
! Input:
! N      = the order of the matrix.
! IA,JA  = sparse structure descriptors of the matrix by rows.
! MAXG   = length of available storage in the IGP array.
!
! Output:
! NGRP   = number of groups.
! JGP    = array of length N containing the column indices by groups.
! IGP    = pointer array of length NGRP + 1 to the locations in JGP
!          of the beginning of each group.
! IER    = error indicator.  IER = 0 if no error occurred, or 1 if
!          MAXG was insufficient.
!
! INCL and JDONE are working arrays of length N.
!-----------------------------------------------------------------------
INTEGER I, J, K, KMIN, KMAX, NCOL, NG
!
IER = 0
DO 10 J = 1,N
10     JDONE(J) = 0
NCOL = 1
DO 60 NG = 1,MAXG
IGP(NG) = NCOL
DO 20 I = 1,N
20       INCL(I) = 0
DO 50 J = 1,N
! Reject column J if it is already in a group.--------------------------
IF (JDONE(J) .EQ. 1) GOTO 50
KMIN = IA(J)
KMAX = IA(J+1) - 1
DO 30 K = KMIN,KMAX
! Reject column J if it overlaps any column already in this group.------
I = JA(K)
IF (INCL(I) .EQ. 1) GOTO 50
30         CONTINUE
! Accept column J into group NG.----------------------------------------
JGP(NCOL) = J
NCOL = NCOL + 1
JDONE(J) = 1
DO 40 K = KMIN,KMAX
I = JA(K)
40         INCL(I) = 1
50       CONTINUE
! Stop if this group is empty (grouping is complete).-------------------
IF (NCOL .EQ. IGP(NG)) GOTO 70
60     CONTINUE
! Error return if not all columns were chosen (MAXG too small).---------
IF (NCOL .LE. N) GOTO 80
NG = MAXG
70   NGRP = NG - 1
RETURN
80   IER = 1
RETURN
!----------------------- End of Subroutine JGROUP ----------------------
END

SUBROUTINE ADJLR (N, ISP, LDIF)
INTEGER N, ISP, LDIF
DIMENSION ISP(*)
!-----------------------------------------------------------------------
! This routine computes an adjustment, LDIF, to the required
! integer storage space in IWK (sparse matrix work space).
! It is called only if the word length ratio is LRAT = 1.
! This is to account for the possibility that the symbolic LU phase
! may require more storage than the numerical LU and solution phases.
!-----------------------------------------------------------------------
INTEGER IP, JLMAX, JUMAX, LNFC, LSFC, NZLU
!
IP = 2*N + 1
! Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
JLMAX = ISP(IP)
JUMAX = ISP(IP+IP)
! NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
NZLU = ISP(N+1) - ISP(1) + ISP(IP+N+1) - ISP(IP+1)
LSFC = 12*N + 3 + 2*MAX(JLMAX,JUMAX)
LNFC = 9*N + 2 + JLMAX + JUMAX + NZLU
LDIF = MAX(0, LSFC - LNFC)
RETURN
!----------------------- End of Subroutine ADJLR -----------------------
END

SUBROUTINE CNTNZU (N, IA, JA, NZSUT)
INTEGER N, IA, JA, NZSUT
DIMENSION IA(*), JA(*)
!-----------------------------------------------------------------------
! This routine counts the number of nonzero elements in the strict
! upper triangle of the matrix M + M(transpose), where the sparsity
! structure of M is given by pointer arrays IA and JA.
! This is needed to compute the storage requirements for the
! sparse matrix reordering operation in ODRV.
!-----------------------------------------------------------------------
INTEGER II, JJ, J, JMIN, JMAX, K, KMIN, KMAX, NUM
!
NUM = 0
DO 50 II = 1,N
JMIN = IA(II)
JMAX = IA(II+1) - 1
IF (JMIN .GT. JMAX) GOTO 50
DO 40 J = JMIN,JMAX
IF (JA(J) - II) 10, 40, 30
10       JJ =JA(J)
KMIN = IA(JJ)
KMAX = IA(JJ+1) - 1
IF (KMIN .GT. KMAX) GOTO 30
DO 20 K = KMIN,KMAX
IF (JA(K) .EQ. II) GOTO 40
20         CONTINUE
30       NUM = NUM + 1
40       CONTINUE
50     CONTINUE
NZSUT = NUM
RETURN
!----------------------- End of Subroutine CNTNZU ----------------------
END

SUBROUTINE DPRJS (NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WK,IWK,F,JAC)
EXTERNAL F,JAC
INTEGER NEQ, NYH, IWK
DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WK
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, &
MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN,&
IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, &
NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMUL, J, JJ, JOK, JMAX, JMIN, K, KMAX, KMIN, NG
DOUBLE PRECISION CON, DI, FAC, HL0, PIJ, R, R0, RCON, RCONT, SRUR, DVNORM
!-----------------------------------------------------------------------
! DPRJS is called to compute and process the matrix
! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
! J is computed by columns, either by the user-supplied routine JAC
! if MITER = 1, or by finite differencing if MITER = 2.
! if MITER = 3, a diagonal approximation to J is used.
! if MITER = 1 or 2, and if the existing value of the Jacobian
! (as contained in P) is considered acceptable, then a new value of
! P is reconstructed from the old value.  In any case, when MITER
! is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
! P and its LU decomposition are stored (separately) in WK.
!
! In addition to variables described previously, communication
! with DPRJS uses the following:
! Y     = array containing predicted values on entry.
! FTEM  = work array of length N (ACOR in DSTODE).
! SAVF  = array containing f evaluated at predicted y.
! WK    = real work space for matrices.  On output it contains the
!         inverse diagonal matrix if MITER = 3, and P and its sparse
!         LU decomposition if MITER is 1 or 2.
!         Storage of matrix elements starts at WK(3).
!         WK also contains the following matrix-related data:
!         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
!         WK(2) = H*EL0, saved for later use if MITER = 3.
! IWK   = integer work space for matrix-related data, assumed to
!         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!         are assumed to have identical locations.
! EL0   = EL(1) (input).
! IERPJ = output error flag (in Common).
!       = 0 if no error.
!       = 1  if zero pivot found in CDRV.
!       = 2  if a singular matrix arose with MITER = 3.
!       = -1 if insufficient storage for CDRV (should not occur here).
!       = -2 if other error found in CDRV (should not occur here).
! JCUR  = output flag showing status of (approximate) Jacobian matrix:
!          = 1 to indicate that the Jacobian is now current, or
!          = 0 to indicate that a saved value was used.
! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
HL0 = H*EL0
CON = -HL0
IF (MITER .EQ. 3) GOTO 300
! See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
JOK = 1
IF (NST .EQ. 0 .OR. NST .GE. NSLJ+MSBJ) JOK = 0
IF (ICF .EQ. 1 .AND. ABS(RC - 1.0D0) .LT. CCMXJ) JOK = 0
IF (ICF .EQ. 2) JOK = 0
IF (JOK .EQ. 1) GOTO 250
!
! MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
20   JCUR = 1
NJE = NJE + 1
NSLJ = NST
IPLOST = 0
CONMIN = ABS(CON)
GOTO (100, 200), MITER
!
! If MITER = 1, call JAC, multiply by scalar, and add identity. --------
100  CONTINUE
KMIN = IWK(IPIAN)
DO 130 J = 1, N
KMAX = IWK(IPIAN+J) - 1
DO 110 I = 1,N
110      FTEM(I) = 0.0D0
CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), FTEM)
DO 120 K = KMIN, KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = FTEM(I)*CON
IF (I .EQ. J) WK(IBA+K) = WK(IBA+K) + 1.0D0
120      CONTINUE
KMIN = KMAX + 1
130    CONTINUE
GOTO 290
!
! If MITER = 2, make NGP calls to F to approximate J and P. ------------
200  CONTINUE
FAC = DVNORM(N, SAVF, EWT)
R0 = 1000.0D0 * ABS(H) * UROUND * N * FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
SRUR = WK(1)
JMIN = IWK(IPIGP)
DO 240 NG = 1,NGP
JMAX = IWK(IPIGP+NG) - 1
DO 210 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
210      Y(JJ) = Y(JJ) + R
CALL F (NEQ, TN, Y, FTEM)
DO 230 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
Y(JJ) = YH(JJ,1)
R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
FAC = -HL0/R
KMIN =IWK(IBIAN+JJ)
KMAX =IWK(IBIAN+JJ+1) - 1
DO 220 K = KMIN,KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = (FTEM(I) - SAVF(I))*FAC
IF (I .EQ. JJ) WK(IBA+K) = WK(IBA+K) + 1.0D0
220        CONTINUE
230      CONTINUE
JMIN = JMAX + 1
240    CONTINUE
NFE = NFE + NGP
GOTO 290
!
! If JOK = 1, reconstruct new P from old P. ----------------------------
250  JCUR = 0
RCON = CON/CON0
RCONT = ABS(CON)/CONMIN
IF (RCONT .GT. RBIG .AND. IPLOST .EQ. 1) GOTO 20
KMIN = IWK(IPIAN)
DO 275 J = 1,N
KMAX = IWK(IPIAN+J) - 1
DO 270 K = KMIN,KMAX
I = IWK(IBJAN+K)
PIJ = WK(IBA+K)
IF (I .NE. J) GOTO 260
PIJ = PIJ - 1.0D0
IF (ABS(PIJ) .GE. PSMALL) GOTO 260
IPLOST = 1
CONMIN = MIN(ABS(CON0),CONMIN)
260      PIJ = PIJ*RCON
IF (I .EQ. J) PIJ = PIJ + 1.0D0
WK(IBA+K) = PIJ
270      CONTINUE
KMIN = KMAX + 1
275    CONTINUE
!
! Do numerical factorization of P matrix. ------------------------------
290  NLU = NLU + 1
CON0 = CON
IERPJ = 0
DO 295 I = 1,N
295    FTEM(I) = 0.0D0
CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),FTEM,FTEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
IF (IYS .EQ. 0) RETURN
IMUL = (IYS - 1)/N
IERPJ = -2
IF (IMUL .EQ. 8) IERPJ = 1
IF (IMUL .EQ. 10) IERPJ = -1
RETURN
!
! If MITER = 3, construct a diagonal approximation to J and P. ---------
300  CONTINUE
JCUR = 1
NJE = NJE + 1
WK(2) = HL0
IERPJ = 0
R = EL0*0.1D0
DO 310 I = 1,N
310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
CALL F (NEQ, TN, Y, WK(3))
NFE = NFE + 1
DO 320 I = 1,N
R0 = H*SAVF(I) - YH(I,2)
DI = 0.1D0*R0 - H*(WK(I+2) - SAVF(I))
WK(I+2) = 1.0D0
IF (ABS(R0) .LT. UROUND/EWT(I)) GOTO 320
IF (ABS(DI) .EQ. 0.0D0) GOTO 330
WK(I+2) = 0.1D0*R0/DI
320    CONTINUE
RETURN
330  IERPJ = 2
RETURN
!----------------------- End of Subroutine DPRJS -----------------------
END

SUBROUTINE DSOLSS (WK, IWK, X, TEM)
INTEGER IWK
DOUBLE PRECISION WK, X, TEM
DIMENSION WK(*), IWK(*), X(*), TEM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I
DOUBLE PRECISION DI, HL0, PHL0, R
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls CDRV to accomplish this.
! If MITER = 3 it updates the coefficient H*EL0 in the diagonal
! matrix, and then computes the solution.
! communication with DSOLSS uses the following variables:
! WK    = real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WK(3).
!         WK also contains the following matrix-related data:
!         WK(1) = SQRT(UROUND) (not used here),
!         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
! IWK   = integer work space for matrix-related data, assumed to
!         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!         are assumed to have identical locations.
! X     = the right-hand side vector on input, and the solution vector
!         on output, of length N.
! TEM   = vector of work space of length N, not used in this version.
! IERSL = output flag (in Common).
!         IERSL = 0  if no trouble occurred.
!         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
!                    This should never occur and is considered fatal.
!         IERSL = 1  if a singular matrix arose with MITER = 3.
! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
IERSL = 0
GOTO (100, 100, 300), MITER
100  CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),WK(IPA),X,X,NSP,IWK(IPISP),WK(IPRSP),IESP,4,IERSL)
IF (IERSL .NE. 0) IERSL = -1
RETURN
!
300  PHL0 = WK(2)
HL0 = H*EL0
WK(2) = HL0
IF (HL0 .EQ. PHL0) GOTO 330
R = HL0/PHL0
DO 320 I = 1,N
DI = 1.0D0 - R*(1.0D0 - 1.0D0/WK(I+2))
IF (ABS(DI) .EQ. 0.0D0) GOTO 390
320    WK(I+2) = 1.0D0/DI
330  DO 340 I = 1,N
340    X(I) = WK(I+2)*X(I)
RETURN
390  IERSL = 1
RETURN
!
!----------------------- End of Subroutine DSOLSS ----------------------
END

SUBROUTINE DSRCMS (RSAV, ISAV, JOB)
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of
! the Common blocks DLS001, DLSS01, which are used
! internally by one or more ODEPACK solvers.
!
! RSAV = real array of length 224 or more.
! ISAV = integer array of length 71 or more.
! JOB  = flag indicating to save or restore the Common blocks:
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
INTEGER ISAV, JOB
INTEGER ILS, ILSS
INTEGER I, LENILS, LENISS, LENRLS, LENRSS
DOUBLE PRECISION RSAV,   RLS, RLSS
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS, LENRSS, LENISS
COMMON /DLS001/ RLS(218), ILS(37)
COMMON /DLSS01/ RLSS(6), ILSS(34)
DATA LENRLS/218/, LENILS/37/, LENRSS/6/, LENISS/34/
!
IF (JOB .EQ. 2) GOTO 100
DO 10 I = 1,LENRLS
10     RSAV(I) = RLS(I)
DO 15 I = 1,LENRSS
15     RSAV(LENRLS+I) = RLSS(I)
!
DO 20 I = 1,LENILS
20     ISAV(I) = ILS(I)
DO 25 I = 1,LENISS
25     ISAV(LENILS+I) = ILSS(I)
!
RETURN
!
100  CONTINUE
DO 110 I = 1,LENRLS
110     RLS(I) = RSAV(I)
DO 115 I = 1,LENRSS
115     RLSS(I) = RSAV(LENRLS+I)
!
DO 120 I = 1,LENILS
120     ILS(I) = ISAV(I)
DO 125 I = 1,LENISS
125     ILSS(I) = ISAV(LENILS+I)
!
RETURN
!----------------------- End of Subroutine DSRCMS ----------------------
END

subroutine odrv(n, ia,ja,a, p,ip, nsp,isp, path, flag)
!                                                                 5/2/83
!***********************************************************************
!  odrv -- driver for sparse matrix reordering routines
!***********************************************************************
!
!  description
!
!    odrv finds a minimum degree ordering of the rows and columns
!    of a matrix m stored in (ia,ja,a) format (see below).  for the
!    reordered matrix, the work and storage required to perform
!    gaussian elimination is (usually) significantly less.
!
!    note.. odrv and its subordinate routines have been modified to
!    compute orderings for general matrices, not necessarily having any
!    symmetry.  the miminum degree ordering is computed for the
!    structure of the symmetric matrix  m + m-transpose.
!    modifications to the original odrv module have been made in
!    the coding in subroutine mdi, and in the initial comments in
!    subroutines odrv and md.
!
!    if only the nonzero entries in the upper triangle of m are being
!    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
!    with the diagonal entries placed first in each row.  this is to
!    ensure that if m(i,j) will be in the upper triangle of m with
!    respect to the new ordering, then m(i,j) is stored in row i (and
!    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
!    strict lower triangle of m, then m(j,i) is stored in row j (and
!    thus m(i,j) is not stored).
!
!
!  storage of sparse matrices
!
!    the nonzero entries of the matrix m are stored row-by-row in the
!    array a.  to identify the individual nonzero entries in each row,
!    we need to know in which column each entry lies.  these column
!    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  to identify the individual rows, we need to know where
!    each row starts.  these row pointers are stored in the array ia.
!    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
!    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
!    the first location following the last element in the last row.
!    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
!    the nonzero entries in the i-th row are stored consecutively in
!
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!
!    and the corresponding column indices are stored consecutively in
!
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!
!    when the coefficient matrix is symmetric, only the nonzero entries
!    in the upper triangle need be stored.  for example, the matrix
!
!             ( 1  0  2  3  0 )
!             ( 0  4  0  0  0 )
!         m = ( 2  0  5  6  0 )
!             ( 3  0  6  7  8 )
!             ( 0  0  0  8  9 )
!
!    could be stored as
!
!            - 1  2  3  4  5  6  7  8  9 10 11 12 13
!         ---+--------------------------------------
!         ia - 1  4  5  8 12 14
!         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
!          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
!
!    or (symmetrically) as
!
!            - 1  2  3  4  5  6  7  8  9
!         ---+--------------------------
!         ia - 1  4  5  7  9 10
!         ja - 1  3  4  2  3  4  4  5  5
!          a - 1  2  3  4  5  6  7  8  9          .
!
!
!  parameters
!
!    n    - order of the matrix
!
!    ia   - integer one-dimensional array containing pointers to delimit
!           rows in ja and a.  dimension = n+1
!
!    ja   - integer one-dimensional array containing the column indices
!           corresponding to the elements of a.  dimension = number of
!           nonzero entries in (the upper triangle of) m
!
!    a    - real one-dimensional array containing the nonzero entries in
!           (the upper triangle of) m, stored by rows.  dimension =
!           number of nonzero entries in (the upper triangle of) m
!
!    p    - integer one-dimensional array used to return the permutation
!           of the rows and columns of m corresponding to the minimum
!           degree ordering.  dimension = n
!
!    ip   - integer one-dimensional array used to return the inverse of
!           the permutation returned in p.  dimension = n
!
!    nsp  - declared dimension of the one-dimensional array isp.  nsp
!           must be at least  3n+4k,  where k is the number of nonzeroes
!           in the strict upper triangle of m
!
!    isp  - integer one-dimensional array used for working storage.
!           dimension = nsp
!
!    path - integer path specification.  values and their meanings are -
!             1  find minimum degree ordering only
!             2  find minimum degree ordering and reorder symmetrically
!                  stored matrix (used when only the nonzero entries in
!                  the upper triangle of m are being stored)
!             3  reorder symmetrically stored matrix as specified by
!                  input permutation (used when an ordering has already
!                  been determined and only the nonzero entries in the
!                  upper triangle of m are being stored)
!             4  same as 2 but put diagonal entries at start of each row
!             5  same as 3 but put diagonal entries at start of each row
!
!    flag - integer error flag.  values and their meanings are -
!               0    no errors detected
!              9n+k  insufficient storage in md
!             10n+1  insufficient storage in odrv
!             11n+1  illegal path specification
!
!
!  conversion from real to double precision
!
!    change the real declarations in odrv and sro to double precision
!    declarations.
!
!-----------------------------------------------------------------------
!
integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,v, l, head,  tmp, q, n, nsp, max, next
!...  real  a(*)
double precision  a(*)
logical  dflag
!
!----initialize error flag and validate path specification
flag = 0
if (path.lt.1 .or. 5.lt.path)  GOTO 111
!
!----allocate storage and find minimum degree ordering
if ((path-1) * (path-2) * (path-4) .ne. 0)  GOTO 1
max = (nsp-n)/2
v    = 1
l    = v     +  max
head = l     +  max
next = head  +  n
if (max.lt.n)  GOTO 110
!
call  md(n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
if (flag.ne.0)  GOTO 100
!
!----allocate storage and symmetrically reorder matrix
1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  GOTO 2
tmp = (nsp+1) -      n
q   = tmp     - (ia(n+1)-1)
if (q.lt.1)  GOTO 110
!
dflag = path.eq.4 .or. path.eq.5
call sro(n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
!
2  return
!
! ** error -- error detected in md
100  return
! ** error -- insufficient storage
110  flag = 10*n + 1
return
! ** error -- illegal path specified
111  flag = 11*n + 1
return
end
subroutine md(n, ia,ja, max, v,l, head,last,next, mark, flag)
!***********************************************************************
!  md -- minimum degree algorithm (based on element model)
!***********************************************************************
!
!  description
!
!    md finds a minimum degree ordering of the rows and columns of a
!    general sparse matrix m stored in (ia,ja,a) format.
!    when the structure of m is nonsymmetric, the ordering is that
!    obtained for the symmetric matrix  m + m-transpose.
!
!
!  additional parameters
!
!    max  - declared dimension of the one-dimensional arrays v and l.
!           max must be at least  n+2k,  where k is the number of
!           nonzeroes in the strict upper triangle of m + m-transpose
!
!    v    - integer one-dimensional work array.  dimension = max
!
!    l    - integer one-dimensional work array.  dimension = max
!
!    head - integer one-dimensional work array.  dimension = n
!
!    last - integer one-dimensional array used to return the permutation
!           of the rows and columns of m corresponding to the minimum
!           degree ordering.  dimension = n
!
!    next - integer one-dimensional array used to return the inverse of
!           the permutation returned in last.  dimension = n
!
!    mark - integer one-dimensional work array (may be the same as v).
!           dimension = n
!
!    flag - integer error flag.  values and their meanings are -
!             0     no errors detected
!             9n+k  insufficient storage in md
!
!
!  definitions of internal parameters
!
!    ---------+---------------------------------------------------------
!    v(s)     - value field of list entry
!    ---------+---------------------------------------------------------
!    l(s)     - link field of list entry  (0 =) end of list)
!    ---------+---------------------------------------------------------
!    l(vi)    - pointer to element list of uneliminated vertex vi
!    ---------+---------------------------------------------------------
!    l(ej)    - pointer to boundary list of active element ej
!    ---------+---------------------------------------------------------
!    head(d)  - vj =) vj head of d-list d
!             -  0 =) no vertex in d-list d
!
!
!             -                  vi uneliminated vertex
!             -          vi in ek           -       vi not in ek
!    ---------+-----------------------------+---------------------------
!    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
!             -                             -  0 =) vi tail of d-list
!    ---------+-----------------------------+---------------------------
!    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
!             --vk =) compute degree        - vj =) vj last in d-list
!             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
!             -  0 =) do not compute degree -
!    ---------+-----------------------------+---------------------------
!    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
!
!
!             -                   vi eliminated vertex
!             -      ei active element      -           otherwise
!    ---------+-----------------------------+---------------------------
!    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
!             -       to be eliminated      -       to be eliminated
!    ---------+-----------------------------+---------------------------
!    last(vi) -  m =) size of ei = m        - undefined
!    ---------+-----------------------------+---------------------------
!    mark(vi) - -m =) overlap count of ei   - undefined
!             -       with ek = m           -
!             - otherwise nonnegative tag   -
!             -       .lt. mark(vk)         -
!
!-----------------------------------------------------------------------
!
integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),mark(*),  flag,  tag, dmin, vk,ek, tail
integer n, max, k
equivalence  (vk,ek)
!
!----initialization
tag = 0
call  mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
if (flag.ne.0)  return
!
k = 0
dmin = 1
!
!----while  k .lt. n  do
1  if (k.ge.n)  GOTO 4
!
!------search for vertex of minimum degree
2    if (head(dmin).gt.0)  GOTO 3
dmin = dmin + 1
GOTO 2
!
!------remove vertex vk of minimum degree from degree list
3    vk = head(dmin)
head(dmin) = next(vk)
if (head(dmin).gt.0)  last(head(dmin)) = -dmin
!
!------number vertex vk, adjust tag, and tag vk
k = k+1
next(vk) = -k
last(ek) = dmin - 1
tag = tag + last(ek)
mark(vk) = tag
!
!------form element ek from uneliminated neighbors of vk
call  mdm(vk,tail, v,l, last,next, mark)
!
!------purge inactive elements and do mass elimination
call  mdp(k,ek,tail, v,l, head,last,next, mark)
!
!------update degrees of uneliminated vertices in ek
call  mdu(ek,dmin, v,l, head,last,next, mark)
!
GOTO 1
!
!----generate inverse permutation from permutation
4  do 5 k=1,n
next(k) = -next(k)
5    last(next(k)) = k
!
return
end
subroutine mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
!***********************************************************************
!  mdi -- initialization
!***********************************************************************
integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*), mark(*), tag,  flag,  sfs, vi,dvi, vj
integer n, max, j, jmax, jmin, k, kmax, lvk, nextvi
!
!----initialize degrees, element lists, and degree lists
do 1 vi=1,n
mark(vi) = 1
l(vi) = 0
1    head(vi) = 0
sfs = n+1
!
!----create nonzero structure
!----for each nonzero entry a(vi,vj)
do 6 vi=1,n
jmin = ia(vi)
jmax = ia(vi+1) - 1
if (jmin.gt.jmax)  GOTO 6
do 5 j=jmin,jmax
vj = ja(j)
if (vj-vi) 2, 5, 4
!
!------if a(vi,vj) is in strict lower triangle
!------check for previous occurrence of a(vj,vi)
2      lvk = vi
kmax = mark(vi) - 1
if (kmax .eq. 0) GOTO 4
do 3 k=1,kmax
lvk = l(lvk)
if (v(lvk).eq.vj) GOTO 5
3        continue
!----for unentered entries a(vi,vj)
4        if (sfs.ge.max)  GOTO 101
!
!------enter vj in element list for vi
mark(vi) = mark(vi) + 1
v(sfs) = vj
l(sfs) = l(vi)
l(vi) = sfs
sfs = sfs+1
!
!------enter vi in element list for vj
mark(vj) = mark(vj) + 1
v(sfs) = vi
l(sfs) = l(vj)
l(vj) = sfs
sfs = sfs+1
5      continue
6    continue
!
!----create degree lists and initialize mark vector
do 7 vi=1,n
dvi = mark(vi)
next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
nextvi = next(vi)
if (nextvi.gt.0)  last(nextvi) = vi
7    mark(vi) = tag
!
return
!
! ** error-  insufficient storage
101  flag = 9*n + vi
return
end
subroutine mdm(vk,tail, v,l, last,next, mark)
!***********************************************************************
!  mdm -- form element from uneliminated neighbors of vk
!***********************************************************************
integer  vk, tail, v(*), l(*), last(*), next(*), mark(*),tag, s,ls,vs,es, b,lb,vb, blp,blpmax
equivalence  (vs, es)
!
!----initialize tag and list of uneliminated neighbors
tag = mark(vk)
tail = vk
!
!----for each vertex/element vs/es in element list of vk
ls = l(vk)
1  s = ls
if (s.eq.0)  GOTO 5
ls = l(s)
vs = v(s)
if (next(vs).lt.0)  GOTO 2
!
!------if vs is uneliminated vertex, then tag and append to list of
!------uneliminated neighbors
mark(vs) = tag
l(tail) = s
tail = s
GOTO 4
!
!------if es is active element, then ...
!--------for each vertex vb in boundary list of element es
2      lb = l(es)
blpmax = last(es)
do 3 blp=1,blpmax
b = lb
lb = l(b)
vb = v(b)
!
!----------if vb is untagged vertex, then tag and append to list of
!----------uneliminated neighbors
if (mark(vb).ge.tag)  GOTO 3
mark(vb) = tag
l(tail) = b
tail = b
3        continue
!
!--------mark es inactive
mark(es) = tag
!
4    GOTO 1
!
!----terminate list of uneliminated neighbors
5  l(tail) = 0
!
return
end
subroutine mdp(k,ek,tail, v,l, head,last,next, mark)
!***********************************************************************
!  mdp -- purge inactive elements and do mass elimination
!***********************************************************************
integer  ek, tail,  v(*), l(*),  head(*), last(*), next(*),mark(*),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax
integer k, i
!
!----initialize tag
tag = mark(ek)
!
!----for each vertex vi in ek
li = ek
ilpmax = last(ek)
if (ilpmax.le.0)  GOTO 12
do 11 ilp=1,ilpmax
i = li
li = l(i)
vi = v(li)
!
!------remove vi from degree list
if (last(vi).eq.0)  GOTO 3
if (last(vi).gt.0)  GOTO 1
head(-last(vi)) = next(vi)
GOTO 2
1        next(last(vi)) = next(vi)
2      if (next(vi).gt.0)  last(next(vi)) = last(vi)
!
!------remove inactive items from element list of vi
3    ls = vi
4    s = ls
ls = l(s)
if (ls.eq.0)  GOTO 6
es = v(ls)
if (mark(es).lt.tag)  GOTO 5
free = ls
l(s) = l(ls)
ls = s
5      GOTO 4
!
!------if vi is interior vertex, then remove from list and eliminate
6    lvi = l(vi)
if (lvi.ne.0)  GOTO 7
l(i) = l(li)
li = i
!
k = k+1
next(vi) = -k
last(ek) = last(ek) - 1
GOTO 11
!
!------else ...
!--------classify vertex vi
7      if (l(lvi).ne.0)  GOTO 9
evi = v(lvi)
if (next(evi).ge.0)  GOTO 9
if (mark(evi).lt.0)  GOTO 8
!
!----------if vi is prototype vertex, then mark as such, initialize
!----------overlap count for corresponding element, and move vi to end
!----------of boundary list
last(vi) = evi
mark(evi) = -1
l(tail) = li
tail = li
l(i) = l(li)
li = i
GOTO 10
!
!----------else if vi is duplicate vertex, then mark as such and adjust
!----------overlap count for corresponding element
8            last(vi) = 0
mark(evi) = mark(evi) - 1
GOTO 10
!
!----------else mark vi to compute degree
9            last(vi) = -ek
!
!--------insert ek in element list of vi
10      v(free) = ek
l(free) = l(vi)
l(vi) = free
11    continue
!
!----terminate boundary list
12  l(tail) = 0
!
return
end
subroutine mdu(ek,dmin, v,l, head,last,next, mark)
!***********************************************************************
!  mdu -- update degrees of uneliminated vertices in ek
!***********************************************************************
integer  ek, dmin,  v(*), l(*),  head(*), last(*), next(*),mark(*),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,blp,blpmax
integer i
equivalence  (vs, es)
!
!----initialize tag
tag = mark(ek) - last(ek)
!
!----for each vertex vi in ek
i = ek
ilpmax = last(ek)
if (ilpmax.le.0)  GOTO 11
do 10 ilp=1,ilpmax
i = l(i)
vi = v(i)
if (last(vi))  1, 10, 8
!
!------if vi neither prototype nor duplicate vertex, then merge elements
!------to compute degree
1      tag = tag + 1
dvi = last(ek)
!
!--------for each vertex/element vs/es in element list of vi
s = l(vi)
2      s = l(s)
if (s.eq.0)  GOTO 9
vs = v(s)
if (next(vs).lt.0)  GOTO 3
!
!----------if vs is uneliminated vertex, then tag and adjust degree
mark(vs) = tag
dvi = dvi + 1
GOTO 5
!
!----------if es is active element, then expand
!------------check for outmatched vertex
3          if (mark(es).lt.0)  GOTO 6
!
!------------for each vertex vb in es
b = es
blpmax = last(es)
do 4 blp=1,blpmax
b = l(b)
vb = v(b)
!
!--------------if vb is untagged, then tag and adjust degree
if (mark(vb).ge.tag)  GOTO 4
mark(vb) = tag
dvi = dvi + 1
4            continue
!
5        GOTO 2
!
!------else if vi is outmatched vertex, then adjust overlaps but do not
!------compute degree
6      last(vi) = 0
mark(es) = mark(es) - 1
7      s = l(s)
if (s.eq.0)  GOTO 10
es = v(s)
if (mark(es).lt.0)  mark(es) = mark(es) - 1
GOTO 7
!
!------else if vi is prototype vertex, then calculate degree by
!------inclusion/exclusion and reset overlap count
8      evi = last(vi)
dvi = last(ek) + last(evi) + mark(evi)
mark(evi) = 0
!
!------insert vi in appropriate degree list
9    next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
if (next(vi).gt.0)  last(next(vi)) = vi
if (dvi.lt.dmin)  dmin = dvi
!
10    continue
!
11  return
end
subroutine sro(n, ip, ia,ja,a, q, r, dflag)
!***********************************************************************
!  sro -- symmetric reordering of sparse symmetric matrix
!***********************************************************************
!
!  description
!
!    the nonzero entries of the matrix m are assumed to be stored
!    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
!    are stored if i ne j).
!
!    sro does not rearrange the order of the rows, but does move
!    nonzeroes from one row to another to ensure that if m(i,j) will be
!    in the upper triangle of m with respect to the new ordering, then
!    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
!    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
!    stored in row j (and thus m(i,j) is not stored).
!
!
!  additional parameters
!
!    q     - integer one-dimensional work array.  dimension = n
!
!    r     - integer one-dimensional work array.  dimension = number of
!            nonzero entries in the upper triangle of m
!
!    dflag - logical variable.  if dflag = .true., then store nonzero
!            diagonal elements at the beginning of the row
!
!-----------------------------------------------------------------------
!
integer  ip(*),  ia(*), ja(*),  q(*), r(*)
integer n, i, ilast, j, jak, jdummy, jmin, jmax, k
!...  real  a(*),  ak
double precision  a(*),  ak
logical  dflag
!
!
!--phase 1 -- find row in which to store each nonzero
!----initialize count of nonzeroes to be stored in each row
do 1 i=1,n
1     q(i) = 0
!
!----for each nonzero element a(j)
do 3 i=1,n
jmin = ia(i)
jmax = ia(i+1) - 1
if (jmin.gt.jmax)  GOTO 3
do 2 j=jmin,jmax
!
!--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
k = ja(j)
if (ip(k).lt.ip(i))  ja(j) = i
if (ip(k).ge.ip(i))  k = i
r(j) = k
!
!--------... and increment count of nonzeroes (=q(r(j)) in that row
2       q(k) = q(k) + 1
3     continue
!
!
!--phase 2 -- find new ia and permutation to apply to (ja,a)
!----determine pointers to delimit rows in permuted (ja,a)
do 4 i=1,n
ia(i+1) = ia(i) + q(i)
4     q(i) = ia(i+1)
!
!----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
!----for each nonzero element (in reverse order)
ilast = 0
jmin = ia(1)
jmax = ia(n+1) - 1
j = jmax
do 6 jdummy=jmin,jmax
i = r(j)
if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  GOTO 5
!
!------if dflag, then put diagonal nonzero at beginning of row
r(j) = ia(i)
ilast = i
GOTO 6
!
!------put (off-diagonal) nonzero in last unused location in row
5       q(i) = q(i) - 1
r(j) = q(i)
!
6     j = j-1
!
!
!--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
do 8 j=jmin,jmax
7     if (r(j).eq.j)  GOTO 8
k = r(j)
r(j) = r(k)
r(k) = k
jak = ja(k)
ja(k) = ja(j)
ja(j) = jak
ak = a(k)
a(k) = a(j)
a(j) = ak
GOTO 7
8     continue
!
return
end

subroutine cdrv(n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
!*** subroutine cdrv
!*** driver for subroutines for solving sparse nonsymmetric systems of
!       linear equations (compressed pointer storage)
!
!
!    parameters
!    class abbreviations are--
!       n - integer variable
!       f - real variable
!       v - supplies a value to the driver
!       r - returns a result from the driver
!       i - used internally by the driver
!       a - array
!
! class - parameter
! ------+----------
!       -
!         the nonzero entries of the coefficient matrix m are stored
!    row-by-row in the array a.  to identify the individual nonzero
!    entries in each row, we need to know in which column each entry
!    lies.  the column indices which correspond to the nonzero entries
!    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  in addition, we need to know where each row starts and
!    how long it is.  the index positions in ja and a where the rows of
!    m begin are stored in the array ia.  i.e., if m(i,j) is the first
!    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
!    ia(i) = k.  moreover, the index in ja and a of the first location
!    following the last element in the last row is stored in ia(n+1).
!    thus, the number of entries in the i-th row is given by
!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!    consecutively in
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!    and the corresponding column indices are stored consecutively in
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!    for example, the 5 by 5 matrix
!                ( 1. 0. 2. 0. 0.)
!                ( 0. 3. 0. 0. 0.)
!            m = ( 0. 4. 5. 6. 0.)
!                ( 0. 0. 0. 7. 0.)
!                ( 0. 0. 0. 8. 9.)
!    would be stored as
!               - 1  2  3  4  5  6  7  8  9
!            ---+--------------------------
!            ia - 1  3  4  7  8 10
!            ja - 1  3  2  2  3  4  4  4  5
!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
! nv    - n     - number of variables/equations.
! fva   - a     - nonzero entries of the coefficient matrix m, stored
!       -           by rows.
!       -           size = number of nonzero entries in m.
! nva   - ia    - pointers to delimit the rows in a.
!       -           size = n+1.
! nva   - ja    - column numbers corresponding to the elements of a.
!       -           size = size of a.
! fva   - b     - right-hand side b.  b and z can the same array.
!       -           size = n.
! fra   - z     - solution x.  b and z can be the same array.
!       -           size = n.
!
!         the rows and columns of the original matrix m can be
!    reordered (e.g., to reduce fillin or ensure numerical stability)
!    before calling the driver.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
!    in the original order.
!         if the columns have been reordered (i.e.,  c(i).ne.i  for some
!    i), then the driver will call a subroutine (nroc) which rearranges
!    each row of ja and a, leaving the rows in the original order, but
!    placing the elements of each row in increasing order with respect
!    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
!    been called already.
!
! nva   - r     - ordering of the rows of m.
!       -           size = n.
! nva   - c     - ordering of the columns of m.
!       -           size = n.
! nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
!       -           ic(c(i)) = i  for i=1,...,n.
!       -           size = n.
!
!         the solution of the system of linear equations is divided into
!    three stages --
!      nsfc -- the matrix m is processed symbolically to determine where
!               fillin will occur during the numeric factorization.
!      nnfc -- the matrix m is factored numerically into the product ldu
!               of a unit lower triangular matrix l, a diagonal matrix
!               d, and a unit upper triangular matrix u, and the system
!               mx = b  is solved.
!      nnsc -- the linear system  mx = b  is solved using the ldu
!  or           factorization from nnfc.
!      nntc -- the transposed linear system  mt x = b  is solved using
!               the ldu factorization from nnf.
!    for several systems whose coefficient matrices have the same
!    nonzero structure, nsfc need be done only once (for the first
!    system).  then nnfc is done once for each additional system.  for
!    several systems with the same coefficient matrix, nsfc and nnfc
!    need be done only once (for the first system).  then nnsc or nntc
!    is done once for each additional right-hand side.
!
! nv    - path  - path specification.  values and their meanings are --
!       -           1  perform nroc, nsfc, and nnfc.
!       -           2  perform nnfc only  (nsfc is assumed to have been
!       -               done in a manner compatible with the storage
!       -               allocation used in the driver).
!       -           3  perform nnsc only  (nsfc and nnfc are assumed to
!       -               have been done in a manner compatible with the
!       -               storage allocation used in the driver).
!       -           4  perform nntc only  (nsfc and nnfc are assumed to
!       -               have been done in a manner compatible with the
!       -               storage allocation used in the driver).
!       -           5  perform nroc and nsfc.
!
!         various errors are detected by the driver and the individual
!    subroutines.
!
! nr    - flag  - error flag.  values and their meanings are --
!       -             0     no errors detected
!       -             n+k   null row in a  --  row = k
!       -            2n+k   duplicate entry in a  --  row = k
!       -            3n+k   insufficient storage in nsfc  --  row = k
!       -            4n+1   insufficient storage in nnfc
!       -            5n+k   null pivot  --  row = k
!       -            6n+k   insufficient storage in nsfc  --  row = k
!       -            7n+1   insufficient storage in nnfc
!       -            8n+k   zero pivot  --  row = k
!       -           10n+1   insufficient storage in cdrv
!       -           11n+1   illegal path specification
!
!         working storage is needed for the factored form of the matrix
!    m plus various temporary vectors.  the arrays isp and rsp should be
!    equivalenced.  integer storage is allocated from the beginning of
!    isp and real storage from the end of rsp.
!
! nv    - nsp   - declared dimension of rsp.  nsp generally must
!       -           be larger than  8n+2 + 2k  (where  k = (number of
!       -           nonzero entries in m)).
! nvira - isp   - integer working storage divided up into various arrays
!       -           needed by the subroutines.  isp and rsp should be
!       -           equivalenced.
!       -           size = lratio*nsp.
! fvira - rsp   - real working storage divided up into various arrays
!       -           needed by the subroutines.  isp and rsp should be
!       -           equivalenced.
!       -           size = nsp.
! nr    - esp   - if sufficient storage was available to perform the
!       -           symbolic factorization (nsfc), then esp is set to
!       -           the amount of excess storage provided (negative if
!       -           insufficient storage was available to perform the
!       -           numeric factorization (nnfc)).
!
!
!  conversion to double precision
!
!    to convert these routines for double precision arrays..
!    (1) use the double precision declarations in place of the real
!    declarations in each subprogram, as given in comment cards.
!    (2) change the data-loaded value of the integer  lratio
!    in subroutine cdrv, as indicated below.
!    (3) change e0 to d0 in the constants in statement number 10
!    in subroutine nnfc and the line following that.
!
integer  r(*), c(*), ic(*),  ia(*), ja(*),  isp(*), esp,  path,flag,  d, u, q, row, tmp, ar,  umax
integer n, nsp, i, ijl, iju, il, ira, irac, irl, iru, iu, j, jl, jlmax, jra, jrl, jru, ju, jumax, jutmp, l, lmax, lratio, max
!     real  a(*), b(*), z(*), rsp(*)
double precision  a(*), b(*), z(*), rsp(*)
!
!  set lratio equal to the ratio between the length of floating point
!  and integer array data.  e. g., lratio = 1 for (real, integer),
!  lratio = 2 for (double precision, integer)
!
data lratio/2/
!
if (path.lt.1 .or. 5.lt.path)  GOTO 111
!******initialize and divide up temporary storage  *******************
il   = 1
ijl  = il  + (n+1)
iu   = ijl +   n
iju  = iu  + (n+1)
irl  = iju +   n
jrl  = irl +   n
jl   = jrl +   n
!
!  ******  reorder a if necessary, call nsfc if flag is set  ***********
if ((path-1) * (path-5) .ne. 0)  GOTO 5
max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
jlmax = max/2
q     = jl   + jlmax
ira   = q    + (n+1)
jra   = ira  +   n
irac  = jra  +   n
iru   = irac +   n
jru   = iru  +   n
jutmp = jru  +   n
jumax = lratio*nsp  + 1 - jutmp
esp = max/lratio
if (jlmax.le.0 .or. jumax.le.0)  GOTO 110
!
do 1 i=1,n
if (c(i).ne.i)  GOTO 2
1      continue
GOTO 3
2    ar = nsp + 1 - n
call  nroc(n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
if (flag.ne.0)  GOTO 100
!
3    call  nsfc(n, r, ic, ia,ja,jlmax, isp(il), isp(jl), isp(ijl),jumax, isp(iu), isp(jutmp), isp(iju),isp(q), isp(ira), &
isp(jra), isp(irac),isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
if(flag .ne. 0)  GOTO 100
!  ******  move ju next to jl  *****************************************
jlmax = isp(ijl+n-1)
ju    = jl + jlmax
jumax = isp(iju+n-1)
if (jumax.le.0)  GOTO 5
do 4 j=1,jumax
4      isp(ju+j-1) = isp(jutmp+j-1)
!
!  ******  call remaining subroutines  *********************************
5  jlmax = isp(ijl+n-1)
ju    = jl  + jlmax
jumax = isp(iju+n-1)
l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
lmax  = isp(il+n) - 1
d     = l   + lmax
u     = d   + n
row   = nsp + 1 - n
tmp   = row - n
umax  = tmp - u
esp   = umax - (isp(iu+n) - 1)
!
if ((path-1) * (path-2) .ne. 0)  GOTO 6
if (umax.lt.0)  GOTO 110
call nnfc(n,  r, c, ic,  ia, ja, a, z, b,lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),umax, isp(iu), &
isp(ju), isp(iju), rsp(u),rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
if(flag .ne. 0)  GOTO 100
!
6  if ((path-3) .ne. 0)  GOTO 7
call nnsc(n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),z, b,  rsp(tmp))
!
7  if ((path-4) .ne. 0)  GOTO 8
call nntc(n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),z, b,  rsp(tmp))
8  return
!
! ** error.. error detected in nroc, nsfc, nnfc, or nnsc
100  return
! ** error.. insufficient storage
110  flag = 10*n + 1
return
! ** error.. illegal path specification
111  flag = 11*n + 1
return
end
subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
!
!       ----------------------------------------------------------------
!
!               yale sparse matrix package - nonsymmetric codes
!                    solving the system of equations mx = b
!
!    i.   calling sequences
!         the coefficient matrix can be processed by an ordering routine
!    (e.g., to reduce fillin or ensure numerical stability) before using
!    the remaining subroutines.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
!    is used, then nroc should be used to reorder the coefficient matrix
!    the calling sequence is --
!        (       (matrix ordering))
!        (nroc   (matrix reordering))
!         nsfc   (symbolic factorization to determine where fillin will
!                  occur during numeric factorization)
!         nnfc   (numeric factorization into product ldu of unit lower
!                  triangular matrix l, diagonal matrix d, and unit
!                  upper triangular matrix u, and solution of linear
!                  system)
!         nnsc   (solution of linear system for additional right-hand
!                  side using ldu factorization from nnfc)
!    (if only one system of equations is to be solved, then the
!    subroutine trk should be used.)
!
!    ii.  storage of sparse matrices
!         the nonzero entries of the coefficient matrix m are stored
!    row-by-row in the array a.  to identify the individual nonzero
!    entries in each row, we need to know in which column each entry
!    lies.  the column indices which correspond to the nonzero entries
!    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  in addition, we need to know where each row starts and
!    how long it is.  the index positions in ja and a where the rows of
!    m begin are stored in the array ia.  i.e., if m(i,j) is the first
!    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
!    ia(i) = k.  moreover, the index in ja and a of the first location
!    following the last element in the last row is stored in ia(n+1).
!    thus, the number of entries in the i-th row is given by
!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!    consecutively in
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!    and the corresponding column indices are stored consecutively in
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!    for example, the 5 by 5 matrix
!                ( 1. 0. 2. 0. 0.)
!                ( 0. 3. 0. 0. 0.)
!            m = ( 0. 4. 5. 6. 0.)
!                ( 0. 0. 0. 7. 0.)
!                ( 0. 0. 0. 8. 9.)
!    would be stored as
!               - 1  2  3  4  5  6  7  8  9
!            ---+--------------------------
!            ia - 1  3  4  7  8 10
!            ja - 1  3  2  2  3  4  4  4  5
!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
!         the strict upper (lower) triangular portion of the matrix
!    u (l) is stored in a similar fashion using the arrays  iu, ju, u
!    (il, jl, l)  except that an additional array iju (ijl) is used to
!    compress storage of ju (jl) by allowing some sequences of column
!    (row) indices to used for more than one row (column)  (n.b., l is
!    stored by columns).  iju(k) (ijl(k)) points to the starting
!    location in ju (jl) of entries for the kth row (column).
!    compression in ju (jl) occurs in two ways.  first, if a row
!    (column) i was merged into the current row (column) k, and the
!    number of elements merged in from (the tail portion of) row
!    (column) i is the same as the final length of row (column) k, then
!    the kth row (column) and the tail of row (column) i are identical
!    and iju(k) (ijl(k)) points to the start of the tail.  second, if
!    some tail portion of the (k-1)st row (column) is identical to the
!    head of the kth row (column), then iju(k) (ijl(k)) points to the
!    start of that tail portion.  for example, the nonzero structure of
!    the strict upper triangular part of the matrix
!            d 0 x x x
!            0 d 0 x x
!            0 0 d x 0
!            0 0 0 d x
!            0 0 0 0 d
!    would be represented as
!                - 1 2 3 4 5 6
!            ----+------------
!             iu - 1 4 6 7 8 8
!             ju - 3 4 5 4
!            iju - 1 2 4 3           .
!    the diagonal entries of l and u are assumed to be equal to one and
!    are not stored.  the array d contains the reciprocals of the
!    diagonal entries of the matrix d.
!
!    iii. additional storage savings
!         in nsfc, r and ic can be the same array in the calling
!    sequence if no reordering of the coefficient matrix has been done.
!         in nnfc, r, c, and ic can all be the same array if no
!    reordering has been done.  if only the rows have been reordered,
!    then c and ic can be the same array.  if the row and column
!    orderings are the same, then r and c can be the same array.  z and
!    row can be the same array.
!         in nnsc or nntc, r and c can be the same array if no
!    reordering has been done or if the row and column orderings are the
!    same.  z and b can be the same array.  however, then b will be
!    destroyed.
!
!    iv.  parameters
!         following is a list of parameters to the programs.  names are
!    uniform among the various subroutines.  class abbreviations are --
!       n - integer variable
!       f - real variable
!       v - supplies a value to a subroutine
!       r - returns a result from a subroutine
!       i - used internally by a subroutine
!       a - array
!
! class - parameter
! ------+----------
! fva   - a     - nonzero entries of the coefficient matrix m, stored
!       -           by rows.
!       -           size = number of nonzero entries in m.
! fva   - b     - right-hand side b.
!       -           size = n.
! nva   - c     - ordering of the columns of m.
!       -           size = n.
! fvra  - d     - reciprocals of the diagonal entries of the matrix d.
!       -           size = n.
! nr    - flag  - error flag.  values and their meanings are --
!       -            0     no errors detected
!       -            n+k   null row in a  --  row = k
!       -           2n+k   duplicate entry in a  --  row = k
!       -           3n+k   insufficient storage for jl  --  row = k
!       -           4n+1   insufficient storage for l
!       -           5n+k   null pivot  --  row = k
!       -           6n+k   insufficient storage for ju  --  row = k
!       -           7n+1   insufficient storage for u
!       -           8n+k   zero pivot  --  row = k
! nva   - ia    - pointers to delimit the rows of a.
!       -           size = n+1.
! nvra  - ijl   - pointers to the first element in each column in jl,
!       -           used to compress storage in jl.
!       -           size = n.
! nvra  - iju   - pointers to the first element in each row in ju, used
!       -           to compress storage in ju.
!       -           size = n.
! nvra  - il    - pointers to delimit the columns of l.
!       -           size = n+1.
! nvra  - iu    - pointers to delimit the rows of u.
!       -           size = n+1.
! nva   - ja    - column numbers corresponding to the elements of a.
!       -           size = size of a.
! nvra  - jl    - row numbers corresponding to the elements of l.
!       -           size = jlmax.
! nv    - jlmax - declared dimension of jl.  jlmax must be larger than
!       -           the number of nonzeros in the strict lower triangle
!       -           of m plus fillin minus compression.
! nvra  - ju    - column numbers corresponding to the elements of u.
!       -           size = jumax.
! nv    - jumax - declared dimension of ju.  jumax must be larger than
!       -           the number of nonzeros in the strict upper triangle
!       -           of m plus fillin minus compression.
! fvra  - l     - nonzero entries in the strict lower triangular portion
!       -           of the matrix l, stored by columns.
!       -           size = lmax.
! nv    - lmax  - declared dimension of l.  lmax must be larger than
!       -           the number of nonzeros in the strict lower triangle
!       -           of m plus fillin  (il(n+1)-1 after nsfc).
! nv    - n     - number of variables/equations.
! nva   - r     - ordering of the rows of m.
!       -           size = n.
! fvra  - u     - nonzero entries in the strict upper triangular portion
!       -           of the matrix u, stored by rows.
!       -           size = umax.
! nv    - umax  - declared dimension of u.  umax must be larger than
!       -           the number of nonzeros in the strict upper triangle
!       -           of m plus fillin  (iu(n+1)-1 after nsfc).
! fra   - z     - solution x.
!       -           size = n.
!
!       ----------------------------------------------------------------
!
!*** subroutine nroc
!*** reorders rows of a, leaving row order unchanged
!
!
!       input parameters.. n, ic, ia, ja, a
!       output parameters.. ja, a, flag
!
!       parameters used internally..
! nia   - p     - at the kth step, p is a linked list of the reordered
!       -           column indices of the kth row of a.  p(n+1) points
!       -           to the first entry in the list.
!       -           size = n+1.
! nia   - jar   - at the kth step,jar contains the elements of the
!       -           reordered column indices of a.
!       -           size = n.
! fia   - ar    - at the kth step, ar contains the elements of the
!       -           reordered row of a.
!       -           size = n.
!
integer  ic(*), ia(*), ja(*), jar(*), p(*), flag
integer n, i, j, jmax, jmin, k, newj
!     real  a(*), ar(*)
double precision  a(*), ar(*)
!
!  ******  for each nonempty row  *******************************
do 5 k=1,n
jmin = ia(k)
jmax = ia(k+1) - 1
if(jmin .gt. jmax) GOTO 5
p(n+1) = n + 1
!  ******  insert each element in the list  *********************
do 3 j=jmin,jmax
newj = ic(ja(j))
i = n + 1
1      if(p(i) .ge. newj) GOTO 2
i = p(i)
GOTO 1
2      if(p(i) .eq. newj) GOTO 102
p(newj) = p(i)
p(i) = newj
jar(newj) = ja(j)
ar(newj) = a(j)
3      continue
!  ******  replace old row in ja and a  *************************
i = n + 1
do 4 j=jmin,jmax
i = p(i)
ja(j) = jar(i)
4      a(j) = ar(i)
5    continue
flag = 0
return
!
! ** error.. duplicate entry in a
102  flag = n + k
return
end
subroutine nsfc(n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,q, ira,jra, irac, irl,jrl, iru,jru, flag)
!*** subroutine nsfc
!*** symbolic ldu-factorization of nonsymmetric sparse matrix
!      (compressed pointer storage)
!
!
!       input variables.. n, r, ic, ia, ja, jlmax, jumax.
!       output variables.. il, jl, ijl, iu, ju, iju, flag.
!
!       parameters used internally..
! nia   - q     - suppose  m*  is the result of reordering  m.  if
!       -           processing of the ith row of  m*  (hence the ith
!       -           row of  u) is being done,  q(j)  is initially
!       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
!       -           values need not be stored, each entry points to the
!       -           next nonzero and  q(n+1)  points to the first.  n+1
!       -           indicates the end of the list.  for example, if n=9
!       -           and the 5th row of  m*  is
!       -              0 x x 0 x 0 0 x 0
!       -           then  q  will initially be
!       -              a a a a 8 a a 10 5           (a - arbitrary).
!       -           as the algorithm proceeds, other elements of  q
!       -           are inserted in the list because of fillin.
!       -           q  is used in an analogous manner to compute the
!       -           ith column of  l.
!       -           size = n+1.
! nia   - ira,  - vectors used to find the columns of  m.  at the kth
! nia   - jra,      step of the factorization,  irac(k)  points to the
! nia   - irac      head of a linked list in  jra  of row indices i
!       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
!       -           indicates the end of the list.  ira(i)  (i.ge.k)
!       -           points to the smallest j such that j .ge. k and
!       -           m(i,j)  is nonzero.
!       -           size of each = n.
! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
! nia   - jrl       of the factorization,  jrl(k)  points to the head
!       -           of a linked list in  jrl  of column indices j
!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!       -           points to the smallest i such that i .ge. k and
!       -           l(i,j)  is nonzero.
!       -           size of each = n.
! nia   - iru,  - vectors used in a manner analogous to  irl and jrl
! nia   - jru       to find the columns of  u.
!       -           size of each = n.
!
!  internal variables..
!    jlptr - points to the last position used in  jl.
!    juptr - points to the last position used in  ju.
!    jmin,jmax - are the indices in  a or u  of the first and last
!                elements to be examined in a given row.
!                for example,  jmin=ia(k), jmax=ia(k+1)-1.
!
integer cend, qm, rend, rk, vj
integer ia(*), ja(*), ira(*), jra(*), il(*), jl(*), ijl(*)
integer iu(*), ju(*), iju(*), irl(*), jrl(*), iru(*), jru(*)
integer r(*), ic(*), q(*), irac(*), flag
integer n, jlmax, jumax, i, i1, iak, irai, irll, irul, j, jaiak, jairai, jlmin, jtmp, jumin, juptr
integer k, lasti, lastid, long, luk, m, jlptr, jmax, jmin, np1
!
!  ******  initialize pointers  ****************************************
np1 = n + 1
jlmin = 1
jlptr = 0
il(1) = 1
jumin = 1
juptr = 0
iu(1) = 1
do 1 k=1,n
irac(k) = 0
jra(k) = 0
jrl(k) = 0
1    jru(k) = 0
!  ******  initialize column pointers for a  ***************************
do 2 k=1,n
rk = r(k)
iak = ia(rk)
if (iak .ge. ia(rk+1))  GOTO 101
jaiak = ic(ja(iak))
if (jaiak .gt. k)  GOTO 105
jra(k) = irac(jaiak)
irac(jaiak) = k
2    ira(k) = iak
!
!  ******  for each column of l and row of u  **************************
do 41 k=1,n
!
!  ******  initialize q for computing kth column of l  *****************
q(np1) = np1
luk = -1
!  ******  by filling in kth column of a  ******************************
vj = irac(k)
if (vj .eq. 0)  GOTO 5
3      qm = np1
4      m = qm
qm =  q(m)
if (qm .lt. vj)  GOTO 4
if (qm .eq. vj)  GOTO 102
luk = luk + 1
q(m) = vj
q(vj) = qm
vj = jra(vj)
if (vj .ne. 0)  GOTO 3
!  ******  link through jru  *******************************************
5    lastid = 0
lasti = 0
ijl(k) = jlptr
i = k
6      i = jru(i)
if (i .eq. 0)  GOTO 10
qm = np1
jmin = irl(i)
jmax = ijl(i) + il(i+1) - il(i) - 1
long = jmax - jmin
if (long .lt. 0)  GOTO 6
jtmp = jl(jmin)
if (jtmp .ne. k)  long = long + 1
if (jtmp .eq. k)  r(i) = -r(i)
if (lastid .ge. long)  GOTO 7
lasti = i
lastid = long
!  ******  and merge the corresponding columns into the kth column  ****
7      do 9 j=jmin,jmax
vj = jl(j)
8        m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 8
if (qm .eq. vj)  GOTO 9
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
9        continue
GOTO 6
!  ******  lasti is the longest column merged into the kth  ************
!  ******  see if it equals the entire kth column  *********************
10    qm = q(np1)
if (qm .ne. k)  GOTO 105
if (luk .eq. 0)  GOTO 17
if (lastid .ne. luk)  GOTO 11
!  ******  if so, jl can be compressed  ********************************
irll = irl(lasti)
ijl(k) = irll + 1
if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
GOTO 17
!  ******  if not, see if kth column can overlap the previous one  *****
11    if (jlmin .gt. jlptr)  GOTO 15
qm = q(qm)
do 12 j=jlmin,jlptr
if (jl(j) - qm)  12, 13, 15
12      continue
GOTO 15
13    ijl(k) = j
do 14 i=j,jlptr
if (jl(i) .ne. qm)  GOTO 15
qm = q(qm)
if (qm .gt. n)  GOTO 17
14      continue
jlptr = j - 1
!  ******  move column indices from q to jl, update vectors  ***********
15    jlmin = jlptr + 1
ijl(k) = jlmin
if (luk .eq. 0)  GOTO 17
jlptr = jlptr + luk
if (jlptr .gt. jlmax)  GOTO 103
qm = q(np1)
do 16 j=jlmin,jlptr
qm = q(qm)
16        jl(j) = qm
17    irl(k) = ijl(k)
il(k+1) = il(k) + luk
!
!  ******  initialize q for computing kth row of u  ********************
q(np1) = np1
luk = -1
!  ******  by filling in kth row of reordered a  ***********************
rk = r(k)
jmin = ira(k)
jmax = ia(rk+1) - 1
if (jmin .gt. jmax)  GOTO 20
do 19 j=jmin,jmax
vj = ic(ja(j))
qm = np1
18      m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 18
if (qm .eq. vj)  GOTO 102
luk = luk + 1
q(m) = vj
q(vj) = qm
19      continue
!  ******  link through jrl,  ******************************************
20    lastid = 0
lasti = 0
iju(k) = juptr
i = k
i1 = jrl(k)
21      i = i1
if (i .eq. 0)  GOTO 26
i1 = jrl(i)
qm = np1
jmin = iru(i)
jmax = iju(i) + iu(i+1) - iu(i) - 1
long = jmax - jmin
if (long .lt. 0)  GOTO 21
jtmp = ju(jmin)
if (jtmp .eq. k)  GOTO 22
!  ******  update irl and jrl, *****************************************
long = long + 1
cend = ijl(i) + il(i+1) - il(i)
irl(i) = irl(i) + 1
if (irl(i) .ge. cend)  GOTO 22
j = jl(irl(i))
jrl(i) = jrl(j)
jrl(j) = i
22      if (lastid .ge. long)  GOTO 23
lasti = i
lastid = long
!  ******  and merge the corresponding rows into the kth row  **********
23      do 25 j=jmin,jmax
vj = ju(j)
24        m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 24
if (qm .eq. vj)  GOTO 25
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
25        continue
GOTO 21
!  ******  update jrl(k) and irl(k)  ***********************************
26    if (il(k+1) .le. il(k))  GOTO 27
j = jl(irl(k))
jrl(k) = jrl(j)
jrl(j) = k
!  ******  lasti is the longest row merged into the kth  ***************
!  ******  see if it equals the entire kth row  ************************
27    qm = q(np1)
if (qm .ne. k)  GOTO 105
if (luk .eq. 0)  GOTO 34
if (lastid .ne. luk)  GOTO 28
!  ******  if so, ju can be compressed  ********************************
irul = iru(lasti)
iju(k) = irul + 1
if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
GOTO 34
!  ******  if not, see if kth row can overlap the previous one  ********
28    if (jumin .gt. juptr)  GOTO 32
qm = q(qm)
do 29 j=jumin,juptr
if (ju(j) - qm)  29, 30, 32
29      continue
GOTO 32
30    iju(k) = j
do 31 i=j,juptr
if (ju(i) .ne. qm)  GOTO 32
qm = q(qm)
if (qm .gt. n)  GOTO 34
31      continue
juptr = j - 1
!  ******  move row indices from q to ju, update vectors  **************
32    jumin = juptr + 1
iju(k) = jumin
if (luk .eq. 0)  GOTO 34
juptr = juptr + luk
if (juptr .gt. jumax)  GOTO 106
qm = q(np1)
do 33 j=jumin,juptr
qm = q(qm)
33        ju(j) = qm
34    iru(k) = iju(k)
iu(k+1) = iu(k) + luk
!
!  ******  update iru, jru  ********************************************
i = k
35      i1 = jru(i)
if (r(i) .lt. 0)  GOTO 36
rend = iju(i) + iu(i+1) - iu(i)
if (iru(i) .ge. rend)  GOTO 37
j = ju(iru(i))
jru(i) = jru(j)
jru(j) = i
GOTO 37
36      r(i) = -r(i)
37      i = i1
if (i .eq. 0)  GOTO 38
iru(i) = iru(i) + 1
GOTO 35
!
!  ******  update ira, jra, irac  **************************************
38    i = irac(k)
if (i .eq. 0)  GOTO 41
39      i1 = jra(i)
ira(i) = ira(i) + 1
if (ira(i) .ge. ia(r(i)+1))  GOTO 40
irai = ira(i)
jairai = ic(ja(irai))
if (jairai .gt. i)  GOTO 40
jra(i) = irac(jairai)
irac(jairai) = i
40      i = i1
if (i .ne. 0)  GOTO 39
41    continue
!
ijl(n) = jlptr
iju(n) = juptr
flag = 0
return
!
! ** error.. null row in a
101  flag = n + rk
return
! ** error.. duplicate entry in a
102  flag = 2*n + rk
return
! ** error.. insufficient storage for jl
103  flag = 3*n + k
return
! ** error.. null pivot
105  flag = 5*n + k
return
! ** error.. insufficient storage for ju
106  flag = 6*n + k
return
end
subroutine nnfc(n, r,c,ic, ia,ja,a, z, b,lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,row, tmp, irl,jrl, flag)
!*** subroutine nnfc
!*** numerical ldu-factorization of sparse nonsymmetric matrix and
!      solution of system of linear equations (compressed pointer
!      storage)
!
!
!       input variables..  n, r, c, ic, ia, ja, a, b,
!                          il, jl, ijl, lmax, iu, ju, iju, umax
!       output variables.. z, l, d, u, flag
!
!       parameters used internally..
! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
! nia   - jrl       of the factorization,  jrl(k)  points to the head
!       -           of a linked list in  jrl  of column indices j
!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!       -           points to the smallest i such that i .ge. k and
!       -           l(i,j)  is nonzero.
!       -           size of each = n.
! fia   - row   - holds intermediate values in calculation of  u and l.
!       -           size = n.
! fia   - tmp   - holds new right-hand side  b*  for solution of the
!       -           equation ux = b*.
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row to
!      be examined.
!    sum - used in calculating  tmp.
!
integer rk,umax
integer  r(*), c(*), ic(*), ia(*), ja(*), il(*), jl(*), ijl(*)
integer  iu(*), ju(*), iju(*), irl(*), jrl(*), flag
integer n, lmax, i, i1, i2, ijlb, j, jmax, jmin, k, mu
!     real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
!     real tmp(*), lki, sum, dk
double precision  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
double precision  tmp(*), lki, sum, dk
!
!  ******  initialize pointers and test storage  ***********************
if(il(n+1)-1 .gt. lmax) GOTO 104
if(iu(n+1)-1 .gt. umax) GOTO 107
do 1 k=1,n
irl(k) = il(k)
jrl(k) = 0
1    continue
!
!  ******  for each row  ***********************************************
do 19 k=1,n
!  ******  reverse jrl and zero row where kth row of l will fill in  ***
row(k) = 0
i1 = 0
if (jrl(k) .eq. 0) GOTO 3
i = jrl(k)
2    i2 = jrl(i)
jrl(i) = i1
i1 = i
row(i) = 0
i = i2
if (i .ne. 0) GOTO 2
!  ******  set row to zero where u will fill in  ***********************
3    jmin = iju(k)
jmax = jmin + iu(k+1) - iu(k) - 1
if (jmin .gt. jmax) GOTO 5
do 4 j=jmin,jmax
4      row(ju(j)) = 0
!  ******  place kth row of a in row  **********************************
5    rk = r(k)
jmin = ia(rk)
jmax = ia(rk+1) - 1
do 6 j=jmin,jmax
row(ic(ja(j))) = a(j)
6      continue
!  ******  initialize sum, and link through jrl  ***********************
sum = b(rk)
i = i1
if (i .eq. 0) GOTO 10
!  ******  assign the kth row of l and adjust row, sum  ****************
7      lki = -row(i)
!  ******  if l is not required, then comment out the following line  **
l(irl(i)) = -lki
sum = sum + lki * tmp(i)
jmin = iu(i)
jmax = iu(i+1) - 1
if (jmin .gt. jmax) GOTO 9
mu = iju(i) - jmin
do 8 j=jmin,jmax
8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
9      i = jrl(i)
if (i .ne. 0) GOTO 7
!
!  ******  assign kth row of u and diagonal d, set tmp(k)  *************
10    if (row(k) .eq. 0.0d0) GOTO 108
dk = 1.0d0 / row(k)
d(k) = dk
tmp(k) = sum * dk
if (k .eq. n) GOTO 19
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax)  GOTO 12
mu = iju(k) - jmin
do 11 j=jmin,jmax
11      u(j) = row(ju(mu+j)) * dk
12    continue
!
!  ******  update irl and jrl, keeping jrl in decreasing order  ********
i = i1
if (i .eq. 0) GOTO 18
14    irl(i) = irl(i) + 1
i1 = jrl(i)
if (irl(i) .ge. il(i+1)) GOTO 17
ijlb = irl(i) - il(i) + ijl(i)
j = jl(ijlb)
15    if (i .gt. jrl(j)) GOTO 16
j = jrl(j)
GOTO 15
16    jrl(i) = jrl(j)
jrl(j) = i
17    i = i1
if (i .ne. 0) GOTO 14
18    if (irl(k) .ge. il(k+1)) GOTO 19
j = jl(ijl(k))
jrl(k) = jrl(j)
jrl(j) = k
19    continue
!
!  ******  solve  ux = tmp  by back substitution  **********************
k = n
do 22 i=1,n
sum =  tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax)  GOTO 21
mu = iju(k) - jmin
do 20 j=jmin,jmax
20      sum = sum - u(j) * tmp(ju(mu+j))
21    tmp(k) =  sum
z(c(k)) =  sum
22    k = k-1
flag = 0
return
!
! ** error.. insufficient storage for l
104  flag = 4*n + 1
return
! ** error.. insufficient storage for u
107  flag = 7*n + 1
return
! ** error.. zero pivot
108  flag = 8*n + k
return
end
subroutine nnsc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
!*** subroutine nnsc
!*** numerical solution of sparse nonsymmetric system of linear
!      equations given ldu-factorization (compressed pointer storage)
!
!
!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!       output variables.. z
!
!       parameters used internally..
! fia   - tmp   - temporary vector which gets result of solving  ly = b.
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row of
!      u or l  to be used.
!
integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
integer n, i, j, jmax, jmin, k, ml, mu
!     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
double precision  l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
!
!  ******  set tmp to reordered b  *************************************
do 1 k=1,n
1    tmp(k) = b(r(k))
!  ******  solve  ly = b  by forward substitution  *********************
do 3 k=1,n
jmin = il(k)
jmax = il(k+1) - 1
tmpk = -d(k) * tmp(k)
tmp(k) = -tmpk
if (jmin .gt. jmax) GOTO 3
ml = ijl(k) - jmin
do 2 j=jmin,jmax
2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
3    continue
!  ******  solve  ux = y  by back substitution  ************************
k = n
do 6 i=1,n
sum = -tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax) GOTO 5
mu = iju(k) - jmin
do 4 j=jmin,jmax
4      sum = sum + u(j) * tmp(ju(mu+j))
5    tmp(k) = -sum
z(c(k)) = -sum
k = k - 1
6    continue
return
end
subroutine nntc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
!*** subroutine nntc
!*** numeric solution of the transpose of a sparse nonsymmetric system
!      of linear equations given lu-factorization (compressed pointer
!      storage)
!
!
!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!       output variables.. z
!
!       parameters used internally..
! fia   - tmp   - temporary vector which gets result of solving ut y = b
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row of
!      u or l  to be used.
!
integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
integer n, i, j, jmax, jmin, k, ml, mu
!     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
double precision l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
!
!  ******  set tmp to reordered b  *************************************
do 1 k=1,n
1    tmp(k) = b(c(k))
!  ******  solve  ut y = b  by forward substitution  *******************
do 3 k=1,n
jmin = iu(k)
jmax = iu(k+1) - 1
tmpk = -tmp(k)
if (jmin .gt. jmax) GOTO 3
mu = iju(k) - jmin
do 2 j=jmin,jmax
2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
3    continue
!  ******  solve  lt x = y  by back substitution  **********************
k = n
do 6 i=1,n
sum = -tmp(k)
jmin = il(k)
jmax = il(k+1) - 1
if (jmin .gt. jmax) GOTO 5
ml = ijl(k) - jmin
do 4 j=jmin,jmax
4      sum = sum + l(j) * tmp(jl(ml+j))
5    tmp(k) = -sum * d(k)
z(r(k)) = tmp(k)
k = k - 1
6    continue
return
end

SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM, F, JAC, PJAC, SLVS)
EXTERNAL F, JAC, PJAC, SLVS
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM
INTEGER LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IOWND2, ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION ROWND2, CM1, CM2, PDEST, PDLAST, RATIO, PDNORM
COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO, PDNORM, IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM
DOUBLE PRECISION ALPHA, DM1,DM2, EXM1,EXM2, PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)
SAVE SM1
DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0, 0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
!-----------------------------------------------------------------------
! DSTODA performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! Note: DSTODA is independent of the value of the iteration method
! indicator MITER, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTODA is done with the following variables:
!
! Y      = an array of length .ge. N used as the Y argument in
!          all calls to F and JAC.
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to F and JAC.
! YH     = an NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = a constant integer .ge. N, the first dimension of YH.
! YH1    = a one-dimensional array occupying the same space as YH.
! EWT    = an array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = an array of working storage, of length N.
! ACOR   = a work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = real and integer work arrays associated with matrix
!          operations in chord iteration (MITER .ne. 0).
! PJAC   = name of routine to evaluate and preprocess Jacobian matrix
!          and P = I - H*EL0*Jac, if a chord method is being used.
!          It also returns an estimate of norm(Jac) in PDNORM.
! SLVS   = name of routine to solve linear system in chord iteration.
! CCMAX  = maximum relative change in H*EL0 before PJAC is called.
! H      = the step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = the minimum absolute value of the step size H to be used.
! HMXI   = inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = the independent variable. TN is updated on each step taken.
! JSTART = an integer used for input only, with the following
!          values and meanings:
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of H,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings:
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in PJAC or SLVS.
!          A return with KFLAG = -1 or -2 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = the maximum order of integration method to be allowed.
! MAXCOR = the maximum number of corrector iterations allowed.
! MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
! MXNCF  = maximum number of convergence failures allowed.
! METH   = current method.
!          METH = 1 means Adams method (nonstiff)
!          METH = 2 means BDF method (stiff)
!          METH may be reset by DSTODA.
! MITER  = corrector iteration method.
!          MITER = 0 means functional iteration.
!          MITER = JT .gt. 0 means a chord iteration corresponding
!          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
!          communicated here as JTYP, but is not used in DSTODA
!          except to load MITER following a method switch.)
!          MITER may be reset by DSTODA.
! N      = the number of first-order differential equations.
!-----------------------------------------------------------------------
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set at 2
! for the next increase.
! DCFODE is called to get the needed coefficients for both methods.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
NSLP = 0
IPUP = MITER
IRET = 3
! Initialize switching parameters.  METH = 1 is assumed initially. -----
ICOUNT = 20
IRFLAG = 0
PDEST = 0.0D0
PDLAST = 0.0D0
RATIO = 5.0D0
CALL DCFODE (2, ELCO, TESCO)
DO 10 I = 1,5
10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
CALL DCFODE (1, ELCO, TESCO)
DO 20 I = 1,12
20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
GOTO 150
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MUSED) GOTO 160
CALL DCFODE (METH, ELCO, TESCO)
IALTH = L
IRET = 1
!-----------------------------------------------------------------------
! The el vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
!-----------------------------------------------------------------------
! If METH = 1, also restrict the new step size by the stability region.
! If this reduces H, set IRFLAG to 1 so that if there are roundoff
! problems later, we can assume that is the cause of the trouble.
!-----------------------------------------------------------------------
IF (METH .EQ. 2) GOTO 178
IRFLAG = 0
PDH = MAX(ABS(H)*PDLAST,0.000001D0)
IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GOTO 178
RH = SM1(NQ)/PDH
IRFLAG = 1
178  CONTINUE
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called, if a Jacobian is involved.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
PNORM = DMNORM (N, YH1, EWT)
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS-norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220  M = 0
RATE = 0.0D0
DEL = 0.0D0
DO 230 I = 1,N
230    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
IF (IPUP .LE. 0) GOTO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
IPUP = 0
RC = 1.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .NE. 0) GOTO 430
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
270  IF (MITER .NE. 0) GOTO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
DO 290 I = 1,N
SAVF(I) = H*SAVF(I) - YH(I,2)
290    Y(I) = SAVF(I) - ACOR(I)
DEL = DMNORM (N, Y, EWT)
DO 300 I = 1,N
Y(I) = YH(I,1) + EL(1)*SAVF(I)
300    ACOR(I) = SAVF(I)
GOTO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
350  DO 360 I = 1,N
360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
CALL SLVS (WM, IWM, Y, SAVF)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
DEL = DMNORM (N, Y, EWT)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + Y(I)
380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!
! We first check for a change of iterates that is the size of
! roundoff error.  If this occurs, the iteration has converged, and a
! new rate estimate is not formed.
! In all other cases, force at least two iterations to estimate a
! local Lipschitz constant estimate for Adams methods.
! On convergence, form PDEST = local maximum Lipschitz constant
! estimate.  PDLAST is the most recent nonzero estimate.
!-----------------------------------------------------------------------
400  CONTINUE
IF (DEL .LE. 100.0D0*PNORM*UROUND) GOTO 450
IF (M .EQ. 0 .AND. METH .EQ. 1) GOTO 405
IF (M .EQ. 0) GOTO 402
RM = 1024.0D0
IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
RATE = MAX(RATE,RM)
CRATE = MAX(0.2D0*CRATE,RM)
402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
IF (DCON .GT. 1.0D0) GOTO 405
PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
IF (PDEST .NE. 0.0D0) PDLAST = PDEST
GOTO 450
405  CONTINUE
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
DELP = DEL
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
GOTO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GOTO 430
ICF = 1
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
RMAX = 2.0D0
TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 680
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 670
IF (NCF .EQ. MXNCF) GOTO 670
RH = 0.25D0
IPUP = MITER
IREDO = 1
GOTO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450  JCUR = 0
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Decrease ICOUNT by 1, and if it is -1, consider switching methods.
! If a method switch is made, reset various parameters,
! rescale the YH array, and exit.  If there is no switch,
! consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
HU = H
NQU = NQ
MUSED = METH
DO 460 J = 1,L
DO 460 I = 1,N
460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
ICOUNT = ICOUNT - 1
IF (ICOUNT .GE. 0) GOTO 488
IF (METH .EQ. 2) GOTO 480
!-----------------------------------------------------------------------
! We are currently using an Adams method.  Consider switching to BDF.
! If the current order is greater than 5, assume the problem is
! not stiff, and skip this section.
! If the Lipschitz constant and error estimate are not polluted
! by roundoff, GOTO 470 and perform the usual test.
! Otherwise, switch to the BDF methods if the last step was
! restricted to insure stability (irflag = 1), and stay with Adams
! method if not.  When switching to BDF with polluted error estimates,
! in the absence of other information, double the step size.
!
! When the estimates are OK, we make the usual test by computing
! the step size we could have (ideally) used on this step,
! with the current (Adams) method, and also that for the BDF.
! If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
! Compare the two step sizes to decide whether to switch.
! The step size advantage must be at least RATIO = 5 to switch.
!-----------------------------------------------------------------------
IF (NQ .GT. 5) GOTO 488
IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0) GOTO 470
IF (IRFLAG .EQ. 0) GOTO 488
RH2 = 2.0D0
NQM2 = MIN(NQ,MXORDS)
GOTO 478
470  CONTINUE
EXSM = 1.0D0/L
RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RH1IT = 2.0D0*RH1
PDH = PDLAST*ABS(H)
IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
RH1 = MIN(RH1,RH1IT)
IF (NQ .LE. MXORDS) GOTO 474
NQM2 = MXORDS
LM2 = MXORDS + 1
EXM2 = 1.0D0/LM2
LM2P1 = LM2 + 1
DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
GOTO 476
474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
NQM2 = NQ
476  CONTINUE
IF (RH2 .LT. RATIO*RH1) GOTO 488
! THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
478  RH = RH2
ICOUNT = 20
METH = 2
MITER = JTYP
PDLAST = 0.0D0
NQ = NQM2
L = NQ + 1
GOTO 170
!-----------------------------------------------------------------------
! We are currently using a BDF method.  Consider switching to Adams.
! Compute the step size we could have (ideally) used on this step,
! with the current (BDF) method, and also that for the Adams.
! If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
! Compare the two step sizes to decide whether to switch.
! The step size advantage must be at least 5/RATIO = 1 to switch.
! If the step size for Adams would be so small as to cause
! roundoff pollution, we stay with BDF.
!-----------------------------------------------------------------------
480  CONTINUE
EXSM = 1.0D0/L
IF (MXORDN .GE. NQ) GOTO 484
NQM1 = MXORDN
LM1 = MXORDN + 1
EXM1 = 1.0D0/LM1
LM1P1 = LM1 + 1
DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
GOTO 486
484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
NQM1 = NQ
EXM1 = EXSM
486  RH1IT = 2.0D0*RH1
PDH = PDNORM*ABS(H)
IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
RH1 = MIN(RH1,RH1IT)
RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
IF (RH1*RATIO .LT. 5.0D0*RH2) GOTO 488
ALPHA = MAX(0.001D0,RH1)
DM1 = (ALPHA**EXM1)*DM1
IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GOTO 488
! The switch test passed.  Reset relevant quantities for Adams. --------
RH = RH1
ICOUNT = 20
METH = 1
MITER = 0
PDLAST = 0.0D0
NQ = NQM1
L = NQ + 1
GOTO 170
!
! No method switch is being made.  Do the usual step/order selection. --
488  CONTINUE
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -3) GOTO 640
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 550
DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
! If METH = 1, limit RH according to the stability region also. --------
550  IF (METH .EQ. 2) GOTO 560
PDH = MAX(ABS(H)*PDLAST,0.000001D0)
IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
RHSM = MIN(RHSM,SM1(NQ)/PDH)
IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
PDEST = 0.0D0
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
! If METH = 1 and H is restricted by stability, bypass 10 percent test.
620  IF (METH .EQ. 2) GOTO 622
IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GOTO 625
622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GOTO 610
625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640  IF (KFLAG .EQ. -10) GOTO 660
RH = 0.1D0
RH = MAX(HMIN/ABS(H),RH)
H = H*RH
DO 645 I = 1,N
645    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
DO 650 I = 1,N
650    YH(I,2) = H*SAVF(I)
IPUP = MITER
IALTH = 5
IF (NQ .EQ. 1) GOTO 200
NQ = 1
L = 2
IRET = 3
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -3
GOTO 720
690  RMAX = 10.0D0
700  R = 1.0D0/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- End of Subroutine DSTODA ----------------------
END

SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,F, JAC)
EXTERNAL F, JAC
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IOWND2, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION ROWND2, ROWNS2, PDNORM
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM, IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP, MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ, DMNORM, DFNORM, DBNORM
!-----------------------------------------------------------------------
! DPRJA is called by DSTODA to compute and process the matrix
! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
! J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
! matrix norm consistent with the weighted max-norm on vectors given
! by DMNORM) is computed, and J is overwritten by P.  P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix.  This is done
! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!
! In addition to variables described previously, communication
! with DPRJA uses the following:
! Y     = array containing predicted values on entry.
! FTEM  = work array of length N (ACOR in DSTODA).
! SAVF  = array containing f evaluated at predicted y.
! WM    = real work space for matrices.  On output it contains the
!         LU decomposition of P.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data:
!         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
! IWM   = integer work space containing pivot information, starting at
!         IWM(21).   IWM also contains the band parameters
!         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! EL0   = EL(1) (input).
! PDNORM= norm of Jacobian matrix. (Output).
! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!         P matrix found to be singular.
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses the Common variables EL0, H, TN, UROUND,
! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
NJE = NJE + 1
IERPJ = 0
JCUR = 1
HL0 = H*EL0
GOTO (100, 200, 300, 400, 500), MITER
! If MITER = 1, call JAC and multiply by scalar. -----------------------
100  LENP = N*N
DO 110 I = 1,LENP
110    WM(I+2) = 0.0D0
CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
CON = -HL0
DO 120 I = 1,LENP
120    WM(I+2) = WM(I+2)*CON
GOTO 240
! If MITER = 2, make N calls to F to approximate J. --------------------
200  FAC = DMNORM (N, SAVF, EWT)
R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
SRUR = WM(1)
J1 = 2
DO 230 J = 1,N
YJ = Y(J)
R = MAX(SRUR*ABS(YJ),R0/EWT(J))
Y(J) = Y(J) + R
FAC = -HL0/R
CALL F (NEQ, TN, Y, FTEM)
DO 220 I = 1,N
220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
Y(J) = YJ
J1 = J1 + N
230    CONTINUE
NFE = NFE + N
240  CONTINUE
! Compute norm of Jacobian. --------------------------------------------
PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
! Add identity matrix. -------------------------------------------------
J = 3
NP1 = N + 1
DO 250 I = 1,N
WM(J) = WM(J) + 1.0D0
250    J = J + NP1
! Do LU decomposition on P. --------------------------------------------
CALL DGEFA (WM(3), N, N, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! Dummy block only, since MITER is never 3 in this routine. ------------
300  RETURN
! If MITER = 4, call JAC and multiply by scalar. -----------------------
400  ML = IWM(1)
MU = IWM(2)
ML3 = ML + 3
MBAND = ML + MU + 1
MEBAND = MBAND + ML
LENP = MEBAND*N
DO 410 I = 1,LENP
410    WM(I+2) = 0.0D0
CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
CON = -HL0
DO 420 I = 1,LENP
420    WM(I+2) = WM(I+2)*CON
GOTO 570
! If MITER = 5, make MBAND calls to F to approximate J. ----------------
500  ML = IWM(1)
MU = IWM(2)
MBAND = ML + MU + 1
MBA = MIN(MBAND,N)
MEBAND = MBAND + ML
MEB1 = MEBAND - 1
SRUR = WM(1)
FAC = DMNORM (N, SAVF, EWT)
R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
DO 560 J = 1,MBA
DO 530 I = J,N,MBAND
YI = Y(I)
R = MAX(SRUR*ABS(YI),R0/EWT(I))
530      Y(I) = Y(I) + R
CALL F (NEQ, TN, Y, FTEM)
DO 550 JJ = J,N,MBAND
Y(JJ) = YH(JJ,1)
YJJ = Y(JJ)
R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
FAC = -HL0/R
I1 = MAX(JJ-MU,1)
I2 = MIN(JJ+ML,N)
II = JJ*MEB1 - ML + 2
DO 540 I = I1,I2
540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
550      CONTINUE
560    CONTINUE
NFE = NFE + MBA
570  CONTINUE
! Compute norm of Jacobian. --------------------------------------------
PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
! Add identity matrix. -------------------------------------------------
II = MBAND + 2
DO 580 I = 1,N
WM(II) = WM(II) + 1.0D0
580    II = II + MEBAND
! Do LU decomposition of P. --------------------------------------------
CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
!----------------------- End of Subroutine DPRJA -----------------------
END

DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
!-----------------------------------------------------------------------
! This function routine computes the weighted max-norm
! of the vector of length N contained in the array V, with weights
! contained in the array w of length N:
!   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
!-----------------------------------------------------------------------
INTEGER N,   I
DOUBLE PRECISION V, W,   VM
DIMENSION V(N), W(N)
VM = 0.0D0
DO 10 I = 1,N
10     VM = MAX(VM,ABS(V(I))*W(I))
DMNORM = VM
RETURN
!----------------------- End of Function DMNORM ------------------------
END

DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
!-----------------------------------------------------------------------
! This function computes the norm of a full N by N matrix,
! stored in the array A, that is consistent with the weighted max-norm
! on vectors, with weights stored in the array W:
!   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!-----------------------------------------------------------------------
INTEGER N,   I, J
DOUBLE PRECISION A,   W, AN, SUM
DIMENSION A(N,N), W(N)
AN = 0.0D0
DO 20 I = 1,N
SUM = 0.0D0
DO 10 J = 1,N
10       SUM = SUM + ABS(A(I,J))/W(J)
AN = MAX(AN,SUM*W(I))
20     CONTINUE
DFNORM = AN
RETURN
!----------------------- End of Function DFNORM ------------------------
END

DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
!-----------------------------------------------------------------------
! This function computes the norm of a banded N by N matrix,
! stored in the array A, that is consistent with the weighted max-norm
! on vectors, with weights stored in the array W.
! ML and MU are the lower and upper half-bandwidths of the matrix.
! NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
! In terms of the matrix elements a(i,j), the norm is given by:
!   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!-----------------------------------------------------------------------
INTEGER N, NRA, ML, MU
INTEGER I, I1, JLO, JHI, J
DOUBLE PRECISION A, W
DOUBLE PRECISION AN, SUM
DIMENSION A(NRA,N), W(N)
AN = 0.0D0
DO 20 I = 1,N
SUM = 0.0D0
I1 = I + MU + 1
JLO = MAX(I-ML,1)
JHI = MIN(I+MU,N)
DO 10 J = JLO,JHI
10       SUM = SUM + ABS(A(I1-J,J))/W(J)
AN = MAX(AN,SUM*W(I))
20     CONTINUE
DBNORM = AN
RETURN
!----------------------- End of Function DBNORM ------------------------
END

SUBROUTINE DSRCMA (RSAV, ISAV, JOB)
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of
! the Common blocks DLS001, DLSA01, which are used
! internally by one or more ODEPACK solvers.
!
! RSAV = real array of length 240 or more.
! ISAV = integer array of length 46 or more.
! JOB  = flag indicating to save or restore the Common blocks:
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
INTEGER ISAV, JOB
INTEGER ILS, ILSA
INTEGER I, LENRLS, LENILS, LENRLA, LENILA
DOUBLE PRECISION RSAV
DOUBLE PRECISION RLS, RLSA
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS, LENRLA, LENILA
COMMON /DLS001/ RLS(218), ILS(37)
COMMON /DLSA01/ RLSA(22), ILSA(9)
DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
!
IF (JOB .EQ. 2) GOTO 100
DO 10 I = 1,LENRLS
10     RSAV(I) = RLS(I)
DO 15 I = 1,LENRLA
15     RSAV(LENRLS+I) = RLSA(I)
!
DO 20 I = 1,LENILS
20     ISAV(I) = ILS(I)
DO 25 I = 1,LENILA
25     ISAV(LENILS+I) = ILSA(I)
!
RETURN
!
100  CONTINUE
DO 110 I = 1,LENRLS
110     RLS(I) = RSAV(I)
DO 115 I = 1,LENRLA
115     RLSA(I) = RSAV(LENRLS+I)
!
DO 120 I = 1,LENILS
120     ILS(I) = ISAV(I)
DO 125 I = 1,LENILA
125     ILSA(I) = ISAV(LENILS+I)
!
RETURN
!----------------------- End of Subroutine DSRCMA ----------------------
END

SUBROUTINE DRCHEK (JOB, G, NEQ, Y, YH,NYH, G0, G1, GX, JROOT, IRT)
EXTERNAL G
INTEGER JOB, NEQ, NYH, JROOT, IRT
DOUBLE PRECISION Y, YH, G0, G1, GX
DIMENSION NEQ(*), Y(*), YH(NYH,*), G0(*), G1(*), GX(*), JROOT(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IOWND3, IOWNR3, IRFND, ITASKC, NGC, NGE
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION ROWNR3, T0, TLAST, TOUTC
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC, IOWND3(3), IOWNR3(2), IRFND, ITASKC, NGC, NGE
INTEGER I, IFLAG, JFLAG
DOUBLE PRECISION HMING, T1, TEMP1, TEMP2, X
LOGICAL ZROOT
!-----------------------------------------------------------------------
! This routine checks for the presence of a root in the vicinity of
! the current T, in a manner depending on the input flag JOB.  It calls
! Subroutine DROOTS to locate the root as precisely as possible.
!
! In addition to variables described previously, DRCHEK
! uses the following for communication:
! JOB    = integer flag indicating type of call:
!          JOB = 1 means the problem is being initialized, and DRCHEK
!                  is to look for a root at or very near the initial T.
!          JOB = 2 means a continuation call to the solver was just
!                  made, and DRCHEK is to check for a root in the
!                  relevant part of the step last taken.
!          JOB = 3 means a successful step was just taken, and DRCHEK
!                  is to look for a root in the interval of the step.
! G0     = array of length NG, containing the value of g at T = T0.
!          G0 is input for JOB .ge. 2, and output in all cases.
! G1,GX  = arrays of length NG for work space.
! IRT    = completion flag:
!          IRT = 0  means no root was found.
!          IRT = -1 means JOB = 1 and a root was found too near to T.
!          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
!                   On return, T0 is the root location, and Y is the
!                   corresponding solution vector.
! T0     = value of T at one endpoint of interval of interest.  Only
!          roots beyond T0 in the direction of integration are sought.
!          T0 is input if JOB .ge. 2, and output in all cases.
!          T0 is updated by DRCHEK, whether a root is found or not.
! TLAST  = last value of T returned by the solver (input only).
! TOUTC  = copy of TOUT (input only).
! IRFND  = input flag showing whether the last step taken had a root.
!          IRFND = 1 if it did, = 0 if not.
! ITASKC = copy of ITASK (input only).
! NGC    = copy of NG (input only).
!-----------------------------------------------------------------------
IRT = 0
DO 10 I = 1,NGC
10     JROOT(I) = 0
HMING = (ABS(TN) + ABS(H))*UROUND*100.0D0
!
GOTO (100, 200, 300), JOB
!
! Evaluate g at initial T, and check for zero values. ------------------
100  CONTINUE
T0 = TN
CALL G (NEQ, T0, Y, NGC, G0)
NGE = 1
ZROOT = .FALSE.
DO 110 I = 1,NGC
110    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
IF (.NOT. ZROOT) GOTO 190
! g has a zero at T.  Look at g at T + (small increment). --------------
TEMP2 = MAX(HMING/ABS(H), 0.1D0)
TEMP1 = TEMP2*H
T0 = T0 + TEMP1
DO 120 I = 1,N
120    Y(I) = Y(I) + TEMP2*YH(I,2)
CALL G (NEQ, T0, Y, NGC, G0)
NGE = NGE + 1
ZROOT = .FALSE.
DO 130 I = 1,NGC
130    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
IF (.NOT. ZROOT) GOTO 190
! g has a zero at T and also close to T.  Take error return. -----------
IRT = -1
RETURN
!
190  CONTINUE
RETURN
!
!
200  CONTINUE
IF (IRFND .EQ. 0) GOTO 260
! If a root was found on the previous step, evaluate G0 = g(T0). -------
CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG)
CALL G (NEQ, T0, Y, NGC, G0)
NGE = NGE + 1
ZROOT = .FALSE.
DO 210 I = 1,NGC
210    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
IF (.NOT. ZROOT) GOTO 260
! g has a zero at T0.  Look at g at T + (small increment). -------------
TEMP1 = SIGN(HMING,H)
T0 = T0 + TEMP1
IF ((T0 - TN)*H .LT. 0.0D0) GOTO 230
TEMP2 = TEMP1/H
DO 220 I = 1,N
220    Y(I) = Y(I) + TEMP2*YH(I,2)
GOTO 240
230  CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG)
240  CALL G (NEQ, T0, Y, NGC, G0)
NGE = NGE + 1
ZROOT = .FALSE.
DO 250 I = 1,NGC
IF (ABS(G0(I)) .GT. 0.0D0) GOTO 250
JROOT(I) = 1
ZROOT = .TRUE.
250    CONTINUE
IF (.NOT. ZROOT) GOTO 260
! g has a zero at T0 and also close to T0.  Return root. ---------------
IRT = 1
RETURN
! G0 has no zero components.  Proceed to check relevant interval. ------
260  IF (TN .EQ. TLAST) GOTO 390
!
300  CONTINUE
! Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
IF (ITASKC.EQ.2 .OR. ITASKC.EQ.3 .OR. ITASKC.EQ.5) GOTO 310
IF ((TOUTC - TN)*H .GE. 0.0D0) GOTO 310
T1 = TOUTC
IF ((T1 - T0)*H .LE. 0.0D0) GOTO 390
CALL DINTDY (T1, 0, YH, NYH, Y, IFLAG)
GOTO 330
310  T1 = TN
DO 320 I = 1,N
320    Y(I) = YH(I,1)
330  CALL G (NEQ, T1, Y, NGC, G1)
NGE = NGE + 1
! Call DROOTS to search for root in interval from T0 to T1. ------------
JFLAG = 0
350  CONTINUE
CALL DROOTS (NGC, HMING, JFLAG, T0, T1, G0, G1, GX, X, JROOT)
IF (JFLAG .GT. 1) GOTO 360
CALL DINTDY (X, 0, YH, NYH, Y, IFLAG)
CALL G (NEQ, X, Y, NGC, GX)
NGE = NGE + 1
GOTO 350
360  T0 = X
CALL DCOPY (NGC, GX, 1, G0, 1)
IF (JFLAG .EQ. 4) GOTO 390
! Found a root.  Interpolate to X and return. --------------------------
CALL DINTDY (X, 0, YH, NYH, Y, IFLAG)
IRT = 1
RETURN
!
390  CONTINUE
RETURN
!----------------------- End of Subroutine DRCHEK ----------------------
END

SUBROUTINE DROOTS (NG, HMIN, JFLAG, X0, X1, G0, G1, GX, X, JROOT)
INTEGER NG, JFLAG, JROOT
DOUBLE PRECISION HMIN, X0, X1, G0, G1, GX, X
DIMENSION G0(NG), G1(NG), GX(NG), JROOT(NG)
INTEGER IOWND3, IMAX, LAST, IDUM3
DOUBLE PRECISION ALPHA, X2, RDUM3
COMMON /DLSR01/ ALPHA, X2, RDUM3(3), IOWND3(3), IMAX, LAST, IDUM3(4)
!-----------------------------------------------------------------------
! This subroutine finds the leftmost root of a set of arbitrary
! functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
! of odd multiplicity (i.e. changes of sign of the gi) are found.
! Here the sign of X1 - X0 is arbitrary, but is constant for a given
! problem, and -leftmost- means nearest to X0.
! The values of the vector-valued function g(x) = (gi, i=1...NG)
! are communicated through the call sequence of DROOTS.
! The method used is the Illinois algorithm.
!
! Reference:
! Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
! Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
! February 1980.
!
! Description of parameters.
!
! NG     = number of functions gi, or the number of components of
!          the vector valued function g(x).  Input only.
!
! HMIN   = resolution parameter in X.  Input only.  When a root is
!          found, it is located only to within an error of HMIN in X.
!          Typically, HMIN should be set to something on the order of
!               100 * UROUND * MAX(ABS(X0),ABS(X1)),
!          where UROUND is the unit roundoff of the machine.
!
! JFLAG  = integer flag for input and output communication.
!
!          On input, set JFLAG = 0 on the first call for the problem,
!          and leave it unchanged until the problem is completed.
!          (The problem is completed when JFLAG .ge. 2 on return.)
!
!          On output, JFLAG has the following values and meanings:
!          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
!                    and call DROOTS again.
!          JFLAG = 2 means a root has been found.  The root is
!                    at X, and GX contains g(X).  (Actually, X is the
!                    rightmost approximation to the root on an interval
!                    (X0,X1) of size HMIN or less.)
!          JFLAG = 3 means X = X1 is a root, with one or more of the gi
!                    being zero at X1 and no sign changes in (X0,X1).
!                    GX contains g(X) on output.
!          JFLAG = 4 means no roots (of odd multiplicity) were
!                    found in (X0,X1) (no sign changes).
!
! X0,X1  = endpoints of the interval where roots are sought.
!          X1 and X0 are input when JFLAG = 0 (first call), and
!          must be left unchanged between calls until the problem is
!          completed.  X0 and X1 must be distinct, but X1 - X0 may be
!          of either sign.  However, the notion of -left- and -right-
!          will be used to mean nearer to X0 or X1, respectively.
!          When JFLAG .ge. 2 on return, X0 and X1 are output, and
!          are the endpoints of the relevant interval.
!
! G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
!          respectively.  When JFLAG = 0, G0 and G1 are input and
!          none of the G0(i) should be zero.
!          When JFLAG .ge. 2 on return, G0 and G1 are output.
!
! GX     = array of length NG containing g(X).  GX is input
!          when JFLAG = 1, and output when JFLAG .ge. 2.
!
! X      = independent variable value.  Output only.
!          When JFLAG = 1 on output, X is the point at which g(x)
!          is to be evaluated and loaded into GX.
!          When JFLAG = 2 or 3, X is the root.
!          When JFLAG = 4, X is the right endpoint of the interval, X1.
!
! JROOT  = integer array of length NG.  Output only.
!          When JFLAG = 2 or 3, JROOT indicates which components
!          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
!          component has a root, and JROOT(i) = 0 otherwise.
!-----------------------------------------------------------------------
INTEGER I, IMXOLD, NXLAST
DOUBLE PRECISION T2, TMAX, FRACINT, FRACSUB, ZERO,HALF,TENTH,FIVE
LOGICAL ZROOT, SGNCHG, XROOT
SAVE ZERO, HALF, TENTH, FIVE
DATA ZERO/0.0D0/, HALF/0.5D0/, TENTH/0.1D0/, FIVE/5.0D0/
!
IF (JFLAG .EQ. 1) GOTO 200
! JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
IMAX = 0
TMAX = ZERO
ZROOT = .FALSE.
DO 120 I = 1,NG
IF (ABS(G1(I)) .GT. ZERO) GOTO 110
ZROOT = .TRUE.
GOTO 120
! At this point, G0(i) has been checked and cannot be zero. ------------
110    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,G1(I))) GOTO 120
T2 = ABS(G1(I)/(G1(I)-G0(I)))
IF (T2 .LE. TMAX) GOTO 120
TMAX = T2
IMAX = I
120    CONTINUE
IF (IMAX .GT. 0) GOTO 130
SGNCHG = .FALSE.
GOTO 140
130  SGNCHG = .TRUE.
140  IF (.NOT. SGNCHG) GOTO 400
! There is a sign change.  Find the first root in the interval. --------
XROOT = .FALSE.
NXLAST = 0
LAST = 1
!
! Repeat until the first root in the interval is found.  Loop point. ---
150  CONTINUE
IF (XROOT) GOTO 300
IF (NXLAST .EQ. LAST) GOTO 160
ALPHA = 1.0D0
GOTO 180
160  IF (LAST .EQ. 0) GOTO 170
ALPHA = 0.5D0*ALPHA
GOTO 180
170  ALPHA = 2.0D0*ALPHA
180  X2 = X1 - (X1 - X0)*G1(IMAX) / (G1(IMAX) - ALPHA*G0(IMAX))
! If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
! distance that is between 0.1 and 0.5. --------------------------------
IF (ABS(X2 - X0) < HALF*HMIN) THEN
FRACINT = ABS(X1 - X0)/HMIN
FRACSUB = TENTH
IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT
X2 = X0 + FRACSUB*(X1 - X0)
ENDIF
IF (ABS(X1 - X2) < HALF*HMIN) THEN
FRACINT = ABS(X1 - X0)/HMIN
FRACSUB = TENTH
IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT
X2 = X1 - FRACSUB*(X1 - X0)
ENDIF
JFLAG = 1
X = X2
! Return to the calling routine to get a value of GX = g(X). -----------
RETURN
! Check to see in which interval g changes sign. -----------------------
200  IMXOLD = IMAX
IMAX = 0
TMAX = ZERO
ZROOT = .FALSE.
DO 220 I = 1,NG
IF (ABS(GX(I)) .GT. ZERO) GOTO 210
ZROOT = .TRUE.
GOTO 220
! Neither G0(i) nor GX(i) can be zero at this point. -------------------
210    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,GX(I))) GOTO 220
T2 = ABS(GX(I)/(GX(I) - G0(I)))
IF (T2 .LE. TMAX) GOTO 220
TMAX = T2
IMAX = I
220    CONTINUE
IF (IMAX .GT. 0) GOTO 230
SGNCHG = .FALSE.
IMAX = IMXOLD
GOTO 240
230  SGNCHG = .TRUE.
240  NXLAST = LAST
IF (.NOT. SGNCHG) GOTO 250
! Sign change between X0 and X2, so replace X1 with X2. ----------------
X1 = X2
CALL DCOPY (NG, GX, 1, G1, 1)
LAST = 1
XROOT = .FALSE.
GOTO 270
250  IF (.NOT. ZROOT) GOTO 260
! Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
X1 = X2
CALL DCOPY (NG, GX, 1, G1, 1)
XROOT = .TRUE.
GOTO 270
! No sign change between X0 and X2.  Replace X0 with X2. ---------------
260  CONTINUE
CALL DCOPY (NG, GX, 1, G0, 1)
X0 = X2
LAST = 0
XROOT = .FALSE.
270  IF (ABS(X1-X0) .LE. HMIN) XROOT = .TRUE.
GOTO 150
!
! Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
300  JFLAG = 2
X = X1
CALL DCOPY (NG, G1, 1, GX, 1)
DO 320 I = 1,NG
JROOT(I) = 0
IF (ABS(G1(I)) .GT. ZERO) GOTO 310
JROOT(I) = 1
GOTO 320
310    IF (SIGN(1.0D0,G0(I)) .NE. SIGN(1.0D0,G1(I))) JROOT(I) = 1
320    CONTINUE
RETURN
!
! No sign change in the interval.  Check for zero at right endpoint. ---
400  IF (.NOT. ZROOT) GOTO 420
!
! Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
X = X1
CALL DCOPY (NG, G1, 1, GX, 1)
DO 410 I = 1,NG
JROOT(I) = 0
IF (ABS(G1(I)) .LE. ZERO) JROOT (I) = 1
410  CONTINUE
JFLAG = 3
RETURN
!
! No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
420  CALL DCOPY (NG, G1, 1, GX, 1)
X = X1
JFLAG = 4
RETURN
!----------------------- End of Subroutine DROOTS ----------------------
END

SUBROUTINE DSRCAR (RSAV, ISAV, JOB)
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of
! the Common blocks DLS001, DLSA01, DLSR01, which are used
! internally by one or more ODEPACK solvers.
!
! RSAV = real array of length 245 or more.
! ISAV = integer array of length 55 or more.
! JOB  = flag indicating to save or restore the Common blocks:
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
INTEGER ISAV, JOB
INTEGER ILS, ILSA, ILSR
INTEGER I, IOFF, LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
DOUBLE PRECISION RSAV
DOUBLE PRECISION RLS, RLSA, RLSR
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
COMMON /DLS001/ RLS(218), ILS(37)
COMMON /DLSA01/ RLSA(22), ILSA(9)
COMMON /DLSR01/ RLSR(5), ILSR(9)
DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
DATA LENRLR/5/, LENILR/9/
!
IF (JOB .EQ. 2) GOTO 100
DO 10 I = 1,LENRLS
10     RSAV(I) = RLS(I)
DO 15 I = 1,LENRLA
15     RSAV(LENRLS+I) = RLSA(I)
IOFF = LENRLS + LENRLA
DO 20 I = 1,LENRLR
20     RSAV(IOFF+I) = RLSR(I)
!
DO 30 I = 1,LENILS
30     ISAV(I) = ILS(I)
DO 35 I = 1,LENILA
35     ISAV(LENILS+I) = ILSA(I)
IOFF = LENILS + LENILA
DO 40 I = 1,LENILR
40     ISAV(IOFF+I) = ILSR(I)
!
RETURN
!
100  CONTINUE
DO 110 I = 1,LENRLS
110     RLS(I) = RSAV(I)
DO 115 I = 1,LENRLA
115     RLSA(I) = RSAV(LENRLS+I)
IOFF = LENRLS + LENRLA
DO 120 I = 1,LENRLR
120     RLSR(I) = RSAV(IOFF+I)
!
DO 130 I = 1,LENILS
130     ILS(I) = ISAV(I)
DO 135 I = 1,LENILA
135     ILSA(I) = ISAV(LENILS+I)
IOFF = LENILS + LENILA
DO 140 I = 1,LENILR
140     ILSR(I) = ISAV(IOFF+I)
!
RETURN
!----------------------- End of Subroutine DSRCAR ----------------------
END

SUBROUTINE DSTODPK (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVX, ACOR, WM, IWM, F, JAC, PSOL)
EXTERNAL F, JAC, PSOL
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVX, ACOR, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), SAVX(*), ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM
INTEGER LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN, JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
!-----------------------------------------------------------------------
! DSTODPK performs one step of the integration of an initial value
! problem for a system of Ordinary Differential Equations.
!-----------------------------------------------------------------------
! The following changes were made to generate Subroutine DSTODPK
! from Subroutine DSTODE:
! 1. The array SAVX was added to the call sequence.
! 2. PJAC and SLVS were replaced by PSOL in the call sequence.
! 3. The Common block /DLPK01/ was added for communication.
! 4. The test constant EPCON is loaded into Common below statement
!    numbers 125 and 155, and used below statement 400.
! 5. The Newton iteration counter MNEWT is set below 220 and 400.
! 6. The call to PJAC was replaced with a call to DPKSET (fixed name),
!    with a longer call sequence, called depending on JACFLG.
! 7. The corrector residual is stored in SAVX (not Y) at 360,
!    and the solution vector is in SAVX in the 380 loop.
! 8. SLVS was renamed DSOLPK and includes NEQ, SAVX, EWT, F, and JAC.
!    SAVX was added because DSOLPK now needs Y and SAVF undisturbed.
! 9. The nonlinear convergence failure count NCFN is set at 430.
!-----------------------------------------------------------------------
! Note: DSTODPK is independent of the value of the iteration method
! indicator MITER, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTODPK is done with the following variables:
!
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to F and JAC.
! Y      = an array of length .ge. N used as the Y argument in
!          all calls to F and JAC.
! YH     = an NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = a constant integer .ge. N, the first dimension of YH.
! YH1    = a one-dimensional array occupying the same space as YH.
! EWT    = an array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = an array of working storage, of length N.
!          Also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! SAVX   = an array of working storage, of length N.
! ACOR   = a work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = real and integer work arrays associated with matrix
!          operations in chord iteration (MITER .ne. 0).
! CCMAX  = maximum relative change in H*EL0 before DPKSET is called.
! H      = the step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = the minimum absolute value of the step size H to be used.
! HMXI   = inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = the independent variable. TN is updated on each step taken.
! JSTART = an integer used for input only, with the following
!          values and meanings:
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of H, MAXORD,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings:
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in DPKSET or DSOLPK.
!          A return with KFLAG = -1 or -2 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = the maximum order of integration method to be allowed.
! MAXCOR = the maximum number of corrector iterations allowed.
! MSBP   = maximum number of steps between DPKSET calls (MITER .gt. 0).
! MXNCF  = maximum number of convergence failures allowed.
! METH/MITER = the method flags.  See description in driver.
! N      = the number of first-order differential equations.
!-----------------------------------------------------------------------
INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
!
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set at 2
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
MEO = METH
NSLP = 0
IPUP = MITER
IRET = 3
GOTO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MEO) GOTO 110
CALL DCFODE (METH, ELCO, TESCO)
MEO = METH
IF (NQ .GT. MAXORD) GOTO 120
IALTH = L
IRET = 1
GOTO 150
110  IF (NQ .LE. MAXORD) GOTO 160
120  NQ = MAXORD
L = LMAX
DO 125 I = 1,L
125    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
EPCON = CONIT*TESCO(2,NQ)
DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
EXDN = 1.0D0/L
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
RH = MIN(RHDN,1.0D0)
IREDO = 3
IF (H .EQ. HOLD) GOTO 170
RH = MIN(RH,ABS(H/HOLD))
H = HOLD
GOTO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140  CALL DCFODE (METH, ELCO, TESCO)
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
EPCON = CONIT*TESCO(2,NQ)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! The flag IPUP is set according to whether matrix data is involved
! (JACFLG .ne. 0) or not (JACFLG = 0), to trigger a call to DPKSET.
! IPUP is set to MITER when RC differs from 1 by more than CCMAX,
! and at least every MSBP steps, when JACFLG = 1.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
!-----------------------------------------------------------------------
200  IF (JACFLG .NE. 0) GOTO 202
IPUP = 0
CRATE = 0.7D0
GOTO 205
202  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
205  TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS-norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220  M = 0
MNEWT = 0
DO 230 I = 1,N
230    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
IF (IPUP .LE. 0) GOTO 250
!-----------------------------------------------------------------------
! If indicated, DPKSET is called to update any matrix data needed,
! before starting the corrector iteration.
! IPUP is set to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL DPKSET (NEQ, Y, YH1, EWT, ACOR, SAVF, WM, IWM, F, JAC)
IPUP = 0
RC = 1.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .NE. 0) GOTO 430
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
270  IF (MITER .NE. 0) GOTO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
DO 290 I = 1,N
SAVF(I) = H*SAVF(I) - YH(I,2)
290    Y(I) = SAVF(I) - ACOR(I)
DEL = DVNORM (N, Y, EWT)
DO 300 I = 1,N
Y(I) = YH(I,1) + EL(1)*SAVF(I)
300    ACOR(I) = SAVF(I)
GOTO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
350  DO 360 I = 1,N
360    SAVX(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
CALL DSOLPK (NEQ, Y, SAVF, SAVX, EWT, WM, IWM, F, PSOL)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
DEL = DVNORM (N, SAVX, EWT)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + SAVX(I)
380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/EPCON
IF (DCON .LE. 1.0D0) GOTO 450
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
MNEWT = M
DELP = DEL
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
GOTO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, DPKSET is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410  IF (MITER.EQ.0 .OR. JCUR.EQ.1 .OR. JACFLG.EQ.0) GOTO 430
ICF = 1
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
NCFN = NCFN + 1
RMAX = 2.0D0
TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 680
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 670
IF (NCF .EQ. MXNCF) GOTO 670
RH = 0.5D0
IPUP = MITER
IREDO = 1
GOTO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450  JCUR = 0
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
HU = H
NQU = NQ
DO 470 J = 1,L
DO 470 I = 1,N
470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -3) GOTO 640
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! the largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 560
DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GOTO 610
IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640  IF (KFLAG .EQ. -10) GOTO 660
RH = 0.1D0
RH = MAX(HMIN/ABS(H),RH)
H = H*RH
DO 645 I = 1,N
645    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
DO 650 I = 1,N
650    YH(I,2) = H*SAVF(I)
IPUP = MITER
IALTH = 5
IF (NQ .EQ. 1) GOTO 200
NQ = 1
L = 2
IRET = 3
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -3
GOTO 720
690  RMAX = 10.0D0
700  R = 1.0D0/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- End of Subroutine DSTODPK ---------------------
END

SUBROUTINE DPKSET (NEQ, Y, YSV, EWT, FTEM, SAVF, WM, IWM, F, JAC)
EXTERNAL F, JAC
INTEGER NEQ, IWM
DOUBLE PRECISION Y, YSV, EWT, FTEM, SAVF, WM
DIMENSION NEQ(*), Y(*), YSV(*), EWT(*), FTEM(*), SAVF(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN, JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
!-----------------------------------------------------------------------
! DPKSET is called by DSTODPK to interface with the user-supplied
! routine JAC, to compute and process relevant parts of
! the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
! as need for preconditioning matrix operations later.
!
! In addition to variables described previously, communication
! with DPKSET uses the following:
! Y     = array containing predicted values on entry.
! YSV   = array containing predicted y, to be saved (YH1 in DSTODPK).
! FTEM  = work array of length N (ACOR in DSTODPK).
! SAVF  = array containing f evaluated at predicted y.
! WM    = real work space for matrices.
!         Space for preconditioning data starts at WM(LOCWP).
! IWM   = integer work space.
!         Space for preconditioning data starts at IWM(LOCIWP).
! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!         JAC returned an error flag.
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
!-----------------------------------------------------------------------
INTEGER IER
DOUBLE PRECISION HL0
!
IERPJ = 0
JCUR = 1
HL0 = EL0*H
CALL JAC (F, NEQ, TN, Y, YSV, EWT, SAVF, FTEM, HL0, WM(LOCWP), IWM(LOCIWP), IER)
NJE = NJE + 1
IF (IER .EQ. 0) RETURN
IERPJ = 1
RETURN
!----------------------- End of Subroutine DPKSET ----------------------
END

SUBROUTINE DSOLPK (NEQ, Y, SAVF, X, EWT, WM, IWM, F, PSOL)
EXTERNAL F, PSOL
INTEGER NEQ, IWM
DOUBLE PRECISION Y, SAVF, X, EWT, WM
DIMENSION NEQ(*), Y(*), SAVF(*), X(*), EWT(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN, JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
!-----------------------------------------------------------------------
! This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
! DUSOL, for the solution of the linear system arising from a Newton
! iteration.  It is called if MITER .ne. 0.
! In addition to variables described elsewhere,
! communication with DSOLPK uses the following variables:
! WM    = real work space containing data for the algorithm
!         (Krylov basis vectors, Hessenberg matrix, etc.)
! IWM   = integer work space containing data for the algorithm
! X     = the right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = output flag (in Common):
!         IERSL =  0 means no trouble occurred.
!         IERSL =  1 means the iterative method failed to converge.
!                    If the preconditioner is out of date, the step
!                    is repeated with a new preconditioner.
!                    Otherwise, the stepsize is reduced (forcing a
!                    new evaluation of the preconditioner) and the
!                    step is repeated.
!         IERSL = -1 means there was a nonrecoverable error in the
!                    iterative solver, and an error exit occurs.
! This routine also uses the Common variables TN, EL0, H, N, MITER,
! DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
! LOCWP, LOCIWP.
!-----------------------------------------------------------------------
INTEGER IFLAG, LB, LDL, LHES, LIOM, LGMR, LPCG, LP, LQ, LR, LV, LW, LWK, LZ, MAXLP1, NPSL
DOUBLE PRECISION DELTA, HL0
!
IERSL = 0
HL0 = H*EL0
DELTA = DELT*EPCON
GOTO (100, 200, 300, 400, 900, 900, 900, 900, 900), MITER
!-----------------------------------------------------------------------
! Use the SPIOM algorithm to solve the linear system P*x = -f.
!-----------------------------------------------------------------------
100  CONTINUE
LV = 1
LB = LV + N*MAXL
LHES = LB + N
LWK = LHES + MAXL*MAXL
CALL DCOPY (N, X, 1, WM(LB), 1)
CALL DSCAL (N, RSQRTN, EWT, 1)
CALL DSPIOM (NEQ, TN, Y, SAVF, WM(LB), EWT, N, MAXL, KMP, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), WM(LHES), IWM,&
LIOM, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
NNI = NNI + 1
NLI = NLI + LIOM
NPS = NPS + NPSL
CALL DSCAL (N, SQRTN, EWT, 1)
IF (IFLAG .NE. 0) NCFL = NCFL + 1
IF (IFLAG .GE. 2) IERSL = 1
IF (IFLAG .LT. 0) IERSL = -1
RETURN
!-----------------------------------------------------------------------
! Use the SPIGMR algorithm to solve the linear system P*x = -f.
!-----------------------------------------------------------------------
200  CONTINUE
MAXLP1 = MAXL + 1
LV = 1
LB = LV + N*MAXL
LHES = LB + N + 1
LQ = LHES + MAXL*MAXLP1
LWK = LQ + 2*MAXL
LDL = LWK + MIN(1,MAXL-KMP)*N
CALL DCOPY (N, X, 1, WM(LB), 1)
CALL DSCAL (N, RSQRTN, EWT, 1)
CALL DSPIGMR (NEQ, TN, Y, SAVF, WM(LB), EWT, N, MAXL, MAXLP1, KMP, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), &
WM(LHES), WM(LQ), LGMR, WM(LOCWP), IWM(LOCIWP), WM(LWK), WM(LDL), IFLAG)
NNI = NNI + 1
NLI = NLI + LGMR
NPS = NPS + NPSL
CALL DSCAL (N, SQRTN, EWT, 1)
IF (IFLAG .NE. 0) NCFL = NCFL + 1
IF (IFLAG .GE. 2) IERSL = 1
IF (IFLAG .LT. 0) IERSL = -1
RETURN
!-----------------------------------------------------------------------
! Use DPCG to solve the linear system P*x = -f
!-----------------------------------------------------------------------
300  CONTINUE
LR = 1
LP = LR + N
LW = LP + N
LZ = LW + N
LWK = LZ + N
CALL DCOPY (N, X, 1, WM(LR), 1)
CALL DPCG (NEQ, TN, Y, SAVF, WM(LR), EWT, N, MAXL, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LP), WM(LW), WM(LZ), LPCG,&
WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
NNI = NNI + 1
NLI = NLI + LPCG
NPS = NPS + NPSL
IF (IFLAG .NE. 0) NCFL = NCFL + 1
IF (IFLAG .GE. 2) IERSL = 1
IF (IFLAG .LT. 0) IERSL = -1
RETURN
!-----------------------------------------------------------------------
! Use DPCGS to solve the linear system P*x = -f
!-----------------------------------------------------------------------
400  CONTINUE
LR = 1
LP = LR + N
LW = LP + N
LZ = LW + N
LWK = LZ + N
CALL DCOPY (N, X, 1, WM(LR), 1)
CALL DPCGS (NEQ, TN, Y, SAVF, WM(LR), EWT, N, MAXL, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LP), WM(LW), WM(LZ), &
LPCG, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
NNI = NNI + 1
NLI = NLI + LPCG
NPS = NPS + NPSL
IF (IFLAG .NE. 0) NCFL = NCFL + 1
IF (IFLAG .GE. 2) IERSL = 1
IF (IFLAG .LT. 0) IERSL = -1
RETURN
!-----------------------------------------------------------------------
! Use DUSOL, which interfaces to PSOL, to solve the linear system
! (no Krylov iteration).
!-----------------------------------------------------------------------
900  CONTINUE
LB = 1
LWK = LB + N
CALL DCOPY (N, X, 1, WM(LB), 1)
CALL DUSOL (NEQ, TN, Y, SAVF, WM(LB), EWT, N, DELTA, HL0, MNEWT, PSOL, NPSL, X, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
NNI = NNI + 1
NPS = NPS + NPSL
IF (IFLAG .NE. 0) NCFL = NCFL + 1
IF (IFLAG .EQ. 3) IERSL = 1
IF (IFLAG .LT. 0) IERSL = -1
RETURN
!----------------------- End of Subroutine DSOLPK ----------------------
END

SUBROUTINE DSPIOM (NEQ, TN, Y, SAVF, B, WGHT, N, MAXL, KMP, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, IPVT, LIOM, &
WP, IWP, WK, IFLAG)
EXTERNAL F, PSOL
INTEGER NEQ,N,MAXL,KMP,JPRE,MNEWT,NPSL,IPVT,LIOM,IWP,IFLAG
DOUBLE PRECISION TN,Y,SAVF,B,WGHT,DELTA,HL0,X,V,HES,WP,WK
DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*), HES(MAXL,MAXL), IPVT(*), WP(*), IWP(*), WK(*)
!-----------------------------------------------------------------------
! This routine solves the linear system A * x = b using a scaled
! preconditioned version of the Incomplete Orthogonalization Method.
! An initial guess of x = 0 is assumed.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!         B    = the right hand side of the system A*x = b.
!                B is also used as work space when computing the
!                final approximation.
!                (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the diagonal
!                scaling matrix D.
!
!         N    = the order of the matrix A, and the lengths
!                of the vectors Y, SAVF, B, WGHT, and X.
!
!         MAXL = the maximum allowable order of the matrix HES.
!
!          KMP = the number of previous vectors the new vector VNEW
!                must be made orthogonal to.  KMP .le. MAXL.
!
!        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array of length N used by DATV and PSOL.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         V    = the N by (LIOM+1) array containing the LIOM
!                orthogonal vectors V(*,1) to V(*,LIOM).
!
!         HES  = the LU factorization of the LIOM by LIOM upper
!                Hessenberg matrix whose entries are the
!                scaled inner products of A*V(*,k) and V(*,i).
!
!         IPVT = an integer array containg pivoting information.
!                It is loaded in DHEFA and used in DHESL.
!
!         LIOM = the number of iterations performed, and current
!                order of the upper Hessenberg matrix HES.
!
!         NPSL = the number of calls to PSOL.
!
!        IFLAG = integer error flag:
!                0 means convergence in LIOM iterations, LIOM.le.MAXL.
!                1 means the convergence test did not pass in MAXL
!                  iterations, but the residual norm is .lt. 1,
!                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
!                2 means the convergence test did not pass in MAXL
!                  iterations, residual .gt. 1, and X is undefined.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
INTEGER I, IER, INFO, J, K, LL, LM1
DOUBLE PRECISION BNRM, BNRM0, PROD, RHO, SNORMW, DNRM2, TEM
!
IFLAG = 0
LIOM = 0
NPSL = 0
!-----------------------------------------------------------------------
! The initial residual is the vector b.  Apply scaling to b, and test
! for an immediate return with X = 0 or X = b.
!-----------------------------------------------------------------------
DO 10 I = 1,N
10     V(I,1) = B(I)*WGHT(I)
BNRM0 = DNRM2 (N, V, 1)
BNRM = BNRM0
IF (BNRM0 .GT. DELTA) GOTO 30
IF (MNEWT .GT. 0) GOTO 20
CALL DCOPY (N, B, 1, X, 1)
RETURN
20   DO 25 I = 1,N
25     X(I) = 0.0D0
RETURN
30   CONTINUE
! Apply inverse of left preconditioner to vector b. --------------------
IER = 0
IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GOTO 55
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 1, IER)
NPSL = 1
IF (IER .NE. 0) GOTO 300
! Calculate norm of scaled vector V(*,1) and normalize it. -------------
DO 50 I = 1,N
50     V(I,1) = B(I)*WGHT(I)
BNRM = DNRM2(N, V, 1)
DELTA = DELTA*(BNRM/BNRM0)
55   TEM = 1.0D0/BNRM
CALL DSCAL (N, TEM, V(1,1), 1)
! Zero out the HES array. ----------------------------------------------
DO 65 J = 1,MAXL
DO 60 I = 1,MAXL
60       HES(I,J) = 0.0D0
65     CONTINUE
!-----------------------------------------------------------------------
! Main loop on LL = l to compute the vectors V(*,2) to V(*,MAXL).
! The running product PROD is needed for the convergence test.
!-----------------------------------------------------------------------
PROD = 1.0D0
DO 90 LL = 1,MAXL
LIOM = LL
!-----------------------------------------------------------------------
! Call routine DATV to compute VNEW = Abar*v(l), where Abar is
! the matrix A with scaling and inverse preconditioner factors applied.
! Call routine DORTHOG to orthogonalize the new vector vnew = V(*,l+1).
! Call routine DHEFA to update the factors of HES.
!-----------------------------------------------------------------------
CALL DATV (NEQ, Y, SAVF, V(1,LL), WGHT, X, F, PSOL, V(1,LL+1), WK, WP, IWP, HL0, JPRE, IER, NPSL)
IF (IER .NE. 0) GOTO 300
CALL DORTHOG (V(1,LL+1), V, HES, N, LL, MAXL, KMP, SNORMW)
CALL DHEFA (HES, MAXL, LL, IPVT, INFO, LL)
LM1 = LL - 1
IF (LL .GT. 1 .AND. IPVT(LM1) .EQ. LM1) PROD = PROD*HES(LL,LM1)
IF (INFO .NE. LL) GOTO 70
!-----------------------------------------------------------------------
! The last pivot in HES was found to be zero.
! If vnew = 0 or l = MAXL, take an error return with IFLAG = 2.
! otherwise, continue the iteration without a convergence test.
!-----------------------------------------------------------------------
IF (SNORMW .EQ. 0.0D0) GOTO 120
IF (LL .EQ. MAXL) GOTO 120
GOTO 80
!-----------------------------------------------------------------------
! Update RHO, the estimate of the norm of the residual b - A*x(l).
! test for convergence.  If passed, compute approximation x(l).
! If failed and l .lt. MAXL, then continue iterating.
!-----------------------------------------------------------------------
70     CONTINUE
RHO = BNRM*SNORMW*ABS(PROD/HES(LL,LL))
IF (RHO .LE. DELTA) GOTO 200
IF (LL .EQ. MAXL) GOTO 100
! If l .lt. MAXL, store HES(l+1,l) and normalize the vector v(*,l+1).
80     CONTINUE
HES(LL+1,LL) = SNORMW
TEM = 1.0D0/SNORMW
CALL DSCAL (N, TEM, V(1,LL+1), 1)
90     CONTINUE
!-----------------------------------------------------------------------
! l has reached MAXL without passing the convergence test:
! If RHO is not too large, compute a solution anyway and return with
! IFLAG = 1.  Otherwise return with IFLAG = 2.
!-----------------------------------------------------------------------
100  CONTINUE
IF (RHO .LE. 1.0D0) GOTO 150
IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GOTO 150
120  CONTINUE
IFLAG = 2
RETURN
150  IFLAG = 1
!-----------------------------------------------------------------------
! Compute the approximation x(l) to the solution.
! Since the vector X was used as work space, and the initial guess
! of the Newton correction is zero, X must be reset to zero.
!-----------------------------------------------------------------------
200  CONTINUE
LL = LIOM
DO 210 K = 1,LL
210    B(K) = 0.0D0
B(1) = BNRM
CALL DHESL (HES, MAXL, LL, IPVT, B)
DO 220 K = 1,N
220    X(K) = 0.0D0
DO 230 I = 1,LL
CALL DAXPY (N, B(I), V(1,I), 1, X, 1)
230    CONTINUE
DO 240 I = 1,N
240    X(I) = X(I)/WGHT(I)
IF (JPRE .LE. 1) RETURN
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, X, 2, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) GOTO 300
RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
300  CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
RETURN
!----------------------- End of Subroutine DSPIOM ----------------------
END

SUBROUTINE DATV (NEQ, Y, SAVF, V, WGHT, FTEM, F, PSOL, Z, VTEM, WP, IWP, HL0, JPRE, IER, NPSL)
EXTERNAL F, PSOL
INTEGER NEQ, IWP, JPRE, IER, NPSL
DOUBLE PRECISION Y, SAVF, V, WGHT, FTEM, Z, VTEM, WP, HL0
DIMENSION NEQ(*), Y(*), SAVF(*), V(*), WGHT(*), FTEM(*), Z(*), VTEM(*), WP(*), IWP(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!-----------------------------------------------------------------------
! This routine computes the product
!
!   (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
!
! where D is a diagonal scaling matrix, and P1 and P2 are the
! left and right preconditioning matrices, respectively.
! v is assumed to have WRMS norm equal to 1.
! The product is stored in z.  This is computed by a
! difference quotient, a call to F, and two calls to PSOL.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            V = real array of length N (can be the same array as Z).
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the matrix D.
!
!         FTEM = work array of length N.
!
!         VTEM = work array of length N used to store the
!                unscaled version of V.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!         JPRE = preconditioner type flag.
!
!
!      On return
!
!            Z = array of length N containing desired scaled
!                matrix-vector product.
!
!          IER = error flag from PSOL.
!
!         NPSL = the number of calls to PSOL.
!
! In addition, this routine uses the Common variables TN, N, NFE.
!-----------------------------------------------------------------------
INTEGER I
DOUBLE PRECISION FAC, RNORM, DNRM2, TEMPN
!
! Set VTEM = D * V.
DO 10 I = 1,N
10     VTEM(I) = V(I)/WGHT(I)
IER = 0
IF (JPRE .GE. 2) GOTO 30
!
! JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
CALL DCOPY (N, Y, 1, Z, 1)
DO 20 I = 1,N
20     Y(I) = Z(I) + VTEM(I)
FAC = HL0
GOTO 60
!
! JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
30   CONTINUE
CALL PSOL (NEQ, TN, Y, SAVF, FTEM, HL0, WP, IWP, VTEM, 2, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) RETURN
! Calculate L-2 norm of (D-inverse) * VTEM.
DO 40 I = 1,N
40     Z(I) = VTEM(I)*WGHT(I)
TEMPN = DNRM2 (N, Z, 1)
RNORM = 1.0D0/TEMPN
! Save Y in Z and increment Y by VTEM/norm.
CALL DCOPY (N, Y, 1, Z, 1)
DO 50 I = 1,N
50     Y(I) = Z(I) + VTEM(I)*RNORM
FAC = HL0*TEMPN
!
! For all JPRE, call F with incremented Y argument, and restore Y.
60   CONTINUE
CALL F (NEQ, TN, Y, FTEM)
NFE = NFE + 1
CALL DCOPY (N, Z, 1, Y, 1)
! Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
DO 70 I = 1,N
70     Z(I) = FTEM(I) - SAVF(I)
DO 80 I = 1,N
80     Z(I) = VTEM(I) - FAC*Z(I)
! Apply inverse of left preconditioner to Z, if nontrivial.
IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GOTO 85
CALL PSOL (NEQ, TN, Y, SAVF, FTEM, HL0, WP, IWP, Z, 1, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) RETURN
85   CONTINUE
! Apply D-inverse to Z and return.
DO 90 I = 1,N
90     Z(I) = Z(I)*WGHT(I)
RETURN
!----------------------- End of Subroutine DATV ------------------------
END

SUBROUTINE DORTHOG (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
INTEGER N, LL, LDHES, KMP
DOUBLE PRECISION VNEW, V, HES, SNORMW
DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
!-----------------------------------------------------------------------
! This routine orthogonalizes the vector VNEW against the previous
! KMP vectors in the V array.  It uses a modified Gram-Schmidt
! orthogonalization procedure with conditional reorthogonalization.
! This is the version of 28 may 1986.
!-----------------------------------------------------------------------
!
!      On entry
!
!         VNEW = the vector of length N containing a scaled product
!                of the Jacobian and the vector V(*,LL).
!
!         V    = the N x l array containing the previous LL
!                orthogonal vectors v(*,1) to v(*,LL).
!
!         HES  = an LL x LL upper Hessenberg matrix containing,
!                in HES(i,k), k.lt.LL, scaled inner products of
!                A*V(*,k) and V(*,i).
!
!        LDHES = the leading dimension of the HES array.
!
!         N    = the order of the matrix A, and the length of VNEW.
!
!         LL   = the current order of the matrix HES.
!
!          KMP = the number of previous vectors the new vector VNEW
!                must be made orthogonal to (KMP .le. MAXL).
!
!
!      On return
!
!         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
!                where i0 = MAX(1, LL-KMP+1).
!
!         HES  = upper Hessenberg matrix with column LL filled in with
!                scaled inner products of A*V(*,LL) and V(*,i).
!
!       SNORMW = L-2 norm of VNEW.
!
!-----------------------------------------------------------------------
INTEGER I, I0
DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
!
! Get norm of unaltered VNEW for later use. ----------------------------
VNRM = DNRM2 (N, VNEW, 1)
!-----------------------------------------------------------------------
! Do modified Gram-Schmidt on VNEW = A*v(LL).
! Scaled inner products give new column of HES.
! Projections of earlier vectors are subtracted from VNEW.
!-----------------------------------------------------------------------
I0 = MAX(1,LL-KMP+1)
DO 10 I = I0,LL
HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
TEM = -HES(I,LL)
CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
10     CONTINUE
!-----------------------------------------------------------------------
! Compute SNORMW = norm of VNEW.
! If VNEW is small compared to its input value (in norm), then
! reorthogonalize VNEW to V(*,1) through V(*,LL).
! Correct if relative correction exceeds 1000*(unit roundoff).
! finally, correct SNORMW using the dot products involved.
!-----------------------------------------------------------------------
SNORMW = DNRM2 (N, VNEW, 1)
IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
SUMDSQ = 0.0D0
DO 30 I = I0,LL
TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GOTO 30
HES(I,LL) = HES(I,LL) - TEM
CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
SUMDSQ = SUMDSQ + TEM**2
30     CONTINUE
IF (SUMDSQ .EQ. 0.0D0) RETURN
ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
SNORMW = SQRT(ARG)
!
RETURN
!----------------------- End of Subroutine DORTHOG ---------------------
END

SUBROUTINE DSPIGMR (NEQ, TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1, KMP, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, Q, &
LGMR, WP, IWP, WK, DL, IFLAG)
EXTERNAL F, PSOL
INTEGER NEQ,N,MAXL,MAXLP1,KMP,JPRE,MNEWT,NPSL,LGMR,IWP,IFLAG
DOUBLE PRECISION TN,Y,SAVF,B,WGHT,DELTA,HL0,X,V,HES,Q,WP,WK,DL
DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*), HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*)
!-----------------------------------------------------------------------
! This routine solves the linear system A * x = b using a scaled
! preconditioned version of the Generalized Minimal Residual method.
! An initial guess of x = 0 is assumed.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            B = the right hand side of the system A*x = b.
!                B is also used as work space when computing
!                the final approximation.
!                (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
!
!         WGHT = the vector of length N containing the nonzero
!                elements of the diagonal scaling matrix.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors WGHT, B and X.
!
!         MAXL = the maximum allowable order of the matrix HES.
!
!       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
!
!          KMP = the number of previous vectors the new vector VNEW
!                must be made orthogonal to.  KMP .le. MAXL.
!
!        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by routine DATV and PSOL.
!
!           DL = real work array used for calculation of the residual
!                norm RHO when the method is incomplete (KMP .lt. MAXL).
!                Not needed or referenced in complete case (KMP = MAXL).
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         LGMR = the number of iterations performed and
!                the current order of the upper Hessenberg
!                matrix HES.
!
!         NPSL = the number of calls to PSOL.
!
!         V    = the N by (LGMR+1) array containing the LGMR
!                orthogonal vectors V(*,1) to V(*,LGMR).
!
!         HES  = the upper triangular factor of the QR decomposition
!                of the (LGMR+1) by lgmr upper Hessenberg matrix whose
!                entries are the scaled inner-products of A*V(*,i)
!                and V(*,k).
!
!         Q    = real array of length 2*MAXL containing the components
!                of the Givens rotations used in the QR decomposition
!                of HES.  It is loaded in DHEQR and used in DHELS.
!
!        IFLAG = integer error flag:
!                0 means convergence in LGMR iterations, LGMR .le. MAXL.
!                1 means the convergence test did not pass in MAXL
!                  iterations, but the residual norm is .lt. 1,
!                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
!                2 means the convergence test did not pass in MAXL
!                  iterations, residual .gt. 1, and X is undefined.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
DOUBLE PRECISION BNRM,BNRM0,C,DLNRM,PROD,RHO,S,SNORMW,DNRM2,TEM
!
IFLAG = 0
LGMR = 0
NPSL = 0
!-----------------------------------------------------------------------
! The initial residual is the vector b.  Apply scaling to b, and test
! for an immediate return with X = 0 or X = b.
!-----------------------------------------------------------------------
DO 10 I = 1,N
10     V(I,1) = B(I)*WGHT(I)
BNRM0 = DNRM2 (N, V, 1)
BNRM = BNRM0
IF (BNRM0 .GT. DELTA) GOTO 30
IF (MNEWT .GT. 0) GOTO 20
CALL DCOPY (N, B, 1, X, 1)
RETURN
20   DO 25 I = 1,N
25     X(I) = 0.0D0
RETURN
30   CONTINUE
! Apply inverse of left preconditioner to vector b. --------------------
IER = 0
IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GOTO 55
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 1, IER)
NPSL = 1
IF (IER .NE. 0) GOTO 300
! Calculate norm of scaled vector V(*,1) and normalize it. -------------
DO 50 I = 1,N
50     V(I,1) = B(I)*WGHT(I)
BNRM = DNRM2 (N, V, 1)
DELTA = DELTA*(BNRM/BNRM0)
55   TEM = 1.0D0/BNRM
CALL DSCAL (N, TEM, V(1,1), 1)
! Zero out the HES array. ----------------------------------------------
DO 65 J = 1,MAXL
DO 60 I = 1,MAXLP1
60       HES(I,J) = 0.0D0
65     CONTINUE
!-----------------------------------------------------------------------
! Main loop to compute the vectors V(*,2) to V(*,MAXL).
! The running product PROD is needed for the convergence test.
!-----------------------------------------------------------------------
PROD = 1.0D0
DO 90 LL = 1,MAXL
LGMR = LL
!-----------------------------------------------------------------------
! Call routine DATV to compute VNEW = Abar*v(ll), where Abar is
! the matrix A with scaling and inverse preconditioner factors applied.
! Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
! Call routine DHEQR to update the factors of HES.
!-----------------------------------------------------------------------
CALL DATV (NEQ, Y, SAVF, V(1,LL), WGHT, X, F, PSOL, V(1,LL+1), WK, WP, IWP, HL0, JPRE, IER, NPSL)
IF (IER .NE. 0) GOTO 300
CALL DORTHOG (V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
HES(LL+1,LL) = SNORMW
CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
IF (INFO .EQ. LL) GOTO 120
!-----------------------------------------------------------------------
! Update RHO, the estimate of the norm of the residual b - A*xl.
! If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
! necessarily orthogonal for LL .gt. KMP.  The vector DL must then
! be computed, and its norm used in the calculation of RHO.
!-----------------------------------------------------------------------
PROD = PROD*Q(2*LL)
RHO = ABS(PROD*BNRM)
IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
IF (LL .EQ. KMP+1) THEN
CALL DCOPY (N, V(1,1), 1, DL, 1)
DO 75 I = 1,KMP
IP1 = I + 1
I2 = I*2
S = Q(I2)
C = Q(I2-1)
DO 70 K = 1,N
70             DL(K) = S*DL(K) + C*V(K,IP1)
75           CONTINUE
ENDIF
S = Q(2*LL)
C = Q(2*LL-1)/SNORMW
LLP1 = LL + 1
DO 80 K = 1,N
80         DL(K) = S*DL(K) + C*V(K,LLP1)
DLNRM = DNRM2 (N, DL, 1)
RHO = RHO*DLNRM
ENDIF
!-----------------------------------------------------------------------
! Test for convergence.  If passed, compute approximation xl.
! if failed and LL .lt. MAXL, then continue iterating.
!-----------------------------------------------------------------------
IF (RHO .LE. DELTA) GOTO 200
IF (LL .EQ. MAXL) GOTO 100
!-----------------------------------------------------------------------
! Rescale so that the norm of V(1,LL+1) is one.
!-----------------------------------------------------------------------
TEM = 1.0D0/SNORMW
CALL DSCAL (N, TEM, V(1,LL+1), 1)
90     CONTINUE
100  CONTINUE
IF (RHO .LE. 1.0D0) GOTO 150
IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GOTO 150
120  CONTINUE
IFLAG = 2
RETURN
150  IFLAG = 1
!-----------------------------------------------------------------------
! Compute the approximation xl to the solution.
! Since the vector X was used as work space, and the initial guess
! of the Newton correction is zero, X must be reset to zero.
!-----------------------------------------------------------------------
200  CONTINUE
LL = LGMR
LLP1 = LL + 1
DO 210 K = 1,LLP1
210    B(K) = 0.0D0
B(1) = BNRM
CALL DHELS (HES, MAXLP1, LL, Q, B)
DO 220 K = 1,N
220    X(K) = 0.0D0
DO 230 I = 1,LL
CALL DAXPY (N, B(I), V(1,I), 1, X, 1)
230    CONTINUE
DO 240 I = 1,N
240    X(I) = X(I)/WGHT(I)
IF (JPRE .LE. 1) RETURN
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, X, 2, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) GOTO 300
RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
300  CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
!
RETURN
!----------------------- End of Subroutine DSPIGMR ---------------------
END

SUBROUTINE DPCG (NEQ, TN, Y, SAVF, R, WGHT, N, MAXL, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, P, W, Z, &
LPCG, WP, IWP, WK, IFLAG)
EXTERNAL F, PSOL
INTEGER NEQ, N, MAXL, JPRE, MNEWT, NPSL, LPCG, IWP, IFLAG
DOUBLE PRECISION TN,Y,SAVF,R,WGHT,DELTA,HL0,X,P,W,Z,WP,WK
DIMENSION NEQ(*), Y(*), SAVF(*), R(*), WGHT(*), X(*), P(*), W(*), Z(*), WP(*), IWP(*), WK(*)
!-----------------------------------------------------------------------
! This routine computes the solution to the system A*x = b using a
! preconditioned version of the Conjugate Gradient algorithm.
! It is assumed here that the matrix A and the preconditioner
! matrix M are symmetric positive definite or nearly so.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            R = the right hand side of the system A*x = b.
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the diagonal
!                scaling matrix D.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
!
!         MAXL = the maximum allowable number of iterates.
!
!        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by routine DATP.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         LPCG = the number of iterations performed, and current
!                order of the upper Hessenberg matrix HES.
!
!         NPSL = the number of calls to PSOL.
!
!        IFLAG = integer error flag:
!                0 means convergence in LPCG iterations, LPCG .le. MAXL.
!                1 means the convergence test did not pass in MAXL
!                  iterations, but the residual norm is .lt. 1,
!                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
!                2 means the convergence test did not pass in MAXL
!                  iterations, residual .gt. 1, and X is undefined.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!                4 means there was a zero denominator in the algorithm.
!                  The system matrix or preconditioner matrix is not
!                  sufficiently close to being symmetric pos. definite.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
INTEGER I, IER
DOUBLE PRECISION ALPHA,BETA,BNRM,PTW,RNRM,DDOT,DVNORM,ZTR,ZTR0
!
IFLAG = 0
NPSL = 0
LPCG = 0
DO 10 I = 1,N
10     X(I) = 0.0D0
BNRM = DVNORM (N, R, WGHT)
! Test for immediate return with X = 0 or X = b. -----------------------
IF (BNRM .GT. DELTA) GOTO 20
IF (MNEWT .GT. 0) RETURN
CALL DCOPY (N, R, 1, X, 1)
RETURN
!
20   ZTR = 0.0D0
! Loop point for PCG iterations. ---------------------------------------
30   CONTINUE
LPCG = LPCG + 1
CALL DCOPY (N, R, 1, Z, 1)
IER = 0
IF (JPRE .EQ. 0) GOTO 40
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, Z, 3, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) GOTO 100
40   CONTINUE
ZTR0 = ZTR
ZTR = DDOT (N, Z, 1, R, 1)
IF (LPCG .NE. 1) GOTO 50
CALL DCOPY (N, Z, 1, P, 1)
GOTO 70
50   CONTINUE
IF (ZTR0 .EQ. 0.0D0) GOTO 200
BETA = ZTR/ZTR0
DO 60 I = 1,N
60     P(I) = Z(I) + BETA*P(I)
70   CONTINUE
!-----------------------------------------------------------------------
!  Call DATP to compute A*p and return the answer in W.
!-----------------------------------------------------------------------
CALL DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
!
PTW = DDOT (N, P, 1, W, 1)
IF (PTW .EQ. 0.0D0) GOTO 200
ALPHA = ZTR/PTW
CALL DAXPY (N, ALPHA, P, 1, X, 1)
ALPHA = -ALPHA
CALL DAXPY (N, ALPHA, W, 1, R, 1)
RNRM = DVNORM (N, R, WGHT)
IF (RNRM .LE. DELTA) RETURN
IF (LPCG .LT. MAXL) GOTO 30
IFLAG = 2
IF (RNRM .LE. 1.0D0) IFLAG = 1
IF (RNRM .LE. BNRM .AND. MNEWT .EQ. 0) IFLAG = 1
RETURN
!-----------------------------------------------------------------------
! This block handles error returns from PSOL.
!-----------------------------------------------------------------------
100  CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
RETURN
!-----------------------------------------------------------------------
! This block handles division by zero errors.
!-----------------------------------------------------------------------
200  CONTINUE
IFLAG = 4
RETURN
!----------------------- End of Subroutine DPCG ------------------------
END

SUBROUTINE DPCGS (NEQ, TN, Y, SAVF, R, WGHT, N, MAXL, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, P, W, Z, LPCG, WP, &
IWP, WK, IFLAG)
EXTERNAL F, PSOL
INTEGER NEQ, N, MAXL, JPRE, MNEWT, NPSL, LPCG, IWP, IFLAG
DOUBLE PRECISION TN,Y,SAVF,R,WGHT,DELTA,HL0,X,P,W,Z,WP,WK
DIMENSION NEQ(*), Y(*), SAVF(*), R(*), WGHT(*), X(*), P(*), W(*), Z(*), WP(*), IWP(*), WK(*)
!-----------------------------------------------------------------------
! This routine computes the solution to the system A*x = b using a
! scaled preconditioned version of the Conjugate Gradient algorithm.
! It is assumed here that the scaled matrix D**-1 * A * D and the
! scaled preconditioner D**-1 * M * D are close to being
! symmetric positive definite.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            R = the right hand side of the system A*x = b.
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the diagonal
!                scaling matrix D.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
!
!         MAXL = the maximum allowable number of iterates.
!
!        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by routine DATP.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         LPCG = the number of iterations performed, and current
!                order of the upper Hessenberg matrix HES.
!
!         NPSL = the number of calls to PSOL.
!
!        IFLAG = integer error flag:
!                0 means convergence in LPCG iterations, LPCG .le. MAXL.
!                1 means the convergence test did not pass in MAXL
!                  iterations, but the residual norm is .lt. 1,
!                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
!                2 means the convergence test did not pass in MAXL
!                  iterations, residual .gt. 1, and X is undefined.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!                4 means there was a zero denominator in the algorithm.
!                  the scaled matrix or scaled preconditioner is not
!                  sufficiently close to being symmetric pos. definite.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
INTEGER I, IER
DOUBLE PRECISION ALPHA, BETA, BNRM, PTW, RNRM, DVNORM, ZTR, ZTR0
!
IFLAG = 0
NPSL = 0
LPCG = 0
DO 10 I = 1,N
10     X(I) = 0.0D0
BNRM = DVNORM (N, R, WGHT)
! Test for immediate return with X = 0 or X = b. -----------------------
IF (BNRM .GT. DELTA) GOTO 20
IF (MNEWT .GT. 0) RETURN
CALL DCOPY (N, R, 1, X, 1)
RETURN
!
20   ZTR = 0.0D0
! Loop point for PCG iterations. ---------------------------------------
30   CONTINUE
LPCG = LPCG + 1
CALL DCOPY (N, R, 1, Z, 1)
IER = 0
IF (JPRE .EQ. 0) GOTO 40
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, Z, 3, IER)
NPSL = NPSL + 1
IF (IER .NE. 0) GOTO 100
40   CONTINUE
ZTR0 = ZTR
ZTR = 0.0D0
DO 45 I = 1,N
45     ZTR = ZTR + Z(I)*R(I)*WGHT(I)**2
IF (LPCG .NE. 1) GOTO 50
CALL DCOPY (N, Z, 1, P, 1)
GOTO 70
50   CONTINUE
IF (ZTR0 .EQ. 0.0D0) GOTO 200
BETA = ZTR/ZTR0
DO 60 I = 1,N
60     P(I) = Z(I) + BETA*P(I)
70   CONTINUE
!-----------------------------------------------------------------------
!  Call DATP to compute A*p and return the answer in W.
!-----------------------------------------------------------------------
CALL DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
!
PTW = 0.0D0
DO 80 I = 1,N
80     PTW = PTW + P(I)*W(I)*WGHT(I)**2
IF (PTW .EQ. 0.0D0) GOTO 200
ALPHA = ZTR/PTW
CALL DAXPY (N, ALPHA, P, 1, X, 1)
ALPHA = -ALPHA
CALL DAXPY (N, ALPHA, W, 1, R, 1)
RNRM = DVNORM (N, R, WGHT)
IF (RNRM .LE. DELTA) RETURN
IF (LPCG .LT. MAXL) GOTO 30
IFLAG = 2
IF (RNRM .LE. 1.0D0) IFLAG = 1
IF (RNRM .LE. BNRM .AND. MNEWT .EQ. 0) IFLAG = 1
RETURN
!-----------------------------------------------------------------------
! This block handles error returns from PSOL.
!-----------------------------------------------------------------------
100  CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
RETURN
!-----------------------------------------------------------------------
! This block handles division by zero errors.
!-----------------------------------------------------------------------
200  CONTINUE
IFLAG = 4
RETURN
!----------------------- End of Subroutine DPCGS -----------------------
END

SUBROUTINE DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
EXTERNAL F
INTEGER NEQ
DOUBLE PRECISION Y, SAVF, P, WGHT, HL0, WK, W
DIMENSION NEQ(*), Y(*), SAVF(*), P(*), WGHT(*), WK(*), W(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, &
MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!-----------------------------------------------------------------------
! This routine computes the product
!
!              w = (I - hl0*df/dy)*p
!
! This is computed by a call to F and a difference quotient.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            P = real array of length N.
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the matrix D.
!
!           WK = work array of length N.
!
!      On return
!
!
!            W = array of length N containing desired
!                matrix-vector product.
!
! In addition, this routine uses the Common variables TN, N, NFE.
!-----------------------------------------------------------------------
INTEGER I
DOUBLE PRECISION FAC, PNRM, RPNRM, DVNORM
!
PNRM = DVNORM (N, P, WGHT)
RPNRM = 1.0D0/PNRM
CALL DCOPY (N, Y, 1, W, 1)
DO 20 I = 1,N
20     Y(I) = W(I) + P(I)*RPNRM
CALL F (NEQ, TN, Y, WK)
NFE = NFE + 1
CALL DCOPY (N, W, 1, Y, 1)
FAC = HL0*PNRM
DO 40 I = 1,N
40     W(I) = P(I) - FAC*(WK(I) - SAVF(I))
RETURN
!----------------------- End of Subroutine DATP ------------------------
END

SUBROUTINE DUSOL (NEQ, TN, Y, SAVF, B, WGHT, N, DELTA, HL0, MNEWT, PSOL, NPSL, X, WP, IWP, WK, IFLAG)
EXTERNAL PSOL
INTEGER NEQ, N, MNEWT, NPSL, IWP, IFLAG
DOUBLE PRECISION TN, Y, SAVF, B, WGHT, DELTA, HL0, X, WP, WK
DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*), WP(*), IWP(*), WK(*)
!-----------------------------------------------------------------------
! This routine solves the linear system A * x = b using only a call
! to the user-supplied routine PSOL (no Krylov iteration).
! If the norm of the right-hand side vector b is smaller than DELTA,
! the vector X returned is X = b (if MNEWT = 0) or X = 0 otherwise.
! PSOL is called with an LR argument of 0.
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            B = the right hand side of the system A*x = b.
!
!         WGHT = the vector of length N containing the nonzero
!                elements of the diagonal scaling matrix.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors WGHT, B and X.
!
!        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
!
!          HL0 = current value of (step size h) * (coefficient l0).
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by PSOL.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         NPSL = the number of calls to PSOL.
!
!        IFLAG = integer error flag:
!                0 means no trouble occurred.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
INTEGER I, IER
DOUBLE PRECISION BNRM, DVNORM
!
IFLAG = 0
NPSL = 0
!-----------------------------------------------------------------------
! Test for an immediate return with X = 0 or X = b.
!-----------------------------------------------------------------------
BNRM = DVNORM (N, B, WGHT)
IF (BNRM .GT. DELTA) GOTO 30
IF (MNEWT .GT. 0) GOTO 10
CALL DCOPY (N, B, 1, X, 1)
RETURN
10   DO 20 I = 1,N
20     X(I) = 0.0D0
RETURN
! Make call to PSOL and copy result from B to X. -----------------------
30   IER = 0
CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 0, IER)
NPSL = 1
IF (IER .NE. 0) GOTO 100
CALL DCOPY (N, B, 1, X, 1)
RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
100  CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
RETURN
!----------------------- End of Subroutine DUSOL -----------------------
END

SUBROUTINE DSRCPK (RSAV, ISAV, JOB)
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of
! the Common blocks DLS001, DLPK01, which are used
! internally by the DLSODPK solver.
!
! RSAV = real array of length 222 or more.
! ISAV = integer array of length 50 or more.
! JOB  = flag indicating to save or restore the Common blocks:
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
INTEGER ISAV, JOB
INTEGER ILS, ILSP
INTEGER I, LENILP, LENRLP, LENILS, LENRLS
DOUBLE PRECISION RSAV,   RLS, RLSP
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS, LENRLP, LENILP
COMMON /DLS001/ RLS(218), ILS(37)
COMMON /DLPK01/ RLSP(4), ILSP(13)
DATA LENRLS/218/, LENILS/37/, LENRLP/4/, LENILP/13/
!
IF (JOB .EQ. 2) GOTO 100
CALL DCOPY (LENRLS, RLS, 1, RSAV, 1)
CALL DCOPY (LENRLP, RLSP, 1, RSAV(LENRLS+1), 1)
DO 20 I = 1,LENILS
20     ISAV(I) = ILS(I)
DO 40 I = 1,LENILP
40     ISAV(LENILS+I) = ILSP(I)
RETURN
!
100  CONTINUE
CALL DCOPY (LENRLS, RSAV, 1, RLS, 1)
CALL DCOPY (LENRLP, RSAV(LENRLS+1), 1, RLSP, 1)
DO 120 I = 1,LENILS
120     ILS(I) = ISAV(I)
DO 140 I = 1,LENILP
140     ILSP(I) = ISAV(LENILS+I)
RETURN
!----------------------- End of Subroutine DSRCPK ----------------------
END

SUBROUTINE DHEFA (A, LDA, N, IPVT, INFO, JOB)
INTEGER LDA, N, IPVT(*), INFO, JOB
DOUBLE PRECISION A(LDA,*)
!-----------------------------------------------------------------------
!     This routine is a modification of the LINPACK routine DGEFA and
!     performs an LU decomposition of an upper Hessenberg matrix A.
!     There are two options available:
!
!          (1)  performing a fresh factorization
!          (2)  updating the LU factors by adding a row and a
!               column to the matrix A.
!-----------------------------------------------------------------------
!     DHEFA factors an upper Hessenberg matrix by elimination.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        JOB     INTEGER
!                JOB = 1    means that a fresh factorization of the
!                           matrix A is desired.
!                JOB .ge. 2 means that the current factorization of A
!                           will be updated by the addition of a row
!                           and a column.
!
!     On return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  U(k,k) .eq. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DHESL will divide by zero if called.
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 7/20/83.  This version dated 6/20/01.
!    
!     BLAS called: DAXPY, IDAMAX
!-----------------------------------------------------------------------
INTEGER IDAMAX, J, K, KM1, KP1, L, NM1
DOUBLE PRECISION T
!
IF (JOB .GT. 1) GOTO 80
!
! A new facorization is desired.  This is essentially the LINPACK
! code with the exception that we know there is only one nonzero
! element below the main diagonal.
!
!     Gaussian elimination with partial pivoting
!
INFO = 0
NM1 = N - 1
IF (NM1 .LT. 1) GOTO 70
DO 60 K = 1, NM1
KP1 = K + 1
!
!        Find L = pivot index
!
L = IDAMAX (2, A(K,K), 1) + K - 1
IPVT(K) = L
!
!        Zero pivot implies this column already triangularized
!
IF (A(L,K) .EQ. 0.0D0) GOTO 40
!
!           Interchange if necessary
!
IF (L .EQ. K) GOTO 10
T = A(L,K)
A(L,K) = A(K,K)
A(K,K) = T
10       CONTINUE
!
!           Compute multipliers
!
T = -1.0D0/A(K,K)
A(K+1,K) = A(K+1,K)*T
!
!           Row elimination with column indexing
!
DO 30 J = KP1, N
T = A(L,J)
IF (L .EQ. K) GOTO 20
A(L,J) = A(K,J)
A(K,J) = T
20          CONTINUE
CALL DAXPY (N-K, T, A(K+1,K), 1, A(K+1,J), 1)
30       CONTINUE
GOTO 50
40    CONTINUE
INFO = K
50    CONTINUE
60 CONTINUE
70 CONTINUE
IPVT(N) = N
IF (A(N,N) .EQ. 0.0D0) INFO = N
RETURN
!
! The old factorization of A will be updated.  A row and a column
! has been added to the matrix A.
! N-1 is now the old order of the matrix.
!
80  CONTINUE
NM1 = N - 1
!
! Perform row interchanges on the elements of the new column, and
! perform elimination operations on the elements using the multipliers.
!
IF (NM1 .LE. 1) GOTO 105
DO 100 K = 2,NM1
KM1 = K - 1
L = IPVT(KM1)
T = A(L,N)
IF (L .EQ. KM1) GOTO 90
A(L,N) = A(KM1,N)
A(KM1,N) = T
90    CONTINUE
A(K,N) = A(K,N) + A(K,KM1)*T
100    CONTINUE
105  CONTINUE
!
! Complete update of factorization by decomposing last 2x2 block.
!
INFO = 0
!
!        Find L = pivot index
!
L = IDAMAX (2, A(NM1,NM1), 1) + NM1 - 1
IPVT(NM1) = L
!
!        Zero pivot implies this column already triangularized
!
IF (A(L,NM1) .EQ. 0.0D0) GOTO 140
!
!           Interchange if necessary
!
IF (L .EQ. NM1) GOTO 110
T = A(L,NM1)
A(L,NM1) = A(NM1,NM1)
A(NM1,NM1) = T
110       CONTINUE
!
!           Compute multipliers
!
T = -1.0D0/A(NM1,NM1)
A(N,NM1) = A(N,NM1)*T
!
!           Row elimination with column indexing
!
T = A(L,N)
IF (L .EQ. NM1) GOTO 120
A(L,N) = A(NM1,N)
A(NM1,N) = T
120          CONTINUE
A(N,N) = A(N,N) + T*A(N,NM1)
GOTO 150
140    CONTINUE
INFO = NM1
150    CONTINUE
IPVT(N) = N
IF (A(N,N) .EQ. 0.0D0) INFO = N
RETURN
!----------------------- End of Subroutine DHEFA -----------------------
END

SUBROUTINE DHESL (A, LDA, N, IPVT, B)
INTEGER LDA, N, IPVT(*)
DOUBLE PRECISION A(LDA,*), B(*)
!-----------------------------------------------------------------------
! This is essentially the LINPACK routine DGESL except for changes
! due to the fact that A is an upper Hessenberg matrix.
!-----------------------------------------------------------------------
!     DHESL solves the real system A * x = b
!     using the factors computed by DHEFA.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DHEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DHEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!     On return
!
!        B       the solution vector  x .
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 7/20/83.  This version dated 6/20/01.
!
!     BLAS called: DAXPY
!-----------------------------------------------------------------------
INTEGER K, KB, L, NM1
DOUBLE PRECISION T
!
NM1 = N - 1
!
!        Solve  A * x = b
!        First solve  L*y = b
!
IF (NM1 .LT. 1) GOTO 30
DO 20 K = 1, NM1
L = IPVT(K)
T = B(L)
IF (L .EQ. K) GOTO 10
B(L) = B(K)
B(K) = T
10       CONTINUE
B(K+1) = B(K+1) + T*A(K+1,K)
20    CONTINUE
30    CONTINUE
!
!        Now solve  U*x = y
!
DO 40 KB = 1, N
K = N + 1 - KB
B(K) = B(K)/A(K,K)
T = -B(K)
CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
40    CONTINUE
RETURN
!----------------------- End of Subroutine DHESL -----------------------
END

SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
INTEGER LDA, N, INFO, IJOB
DOUBLE PRECISION A(LDA,*), Q(*)
!-----------------------------------------------------------------------
!     This routine performs a QR decomposition of an upper
!     Hessenberg matrix A.  There are two options available:
!
!          (1)  performing a fresh decomposition
!          (2)  updating the QR factors by adding a row and a
!               column to the matrix A.
!-----------------------------------------------------------------------
!     DHEQR decomposes an upper Hessenberg matrix by using Givens
!     rotations.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be decomposed.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                A is an (N+1) by N Hessenberg matrix.
!
!        IJOB    INTEGER
!                = 1     means that a fresh decomposition of the
!                        matrix A is desired.
!                .ge. 2  means that the current decomposition of A
!                        will be updated by the addition of a row
!                        and a column.
!     On return
!
!        A       the upper triangular matrix R.
!                The factorization can be written Q*A = R, where
!                Q is a product of Givens rotations and R is upper
!                triangular.
!
!        Q       DOUBLE PRECISION(2*N)
!                the factors c and s of each Givens rotation used
!                in decomposing A.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  A(k,k) .eq. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DHELS will divide by zero
!                     if called.
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 1/13/86.  This version dated 6/20/01.
!-----------------------------------------------------------------------
INTEGER I, IQ, J, K, KM1, KP1, NM1
DOUBLE PRECISION C, S, T, T1, T2
!
IF (IJOB .GT. 1) GOTO 70
!
! A new facorization is desired.
!
!     QR decomposition without pivoting
!
INFO = 0
DO 60 K = 1, N
KM1 = K - 1
KP1 = K + 1
!
!           Compute kth column of R.
!           First, multiply the kth column of A by the previous
!           k-1 Givens rotations.
!
IF (KM1 .LT. 1) GOTO 20
DO 10 J = 1, KM1
I = 2*(J-1) + 1
T1 = A(J,K)
T2 = A(J+1,K)
C = Q(I)
S = Q(I+1)
A(J,K) = C*T1 - S*T2
A(J+1,K) = S*T1 + C*T2
10         CONTINUE
!
!           Compute Givens components c and s
!
20       CONTINUE
IQ = 2*KM1 + 1
T1 = A(K,K)
T2 = A(KP1,K)
IF (T2 .NE. 0.0D0) GOTO 30
C = 1.0D0
S = 0.0D0
GOTO 50
30       CONTINUE
IF (ABS(T2) .LT. ABS(T1)) GOTO 40
T = T1/T2
S = -1.0D0/SQRT(1.0D0+T*T)
C = -S*T
GOTO 50
40       CONTINUE
T = T2/T1
C = 1.0D0/SQRT(1.0D0+T*T)
S = -C*T
50       CONTINUE
Q(IQ) = C
Q(IQ+1) = S
A(K,K) = C*T1 - S*T2
IF (A(K,K) .EQ. 0.0D0) INFO = K
60 CONTINUE
RETURN
!
! The old factorization of A will be updated.  A row and a column
! has been added to the matrix A.
! N by N-1 is now the old size of the matrix.
!
70  CONTINUE
NM1 = N - 1
!
! Multiply the new column by the N previous Givens rotations.
!
DO 100 K = 1,NM1
I = 2*(K-1) + 1
T1 = A(K,N)
T2 = A(K+1,N)
C = Q(I)
S = Q(I+1)
A(K,N) = C*T1 - S*T2
A(K+1,N) = S*T1 + C*T2
100    CONTINUE
!
! Complete update of decomposition by forming last Givens rotation,
! and multiplying it times the column vector (A(N,N), A(N+1,N)).
!
INFO = 0
T1 = A(N,N)
T2 = A(N+1,N)
IF (T2 .NE. 0.0D0) GOTO 110
C = 1.0D0
S = 0.0D0
GOTO 130
110  CONTINUE
IF (ABS(T2) .LT. ABS(T1)) GOTO 120
T = T1/T2
S = -1.0D0/SQRT(1.0D0+T*T)
C = -S*T
GOTO 130
120  CONTINUE
T = T2/T1
C = 1.0D0/SQRT(1.0D0+T*T)
S = -C*T
130  CONTINUE
IQ = 2*N - 1
Q(IQ) = C
Q(IQ+1) = S
A(N,N) = C*T1 - S*T2
IF (A(N,N) .EQ. 0.0D0) INFO = N
RETURN
!----------------------- End of Subroutine DHEQR -----------------------
END

SUBROUTINE DHELS (A, LDA, N, Q, B)
INTEGER LDA, N
DOUBLE PRECISION A(LDA,*), B(*), Q(*)
!-----------------------------------------------------------------------
! This is part of the LINPACK routine DGESL with changes
! due to the fact that A is an upper Hessenberg matrix.
!-----------------------------------------------------------------------
!     DHELS solves the least squares problem
!
!           min (b-A*x, b-A*x)
!
!     using the factors computed by DHEQR.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DHEQR which contains the upper
!                triangular factor R in the QR decomposition of A.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                A is originally an (N+1) by N matrix.
!
!        Q       DOUBLE PRECISION(2*N)
!                The coefficients of the N givens rotations
!                used in the QR factorization of A.
!
!        B       DOUBLE PRECISION(N+1)
!                the right hand side vector.
!
!     On return
!
!        B       the solution vector  x .
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 1/13/86.  This version dated 6/20/01.
!
!     BLAS called: DAXPY
!-----------------------------------------------------------------------
INTEGER IQ, K, KB, KP1
DOUBLE PRECISION C, S, T, T1, T2
!
!        Minimize (b-A*x, b-A*x)
!        First form Q*b.
!
DO 20 K = 1, N
KP1 = K + 1
IQ = 2*(K-1) + 1
C = Q(IQ)
S = Q(IQ+1)
T1 = B(K)
T2 = B(KP1)
B(K) = C*T1 - S*T2
B(KP1) = S*T1 + C*T2
20    CONTINUE
!
!        Now solve  R*x = Q*b.
!
DO 40 KB = 1, N
K = N + 1 - KB
B(K) = B(K)/A(K,K)
T = -B(K)
CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
40    CONTINUE
RETURN
!----------------------- End of Subroutine DHELS -----------------------
END

SUBROUTINE DLHIN (NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND, EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
EXTERNAL F
DOUBLE PRECISION T0, Y0, YDOT, TOUT, UROUND, EWT, ATOL, Y, TEMP, H0
INTEGER NEQ, N, ITOL, NITER, IER
DIMENSION NEQ(*), Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*), TEMP(*)
!-----------------------------------------------------------------------
! Call sequence input -- NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! Common block variables accessed -- None
!
! Subroutines called by DLHIN: F, DCOPY
! Function routines called by DLHIN: DVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define H from WRMS-norm(H**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with DLHIN is done with the following variables:
!
! NEQ    = NEQ array of solver, passed to F.
! N      = size of ODE system, input.
! T0     = initial value of independent variable, input.
! Y0     = vector of initial conditions, input.
! YDOT   = vector of initial first derivatives, input.
! F      = name of subroutine for right-hand side f(t,y), input.
! TOUT   = first output value of independent variable
! UROUND = machine unit roundoff
! EWT, ITOL, ATOL = error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = work arrays of length N.
! H0     = step size to be attempted, output.
! NITER  = number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = the error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and t0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION AFI, ATOLI, DELYI, HALF, HG, HLB, HNEW, HRAT, HUB, HUN, PT1, T1, TDIST, TROUND, TWO, DVNORM, YDDNRM
INTEGER I, ITER
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE HALF, HUN, PT1, TWO
DATA HALF /0.5D0/, HUN /100.0D0/, PT1 /0.1D0/, TWO /2.0D0/
!
NITER = 0
TDIST = ABS(TOUT - T0)
TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
IF (TDIST .LT. TWO*TROUND) GOTO 100
!
! Set a lower bound on H based on the roundoff level in T0 and TOUT. ---
HLB = HUN*TROUND
! Set an upper bound on H based on TOUT-T0 and the initial Y and YDOT. -
HUB = PT1*TDIST
ATOLI = ATOL(1)
DO 10 I = 1,N
IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
DELYI = PT1*ABS(Y0(I)) + ATOLI
AFI = ABS(YDOT(I))
IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
10     CONTINUE
!
! Set initial guess for H as geometric mean of upper and lower bounds. -
ITER = 0
HG = SQRT(HLB*HUB)
! If the bounds have crossed, exit with the mean value. ----------------
IF (HUB .LT. HLB) THEN
H0 = HG
GOTO 90
ENDIF
!
! Looping point for iteration. -----------------------------------------
50   CONTINUE
! Estimate the second derivative as a difference quotient in f. --------
T1 = T0 + HG
DO 60 I = 1,N
60     Y(I) = Y0(I) + HG*YDOT(I)
CALL F (NEQ, T1, Y, TEMP)
DO 70 I = 1,N
70     TEMP(I) = (TEMP(I) - YDOT(I))/HG
YDDNRM = DVNORM (N, TEMP, EWT)
! Get the corresponding new value of H. --------------------------------
IF (YDDNRM*HUB*HUB .GT. TWO) THEN
HNEW = SQRT(TWO/YDDNRM)
ELSE
HNEW = SQRT(HG*HUB)
ENDIF
ITER = ITER + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous H values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous H
! if hnew/hg .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
IF (ITER .GE. 4) GOTO 80
HRAT = HNEW/HG
IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GOTO 80
IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
HNEW = HG
GOTO 80
ENDIF
HG = HNEW
GOTO 50
!
! Iteration done.  Apply bounds, bias factor, and sign. ----------------
80   H0 = HNEW*HALF
IF (H0 .LT. HLB) H0 = HLB
IF (H0 .GT. HUB) H0 = HUB
90   H0 = SIGN(H0, TOUT - T0)
! Restore Y array from Y0, then exit. ----------------------------------
CALL DCOPY (N, Y0, 1, Y, 1)
NITER = ITER
IER = 0
RETURN
! Error return for TOUT - T0 too small. --------------------------------
100  IER = -1
RETURN
!----------------------- End of Subroutine DLHIN -----------------------
END

SUBROUTINE DSTOKA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVX, ACOR, WM, IWM, F, JAC, PSOL)
EXTERNAL F, JAC, PSOL
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVX, ACOR, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), SAVX(*), ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM
INTEGER LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER NEWT, NSFI, NSLJ, NJEV
INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION STIFR
DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,& 
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLS002/ STIFR, NEWT, NSFI, NSLJ, NJEV
COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN, JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
!-----------------------------------------------------------------------
! DSTOKA performs one step of the integration of an initial value
! problem for a system of Ordinary Differential Equations.
!
! This routine was derived from Subroutine DSTODPK in the DLSODPK
! package by the addition of automatic functional/Newton iteration
! switching and logic for re-use of Jacobian data.
!-----------------------------------------------------------------------
! Note: DSTOKA is independent of the value of the iteration method
! indicator MITER, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTOKA is done with the following variables:
!
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to F and JAC.
! Y      = an array of length .ge. N used as the Y argument in
!          all calls to F and JAC.
! YH     = an NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = a constant integer .ge. N, the first dimension of YH.
! YH1    = a one-dimensional array occupying the same space as YH.
! EWT    = an array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = an array of working storage, of length N.
!          Also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! SAVX   = an array of working storage, of length N.
! ACOR   = a work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = real and integer work arrays associated with matrix
!          operations in chord iteration (MITER .ne. 0).
! CCMAX  = maximum relative change in H*EL0 before DSETPK is called.
! H      = the step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = the minimum absolute value of the step size H to be used.
! HMXI   = inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = the independent variable. TN is updated on each step taken.
! JSTART = an integer used for input only, with the following
!          values and meanings:
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of H, MAXORD,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings:
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in DSETPK or DSOLPK.
!          A return with KFLAG = -1 or -2 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = the maximum order of integration method to be allowed.
! MAXCOR = the maximum number of corrector iterations allowed.
! MSBP   = maximum number of steps between DSETPK calls (MITER .gt. 0).
! MXNCF  = maximum number of convergence failures allowed.
! METH/MITER = the method flags.  See description in driver.
! N      = the number of first-order differential equations.
!-----------------------------------------------------------------------
INTEGER I, I1, IREDO, IRET, J, JB, JOK, M, NCF, NEWQ, NSLOW
DOUBLE PRECISION DCON, DDN, DEL, DELP, DRC, DSM, DUP, EXDN, EXSM, EXUP, DFNORM, R, RH
DOUBLE PRECISION RHDN, RHSM, RHUP, ROC, STIFF, TOLD, DVNORM
!
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set at 2
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
MEO = METH
NSLP = 0
NSLJ = 0
IPUP = 0
IRET = 3
NEWT = 0
STIFR = 0.0D0
GOTO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MEO) GOTO 110
CALL DCFODE (METH, ELCO, TESCO)
MEO = METH
IF (NQ .GT. MAXORD) GOTO 120
IALTH = L
IRET = 1
GOTO 150
110  IF (NQ .LE. MAXORD) GOTO 160
120  NQ = MAXORD
L = LMAX
DO 125 I = 1,L
125    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
EPCON = CONIT*TESCO(2,NQ)
DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
EXDN = 1.0D0/L
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
RH = MIN(RHDN,1.0D0)
IREDO = 3
IF (H .EQ. HOLD) GOTO 170
RH = MIN(RH,ABS(H/HOLD))
H = HOLD
GOTO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140  CALL DCFODE (METH, ELCO, TESCO)
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
EPCON = CONIT*TESCO(2,NQ)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! The flag IPUP is set according to whether matrix data is involved
! (NEWT .gt. 0 .and. JACFLG .ne. 0) or not, to trigger a call to DSETPK.
! IPUP is set to MITER when RC differs from 1 by more than CCMAX,
! and at least every MSBP steps, when JACFLG = 1.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
!-----------------------------------------------------------------------
200  IF (NEWT .EQ. 0 .OR. JACFLG .EQ. 0) THEN
DRC = 0.0D0
IPUP = 0
CRATE = 0.7D0
ELSE
DRC = ABS(RC - 1.0D0)
IF (DRC .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
ENDIF
TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS-norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
! Within the corrector loop, an estimated rate of convergence (ROC)
! and a stiffness ratio estimate (STIFF) are kept.  Corresponding
! global estimates are kept as CRATE and stifr.
!-----------------------------------------------------------------------
220  M = 0
MNEWT = 0
STIFF = 0.0D0
ROC = 0.05D0
NSLOW = 0
DO 230 I = 1,N
230    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
IF (NEWT .EQ. 0 .OR. IPUP .LE. 0) GOTO 250
!-----------------------------------------------------------------------
! If indicated, DSETPK is called to update any matrix data needed,
! before starting the corrector iteration.
! JOK is set to indicate if the matrix data need not be recomputed.
! IPUP is set to 0 as an indicator that the matrix data is up to date.
!-----------------------------------------------------------------------
JOK = 1
IF (NST .EQ. 0 .OR. NST .GT. NSLJ+50) JOK = -1
IF (ICF .EQ. 1 .AND. DRC .LT. 0.2D0) JOK = -1
IF (ICF .EQ. 2) JOK = -1
IF (JOK .EQ. -1) THEN
NSLJ = NST
NJEV = NJEV + 1
ENDIF
CALL DSETPK (NEQ, Y, YH1, EWT, ACOR, SAVF, JOK, WM, IWM, F, JAC)
IPUP = 0
RC = 1.0D0
DRC = 0.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .NE. 0) GOTO 430
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
270  IF (NEWT .NE. 0) GOTO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation, and STIFF is set to 1.0.
!-----------------------------------------------------------------------
DO 290 I = 1,N
SAVF(I) = H*SAVF(I) - YH(I,2)
290    Y(I) = SAVF(I) - ACOR(I)
DEL = DVNORM (N, Y, EWT)
DO 300 I = 1,N
Y(I) = YH(I,1) + EL(1)*SAVF(I)
300    ACOR(I) = SAVF(I)
STIFF = 1.0D0
GOTO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.  STIFF is set to the ratio of the norms
! of the residual and the correction vector.
!-----------------------------------------------------------------------
350  DO 360 I = 1,N
360    SAVX(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
DFNORM = DVNORM (N, SAVX, EWT)
CALL DSOLPK (NEQ, Y, SAVF, SAVX, EWT, WM, IWM, F, PSOL)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
DEL = DVNORM (N, SAVX, EWT)
IF (DEL .GT. 1.0D-8) STIFF = MAX(STIFF, DFNORM/DEL)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + SAVX(I)
380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is made for the iteration switch, and is also used
! in the convergence test.   If the iteration seems to be diverging or
! converging at a slow rate (.gt. 0.8 more than once), it is stopped.
!-----------------------------------------------------------------------
400  IF (M .NE. 0) THEN
ROC = MAX(0.05D0, DEL/DELP)
CRATE = MAX(0.2D0*CRATE,ROC)
ENDIF
DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/EPCON
IF (DCON .LE. 1.0D0) GOTO 450
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
IF (ROC .GT. 10.0D0) GOTO 410
IF (ROC .GT. 0.8D0) NSLOW = NSLOW + 1
IF (NSLOW .GE. 2) GOTO 410
MNEWT = M
DELP = DEL
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
GOTO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If functional iteration is being done (NEWT = 0) and MITER .gt. 0
! (and this is not the first step), then switch to Newton
! (NEWT = MITER), and retry the step.  (Setting STIFR = 1023 insures
! that a switch back will not occur for 10 step attempts.)
! If Newton iteration is being done, but using a preconditioner that
! is out of date (JACFLG .ne. 0 .and. JCUR = 0), then signal for a
! re-evalutation of the preconditioner, and retry the step.
! In all other cases, the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410  ICF = 1
IF (NEWT .EQ. 0) THEN
IF (NST .EQ. 0) GOTO 430
IF (MITER .EQ. 0) GOTO 430
NEWT = MITER
STIFR = 1023.0D0
IPUP = MITER
GOTO 220
ENDIF
IF (JCUR.EQ.1 .OR. JACFLG.EQ.0) GOTO 430
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
NCFN = NCFN + 1
RMAX = 2.0D0
TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 680
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 670
IF (NCF .EQ. MXNCF) GOTO 670
RH = 0.5D0
IPUP = MITER
IREDO = 1
GOTO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0 to signal that the
! preconditioner involved may need updating later.
! The stiffness ratio STIFR is updated using the latest STIFF value.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450  JCUR = 0
IF (NEWT .GT. 0) STIFR = 0.5D0*(STIFR + STIFF)
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! If Newton iteration is being done and STIFR is less than 1.5,
! then switch to functional iteration.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
IF (NEWT .EQ. 0) NSFI = NSFI + 1
IF (NEWT .GT. 0 .AND. STIFR .LT. 1.5D0) NEWT = 0
HU = H
NQU = NQ
DO 470 J = 1,L
DO 470 I = 1,N
470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -3) GOTO 640
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! in the case of failure, RHUP = 0.0 to avoid an order increase.
! the largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 560
DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GOTO 610
IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640  IF (KFLAG .EQ. -10) GOTO 660
RH = 0.1D0
RH = MAX(HMIN/ABS(H),RH)
H = H*RH
DO 645 I = 1,N
645    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
DO 650 I = 1,N
650    YH(I,2) = H*SAVF(I)
IPUP = MITER
IALTH = 5
IF (NQ .EQ. 1) GOTO 200
NQ = 1
L = 2
IRET = 3
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -3
GOTO 720
690  RMAX = 10.0D0
700  R = 1.0D0/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- End of Subroutine DSTOKA ----------------------
END

SUBROUTINE DSETPK (NEQ, Y, YSV, EWT, FTEM, SAVF, JOK, WM, IWM, F, JAC)
EXTERNAL F, JAC
INTEGER NEQ, JOK, IWM
DOUBLE PRECISION Y, YSV, EWT, FTEM, SAVF, WM
DIMENSION NEQ(*), Y(*), YSV(*), EWT(*), FTEM(*), SAVF(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN, JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT, NNI, NLI, NPS, NCFN, NCFL
!-----------------------------------------------------------------------
! DSETPK is called by DSTOKA to interface with the user-supplied
! routine JAC, to compute and process relevant parts of
! the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
! as need for preconditioning matrix operations later.
!
! In addition to variables described previously, communication
! with DSETPK uses the following:
! Y     = array containing predicted values on entry.
! YSV   = array containing predicted y, to be saved (YH1 in DSTOKA).
! FTEM  = work array of length N (ACOR in DSTOKA).
! SAVF  = array containing f evaluated at predicted y.
! JOK   = input flag showing whether it was judged that Jacobian matrix
!         data need not be recomputed (JOK = 1) or needs to be
!         (JOK = -1).
! WM    = real work space for matrices.
!         Space for preconditioning data starts at WM(LOCWP).
! IWM   = integer work space.
!         Space for preconditioning data starts at IWM(LOCIWP).
! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!         JAC returned an error flag.
! JCUR  = output flag to indicate whether the matrix data involved
!         is now current (JCUR = 1) or not (JCUR = 0).
! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
!-----------------------------------------------------------------------
INTEGER IER
DOUBLE PRECISION HL0
!
IERPJ = 0
JCUR = 0
IF (JOK .EQ. -1) JCUR = 1
HL0 = EL0*H
CALL JAC (F, NEQ, TN, Y, YSV, EWT, SAVF, FTEM, HL0, JOK, WM(LOCWP), IWM(LOCIWP), IER)
NJE = NJE + 1
IF (IER .EQ. 0) RETURN
IERPJ = 1
RETURN
!----------------------- End of Subroutine DSETPK ----------------------
END

SUBROUTINE DSRCKR (RSAV, ISAV, JOB)
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of
! the Common blocks DLS001, DLS002, DLSR01, DLPK01, which
! are used internally by the DLSODKR solver.
!
! RSAV = real array of length 228 or more.
! ISAV = integer array of length 63 or more.
! JOB  = flag indicating to save or restore the Common blocks:
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
INTEGER ISAV, JOB
INTEGER ILS, ILS2, ILSR, ILSP
INTEGER I, IOFF, LENILP, LENRLP, LENILS, LENRLS, LENILR, LENRLR
DOUBLE PRECISION RSAV,   RLS, RLS2, RLSR, RLSP
DIMENSION RSAV(*), ISAV(*)
SAVE LENRLS, LENILS, LENRLP, LENILP, LENRLR, LENILR
COMMON /DLS001/ RLS(218), ILS(37)
COMMON /DLS002/ RLS2, ILS2(4)
COMMON /DLSR01/ RLSR(5), ILSR(9)
COMMON /DLPK01/ RLSP(4), ILSP(13)
DATA LENRLS/218/, LENILS/37/, LENRLP/4/, LENILP/13/
DATA LENRLR/5/, LENILR/9/
!
IF (JOB .EQ. 2) GOTO 100
CALL DCOPY (LENRLS, RLS, 1, RSAV, 1)
RSAV(LENRLS+1) = RLS2
CALL DCOPY (LENRLR, RLSR, 1, RSAV(LENRLS+2), 1)
CALL DCOPY (LENRLP, RLSP, 1, RSAV(LENRLS+LENRLR+2), 1)
DO 20 I = 1,LENILS
20     ISAV(I) = ILS(I)
ISAV(LENILS+1) = ILS2(1)
ISAV(LENILS+2) = ILS2(2)
ISAV(LENILS+3) = ILS2(3)
ISAV(LENILS+4) = ILS2(4)
IOFF = LENILS + 2
DO 30 I = 1,LENILR
30     ISAV(IOFF+I) = ILSR(I)
IOFF = IOFF + LENILR
DO 40 I = 1,LENILP
40     ISAV(IOFF+I) = ILSP(I)
RETURN
!
100  CONTINUE
CALL DCOPY (LENRLS, RSAV, 1, RLS, 1)
RLS2 = RSAV(LENRLS+1)
CALL DCOPY (LENRLR, RSAV(LENRLS+2), 1, RLSR, 1)
CALL DCOPY (LENRLP, RSAV(LENRLS+LENRLR+2), 1, RLSP, 1)
DO 120 I = 1,LENILS
120    ILS(I) = ISAV(I)
ILS2(1) = ISAV(LENILS+1)
ILS2(2) = ISAV(LENILS+2)
ILS2(3) = ISAV(LENILS+3)
ILS2(4) = ISAV(LENILS+4)
IOFF = LENILS + 2
DO 130 I = 1,LENILR
130    ILSR(I) = ISAV(IOFF+I)
IOFF = IOFF + LENILR
DO 140 I = 1,LENILP
140    ILSP(I) = ISAV(IOFF+I)
RETURN
!----------------------- End of Subroutine DSRCKR ----------------------
END

SUBROUTINE DAINVG (RES, ADDA, NEQ, T, Y, YDOT, MITER, ML, MU, PW, IPVT, IER )
EXTERNAL RES, ADDA
INTEGER NEQ, MITER, ML, MU, IPVT, IER
INTEGER I, LENPW, MLP1, NROWPW
DOUBLE PRECISION T, Y, YDOT, PW
DIMENSION Y(*), YDOT(*), PW(*), IPVT(*)
!-----------------------------------------------------------------------
! This subroutine computes the initial value
! of the vector YDOT satisfying
!     A * YDOT = g(t,y)
! when A is nonsingular.  It is called by DLSODI for
! initialization only, when ISTATE = 0 .
! DAINVG returns an error flag IER:
!   IER  =  0  means DAINVG was successful.
!   IER .ge. 2 means RES returned an error flag IRES = IER.
!   IER .lt. 0 means the a-matrix was found to be singular.
!-----------------------------------------------------------------------
!
IF (MITER .GE. 4)  GOTO 100
!
! Full matrix case -----------------------------------------------------
!
LENPW = NEQ*NEQ
DO 10  I = 1, LENPW
10    PW(I) = 0.0D0
!
IER = 1
CALL RES ( NEQ, T, Y, PW, YDOT, IER )
IF (IER .GT. 1) RETURN
!
CALL ADDA ( NEQ, T, Y, 0, 0, PW, NEQ )
CALL DGEFA ( PW, NEQ, NEQ, IPVT, IER )
IF (IER .EQ. 0) GOTO 20
IER = -IER
RETURN
20 CALL DGESL ( PW, NEQ, NEQ, IPVT, YDOT, 0 )
RETURN
!
! Band matrix case -----------------------------------------------------
!
100 CONTINUE
NROWPW = 2*ML + MU + 1
LENPW = NEQ * NROWPW
DO 110  I = 1, LENPW
110    PW(I) = 0.0D0
!
IER = 1
CALL RES ( NEQ, T, Y, PW, YDOT, IER )
IF (IER .GT. 1) RETURN
!
MLP1 = ML + 1
CALL ADDA ( NEQ, T, Y, ML, MU, PW(MLP1), NROWPW )
CALL DGBFA ( PW, NROWPW, NEQ, ML, MU, IPVT, IER )
IF (IER .EQ. 0) GOTO 120
IER = -IER
RETURN
120 CALL DGBSL ( PW, NROWPW, NEQ, ML, MU, IPVT, YDOT, 0 )
RETURN
!----------------------- End of Subroutine DAINVG ----------------------
END

SUBROUTINE DSTODI (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVR, ACOR, WM, IWM, RES, ADDA, JAC, PJAC, SLVS )
EXTERNAL RES, ADDA, JAC, PJAC, SLVS
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVR, ACOR, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), SAVR(*), ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM
INTEGER LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, I1, IREDO, IRES, IRET, J, JB, KGO, M, NCF, NEWQ
DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, ELJH, EL1H, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
!-----------------------------------------------------------------------
! DSTODI performs one step of the integration of an initial value
! problem for a system of Ordinary Differential Equations.
! Note: DSTODI is independent of the value of the iteration method
! indicator MITER, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTODI is done with the following variables:
!
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to RES, ADDA,
!          and JAC.
! Y      = an array of length .ge. N used as the Y argument in
!          all calls to RES, JAC, and ADDA.
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls tO RES, G, ADDA,
!          and JAC.
! YH     = an NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = a constant integer .ge. N, the first dimension of YH.
! YH1    = a one-dimensional array occupying the same space as YH.
! EWT    = an array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = an array of working storage, of length N. also used for
!          input of YH(*,MAXORD+2) when JSTART = -1 and MAXORD is less
!          than the current order NQ.
!          Same as YDOTI in the driver.
! SAVR   = an array of working storage, of length N.
! ACOR   = a work array of length N used for the accumulated
!          corrections. On a succesful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = real and integer work arrays associated with matrix
!          operations in chord iteration.
! PJAC   = name of routine to evaluate and preprocess Jacobian matrix.
! SLVS   = name of routine to solve linear system in chord iteration.
! CCMAX  = maximum relative change in H*EL0 before PJAC is called.
! H      = the step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = the minimum absolute value of the step size H to be used.
! HMXI   = inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = the independent variable. TN is updated on each step taken.
! JSTART = an integer used for input only, with the following
!          values and meanings:
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of H, MAXORD,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings:
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  RES ordered immediate return.
!              -4  error condition from RES could not be avoided.
!              -5  fatal error in PJAC or SLVS.
!          A return with KFLAG = -1, -2, or -4 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = the maximum order of integration method to be allowed.
! MAXCOR = the maximum number of corrector iterations allowed.
! MSBP   = maximum number of steps between PJAC calls.
! MXNCF  = maximum number of convergence failures allowed.
! METH/MITER = the method flags.  See description in driver.
! N      = the number of first-order differential equations.
!-----------------------------------------------------------------------
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set at 2
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
MEO = METH
NSLP = 0
IPUP = MITER
IRET = 3
GOTO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MEO) GOTO 110
CALL DCFODE (METH, ELCO, TESCO)
MEO = METH
IF (NQ .GT. MAXORD) GOTO 120
IALTH = L
IRET = 1
GOTO 150
110  IF (NQ .LE. MAXORD) GOTO 160
120  NQ = MAXORD
L = LMAX
DO 125 I = 1,L
125    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
EXDN = 1.0D0/L
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
RH = MIN(RHDN,1.0D0)
IREDO = 3
IF (H .EQ. HOLD) GOTO 170
RH = MIN(RH,ABS(H/HOLD))
H = HOLD
GOTO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140  CALL DCFODE (METH, ELCO, TESCO)
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS-norm of each correction, weighted by H and the
! error weight vector EWT.  The sum of the corrections is accumulated
! in ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220  M = 0
DO 230 I = 1,N
SAVF(I) = YH(I,2) / H
230    Y(I) = YH(I,1)
IF (IPUP .LE. 0) GOTO 240
!-----------------------------------------------------------------------
! If indicated, the matrix P = A - H*EL(1)*dr/dy is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVR, SAVF, WM, IWM, RES, JAC, ADDA )
IPUP = 0
RC = 1.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .EQ. 0) GOTO 250
IF (IERPJ .LT. 0) GOTO 435
IRES = IERPJ
GOTO (430, 435, 430), IRES
! Get residual at predicted values, if not already done in PJAC. -------
240  IRES = 1
CALL RES ( NEQ, TN, Y, SAVF, SAVR, IRES )
NFE = NFE + 1
KGO = ABS(IRES)
GOTO ( 250, 435, 430 ) , KGO
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
!-----------------------------------------------------------------------
! Solve the linear system with the current residual as
! right-hand side and P as coefficient matrix.
!-----------------------------------------------------------------------
270  CONTINUE
CALL SLVS (WM, IWM, SAVR, SAVF)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
EL1H = EL(1) * H
DEL = DVNORM (N, SAVR, EWT) * ABS(H)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + SAVR(I)
SAVF(I) = ACOR(I) + YH(I,2)/H
380    Y(I) = YH(I,1) + EL1H*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
IF (DCON .LE. 1.0D0) GOTO 460
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
DELP = DEL
IRES = 1
CALL RES ( NEQ, TN, Y, SAVF, SAVR, IRES )
NFE = NFE + 1
KGO = ABS(IRES)
GOTO ( 270, 435, 410 ) , KGO
!-----------------------------------------------------------------------
! The correctors failed to converge, or RES has returned abnormally.
! on a convergence failure, if the Jacobian is out of date, PJAC is
! called for the next try.  Otherwise the YH array is retracted to its
! values before prediction, and H is reduced, if possible.
! take an error exit if IRES = 2, or H cannot be reduced, or MXNCF
! failures have occurred, or a fatal error occurred in PJAC or SLVS.
!-----------------------------------------------------------------------
410  ICF = 1
IF (JCUR .EQ. 1) GOTO 430
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
RMAX = 2.0D0
435  TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IRES .EQ. 2) GOTO 680
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 685
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 450
IF (NCF .EQ. MXNCF) GOTO 450
RH = 0.25D0
IPUP = MITER
IREDO = 1
GOTO 170
450  IF (IRES .EQ. 3) GOTO 680
GOTO 670
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
460  JCUR = 0
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = ABS(H) * DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
HU = H
NQU = NQ
DO 470 J = 1,L
ELJH = EL(J)*H
DO 470 I = 1,N
470      YH(I,J) = YH(I,J) + ELJH*ACOR(I)
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.1 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!DIR$ IVDEP
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -7) GOTO 660
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = ABS(H) * DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 560
DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = H*EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GOTO 610
IF (KFLAG .LE. -2) RH = MIN(RH,0.1D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -1 - IRES
GOTO 720
685  KFLAG = -5
GOTO 720
690  RMAX = 10.0D0
700  R = H/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- End of Subroutine DSTODI ----------------------
END

SUBROUTINE DPREPJI (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WM, IWM, RES, JAC, ADDA)
EXTERNAL RES, JAC, ADDA
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*), S(*), SAVR(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, I1, I2, IER, II, IRES, J, J1, JJ, LENP, MBA, MBAND, MEB1, MEBAND, ML, ML3, MU
DOUBLE PRECISION CON, FAC, HL0, R, SRUR, YI, YJ, YJJ
!-----------------------------------------------------------------------
! DPREPJI is called by DSTODI to compute and process the matrix
! P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
! where r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
! routine JAC if MITER = 1 or 4, or by finite differencing if MITER =
! 2 or 5.  J is stored in WM, rescaled, and ADDA is called to generate
! P. P is then subjected to LU decomposition in preparation
! for later solution of linear systems with P as coefficient
! matrix.  This is done by DGEFA if MITER = 1 or 2, and by
! DGBFA if MITER = 4 or 5.
!
! In addition to variables described previously, communication
! with DPREPJI uses the following:
! Y     = array containing predicted values on entry.
! RTEM  = work array of length N (ACOR in DSTODI).
! SAVR  = array used for output only.  On output it contains the
!         residual evaluated at current values of t and y.
! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
! WM    = real work space for matrices.  On output it contains the
!         LU decomposition of P.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data:
!         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
! IWM   = integer work space containing pivot information, starting at
!         IWM(21).  IWM also contains the band parameters
!         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! EL0   = el(1) (input).
! IERPJ = output error flag.
!         = 0 if no trouble occurred,
!         = 1 if the P matrix was found to be singular,
!         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses the Common variables EL0, H, TN, UROUND,
! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
NJE = NJE + 1
HL0 = H*EL0
IERPJ = 0
JCUR = 1
GOTO (100, 200, 300, 400, 500), MITER
! If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
100  IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
LENP = N*N
DO 110 I = 1,LENP
110    WM(I+2) = 0.0D0
CALL JAC ( NEQ, TN, Y, S, 0, 0, WM(3), N )
CON = -HL0
DO 120 I = 1,LENP
120    WM(I+2) = WM(I+2)*CON
GOTO 240
! If MITER = 2, make N + 1 calls to RES to approximate J. --------------
200  CONTINUE
IRES = -1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
SRUR = WM(1)
J1 = 2
DO 230 J = 1,N
YJ = Y(J)
R = MAX(SRUR*ABS(YJ),0.01D0/EWT(J))
Y(J) = Y(J) + R
FAC = -HL0/R
CALL RES ( NEQ, TN, Y, S, RTEM, IRES )
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
DO 220 I = 1,N
220      WM(I+J1) = (RTEM(I) - SAVR(I))*FAC
Y(J) = YJ
J1 = J1 + N
230    CONTINUE
IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
! Add matrix A. --------------------------------------------------------
240  CONTINUE
CALL ADDA(NEQ, TN, Y, 0, 0, WM(3), N)
! Do LU decomposition on P. --------------------------------------------
CALL DGEFA (WM(3), N, N, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! Dummy section for MITER = 3
300  RETURN
! If MITER = 4, call RES, then JAC, and multiply by scalar. ------------
400  IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
ML = IWM(1)
MU = IWM(2)
ML3 = ML + 3
MBAND = ML + MU + 1
MEBAND = MBAND + ML
LENP = MEBAND*N
DO 410 I = 1,LENP
410    WM(I+2) = 0.0D0
CALL JAC ( NEQ, TN, Y, S, ML, MU, WM(ML3), MEBAND)
CON = -HL0
DO 420 I = 1,LENP
420    WM(I+2) = WM(I+2)*CON
GOTO 570
! If MITER = 5, make ML + MU + 2 calls to RES to approximate J. --------
500  CONTINUE
IRES = -1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
ML = IWM(1)
MU = IWM(2)
ML3 = ML + 3
MBAND = ML + MU + 1
MBA = MIN(MBAND,N)
MEBAND = MBAND + ML
MEB1 = MEBAND - 1
SRUR = WM(1)
DO 560 J = 1,MBA
DO 530 I = J,N,MBAND
YI = Y(I)
R = MAX(SRUR*ABS(YI),0.01D0/EWT(I))
530      Y(I) = Y(I) + R
CALL RES ( NEQ, TN, Y, S, RTEM, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
DO 550 JJ = J,N,MBAND
Y(JJ) = YH(JJ,1)
YJJ = Y(JJ)
R = MAX(SRUR*ABS(YJJ),0.01D0/EWT(JJ))
FAC = -HL0/R
I1 = MAX(JJ-MU,1)
I2 = MIN(JJ+ML,N)
II = JJ*MEB1 - ML + 2
DO 540 I = I1,I2
540        WM(II+I) = (RTEM(I) - SAVR(I))*FAC
550      CONTINUE
560    CONTINUE
IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
! Add matrix A. --------------------------------------------------------
570 CONTINUE
CALL ADDA(NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
! Do LU decomposition of P. --------------------------------------------
CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
600  IERPJ = IRES
RETURN
!----------------------- End of Subroutine DPREPJI ---------------------
END

SUBROUTINE DAIGBT (RES, ADDA, NEQ, T, Y, YDOT, MB, NB, PW, IPVT, IER )
EXTERNAL RES, ADDA
INTEGER NEQ, MB, NB, IPVT, IER
INTEGER I, LENPW, LBLOX, LPB, LPC
DOUBLE PRECISION T, Y, YDOT, PW
DIMENSION Y(*), YDOT(*), PW(*), IPVT(*), NEQ(*)
!-----------------------------------------------------------------------
! This subroutine computes the initial value
! of the vector YDOT satisfying
!     A * YDOT = g(t,y)
! when A is nonsingular.  It is called by DLSOIBT for
! initialization only, when ISTATE = 0 .
! DAIGBT returns an error flag IER:
!   IER  =  0  means DAIGBT was successful.
!   IER .ge. 2 means RES returned an error flag IRES = IER.
!   IER .lt. 0 means the A matrix was found to have a singular
!              diagonal block (hence YDOT could not be solved for).
!-----------------------------------------------------------------------
LBLOX = MB*MB*NB
LPB = 1 + LBLOX
LPC = LPB + LBLOX
LENPW = 3*LBLOX
DO 10 I = 1,LENPW
10     PW(I) = 0.0D0
IER = 1
CALL RES (NEQ, T, Y, PW, YDOT, IER)
IF (IER .GT. 1) RETURN
CALL ADDA (NEQ, T, Y, MB, NB, PW(1), PW(LPB), PW(LPC) )
CALL DDECBT (MB, NB, PW, PW(LPB), PW(LPC), IPVT, IER)
IF (IER .EQ. 0) GOTO 20
IER = -IER
RETURN
20   CALL DSOLBT (MB, NB, PW, PW(LPB), PW(LPC), YDOT, IPVT)
RETURN
!----------------------- End of Subroutine DAIGBT ----------------------
END

SUBROUTINE DPJIBT (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WM, IWM, RES, JAC, ADDA)
EXTERNAL RES, JAC, ADDA
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WM
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*), S(*), SAVR(*), WM(*), IWM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, IER, IIA, IIB, IIC, IPA, IPB, IPC, IRES, J, J1, J2, K, K1, LENP, LBLOX, LPB, LPC, MB, MBSQ, MWID, NB
DOUBLE PRECISION CON, FAC, HL0, R, SRUR
!-----------------------------------------------------------------------
! DPJIBT is called by DSTODI to compute and process the matrix
! P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
! and r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
! routine JAC if MITER = 1, or by finite differencing if MITER = 2.
! J is stored in WM, rescaled, and ADDA is called to generate P.
! P is then subjected to LU decomposition by DDECBT in preparation
! for later solution of linear systems with P as coefficient matrix.
!
! In addition to variables described previously, communication
! with DPJIBT uses the following:
! Y     = array containing predicted values on entry.
! RTEM  = work array of length N (ACOR in DSTODI).
! SAVR  = array used for output only.  On output it contains the
!         residual evaluated at current values of t and y.
! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
! WM    = real work space for matrices.  On output it contains the
!         LU decomposition of P.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data:
!         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
! IWM   = integer work space containing pivot information, starting at
!         IWM(21).  IWM also contains block structure parameters
!         MB = IWM(1) and NB = IWM(2).
! EL0   = EL(1) (input).
! IERPJ = output error flag.
!         = 0 if no trouble occurred,
!         = 1 if the P matrix was found to be unfactorable,
!         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses the Common variables EL0, H, TN, UROUND,
! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
NJE = NJE + 1
HL0 = H*EL0
IERPJ = 0
JCUR = 1
MB = IWM(1)
NB = IWM(2)
MBSQ = MB*MB
LBLOX = MBSQ*NB
LPB = 3 + LBLOX
LPC = LPB + LBLOX
LENP = 3*LBLOX
GOTO (100, 200), MITER
! If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
100  IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
DO 110 I = 1,LENP
110    WM(I+2) = 0.0D0
CALL JAC (NEQ, TN, Y, S, MB, NB, WM(3), WM(LPB), WM(LPC))
CON = -HL0
DO 120 I = 1,LENP
120    WM(I+2) = WM(I+2)*CON
GOTO 260
!
! If MITER = 2, make 3*MB + 1 calls to RES to approximate J. -----------
200  CONTINUE
IRES = -1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
MWID = 3*MB
SRUR = WM(1)
DO 205 I = 1,LENP
205    WM(2+I) = 0.0D0
DO 250 K = 1,3
DO 240 J = 1,MB
!         Increment Y(I) for group of column indices, and call RES. ----
J1 = J+(K-1)*MB
DO 210 I = J1,N,MWID
R = MAX(SRUR*ABS(Y(I)),0.01D0/EWT(I))
Y(I) = Y(I) + R
210      CONTINUE
CALL RES (NEQ, TN, Y, S, RTEM, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
DO 215 I = 1,N
215        RTEM(I) = RTEM(I) - SAVR(I)
K1 = K
DO 230 I = J1,N,MWID
!           Get Jacobian elements in column I (block-column K1). -------
Y(I) = YH(I,1)
R = MAX(SRUR*ABS(Y(I)),0.01D0/EWT(I))
FAC = -HL0/R
!           Compute and load elements PA(*,J,K1). ----------------------
IIA = I - J
IPA = 2 + (J-1)*MB + (K1-1)*MBSQ
DO 221 J2 = 1,MB
221          WM(IPA+J2) = RTEM(IIA+J2)*FAC
IF (K1 .LE. 1) GOTO 223
!           Compute and load elements PB(*,J,K1-1). --------------------
IIB = IIA - MB
IPB = IPA + LBLOX - MBSQ
DO 222 J2 = 1,MB
222          WM(IPB+J2) = RTEM(IIB+J2)*FAC
223        CONTINUE
IF (K1 .GE. NB) GOTO 225
!           Compute and load elements PC(*,J,K1+1). --------------------
IIC = IIA + MB
IPC = IPA + 2*LBLOX + MBSQ
DO 224 J2 = 1,MB
224          WM(IPC+J2) = RTEM(IIC+J2)*FAC
225        CONTINUE
IF (K1 .NE. 3) GOTO 227
!           Compute and load elements PC(*,J,1). -----------------------
IPC = IPA - 2*MBSQ + 2*LBLOX
DO 226 J2 = 1,MB
226          WM(IPC+J2) = RTEM(J2)*FAC
227        CONTINUE
IF (K1 .NE. NB-2) GOTO 229
!           Compute and load elements PB(*,J,NB). ----------------------
IIB = N - MB
IPB = IPA + 2*MBSQ + LBLOX
DO 228 J2 = 1,MB
228          WM(IPB+J2) = RTEM(IIB+J2)*FAC
229      K1 = K1 + 3
230      CONTINUE
240    CONTINUE
250  CONTINUE
! RES call for first corrector iteration. ------------------------------
IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
! Add matrix A. --------------------------------------------------------
260  CONTINUE
CALL ADDA (NEQ, TN, Y, MB, NB, WM(3), WM(LPB), WM(LPC))
! Do LU decomposition on P. --------------------------------------------
CALL DDECBT (MB, NB, WM(3), WM(LPB), WM(LPC), IWM(21), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
600  IERPJ = IRES
RETURN
!----------------------- End of Subroutine DPJIBT ----------------------
END

SUBROUTINE DSLSBT (WM, IWM, X, TEM)
INTEGER IWM
INTEGER LBLOX, LPB, LPC, MB, NB
DOUBLE PRECISION WM, X, TEM
DIMENSION WM(*), IWM(*), X(*), TEM(*)
!-----------------------------------------------------------------------
! This routine acts as an interface between the core integrator
! routine and the DSOLBT routine for the solution of the linear system
! arising from chord iteration.
! Communication with DSLSBT uses the following variables:
! WM    = real work space containing the LU decomposition,
!         starting at WM(3).
! IWM   = integer work space containing pivot information, starting at
!         IWM(21).  IWM also contains block structure parameters
!         MB = IWM(1) and NB = IWM(2).
! X     = the right-hand side vector on input, and the solution vector
!         on output, of length N.
! TEM   = vector of work space of length N, not used in this version.
!-----------------------------------------------------------------------
MB = IWM(1)
NB = IWM(2)
LBLOX = MB*MB*NB
LPB = 3 + LBLOX
LPC = LPB + LBLOX
CALL DSOLBT (MB, NB, WM(3), WM(LPB), WM(LPC), X, IWM(21))
RETURN
!----------------------- End of Subroutine DSLSBT ----------------------
END

SUBROUTINE DDECBT (M, N, A, B, C, IP, IER)
INTEGER M, N, IP(M,N), IER
DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N)
!-----------------------------------------------------------------------
! Block-tridiagonal matrix decomposition routine.
! Written by A. C. Hindmarsh.
! Latest revision:  November 10, 1983 (ACH)
! Reference:  UCID-30150
!             Solution of Block-Tridiagonal Systems of Linear
!             Algebraic Equations
!             A.C. Hindmarsh
!             February 1977
! The input matrix contains three blocks of elements in each block-row,
! including blocks in the (1,3) and (N,N-2) block positions.
! DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
! for solution of blocks.  Partial pivoting is done within
! block-rows only.
!
! Note: this version uses LINPACK routines DGEFA/DGESL instead of
! of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
! for dot product calculations.
!
! Input:
!     M = order of each block.
!     N = number of blocks in each direction of the matrix.
!         N must be 4 or more.  The complete matrix has order M*N.
!     A = M by M by N array containing diagonal blocks.
!         A(i,j,k) contains the (i,j) element of the k-th block.
!     B = M by M by N array containing the super-diagonal blocks
!         (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
!         block position (in B(*,*,N)).
!     C = M by M by N array containing the subdiagonal blocks
!         (in C(*,*,k) for k = 2,3,...,N) and the block in the
!         (1,3) block position (in C(*,*,1)).
!    IP = integer array of length M*N for working storage.
! Output:
! A,B,C = M by M by N arrays containing the block-LU decomposition
!         of the input matrix.
!    IP = M by N array of pivot information.  IP(*,k) contains
!         information for the k-th digonal block.
!   IER = 0  if no trouble occurred, or
!       = -1 if the input value of M or N was illegal, or
!       = k  if a singular matrix was found in the k-th diagonal block.
! Use DSOLBT to solve the associated linear system.
!
! External routines required: DGEFA and DGESL (from LINPACK) and
! DDOT (from the BLAS, or Basic Linear Algebra package).
!-----------------------------------------------------------------------
INTEGER NM1, NM2, KM1, I, J, K
DOUBLE PRECISION DP, DDOT
IF (M .LT. 1 .OR. N .LT. 4) GOTO 210
NM1 = N - 1
NM2 = N - 2
! Process the first block-row. -----------------------------------------
CALL DGEFA (A, M, M, IP, IER)
K = 1
IF (IER .NE. 0) GOTO 200
DO 10 J = 1,M
CALL DGESL (A, M, M, IP, B(1,J,1), 0)
CALL DGESL (A, M, M, IP, C(1,J,1), 0)
10     CONTINUE
! Adjust B(*,*,2). -----------------------------------------------------
DO 40 J = 1,M
DO 30 I = 1,M
DP = DDOT (M, C(I,1,2), M, C(1,J,1), 1)
B(I,J,2) = B(I,J,2) - DP
30       CONTINUE
40     CONTINUE
! Main loop.  Process block-rows 2 to N-1. -----------------------------
DO 100 K = 2,NM1
KM1 = K - 1
DO 70 J = 1,M
DO 60 I = 1,M
DP = DDOT (M, C(I,1,K), M, B(1,J,KM1), 1)
A(I,J,K) = A(I,J,K) - DP
60         CONTINUE
70       CONTINUE
CALL DGEFA (A(1,1,K), M, M, IP(1,K), IER)
IF (IER .NE. 0) GOTO 200
DO 80 J = 1,M
80       CALL DGESL (A(1,1,K), M, M, IP(1,K), B(1,J,K), 0)
100    CONTINUE
! Process last block-row and return. -----------------------------------
DO 130 J = 1,M
DO 120 I = 1,M
DP = DDOT (M, B(I,1,N), M, B(1,J,NM2), 1)
C(I,J,N) = C(I,J,N) - DP
120      CONTINUE
130    CONTINUE
DO 160 J = 1,M
DO 150 I = 1,M
DP = DDOT (M, C(I,1,N), M, B(1,J,NM1), 1)
A(I,J,N) = A(I,J,N) - DP
150      CONTINUE
160    CONTINUE
CALL DGEFA (A(1,1,N), M, M, IP(1,N), IER)
K = N
IF (IER .NE. 0) GOTO 200
RETURN
! Error returns. -------------------------------------------------------
200  IER = K
RETURN
210  IER = -1
RETURN
!----------------------- End of Subroutine DDECBT ----------------------
END

SUBROUTINE DSOLBT (M, N, A, B, C, Y, IP)
INTEGER M, N, IP(M,N)
DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N), Y(M,N)
!-----------------------------------------------------------------------
! Solution of block-tridiagonal linear system.
! Coefficient matrix must have been previously processed by DDECBT.
! M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
! Written by A. C. Hindmarsh.
! Input:
!     M = order of each block.
!     N = number of blocks in each direction of matrix.
! A,B,C = M by M by N arrays containing block LU decomposition
!         of coefficient matrix from DDECBT.
!    IP = M by N integer array of pivot information from DDECBT.
!     Y = array of length M*N containg the right-hand side vector
!         (treated as an M by N array here).
! Output:
!     Y = solution vector, of length M*N.
!
! External routines required: DGESL (LINPACK) and DDOT (BLAS).
!-----------------------------------------------------------------------
!
INTEGER NM1, NM2, I, K, KB, KM1, KP1
DOUBLE PRECISION DP, DDOT
NM1 = N - 1
NM2 = N - 2
! Forward solution sweep. ----------------------------------------------
CALL DGESL (A, M, M, IP, Y, 0)
DO 30 K = 2,NM1
KM1 = K - 1
DO 20 I = 1,M
DP = DDOT (M, C(I,1,K), M, Y(1,KM1), 1)
Y(I,K) = Y(I,K) - DP
20       CONTINUE
CALL DGESL (A(1,1,K), M, M, IP(1,K), Y(1,K), 0)
30     CONTINUE
DO 50 I = 1,M
DP = DDOT (M, C(I,1,N), M, Y(1,NM1), 1) + DDOT (M, B(I,1,N), M, Y(1,NM2), 1)
Y(I,N) = Y(I,N) - DP
50     CONTINUE
CALL DGESL (A(1,1,N), M, M, IP(1,N), Y(1,N), 0)
! Backward solution sweep. ---------------------------------------------
DO 80 KB = 1,NM1
K = N - KB
KP1 = K + 1
DO 70 I = 1,M
DP = DDOT (M, B(I,1,K), M, Y(1,KP1), 1)
Y(I,K) = Y(I,K) - DP
70       CONTINUE
80     CONTINUE
DO 100 I = 1,M
DP = DDOT (M, C(I,1,1), M, Y(1,3), 1)
Y(I,1) = Y(I,1) - DP
100    CONTINUE
RETURN
!----------------------- End of Subroutine DSOLBT ----------------------
END

SUBROUTINE DIPREPI (NEQ, Y, S, RWORK, IA, JA, IC, JC, IPFLAG, RES, JAC, ADDA)
EXTERNAL RES, JAC, ADDA
INTEGER NEQ, IA, JA, IC, JC, IPFLAG
DOUBLE PRECISION Y, S, RWORK
DIMENSION NEQ(*), Y(*), S(*), RWORK(*), IA(*), JA(*), IC(*), JC(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, &
MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMAX, LEWTN, LYHD, LYHN
!-----------------------------------------------------------------------
! This routine serves as an interface between the driver and
! Subroutine DPREPI.  Tasks performed here are:
!  * call DPREPI,
!  * reset the required WM segment length LENWK,
!  * move YH back to its final location (following WM in RWORK),
!  * reset pointers for YH, SAVR, EWT, and ACOR, and
!  * move EWT to its new position if ISTATE = 0 or 1.
! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
! no trouble, and IPFLAG is the value of the DPREPI error flag IPPER
! if there was trouble in Subroutine DPREPI.
!-----------------------------------------------------------------------
IPFLAG = 0
! Call DPREPI to do matrix preprocessing operations. -------------------
CALL DPREPI (NEQ, Y, S, RWORK(LYH), RWORK(LSAVF), RWORK(LEWT), RWORK(LACOR), IA, JA, IC, JC, RWORK(LWM), RWORK(LWM), IPFLAG,&
RES, JAC, ADDA)
LENWK = MAX(LREQ,LWMIN)
IF (IPFLAG .LT. 0) RETURN
! If DPREPI was successful, move YH to end of required space for WM. ---
LYHN = LWM + LENWK
IF (LYHN .GT. LYH) RETURN
LYHD = LYH - LYHN
IF (LYHD .EQ. 0) GOTO 20
IMAX = LYHN - 1 + LENYHM
DO 10 I=LYHN,IMAX
10     RWORK(I) = RWORK(I+LYHD)
LYH = LYHN
! Reset pointers for SAVR, EWT, and ACOR. ------------------------------
20   LSAVF = LYH + LENYH
LEWTN = LSAVF + N
LACOR = LEWTN + N
IF (ISTATC .EQ. 3) GOTO 40
! If ISTATE = 1, move EWT (left) to its new position. ------------------
IF (LEWTN .GT. LEWT) RETURN
DO 30 I=1,N
30     RWORK(I+LEWTN-1) = RWORK(I+LEWT-1)
40   LEWT = LEWTN
RETURN
!----------------------- End of Subroutine DIPREPI ---------------------
END

SUBROUTINE DPREPI (NEQ, Y, S, YH, SAVR, EWT, RTEM, IA, JA, IC, JC, WK, IWK, IPPER, RES, JAC, ADDA)
EXTERNAL RES, JAC, ADDA
INTEGER NEQ, IA, JA, IC, JC, IWK, IPPER
DOUBLE PRECISION Y, S, YH, SAVR, EWT, RTEM, WK
DIMENSION NEQ(*), Y(*), S(*), YH(*), SAVR(*), EWT(*), RTEM(*), IA(*), JA(*), IC(*), JC(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IBR, IER, IPIL, IPIU, IPTT1, IPTT2, J, K, KNEW, KAMAX, KAMIN, KCMAX, KCMIN, LDIF, LENIGP, LENWK1, LIWK, LJFO
INTEGER MAXG, NP1, NZSUT
DOUBLE PRECISION ERWT, FAC, YJ
!-----------------------------------------------------------------------
! This routine performs preprocessing related to the sparse linear
! systems that must be solved.
! The operations that are performed here are:
!  * compute sparseness structure of the iteration matrix
!      P = A - con*J  according to MOSS,
!  * compute grouping of column indices (MITER = 2),
!  * compute a new ordering of rows and columns of the matrix,
!  * reorder JA corresponding to the new ordering,
!  * perform a symbolic LU factorization of the matrix, and
!  * set pointers for segments of the IWK/WK array.
! In addition to variables described previously, DPREPI uses the
! following for communication:
! YH     = the history array.  Only the first column, containing the
!          current Y vector, is used.  Used only if MOSS .ne. 0.
! S      = array of length NEQ, identical to YDOTI in the driver, used
!          only if MOSS .ne. 0.
! SAVR   = a work array of length NEQ, used only if MOSS .ne. 0.
! EWT    = array of length NEQ containing (inverted) error weights.
!          Used only if MOSS = 2 or 4 or if ISTATE = MOSS = 1.
! RTEM   = a work array of length NEQ, identical to ACOR in the driver,
!          used only if MOSS = 2 or 4.
! WK     = a real work array of length LENWK, identical to WM in
!          the driver.
! IWK    = integer work array, assumed to occupy the same space as WK.
! LENWK  = the length of the work arrays WK and IWK.
! ISTATC = a copy of the driver input argument ISTATE (= 1 on the
!          first call, = 3 on a continuation call).
! IYS    = flag value from ODRV or CDRV.
! IPPER  = output error flag , with the following values and meanings:
!        =   0  no error.
!        =  -1  insufficient storage for internal structure pointers.
!        =  -2  insufficient storage for JGROUP.
!        =  -3  insufficient storage for ODRV.
!        =  -4  other error flag from ODRV (should never occur).
!        =  -5  insufficient storage for CDRV.
!        =  -6  other error flag from CDRV.
!        =  -7  if the RES routine returned error flag IRES = IER = 2.
!        =  -8  if the RES routine returned error flag IRES = IER = 3.
!-----------------------------------------------------------------------
IBIAN = LRAT*2
IPIAN = IBIAN + 1
NP1 = N + 1
IPJAN = IPIAN + NP1
IBJAN = IPJAN - 1
LENWK1 = LENWK - N
LIWK = LENWK*LRAT
IF (MOSS .EQ. 0) LIWK = LIWK - N
IF (MOSS .EQ. 1 .OR. MOSS .EQ. 2) LIWK = LENWK1*LRAT
IF (IPJAN+N-1 .GT. LIWK) GOTO 310
IF (MOSS .EQ. 0) GOTO 30
!
IF (ISTATC .EQ. 3) GOTO 20
! ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination.
! Initialize S with random nonzero elements for structure determination.
DO 10 I=1,N
ERWT = 1.0D0/EWT(I)
FAC = 1.0D0 + 1.0D0/(I + 1.0D0)
Y(I) = Y(I) + FAC*SIGN(ERWT,Y(I))
S(I) = 1.0D0 + FAC*ERWT
10     CONTINUE
GOTO (70, 100, 150, 200), MOSS
!
20   CONTINUE
! ISTATE = 3 and MOSS .ne. 0. Load Y from YH(*,1) and S from YH(*,2). --
DO 25 I = 1,N
Y(I) = YH(I)
25      S(I) = YH(N+I)
GOTO (70, 100, 150, 200), MOSS
!
! MOSS = 0. Process user's IA,JA and IC,JC. ----------------------------
30   KNEW = IPJAN
KAMIN = IA(1)
KCMIN = IC(1)
IWK(IPIAN) = 1
DO 60 J = 1,N
DO 35 I = 1,N
35       IWK(LIWK+I) = 0
KAMAX = IA(J+1) - 1
IF (KAMIN .GT. KAMAX) GOTO 45
DO 40 K = KAMIN,KAMAX
I = JA(K)
IWK(LIWK+I) = 1
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
40       CONTINUE
45     KAMIN = KAMAX + 1
KCMAX = IC(J+1) - 1
IF (KCMIN .GT. KCMAX) GOTO 55
DO 50 K = KCMIN,KCMAX
I = JC(K)
IF (IWK(LIWK+I) .NE. 0) GOTO 50
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
50       CONTINUE
55     IWK(IPIAN+J) = KNEW + 1 - IPJAN
KCMIN = KCMAX + 1
60     CONTINUE
GOTO 240
!
! MOSS = 1. Compute structure from user-supplied Jacobian routine JAC. -
70   CONTINUE
! A dummy call to RES allows user to create temporaries for use in JAC.
IER = 1
CALL RES (NEQ, TN, Y, S, SAVR, IER)
IF (IER .GT. 1) GOTO 370
DO 75 I = 1,N
SAVR(I) = 0.0D0
75     WK(LENWK1+I) = 0.0D0
K = IPJAN
IWK(IPIAN) = 1
DO 95 J = 1,N
CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), WK(LENWK1+1))
CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), SAVR)
DO 90 I = 1,N
LJFO = LENWK1 + I
IF (WK(LJFO) .EQ. 0.0D0) GOTO 80
WK(LJFO) = 0.0D0
SAVR(I) = 0.0D0
GOTO 85
80       IF (SAVR(I) .EQ. 0.0D0) GOTO 90
SAVR(I) = 0.0D0
85       IF (K .GT. LIWK) GOTO 310
IWK(K) = I
K = K+1
90       CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
95     CONTINUE
GOTO 240
!
! MOSS = 2. Compute structure from results of N + 1 calls to RES. ------
100  DO 105 I = 1,N
105    WK(LENWK1+I) = 0.0D0
K = IPJAN
IWK(IPIAN) = 1
IER = -1
IF (MITER .EQ. 1) IER = 1
CALL RES (NEQ, TN, Y, S, SAVR, IER)
IF (IER .GT. 1) GOTO 370
DO 130 J = 1,N
CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), WK(LENWK1+1))
YJ = Y(J)
ERWT = 1.0D0/EWT(J)
Y(J) = YJ + SIGN(ERWT,YJ)
CALL RES (NEQ, TN, Y, S, RTEM, IER)
IF (IER .GT. 1) RETURN
Y(J) = YJ
DO 120 I = 1,N
LJFO = LENWK1 + I
IF (WK(LJFO) .EQ. 0.0D0) GOTO 110
WK(LJFO) = 0.0D0
GOTO 115
110      IF (RTEM(I) .EQ. SAVR(I)) GOTO 120
115      IF (K .GT. LIWK) GOTO 310
IWK(K) = I
K = K + 1
120      CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
130    CONTINUE
GOTO 240
!
! MOSS = 3. Compute structure from the user's IA/JA and JAC routine. ---
150  CONTINUE
! A dummy call to RES allows user to create temporaries for use in JAC.
IER = 1
CALL RES (NEQ, TN, Y, S, SAVR, IER)
IF (IER .GT. 1) GOTO 370
DO 155 I = 1,N
155    SAVR(I) = 0.0D0
KNEW = IPJAN
KAMIN = IA(1)
IWK(IPIAN) = 1
DO 190 J = 1,N
CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), SAVR)
KAMAX = IA(J+1) - 1
IF (KAMIN .GT. KAMAX) GOTO 170
DO 160 K = KAMIN,KAMAX
I = JA(K)
SAVR(I) = 0.0D0
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
160      CONTINUE
170    KAMIN = KAMAX + 1
DO 180 I = 1,N
IF (SAVR(I) .EQ. 0.0D0) GOTO 180
SAVR(I) = 0.0D0
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
180      CONTINUE
IWK(IPIAN+J) = KNEW + 1 - IPJAN
190    CONTINUE
GOTO 240
!
! MOSS = 4. Compute structure from user's IA/JA and N + 1 RES calls. ---
200  KNEW = IPJAN
KAMIN = IA(1)
IWK(IPIAN) = 1
IER = -1
IF (MITER .EQ. 1) IER = 1
CALL RES (NEQ, TN, Y, S, SAVR, IER)
IF (IER .GT. 1) GOTO 370
DO 235 J = 1,N
YJ = Y(J)
ERWT = 1.0D0/EWT(J)
Y(J) = YJ + SIGN(ERWT,YJ)
CALL RES (NEQ, TN, Y, S, RTEM, IER)
IF (IER .GT. 1) RETURN
Y(J) = YJ
KAMAX = IA(J+1) - 1
IF (KAMIN .GT. KAMAX) GOTO 225
DO 220 K = KAMIN,KAMAX
I = JA(K)
RTEM(I) = SAVR(I)
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
220      CONTINUE
225    KAMIN = KAMAX + 1
DO 230 I = 1,N
IF (RTEM(I) .EQ. SAVR(I)) GOTO 230
IF (KNEW .GT. LIWK) GOTO 310
IWK(KNEW) = I
KNEW = KNEW + 1
230      CONTINUE
IWK(IPIAN+J) = KNEW + 1 - IPJAN
235    CONTINUE
!
240  CONTINUE
IF (MOSS .EQ. 0 .OR. ISTATC .EQ. 3) GOTO 250
! If ISTATE = 0 or 1 and MOSS .ne. 0, restore Y from YH. ---------------
DO 245 I = 1,N
245    Y(I) = YH(I)
250  NNZ = IWK(IPIAN+N) - 1
IPPER = 0
NGP = 0
LENIGP = 0
IPIGP = IPJAN + NNZ
IF (MITER .NE. 2) GOTO 260
!
! Compute grouping of column indices (MITER = 2). ----------------------
!
MAXG = NP1
IPJGP = IPJAN + NNZ
IBJGP = IPJGP - 1
IPIGP = IPJGP + N
IPTT1 = IPIGP + NP1
IPTT2 = IPTT1 + N
LREQ = IPTT2 + N - 1
IF (LREQ .GT. LIWK) GOTO 320
CALL JGROUP (N, IWK(IPIAN), IWK(IPJAN), MAXG, NGP, IWK(IPIGP), IWK(IPJGP), IWK(IPTT1), IWK(IPTT2), IER)
IF (IER .NE. 0) GOTO 320
LENIGP = NGP + 1
!
! Compute new ordering of rows/columns of Jacobian. --------------------
260  IPR = IPIGP + LENIGP
IPC = IPR
IPIC = IPC + N
IPISP = IPIC + N
IPRSP = (IPISP-2)/LRAT + 2
IESP = LENWK + 1 - IPRSP
IF (IESP .LT. 0) GOTO 330
IBR = IPR - 1
DO 270 I = 1,N
270    IWK(IBR+I) = I
NSP = LIWK + 1 - IPISP
CALL ODRV(N, IWK(IPIAN), IWK(IPJAN), WK, IWK(IPR), IWK(IPIC), NSP, IWK(IPISP), 1, IYS)
IF (IYS .EQ. 11*N+1) GOTO 340
IF (IYS .NE. 0) GOTO 330
!
! Reorder JAN and do symbolic LU factorization of matrix. --------------
IPA = LENWK + 1 - NNZ
NSP = IPA - IPRSP
LREQ = MAX(12*N/LRAT, 6*N/LRAT+2*N+NNZ) + 3
LREQ = LREQ + IPRSP - 1 + NNZ
IF (LREQ .GT. LENWK) GOTO 350
IBA = IPA - 1
DO 280 I = 1,NNZ
280    WK(IBA+I) = 0.0D0
IPISP = LRAT*(IPRSP - 1) + 1
CALL CDRV(N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)
LREQ = LENWK - IESP
IF (IYS .EQ. 10*N+1) GOTO 350
IF (IYS .NE. 0) GOTO 360
IPIL = IPISP
IPIU = IPIL + 2*N + 1
NZU = IWK(IPIL+N) - IWK(IPIL)
NZL = IWK(IPIU+N) - IWK(IPIU)
IF (LRAT .GT. 1) GOTO 290
CALL ADJLR (N, IWK(IPISP), LDIF)
LREQ = LREQ + LDIF
290  CONTINUE
IF (LRAT .EQ. 2 .AND. NNZ .EQ. N) LREQ = LREQ + 1
NSP = NSP + LREQ - LENWK
IPA = LREQ + 1 - NNZ
IBA = IPA - 1
IPPER = 0
RETURN
!
310  IPPER = -1
LREQ = 2 + (2*N + 1)/LRAT
LREQ = MAX(LENWK+1,LREQ)
RETURN
!
320  IPPER = -2
LREQ = (LREQ - 1)/LRAT + 1
RETURN
!
330  IPPER = -3
CALL CNTNZU (N, IWK(IPIAN), IWK(IPJAN), NZSUT)
LREQ = LENWK - IESP + (3*N + 4*NZSUT - 1)/LRAT + 1
RETURN
!
340  IPPER = -4
RETURN
!
350  IPPER =  -5
RETURN
!
360  IPPER = -6
LREQ = LENWK
RETURN
!
370  IPPER = -IER - 5
LREQ = 2 + (2*N + 1)/LRAT
RETURN
!----------------------- End of Subroutine DPREPI ----------------------
END

SUBROUTINE DAINVGS (NEQ, T, Y, WK, IWK, TEM, YDOT, IER, RES, ADDA)
EXTERNAL RES, ADDA
INTEGER NEQ, IWK, IER
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMUL, J, K, KMIN, KMAX
DOUBLE PRECISION T, Y, WK, TEM, YDOT
DOUBLE PRECISION RLSS
DIMENSION Y(*), WK(*), IWK(*), TEM(*), YDOT(*)
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
!-----------------------------------------------------------------------
! This subroutine computes the initial value of the vector YDOT
! satisfying
!     A * YDOT = g(t,y)
! when A is nonsingular.  It is called by DLSODIS for initialization
! only, when ISTATE = 0.  The matrix A is subjected to LU
! decomposition in CDRV.  Then the system A*YDOT = g(t,y) is solved
! in CDRV.
! In addition to variables described previously, communication
! with DAINVGS uses the following:
! Y     = array of initial values.
! WK    = real work space for matrices.  On output it contains A and
!         its LU decomposition.  The LU decomposition is not entirely
!         sparse unless the structure of the matrix A is identical to
!         the structure of the Jacobian matrix dr/dy.
!         Storage of matrix elements starts at WK(3).
!         WK(1) = SQRT(UROUND), not used here.
! IWK   = integer work space for matrix-related data, assumed to
!         be equivalenced to WK.  In addition, WK(IPRSP) and WK(IPISP)
!         are assumed to have identical locations.
! TEM   = vector of work space of length N (ACOR in DSTODI).
! YDOT  = output vector containing the initial dy/dt. YDOT(i) contains
!         dy(i)/dt when the matrix A is non-singular.
! IER   = output error flag with the following values and meanings:
!       = 0  if DAINVGS was successful.
!       = 1  if the A-matrix was found to be singular.
!       = 2  if RES returned an error flag IRES = IER = 2.
!       = 3  if RES returned an error flag IRES = IER = 3.
!       = 4  if insufficient storage for CDRV (should not occur here).
!       = 5  if other error found in CDRV (should not occur here).
!-----------------------------------------------------------------------
!
DO 10 I = 1,NNZ
10     WK(IBA+I) = 0.0D0
!
IER = 1
CALL RES (NEQ, T, Y, WK(IPA), YDOT, IER)
IF (IER .GT. 1) RETURN
!
KMIN = IWK(IPIAN)
DO 30 J = 1,NEQ
KMAX = IWK(IPIAN+J) - 1
DO 15 K = KMIN,KMAX
I = IWK(IBJAN+K)
15       TEM(I) = 0.0D0
CALL ADDA (NEQ, T, Y, J, IWK(IPIAN), IWK(IPJAN), TEM)
DO 20 K = KMIN,KMAX
I = IWK(IBJAN+K)
20       WK(IBA+K) = TEM(I)
KMIN = KMAX + 1
30   CONTINUE
NLU = NLU + 1
IER = 0
DO 40 I = 1,NEQ
40     TEM(I) = 0.0D0
!
! Numerical factorization of matrix A. ---------------------------------
CALL CDRV (NEQ,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),TEM,TEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
IF (IYS .EQ. 0) GOTO 50
IMUL = (IYS - 1)/NEQ
IER = 5
IF (IMUL .EQ. 8) IER = 1
IF (IMUL .EQ. 10) IER = 4
RETURN
!
! Solution of the linear system. ---------------------------------------
50   CALL CDRV (NEQ,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),YDOT,YDOT,NSP,IWK(IPISP),WK(IPRSP),IESP,4,IYS)
IF (IYS .NE. 0) IER = 5
RETURN
!----------------------- End of Subroutine DAINVGS ---------------------
END

SUBROUTINE DPRJIS (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WK, IWK, RES, JAC, ADDA)
EXTERNAL RES, JAC, ADDA
INTEGER NEQ, NYH, IWK
DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WK
DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*), S(*), SAVR(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM
INTEGER LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP
INTEGER IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMUL, IRES, J, JJ, JMAX, JMIN, K, KMAX, KMIN, NG
DOUBLE PRECISION CON, FAC, HL0, R, SRUR
!-----------------------------------------------------------------------
! DPRJIS is called to compute and process the matrix
! P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
! where r = g(t,y) - A(t,y)*s.  J is computed by columns, either by
! the user-supplied routine JAC if MITER = 1, or by finite differencing
! if MITER = 2.  J is stored in WK, rescaled, and ADDA is called to
! generate P.  The matrix P is subjected to LU decomposition in CDRV.
! P and its LU decomposition are stored separately in WK.
!
! In addition to variables described previously, communication
! with DPRJIS uses the following:
! Y     = array containing predicted values on entry.
! RTEM  = work array of length N (ACOR in DSTODI).
! SAVR  = array containing r evaluated at predicted y. On output it
!         contains the residual evaluated at current values of t and y.
! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
! WK    = real work space for matrices.  On output it contains P and
!         its sparse LU decomposition.  Storage of matrix elements
!         starts at WK(3).
!         WK also contains the following matrix-related data.
!         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
! IWK   = integer work space for matrix-related data, assumed to be
!         equivalenced to WK.  In addition,  WK(IPRSP) and IWK(IPISP)
!         are assumed to have identical locations.
! EL0   = EL(1) (input).
! IERPJ = output error flag (in COMMON).
!         =  0 if no error.
!         =  1 if zero pivot found in CDRV.
!         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
!         = -1 if insufficient storage for CDRV (should not occur).
!         = -2 if other error found in CDRV (should not occur here).
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
HL0 = H*EL0
CON = -HL0
JCUR = 1
NJE = NJE + 1
GOTO (100, 200), MITER
!
! If MITER = 1, call RES, then call JAC and ADDA for each column. ------
100  IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
KMIN = IWK(IPIAN)
DO 130 J = 1,N
KMAX = IWK(IPIAN+J)-1
DO 110 I = 1,N
110      RTEM(I) = 0.0D0
CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), RTEM)
DO 120 I = 1,N
120      RTEM(I) = RTEM(I)*CON
CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), RTEM)
DO 125 K = KMIN,KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = RTEM(I)
125      CONTINUE
KMIN = KMAX + 1
130    CONTINUE
GOTO 290
!
! If MITER = 2, make NGP + 1 calls to RES to approximate J and P. ------
200  CONTINUE
IRES = -1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
SRUR = WK(1)
JMIN = IWK(IPIGP)
DO 240 NG = 1,NGP
JMAX = IWK(IPIGP+NG) - 1
DO 210 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
R = MAX(SRUR*ABS(Y(JJ)),0.01D0/EWT(JJ))
210      Y(JJ) = Y(JJ) + R
CALL RES (NEQ,TN,Y,S,RTEM,IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
DO 230 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
Y(JJ) = YH(JJ,1)
R = MAX(SRUR*ABS(Y(JJ)),0.01D0/EWT(JJ))
FAC = -HL0/R
KMIN = IWK(IBIAN+JJ)
KMAX = IWK(IBIAN+JJ+1) - 1
DO 220 K = KMIN,KMAX
I = IWK(IBJAN+K)
RTEM(I) = (RTEM(I) - SAVR(I))*FAC
220        CONTINUE
CALL ADDA (NEQ, TN, Y, JJ, IWK(IPIAN), IWK(IPJAN), RTEM)
DO 225 K = KMIN,KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = RTEM(I)
225      CONTINUE
230      CONTINUE
JMIN = JMAX + 1
240    CONTINUE
IRES = 1
CALL RES (NEQ, TN, Y, S, SAVR, IRES)
NFE = NFE + 1
IF (IRES .GT. 1) GOTO 600
!
! Do numerical factorization of P matrix. ------------------------------
290  NLU = NLU + 1
IERPJ = 0
DO 295 I = 1,N
295    RTEM(I) = 0.0D0
CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),RTEM,RTEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
IF (IYS .EQ. 0) RETURN
IMUL = (IYS - 1)/N
IERPJ = -2
IF (IMUL .EQ. 8) IERPJ = 1
IF (IMUL .EQ. 10) IERPJ = -1
RETURN
! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
600  IERPJ = IRES
RETURN
!----------------------- End of Subroutine DPRJIS ----------------------
END

