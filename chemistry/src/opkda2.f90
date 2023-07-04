
SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
!***BEGIN PROLOGUE  DGEFA
!***PURPOSE  Factor a matrix using Gaussian elimination.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGEFA factors a double precision matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
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
!     On Return
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
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGEFA
INTEGER LDA,N,IPVT(*),INFO
DOUBLE PRECISION A(LDA,*)
!
DOUBLE PRECISION T
INTEGER IDAMAX,J,K,KP1,L,NM1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  DGEFA
INFO = 0
NM1 = N - 1
IF (NM1 .LT. 1) GOTO 70
DO 60 K = 1, NM1
KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
L = IDAMAX(N-K+1,A(K,K),1) + K - 1
IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
IF (A(L,K) .EQ. 0.0D0) GOTO 40
!
!           INTERCHANGE IF NECESSARY
!
IF (L .EQ. K) GOTO 10
T = A(L,K)
A(L,K) = A(K,K)
A(K,K) = T
10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
T = -1.0D0/A(K,K)
CALL DSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
DO 30 J = KP1, N
T = A(L,J)
IF (L .EQ. K) GOTO 20
A(L,J) = A(K,J)
A(K,J) = T
20          CONTINUE
CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
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
END

SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
!***BEGIN PROLOGUE  DGESL
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
!            factors computed by DGECO or DGEFA.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GOTO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL
INTEGER LDA,N,IPVT(*),JOB
DOUBLE PRECISION A(LDA,*),B(*)
!
DOUBLE PRECISION DDOT,T
INTEGER K,KB,L,NM1
!***FIRST EXECUTABLE STATEMENT  DGESL
NM1 = N - 1
IF (JOB .NE. 0) GOTO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
IF (NM1 .LT. 1) GOTO 30
DO 20 K = 1, NM1
L = IPVT(K)
T = B(L)
IF (L .EQ. K) GOTO 10
B(L) = B(K)
B(K) = T
10       CONTINUE
CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
20    CONTINUE
30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
DO 40 KB = 1, N
K = N + 1 - KB
B(K) = B(K)/A(K,K)
T = -B(K)
CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
40    CONTINUE
GOTO 100
50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
DO 60 K = 1, N
T = DDOT(K-1,A(1,K),1,B(1),1)
B(K) = (B(K) - T)/A(K,K)
60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
IF (NM1 .LT. 1) GOTO 90
DO 80 KB = 1, NM1
K = N - KB
B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
L = IPVT(K)
IF (L .EQ. K) GOTO 70
T = B(L)
B(L) = B(K)
B(K) = T
70       CONTINUE
80    CONTINUE
90    CONTINUE
100 CONTINUE
RETURN
END

SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
!***BEGIN PROLOGUE  DGBFA
!***PURPOSE  Factor a band matrix using Gaussian elimination.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBFA factors a double precision band matrix by elimination.
!
!     DGBFA is usually called by DGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .
!     On Return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL will divide by zero if
!                     called.  Use  RCOND  in DGBCO for a reliable
!                     indication of singularity.
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBFA
INTEGER LDA,N,ML,MU,IPVT(*),INFO
DOUBLE PRECISION ABD(LDA,*)
!
DOUBLE PRECISION T
INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
!
!***FIRST EXECUTABLE STATEMENT  DGBFA
M = ML + MU + 1
INFO = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
J0 = MU + 2
J1 = MIN(N,M) - 1
IF (J1 .LT. J0) GOTO 30
DO 20 JZ = J0, J1
I0 = M + 1 - JZ
DO 10 I = I0, ML
ABD(I,JZ) = 0.0D0
10    CONTINUE
20 CONTINUE
30 CONTINUE
JZ = J1
JU = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
NM1 = N - 1
IF (NM1 .LT. 1) GOTO 130
DO 120 K = 1, NM1
KP1 = K + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
JZ = JZ + 1
IF (JZ .GT. N) GOTO 50
IF (ML .LT. 1) GOTO 50
DO 40 I = 1, ML
ABD(I,JZ) = 0.0D0
40       CONTINUE
50    CONTINUE
!
!        FIND L = PIVOT INDEX
!
LM = MIN(ML,N-K)
L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
IPVT(K) = L + K - M
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
IF (ABD(L,K) .EQ. 0.0D0) GOTO 100
!
!           INTERCHANGE IF NECESSARY
!
IF (L .EQ. M) GOTO 60
T = ABD(L,K)
ABD(L,K) = ABD(M,K)
ABD(M,K) = T
60       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
T = -1.0D0/ABD(M,K)
CALL DSCAL(LM,T,ABD(M+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
JU = MIN(MAX(JU,MU+IPVT(K)),N)
MM = M
IF (JU .LT. KP1) GOTO 90
DO 80 J = KP1, JU
L = L - 1
MM = MM - 1
T = ABD(L,J)
IF (L .EQ. MM) GOTO 70
ABD(L,J) = ABD(MM,J)
ABD(MM,J) = T
70          CONTINUE
CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
80       CONTINUE
90       CONTINUE
GOTO 110
100    CONTINUE
INFO = K
110    CONTINUE
120 CONTINUE
130 CONTINUE
IPVT(N) = N
IF (ABD(M,N) .EQ. 0.0D0) INFO = N
RETURN
END

SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
!***BEGIN PROLOGUE  DGBSL
!***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
!            the factors computed by DGBCO or DGBFA.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DGBCO or DGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND .GT. 0.0
!        or DGBFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GOTO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBSL
INTEGER LDA,N,ML,MU,IPVT(*),JOB
DOUBLE PRECISION ABD(LDA,*),B(*)
!
DOUBLE PRECISION DDOT,T
INTEGER K,KB,L,LA,LB,LM,M,NM1
!***FIRST EXECUTABLE STATEMENT  DGBSL
M = MU + ML + 1
NM1 = N - 1
IF (JOB .NE. 0) GOTO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
IF (ML .EQ. 0) GOTO 30
IF (NM1 .LT. 1) GOTO 30
DO 20 K = 1, NM1
LM = MIN(ML,N-K)
L = IPVT(K)
T = B(L)
IF (L .EQ. K) GOTO 10
B(L) = B(K)
B(K) = T
10          CONTINUE
CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
20       CONTINUE
30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
DO 40 KB = 1, N
K = N + 1 - KB
B(K) = B(K)/ABD(M,K)
LM = MIN(K,M) - 1
LA = M - LM
LB = K - LM
T = -B(K)
CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
40    CONTINUE
GOTO 100
50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
DO 60 K = 1, N
LM = MIN(K,M) - 1
LA = M - LM
LB = K - LM
T = DDOT(LM,ABD(LA,K),1,B(LB),1)
B(K) = (B(K) - T)/ABD(M,K)
60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
IF (ML .EQ. 0) GOTO 90
IF (NM1 .LT. 1) GOTO 90
DO 80 KB = 1, NM1
K = N - KB
LM = MIN(ML,N-K)
B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
L = IPVT(K)
IF (L .EQ. K) GOTO 70
T = B(L)
B(L) = B(K)
B(K) = T
70          CONTINUE
80       CONTINUE
90    CONTINUE
100 CONTINUE
RETURN
END

SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)
!
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DAXPY
integer n, incx, incy, i, ix, iy, m, mp1, ns
DOUBLE PRECISION DX(*), DY(*), DA
!***FIRST EXECUTABLE STATEMENT  DAXPY
IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
5 IX = 1
IY = 1
IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
DO 10 I = 1,N
DY(IY) = DY(IY) + DA*DX(IX)
IX = IX + INCX
IY = IY + INCY
10 CONTINUE
RETURN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 4.
!
20 M = MOD(N,4)
IF (M .EQ. 0) GOTO 40
DO 30 I = 1,M
DY(I) = DY(I) + DA*DX(I)
30 CONTINUE
IF (N .LT. 4) RETURN
40 MP1 = M + 1
DO 50 I = MP1,N,4
DY(I) = DY(I) + DA*DX(I)
DY(I+1) = DY(I+1) + DA*DX(I+1)
DY(I+2) = DY(I+2) + DA*DX(I+2)
DY(I+3) = DY(I+3) + DA*DX(I+3)
50 CONTINUE
RETURN
!
!     Code for equal, positive, non-unit increments.
!
60 NS = N*INCX
DO 70 I = 1,NS,INCX
DY(I) = DA*DX(I) + DY(I)
70 CONTINUE
RETURN
END

SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DCOPY
!***PURPOSE  Copy a vector.
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
!***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  copy of vector DX (unchanged if N .LE. 0)
!
!     Copy double precision DX to double precision DY.
!     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCOPY
integer n, incx, incy, i, ix, iy, m, mp1, ns
DOUBLE PRECISION DX(*), DY(*)
!***FIRST EXECUTABLE STATEMENT  DCOPY
IF (N .LE. 0) RETURN
IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
5 IX = 1
IY = 1
IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
DO 10 I = 1,N
DY(IY) = DX(IX)
IX = IX + INCX
IY = IY + INCY
10 CONTINUE
RETURN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 7.
!
20 M = MOD(N,7)
IF (M .EQ. 0) GOTO 40
DO 30 I = 1,M
DY(I) = DX(I)
30 CONTINUE
IF (N .LT. 7) RETURN
40 MP1 = M + 1
DO 50 I = MP1,N,7
DY(I) = DX(I)
DY(I+1) = DX(I+1)
DY(I+2) = DX(I+2)
DY(I+3) = DX(I+3)
DY(I+4) = DX(I+4)
DY(I+5) = DX(I+5)
DY(I+6) = DX(I+6)
50 CONTINUE
RETURN
!
!     Code for equal, positive, non-unit increments.
!
60 NS = N*INCX
DO 70 I = 1,NS,INCX
DY(I) = DX(I)
70 CONTINUE
RETURN
END

DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DDOT
!***PURPOSE  Compute the inner product of two vectors.
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!     DDOT  double precision dot product (zero if N .LE. 0)
!
!     Returns the dot product of double precision DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DDOT
integer n, incx, incy, i, ix, iy, m, mp1, ns
DOUBLE PRECISION DX(*), DY(*)
!***FIRST EXECUTABLE STATEMENT  DDOT
DDOT = 0.0D0
IF (N .LE. 0) RETURN
IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
5 IX = 1
IY = 1
IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
DO 10 I = 1,N
DDOT = DDOT + DX(IX)*DY(IY)
IX = IX + INCX
IY = IY + INCY
10 CONTINUE
RETURN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
20 M = MOD(N,5)
IF (M .EQ. 0) GOTO 40
DO 30 I = 1,M
DDOT = DDOT + DX(I)*DY(I)
30 CONTINUE
IF (N .LT. 5) RETURN
40 MP1 = M + 1
DO 50 I = MP1,N,5
DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
50 CONTINUE
RETURN
!
!     Code for equal, positive, non-unit increments.
!
60 NS = N*INCX
DO 70 I = 1,NS,INCX
DDOT = DDOT + DX(I)*DY(I)
70 CONTINUE
RETURN
END

DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)
!***BEGIN PROLOGUE  DNRM2
!***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
!***CATEGORY  D1A3B
!***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
!***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!             LINEAR ALGEBRA, UNITARY, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!    DNRM2  double precision result (zero if N .LE. 0)
!
!     Euclidean norm of the N-vector stored in DX with storage
!     increment INCX.
!     If N .LE. 0, return with result = 0.
!     If N .GE. 1, then INCX must be .GE. 1
!
!     Four phase method using two built-in constants that are
!     hopefully applicable to all machines.
!         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
!         CUTHI = minimum of  SQRT(V)      over all known machines.
!     where
!         EPS = smallest no. such that EPS + 1. .GT. 1.
!         U   = smallest positive no.   (underflow limit)
!         V   = largest  no.            (overflow  limit)
!
!     Brief outline of algorithm.
!
!     Phase 1 scans zero components.
!     move to phase 2 when a component is nonzero and .LE. CUTLO
!     move to phase 3 when a component is .GT. CUTLO
!     move to phase 4 when a component is .GE. CUTHI/M
!     where M = N for X() real and M = 2*N for complex.
!
!     Values for CUTLO and CUTHI.
!     From the environmental parameters listed in the IMSL converter
!     document the limiting values are as follows:
!     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
!                   Univac and DEC at 2**(-103)
!                   Thus CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
!                   Thus CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
!                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNRM2
INTEGER NEXT
integer n, incx, i, j, nn
DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
SAVE CUTLO, CUTHI, ZERO, ONE
DATA ZERO, ONE /0.0D0, 1.0D0/
!
DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!***FIRST EXECUTABLE STATEMENT  DNRM2
IF (N .GT. 0) GOTO 10
DNRM2  = ZERO
GOTO 300
!
10 NEXT = 30
SUM = ZERO
NN = N * INCX
!
!                                                 BEGIN MAIN LOOP
!
I = 1
20    GOTO (30, 50, 70, 110) NEXT
30 IF (ABS(DX(I)) .GT. CUTLO) GOTO 85
NEXT = 50
XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
50 IF (DX(I) .EQ. ZERO) GOTO 200
IF (ABS(DX(I)) .GT. CUTLO) GOTO 85
!
!                                PREPARE FOR PHASE 2.
!
NEXT = 70
GOTO 105
!
!                                PREPARE FOR PHASE 4.
!
100 I = J
NEXT = 110
SUM = (SUM / DX(I)) / DX(I)
105 XMAX = ABS(DX(I))
GOTO 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
70 IF (ABS(DX(I)) .GT. CUTLO) GOTO 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
110 IF (ABS(DX(I)) .LE. XMAX) GOTO 115
SUM = ONE + SUM * (XMAX / DX(I))**2
XMAX = ABS(DX(I))
GOTO 200
!
115 SUM = SUM + (DX(I)/XMAX)**2
GOTO 200
!
!                  PREPARE FOR PHASE 3.
!
75 SUM = (SUM * XMAX) * XMAX
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
85 HITEST = CUTHI / N
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
DO 95 J = I,NN,INCX
IF (ABS(DX(J)) .GE. HITEST) GOTO 100
95    SUM = SUM + DX(J)**2
DNRM2 = SQRT(SUM)
GOTO 300
!
200 CONTINUE
I = I + INCX
IF (I .LE. NN) GOTO 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
DNRM2 = XMAX * SQRT(SUM)
300 CONTINUE
RETURN
END

SUBROUTINE DSCAL (N, DA, DX, INCX)
!***BEGIN PROLOGUE  DSCAL
!***PURPOSE  Multiply a vector by a constant.
!***CATEGORY  D1A6
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  double precision result (unchanged if N.LE.0)
!
!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSCAL
DOUBLE PRECISION DA, DX(*)
INTEGER I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  DSCAL
IF (N .LE. 0) RETURN
IF (INCX .EQ. 1) GOTO 20
!
!     Code for increment not equal to 1.
!
IX = 1
IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
DO 10 I = 1,N
DX(IX) = DA*DX(IX)
IX = IX + INCX
10 CONTINUE
RETURN
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
20 M = MOD(N,5)
IF (M .EQ. 0) GOTO 40
DO 30 I = 1,M
DX(I) = DA*DX(I)
30 CONTINUE
IF (N .LT. 5) RETURN
40 MP1 = M + 1
DO 50 I = MP1,N,5
DX(I) = DA*DX(I)
DX(I+1) = DA*DX(I+1)
DX(I+2) = DA*DX(I+2)
DX(I+3) = DA*DX(I+3)
DX(I+4) = DA*DX(I+4)
50 CONTINUE
RETURN
END

INTEGER FUNCTION IDAMAX (N, DX, INCX)
!***BEGIN PROLOGUE  IDAMAX
!***PURPOSE  Find the smallest index of that component of a vector
!            having the maximum magnitude.
!***CATEGORY  D1A2
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IDAMAX
DOUBLE PRECISION DX(*), DMAX, XMAG
INTEGER I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  IDAMAX
IDAMAX = 0
IF (N .LE. 0) RETURN
IDAMAX = 1
IF (N .EQ. 1) RETURN
!
IF (INCX .EQ. 1) GOTO 20
!
!     Code for increments not equal to 1.
!
IX = 1
IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
DMAX = ABS(DX(IX))
IX = IX + INCX
DO 10 I = 2,N
XMAG = ABS(DX(IX))
IF (XMAG .GT. DMAX) THEN
IDAMAX = I
DMAX = XMAG
ENDIF
IX = IX + INCX
10 CONTINUE
RETURN
!
!     Code for increments equal to 1.
!
20 DMAX = ABS(DX(1))
DO 30 I = 2,N
XMAG = ABS(DX(I))
IF (XMAG .GT. DMAX) THEN
IDAMAX = I
DMAX = XMAG
ENDIF
30 CONTINUE
RETURN
END

SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.
!
DOUBLE PRECISION R1, R2
INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
CHARACTER*(*) MSG
!
!  Declare local variables.
!
INTEGER LUNIT, IXSAV, MESFLG
!
!  Get logical unit number and message print flag.
!
!***FIRST EXECUTABLE STATEMENT  XERRWD
LUNIT = IXSAV (1, 0, .FALSE.)
MESFLG = IXSAV (2, 0, .FALSE.)
IF (MESFLG .EQ. 0) GOTO 100
!
!  Write the message.
!
WRITE (LUNIT,10)  MSG
10   FORMAT(1X,A)
IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
20   FORMAT(6X,'In above message,  I1 =',I10)
IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
40   FORMAT(6X,'In above message,  R1 =',D21.13)
IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
!
!  Abort the run if LEVEL = 2.
!
100  IF (LEVEL .NE. 2) RETURN
STOP
!----------------------- End of Subroutine XERRWD ----------------------
END

SUBROUTINE XSETF (MFLAG)
!***BEGIN PROLOGUE  XSETF
!***PURPOSE  Reset the error print control flag.
!***CATEGORY  R3A
!***TYPE      ALL (XSETF-A)
!***KEYWORDS  ERROR CONTROL
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!   XSETF sets the error print control flag to MFLAG:
!      MFLAG=1 means print all messages (the default).
!      MFLAG=0 means no printing.
!
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETF
!
! Subroutines called by XSETF.. None
! Function routine called by XSETF.. IXSAV
!-----------------------------------------------------------------------
!**End
INTEGER MFLAG, JUNK, IXSAV
!
!***FIRST EXECUTABLE STATEMENT  XSETF
IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
RETURN
!----------------------- End of Subroutine XSETF -----------------------
END

SUBROUTINE XSETUN (LUN)
!***BEGIN PROLOGUE  XSETUN
!***PURPOSE  Reset the logical unit number for error messages.
!***CATEGORY  R3B
!***TYPE      ALL (XSETUN-A)
!***KEYWORDS  ERROR CONTROL
!***DESCRIPTION
!
!   XSETUN sets the logical unit number for error messages to LUN.
!
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETUN
!
! Subroutines called by XSETUN.. None
! Function routine called by XSETUN.. IXSAV
!-----------------------------------------------------------------------
!**End
INTEGER LUN, JUNK, IXSAV
!
!***FIRST EXECUTABLE STATEMENT  XSETUN
IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
RETURN
!----------------------- End of Subroutine XSETUN ----------------------
END

INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!    LUNIT, the logical unit number to which messages are printed, and
!    MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!   LUNIT  = Logical unit number for messages.  The default is obtained
!            by a call to IUMACH (may be machine-dependent).
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  On input..
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!    ISET   = Logical flag to indicate whether to read or write.
!             If ISET = .TRUE., the parameter will be given
!             the value IVALUE.  If ISET = .FALSE., the parameter
!             will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!    IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Modified prologue to SLATEC format. (FNF)
!   930915  Added IUMACH call to get default output unit.  (ACH)
!   930922  Minor cosmetic changes. (FNF)
!   010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
LOGICAL ISET
INTEGER IPAR, IVALUE
!-----------------------------------------------------------------------
INTEGER IUMACH, LUNIT, MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
SAVE LUNIT, MESFLG
DATA LUNIT/-1/, MESFLG/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
IF (IPAR .EQ. 1) THEN
IF (LUNIT .EQ. -1) LUNIT = IUMACH()
IXSAV = LUNIT
IF (ISET) LUNIT = IVALUE
ENDIF
!
IF (IPAR .EQ. 2) THEN
IXSAV = MESFLG
IF (ISET) MESFLG = IVALUE
ENDIF
!
RETURN
!----------------------- End of Function IXSAV -------------------------
END

INTEGER FUNCTION IUMACH()
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
IUMACH = 6
!
RETURN
!----------------------- End of Function IUMACH ------------------------
END
