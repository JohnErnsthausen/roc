      PROGRAM SAFMIN
      DOUBLE PRECISION EPMACH, SAFEMIN

      EPMACH = 0.0D0
      SAFEMIN = 0.0D0

      CALL DMACH(EPMACH,SAFEMIN)

      WRITE(6,50) EPMACH,SAFEMIN
   50 FORMAT(' EPMACH= ',G23.16,' SAFEMIN= ',G23.16)
      STOP
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE DMACH( EPMACH, SAFMIN )
C
      DOUBLE PRECISION EPMACH, SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  To determine the relative machine precision EPMACH and the 
C  safe minimum SAFMIN.
C
C  Variables in the calling sequence
C  ---------------------------------
C  EPMACH D   OUT   The smallest floating point number such that
C                   1.0D0+EPS .NE. 1.0D0
C  SAFMIN D   OUT   SAFMIN = BETA**M where M is the largest in 
C                   magnitude negative integer M such that BETA**M 
C                   is positive and normalized
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External function
C
      DOUBLE PRECISION ABS
C
C.....Parameters
C
      DOUBLE PRECISION ZER, ONE
      PARAMETER( ZER=0.0D0, ONE=1.0D0 )
C
C.....Local variables
C
      INTEGER I,IBETA,IT,N
      DOUBLE PRECISION A,B,BETA,BETAIN,C,D,T,Y,Z
C
C.....Variables to be saved
C
      LOGICAL DONE
      DOUBLE PRECISION EPMCH,SFMIN
      SAVE DONE,EPMCH,SFMIN
C
C.....Data statements
C
      DATA DONE/.FALSE./, EPMCH/0.0D0/, SFMIN/0.0D0/
C
C.......................Executable statements...........................
C
C.....If this is not the first call, return the previously
C.....computed values
C
      IF( DONE ) THEN
         EPMACH = EPMCH
         SAFMIN = SFMIN
         RETURN
      ENDIF
      DONE = .TRUE.
C
C.....Determine the base of the floating point representation
C
      A = ONE
   10 CONTINUE
         A = A + A
         B = A + ONE
         C = B - A
      IF (C-ONE .EQ. ZER) GOTO 10
      B = ONE
   20 CONTINUE
         B = B + B
         C = A + B
         IBETA = INT(C - A)
      IF (IBETA .EQ. 0) GOTO 20
      BETA = DBLE(IBETA)
C
C.....Determine the number of digits in the representation
C
      IT = 0
      B = ONE
   30 CONTINUE
         IT = IT + 1
         B = B * BETA
         C = B + ONE
         D = C - B
      IF (D-ONE .EQ. ZER) GOTO 30
C
C.....Determine the smallest positive floating-point number 
C.....such that 1.0D0 + EPS .NE. 1.0D0
C
      N = IT + 3
      BETAIN = ONE/BETA
      A = ONE
      DO 40 I = 1, N
         A = A * BETAIN
   40 CONTINUE
C
   50 CONTINUE
         B = ONE + A
         IF (B - ONE .NE. ZER) GOTO 60
         A = A * BETA
         GOTO 50
   60 CONTINUE 
      EPMCH = A
C
C.....Determine the largest in magnitude negative integer M 
C.....such that BETA**M is positive and normalized and set 
C.....SAFMIN = BETA**M.
C
      Z = BETAIN
      T = ONE + EPMCH
   70 CONTINUE
         Y = Z
         Z = Y * Y
         A = Z * ONE
         C = Z * T
         IF ((A+A .EQ. ZER) .OR. (ABS(Z) .GE. Y)) GOTO 80
         D = C * BETAIN
         IF (D*BETA .EQ. Z) GOTO 80
         GOTO 70
   80 CONTINUE
         SFMIN = Y
         Y = Y * BETAIN
         A = Y * ONE
         C = Y * T
         IF (((A+A) .EQ. ZER) .OR. (ABS(Y) .GE. SFMIN)) GOTO 90
         D = C * BETAIN
         IF ((D*BETA .NE. Y) .OR. (C .EQ. Y)) GOTO 80
         SFMIN = Y
   90 CONTINUE
C
      EPMACH = EPMCH
      SAFMIN = SFMIN
C
      RETURN
C
C.....End of DMACH
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
