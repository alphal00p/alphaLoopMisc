ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        MP__G = 2 * SQRT(MP__AS*MP__PI)  ! for the first init

        MP__MDL_CONJG__CKM3X3 = 1.000000E+00_16

        MP__MDL_CONJG__CKM1X1 = 1.000000E+00_16

        MP__MDL_CKM3X3 = 1.000000E+00_16

        MP__MDL_COMPLEXI = CMPLX(0.000000E+00_16,1.000000E+00_16
     $   ,KIND=16)

        MP__MDL_MZ__EXP__2 = MP__MDL_MZ**2

        MP__MDL_MZ__EXP__4 = MP__MDL_MZ**4

        MP__MDL_SQRT__2 = SQRT(CMPLX((2.000000E+00_16),KIND=16))

        MP__MDL_MH__EXP__2 = MP__MDL_MH**2

        MP__MDL_AEW = 1.000000E+00_16/MP__AEWM1

        MP__MDL_MW = SQRT(CMPLX((MP__MDL_MZ__EXP__2/2.000000E+00_16
     $   +SQRT(CMPLX((MP__MDL_MZ__EXP__4/4.000000E+00_16-(MP__MDL_AEW
     $   *MP__PI*MP__MDL_MZ__EXP__2)/(MP__MDL_GF*MP__MDL_SQRT__2))
     $   ,KIND=16))),KIND=16))

        MP__MDL_SQRT__AEW = SQRT(CMPLX((MP__MDL_AEW),KIND=16))

        MP__MDL_EE = 2.000000E+00_16*MP__MDL_SQRT__AEW
     $   *SQRT(CMPLX((MP__PI),KIND=16))

        MP__MDL_MW__EXP__2 = MP__MDL_MW**2

        MP__MDL_SW2 = 1.000000E+00_16-MP__MDL_MW__EXP__2
     $   /MP__MDL_MZ__EXP__2

        MP__MDL_CW = SQRT(CMPLX((1.000000E+00_16-MP__MDL_SW2),KIND=16))

        MP__MDL_SQRT__SW2 = SQRT(CMPLX((MP__MDL_SW2),KIND=16))

        MP__MDL_SW = MP__MDL_SQRT__SW2

        MP__MDL_G1 = MP__MDL_EE/MP__MDL_CW

        MP__MDL_GW = MP__MDL_EE/MP__MDL_SW

        MP__MDL_VEV = (2.000000E+00_16*MP__MDL_MW*MP__MDL_SW)
     $   /MP__MDL_EE

        MP__MDL_VEV__EXP__2 = MP__MDL_VEV**2

        MP__MDL_LAM = MP__MDL_MH__EXP__2/(2.000000E+00_16
     $   *MP__MDL_VEV__EXP__2)

        MP__MDL_YB = (MP__MDL_YMB*MP__MDL_SQRT__2)/MP__MDL_VEV

        MP__MDL_YT = (MP__MDL_YMT*MP__MDL_SQRT__2)/MP__MDL_VEV

        MP__MDL_YTAU = (MP__MDL_YMTAU*MP__MDL_SQRT__2)/MP__MDL_VEV

        MP__MDL_MUH = SQRT(CMPLX((MP__MDL_LAM*MP__MDL_VEV__EXP__2)
     $   ,KIND=16))

        MP__MDL_I1X33 = MP__MDL_YB*MP__MDL_CONJG__CKM3X3

        MP__MDL_I2X33 = MP__MDL_YT*MP__MDL_CONJG__CKM3X3

        MP__MDL_I3X33 = MP__MDL_CKM3X3*MP__MDL_YT

        MP__MDL_I4X33 = MP__MDL_CKM3X3*MP__MDL_YB

        MP__MDL_EE__EXP__2 = MP__MDL_EE**2

        MP__MDL_SW__EXP__2 = MP__MDL_SW**2

        MP__MDL_CW__EXP__2 = MP__MDL_CW**2

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      MP__AS = MP__G**2/4/MP__PI

      MP__MDL_SQRT__AS = SQRT(CMPLX((MP__AS),KIND=16))

      MP__MDL_G__EXP__2 = MP__G**2

C     
C     Parameters that should be updated for the loops.
C     
C     
C     Definition of the EW coupling used in the write out of aqed
C     
      MP__GAL(1) = 2 * SQRT(MP__PI/ABS(MP__AEWM1))
      MP__GAL(2) = 1D0

