	MODULE MOD_XRAY_FLUXES
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: X_NU(:)
	REAL*8, ALLOCATABLE :: LOG_X_TEMP(:)
	REAL*8, ALLOCATABLE :: LOG_X_ED(:)
	REAL*8, ALLOCATABLE :: XRAY_FLUXES(:,:,:)
!
	REAL*8, ALLOCATABLE :: X_EMISS1(:)
	REAL*8, ALLOCATABLE :: X_EMISS2(:)
!               
	REAL*8 BIN_MIN
	REAL*8 BIN_SIZE
	REAL*8 LOG_T_MIN
	REAL*8 DEL_LOG_T
	REAL*8 LOG_ED_MIN
	REAL*8 DEL_LOG_ED
!
	REAL*8 T_SHOCK1_SAV
	REAL*8 T_SHOCK2_SAV
!
	INTEGER N_BINS
	INTEGER N_TEMP
	INTEGER N_ED
!
	END MODULE MOD_XRAY_FLUXES
!
!
! Routine to return X-ray EMISSIVITIES for a set of NFREQ frequencies,
! and for 2 different shock temperatures (in units of 10^4 K).
! The returned emissivities have units ergs/cm^3/s/Hz/steradian.
! At present, we assume that the X-ray emissivity is independent of density.
!
! On the first call the tabulated RS data is read in from a file 
! RS_XRAY_FLUXES. FREQ may be monotonically increasing, or decreasing.
!
	PROGRAM RD_XRAY_SPEC
	USE MOD_XRAY_FLUXES
	IMPLICIT NONE
!
	INTEGER LU_IN
	INTEGER LU_OUT
!
! Local variables:
!
	REAL*8, PARAMETER :: EV_TO_HZ=0.241838D0
	REAL*8 T1
!
	REAL*8 LOG_T_SHOCK1
	REAL*8 LOG_T_SHOCK2
!
	INTEGER I,J
	INTEGER IOS
	INTEGER T_INDX
	INTEGER ED_INDX
!
	CHARACTER*132 STRING
!
	REAL*8 FUN_PI,PI
	INTEGER ERROR_LU,LU_ER
	EXTERNAL ERROR_LU,FUN_PI
!
	LU_IN=7
	LU_OUT=10
	LU_ER=ERROR_LU()
!
! Read in RS data table. The fluxes are assumed to be units of
! 10^{-23} ergs/cm3/sec. They represent the cooling function per
! unit density of electrons and H ions.
!                      
	  OPEN(UNIT=LU_IN,FILE='RS_XRAY_FLUXES',ACTION='READ',
	1       STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error opening RS_XRAY_FLUXES'
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    STOP
	  END IF
!
! RS data assumed to be tabulated in bins equally spaced in eV, and ordered
! in increasing eV.
!
	  CALL RD_INT(N_BINS,'N_BINS',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(BIN_MIN,'BIN_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(BIN_SIZE,'BIN_SIZE',LU_IN,LU_ER,' ')
	  BIN_MIN=BIN_MIN*EV_TO_HZ	!Convert to units of 10^15 HZ
	  BIN_SIZE=BIN_SIZE*EV_TO_HZ
! 
! Temperature tabulated in equal increments of Log T.
!
	  CALL RD_INT(N_TEMP,'N_TEMP',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(LOG_T_MIN,'LOG_T_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(DEL_LOG_T,'DEL_LOG_T',LU_IN,LU_ER,' ')
	  LOG_T_MIN=LOG_T_MIN-4.0		!Convert from K to units of 10^4 K
! 
! Temperature tabulated in equal increments of Log Ne.
!
	  CALL RD_INT(N_ED,'N_ED',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(LOG_ED_MIN,'LOG_ED_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(DEL_LOG_ED,'DEL_LOG_ED',LU_IN,LU_ER,' ')
!
! Now that we have the vector sizes, we cab allocate memory for the X-ray 
! table and vectors.
!
	  ALLOCATE (X_NU(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (X_EMISS1(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (X_EMISS2(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LOG_X_TEMP(N_TEMP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LOG_X_ED(N_ED),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XRAY_FLUXES(N_BINS,N_TEMP,N_ED),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error allocating memory (1)'
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    WRITE(LU_ER,*)'STAT=',IOS
	    STOP
	  END IF
!
! Read in fluxes for each parameter set. Blank lines, comments,
! and a single header record are ignored.
!
	  DO J=1,N_ED
	    DO I=1,N_TEMP
	      STRING='!'
	      DO WHILE (STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!'
	1           .OR. INDEX(STRING,'Temperature(K)') .NE. 0)
	        READ(LU_IN,'(A)')STRING
		  IF(INDEX(STRING,'Temperature(K)') .NE. 0)WRITE(LU_OUT,'(A)')TRIM(STRING)
	       END DO
	       BACKSPACE(LU_IN)
	       READ(LU_IN,*)XRAY_FLUXES(:,I,J)
	       XRAY_FLUXES(:,I,J)=XRAY_FLUXES(:,I,J)*1000.0D0
	       WRITE(LU_OUT,'(10ES12.3)')XRAY_FLUXES(:,I,J)
	    END DO
	 END DO
!
	STOP
	END
