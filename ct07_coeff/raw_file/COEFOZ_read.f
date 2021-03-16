      PROGRAM READ_OZONE 
C
      PARAMETER (NLAT=64,NLEV=91,NCOEFF=9,NMONTH=12)
C
      PARAMETER (NIN=1)
C
      REAL PHI(NLAT), PRES(NLEV)
C
      REAL COEFF(NLAT,NLEV,NCOEFF,NMONTH)
C
      CHARACTER*5 TYEAR
      CHARACTER*6 TMONTH
C
C  READ_OZONE reads the 8 coefficients for the ozone photochemistry
C             linear parametrisation for each month 
C             (see O3_readme.txt file)
C
C    input :  NIN  logical unit of the file to be read (COEFOZ_v2.x.ascii)
C
C   output :  PHI  (NLAT)   array containing latitudes of the grid (units:degree, negative in the SH)
C             PRES (NLEV)   array containing pressures of the grid (units: hPa)
C             COEFF(NLATxNLEVxNCOEFFxNMONTH) array containing the coefficients 
C             on the referenced grid (see O3_readme.txt file for units)
C
C
C
      DO 1 IM=1,NMONTH
      READ (NIN,*)TYEAR, IYEAR
      READ (NIN,*)TMONTH, IMONTH
      PRINT *,TYEAR, IYEAR
      PRINT *,TMONTH, IMONTH
      DO 2 K=1,NLAT
      DO 2 L=1,NLEV
      READ (NIN,*) PHI(K),PRES(L),(COEFF(K,L,JC,IM),JC=1,NCOEFF)
 2    CONTINUE
C   
C     PRINTING FOR CONTROL
C
      PRINT *, PHI
      PRINT *, PRES
      DO 3 JC=1,NCOEFF
      PRINT *, COEFF (1,1,JC,IM), COEFF(NLAT,NLEV,JC,IM)
 3    CONTINUE
 1    CONTINUE
      STOP
      END
