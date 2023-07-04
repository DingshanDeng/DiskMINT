!########################################################
PROGRAM DISK_EXTRACT
!########################################################

   
   USE global_variables
   USE init_chem
   USE init_struc
   USE CHEMISTRY

   IMPLICIT NONE
   INTEGER :: i,j,k,l
   CHARACTER(len=3) :: usrans
   INTEGER :: output
   LOGICAL :: isDefined
   CHARACTER(len=80) :: filename_output
   CHARACTER(len=80) :: output_format
   REAL(KIND=dp) :: time
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: timetmp
   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: abtmp
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: abin
   REAL(KIND=dp) :: fuv11
   REAL(KIND=dp) :: fuv10
   REAL(KIND=dp) :: fuvLy
   REAL(KIND=dp) :: fuv07
   REAL(KIND=dp) :: fuv06
   REAL(KIND=dp) :: fuv04
   REAL(KIND=dp) :: fuv03
   REAL(KIND=dp) :: fuv02
   REAL(KIND=dp) :: fuv01
   
   REAL(KIND=dp), DIMENSION(2) :: tmp1, tmp2
   
   REAL(kind=dp) :: wvl, dum1, dum2, dum3
   REAL(KIND=dp), DIMENSION(:,:,:),  ALLOCATABLE :: H2_emiss
   REAL(KIND=dp), DIMENSION(:),  ALLOCATABLE :: H2_emiss_tot
   
   REAL(KIND=dp) :: G0star, G0ismup, G0ismdn
   REAL(kind=dp) :: ratio 
   INTEGER :: indkco = 841

   real(kind=dp) :: tetaco,tetacoup,tetacodn
   
      
   ! Initialise the chemistry 
   CALL initial_chem
   
   ! Initialise the structure
   CALL initial_struc

   WRITE(filename_output, '(a,i0.6,a)') 'out/abundances.000010.out'
   print*, filename_output
   OPEN(UNIT=35, file=filename_output, form='unformatted',status='old')
   READ(35) time
   nz = maxval(nztmp)
   READ(35) abdisk(1:nspe,1:nr,1:nz)
   CLOSE(35)
   
   isDefined = .true.
   INQUIRE(file='ab', exist=isDefined)
   
   ! We create the folder 'ab' if he doesn't exists.
   IF (.not.isDefined) THEN
     CALL system("mkdir out/ab")
   END IF
   
   DO i = 1, nspe
      WRITE(filename_output, '(a,a,a)') 'out/ab/', trim(spe(i)%name),'.txt'
      !print*, filename_output
      OPEN(10, file=filename_output)
      WRITE(10,'(a,ES10.2,a)') '# Time = ', time, '[yr]'
      WRITE(output_format, *) '(3ES16.5)'
      DO j=1,nr
         nz = nztmp(j)
         DO k = 1, nz
            WRITE(10,output_format) rdisk(j)/au, zdisk(j,k)/au,log10(abdisk(i,j,k))
         ENDDO
     ENDDO
     CLOSE(10)
   ENDDO      

END PROGRAM DISK_EXTRACT
!########################################################
