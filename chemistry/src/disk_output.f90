MODULE disk_output
   
   USE global_variables
   
   PRIVATE
   PUBLIC :: write_output1, write_output2
   
   CONTAINS
   
   !########################################################
   SUBROUTINE write_output1(i,time)
   !########################################################
      
   IMPLICIT NONE
   INTEGER, INTENT(in) :: i
   LOGICAL :: isDefined
   REAL(KIND=dp), INTENT(in) :: time
   CHARACTER(len=80) :: filename_output
   
   INQUIRE(file='out', exist=isDefined)
   ! We create the folder 'out' if he doesn't exists.
   IF (.not.isDefined) then
      CALL system("mkdir out")
   END IF
   
   nz = maxval(nztmp)
   
   ! Write abundances
   WRITE(filename_output, '(a,i0.6,a)') 'out/abundances.',i,'.out'
   OPEN(UNIT=35, file=filename_output, form='unformatted')
   WRITE(35) time/year
   WRITE(35) abdisk(1:nspe,1:nr,1:nz)
   CLOSE(35)
   
   ! Write rates
   WRITE(filename_output, '(a,i0.6,a)') 'out/rates.',i,'.out'
   OPEN(45, file=filename_output, form='unformatted')
   WRITE(45) ratesdisk(1:nreac,1:nr,1:nz)
   CLOSE(45)
      
   END SUBROUTINE write_output1
   !########################################################
   
   !########################################################
   SUBROUTINE write_output2(i,time)
   !########################################################
      
   IMPLICIT NONE
   INTEGER, INTENT(in) :: i
   REAL(KIND=dp), INTENT(in) :: time
   CHARACTER(len=80) :: filename_output
   
   nz = maxval(nztmp)
   
   ! Write column densities
   WRITE(filename_output, '(a,i0.6,a)') 'out/columns.',i,'.out'
   OPEN(UNIT=35, file=filename_output, form='unformatted')
   WRITE(35) time/year
   WRITE(35) ColDensdisk(1:nr,1:nz,1:5)
   WRITE(35) ColDensupdisk(1:nr,1:nz,1:5)
   WRITE(35) ColDensdndisk(1:nr,1:nz,1:5)
   CLOSE(35)
         
   END SUBROUTINE write_output2
   !########################################################
   
END MODULE disk_output