MODULE init_struc
   
   USE global_variables
   
   PRIVATE
   PUBLIC :: initial_struc, get_disk_columns, get_taus, save_taus, save_disk_columns, get_abun

   CONTAINS
   
   !########################################################
   SUBROUTINE initial_struc
   !########################################################

   IMPLICIT NONE
   INTEGER :: i

   ! Read structure files
   CALL read_struc

   ! Put everything in cgs units
   rdisk = au * 10.0_dp**rdisk
   zdisk = au * 10.0_dp**zdisk
   ndisk = 10.0_dp**ndisk
   Tgdisk = 10.0_dp**Tgdisk
   where(Tgdisk.ge.8.0e3_dp) Tgdisk = 8.0e3_dp
   Tddisk = 10.0_dp**Tddisk
   where(Tddisk.ge.2.5e2_dp) Tddisk = 2.5e2_dp
   where(Avdisk.ge.250.0_dp) Avdisk = 250.0_dp
   
   END SUBROUTINE initial_struc
   !########################################################
   
   !########################################################
   SUBROUTINE read_struc
   !########################################################
   
   IMPLICIT NONE
   INTEGER            :: i, j,k
   INTEGER            :: ntmp
   CHARACTER(len=80)  :: line
   CHARACTER(LEN=80)  :: input_file
   REAL(kind=dp)      :: dummy, tol
   
   REAL(kind=dp),DIMENSION(:), ALLOCATABLE :: rtmp
   
   input_file = "data/struc/diskstructure"
   OPEN(10,file=input_file,status="old")
   READ(10,*) nr
   ALLOCATE(nztmp(nr))
   do i = 1, nr
      READ(10,*) ntmp
      nztmp(i) = ntmp
      DO j = 1, ntmp
         READ(10,*)
      ENDDO
   ENDDO
   CLOSE(10)

   nz = maxval(nztmp)
   
   ALLOCATE(rdisk(nr))
   ALLOCATE(zdisk(nr,nz))
   ALLOCATE(ndisk(nr,nz))
   ALLOCATE(Avdisk(nr,nz))
   ALLOCATE(Avupdisk(nr,nz))
   ALLOCATE(Avdndisk(nr,nz))
   ALLOCATE(Tgdisk(nr,nz))
   ALLOCATE(Tddisk(nr,nz))
   ALLOCATE(ColDensdisk(nr,nz,5))
   ALLOCATE(ColDensupdisk(nr,nz,5))
   ALLOCATE(ColDensdndisk(nr,nz,5))
   !ALLOCATE(Xray_table(nr,nz,9))
   !ALLOCATE(pah_rate_table(nr,nz,5))
   !ALLOCATE(in_ndusta2(nr,nz))
   !ALLOCATE(in_ndusta(nr,nz))
   !ALLOCATE(in_ndust(nr,nz))
   ALLOCATE(abdisk(nspe,nr,nz))
   ALLOCATE(ratesdisk(nreac,nr,nz))
   
   ALLOCATE(rtmp(nr))
   
   rdisk = 1.0e-99_dp
   zdisk = 1.0e-99_dp
   ndisk = 1.0e-99_dp
   G0disk = 1.0e-99_dp
   Avdisk = 1.0e-99_dp
   Avupdisk = 1.0e-99_dp
   Avdndisk = 1.0e-99_dp
   Tgdisk = 1.0e-99_dp
   Tddisk = 1.0e-99_dp
   ColDensdisk = 1.0e-99_dp
   ColDensupdisk = 1.0e-99_dp
   ColDensdndisk = 1.0e-99_dp
   in_ndust = 1.0e-99_dp
   in_ndusta = 1.0e-99_dp
   in_ndusta2 = 1.0e-99_dp
   abdisk = 1.0e-99_dp
   ratesdisk = 1.0e-99_dp
   !Xray_table = 1.0e-99_dp
   !pah_rate_table = 1.0e-99_dp

   input_file = "data/struc/diskparameters"
   OPEN(10,file=input_file,status="old")
   READ(10,*) Rstar
   READ(10,*) G0disk
   READ(10,*) in_ndust
   READ(10,*) in_ndusta
   READ(10,*) in_ndusta2
   CLOSE(10)

   input_file = "data/struc/diskstructure"
   OPEN(10,file=input_file,status="old")
   READ(10,*) nr
   do i = 1, nr
      READ(10,*)
      nz = nztmp(i)
      DO j = 1, nz
         READ(10,*) rdisk(i), zdisk(i,j), ndisk(i,j), Tgdisk(i,j),&
                    Tddisk(i,j),Avdisk(i,j), Avupdisk(i,j), Avdndisk(i,j)
      ENDDO
   ENDDO
   CLOSE(10)
         
  END SUBROUTINE read_struc
   !########################################################
   
   !########################################################
   SUBROUTINE get_disk_columns(i,j)
   !########################################################

   USE global_variables

   IMPLICIT NONE
   ! Inputs
   INTEGER, INTENT(IN) :: i,j
   
   ! Outputs
   !REAL(KIND=dp), INTENT(OUT) :: col, colup,coldown
   
   ! Locals
   INTEGER :: k,l,index
   REAL(KIND=dp) :: dr, drk, zray
   REAL(KIND=dp) :: alph, alphtmp, dummy
   
   ! Compute the colum density from the star
   dr = rdisk(i+1) - rdisk(i)
   Ntot = 0.5 * dr * ndisk(i,j)
   NH   = 0.5 * dr * ndisk(i,j) * abdisk(indH,i,j)
   NH2  = 0.5 * dr * ndisk(i,j) * abdisk(indH2,i,j)
   NH2v = 0.5 * dr * ndisk(i,j) * abdisk(indH2v,i,j)
   NCO  = 0.5 * dr * ndisk(i,j) * abdisk(indCO,i,j)
   if (ind13CO>0) N13CO  = 0.5 * dr * ndisk(i,j) * abdisk(ind13CO,i,j)
   if (indN2>0) NN2  = 0.5 * dr * ndisk(i,j) * abdisk(indN2,i,j)
   IF(i.eq.nr) THEN 
      Ntot =1.0e-99_dp  
      NH   =1.0e-99_dp
      NH2  =1.0e-99_dp
      NH2v =1.0e-99_dp
      NCO  =1.0e-99_dp
      N13CO = 1.0e-99_dp
      NN2  =1.0e-99_dp
   ELSE
   alph = atan(zdisk(i,j)/rdisk(i))        ! Angle from the center of the disk to the (i,j) point
   DO k = 1,i-1                            ! Integrate from the center to the desired (i,j) point
      drk = rdisk(k+1) - rdisk(k)
      nz = nztmp(k)
      zray = rdisk(k)*tan(alph)
      index = 0
      IF (zray.le.zdisk(k,1)) THEN
         index = 1
      ELSE
         DO l = 1, nz-1
            IF(zdisk(k,l).le.zray.and.zdisk(k,l+1).gt.zray) THEN
               index = l+1
               exit
            END IF
         END DO
      END IF
      IF (index.gt.0) THEN
         Ntot = Ntot + drk * ndisk(k,index)
         NH   = NH   + drk * ndisk(k,index) * abdisk(indH,k,index)
         NH2  = NH2  + drk * ndisk(k,index) * abdisk(indH2,k,index)
         NH2v = NH2v + drk * ndisk(k,index) * abdisk(indH2v,k,index)
         NCO  = NCO  + drk * ndisk(k,index) * abdisk(indCO,k,index)
         if (ind13CO>0) N13CO  = N13CO  + drk * ndisk(k,index) * abdisk(ind13CO,k,index)
         if (indN2>0) NN2  = NN2  + drk * ndisk(k,index) * abdisk(indN2,k,index)
      END IF
   END DO
   END IF
   
   ! Compute the colum density up
   NH_up   =1.0e-99_dp
   NH2_up  =1.0e-99_dp
   NH2v_up =1.0e-99_dp
   NCO_up  =1.0e-99_dp
   N13CO_up = 1.0e-99_dp
   NN2_up  =1.0e-99_dp
   nz = nztmp(i)
   DO k = nz-1,j,-1
      NH_up   = NH_up   + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH,i,k+1)
      NH2_up  = NH2_up  + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH2,i,k+1)
      NH2v_up = NH2v_up + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH2v,i,k+1)
      NCO_up  = NCO_up  + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indCO,i,k+1)
      if (ind13CO>0) N13CO_up  = N13CO_up  + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(ind13CO,i,k+1)
      if (indN2>0) NN2_up  = NN2_up  + (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indN2,i,k+1)
   ENDDO
   
   ! Compute the colum density down
   NH_dn   = NH_up 
   NH2_dn  = NH2_up
   NH2v_dn = NH2v_up
   NCO_dn  = NCO_up
   N13CO_dn  = N13CO_up
   NN2_dn  = NN2_up
   DO k = 1,j-1
      NH_dn   = NH_dn   + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH,i,k+1)
      NH2_dn  = NH2_dn  + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH2,i,k+1)
      NH2v_dn = NH2v_dn + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indH2v,i,k+1)
      NCO_dn  = NCO_dn  + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indCO,i,k+1)
      if (ind13CO>0) N13CO_dn  = N13CO_dn  + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(ind13CO,i,k+1)
      if (indN2>0) NN2_dn  = NN2_dn  + 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)*abdisk(indN2,i,k+1)
   END DO

   END SUBROUTINE get_disk_columns
   !########################################################
   
   !########################################################
   SUBROUTINE get_taus(i,j)
   !########################################################
   
   USE global_variables
   
   IMPLICIT NONE
   ! Inputs
   INTEGER, INTENT(IN) :: i,j
   
   ! Locals
   INTEGER :: k,l,index, idx
   REAL(KIND=dp) :: dr, drk, zray
   REAL(KIND=dp) :: alph, alphtmp, dummy
   REAL(KIND=dp) :: col
   
   ! Compute the colum density from the star
   dr = rdisk(i+1) - rdisk(i)
   col   = 0.5 * dr * ndisk(i,j)
   tauh2  = 1.2e-14_dp*abdisk( indh2,i,j)*col/(0.13_dp*sqrt(Tgdisk(i,j)/2.02_dp))
   tauh2v = 1.0e-15_dp*abdisk(indh2v,i,j)*col/(0.13_dp*sqrt(Tgdisk(i,j)/2.02_dp))
   IF(i.eq.nr) THEN
      tauh2 =1.0e-99_dp
      tauh2v =1.0e-99_dp
   ELSE
   alph = atan(zdisk(i,j)/rdisk(i))        ! Angle from the center of the disk to the (i,j) point
   DO k = 1,i-1                            ! Integrate from the center to the desired (i,j) point
      drk = rdisk(k+1) - rdisk(k)
      nz = nztmp(k)
      zray = rdisk(k)*tan(alph)
      index = 0
      IF (zray.le.zdisk(k,1)) THEN
         index = 1
      ELSE
         DO l = 1, nz-1
            IF(zdisk(k,l).le.zray.and.zdisk(k,l+1).gt.zray) THEN
               index = l+1
               exit
            END IF
         END DO
      END IF
      IF (index.gt.0) THEN
         col = drk * ndisk(k,index)
         tauh2  = tauh2  + 1.2e-14_dp*abdisk( indh2,k,index)*col/(0.13_dp*sqrt(Tgdisk(k,index)/2.02_dp))
         tauh2v = tauh2v + 1.0e-15_dp*abdisk(indh2v,k,index)*col/(0.13_dp*sqrt(Tgdisk(k,index)/2.02_dp))
      ENDIF
   END DO
   END IF
   
   ! Compute the colum density up
   tauuph2  =1.0e-99_dp
   tauuph2v =1.0e-99_dp
   nz = nztmp(i)
   DO k = nz-1,j,-1
      col   = (zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)
      tauuph2  = tauuph2  + 1.2e-14_dp*abdisk( indh2,i,k+1)*col/(0.13_dp*sqrt(Tgdisk(i,k+1)/2.02_dp))
      tauuph2v = tauuph2v + 1.0e-15_dp*abdisk(indh2v,i,k+1)*col/(0.13_dp*sqrt(Tgdisk(i,k+1)/2.02_dp))
   ENDDO
   
   ! Compute the colum density down
   taudnh2 = tauuph2
   taudnh2v = tauuph2v
   DO k = 1,j-1
      col = 2.0_dp*(zdisk(i,k+1)-zdisk(i,k))*ndisk(i,k+1)
      taudnh2  = taudnh2  + 1.2e-14_dp*abdisk( indh2,i,k+1)*col/(0.13_dp*sqrt(Tgdisk(i,k+1)/2.02_dp))
      taudnh2v = taudnh2v + 1.0e-15_dp*abdisk(indh2v,i,k+1)*col/(0.13_dp*sqrt(Tgdisk(i,k+1)/2.02_dp))
   END DO
      
   END SUBROUTINE get_taus
   !########################################################
   
   !########################################################
   SUBROUTINE get_abun(i,j)
   !########################################################

   USE global_variables

   IMPLICIT NONE
   ! Inputs
   INTEGER, INTENT(IN) :: i,j
   
   ! Locals
   INTEGER :: k,l,index
   REAL(KIND=dp) :: dr, drk, zray
   REAL(KIND=dp) :: alph, alphtmp, dummy
   
   if (i.eq.1) return
   
   alph = atan(zdisk(i,j)/rdisk(i))        ! Angle from the center of the disk to the (i,j) point
   k = i-1
   nz = nztmp(k)
   zray = rdisk(k)*tan(alph)
   index = 0
   IF (zray.le.zdisk(k,1)) THEN
      index = 1
   ELSE
      DO l = 1, nz-1
         IF(zdisk(k,l).le.zray.and.zdisk(k,l+1).gt.zray) THEN
            index = l+1
            exit
         END IF
      END DO
   END IF
   IF (index.gt.0) THEN
      abdisk(1:nspe,i,j) = abdisk(1:nspe,k,index)
   END IF
   
   END SUBROUTINE get_abun
   !########################################################

   !########################################################
   SUBROUTINE save_disk_columns(i,j)
   !########################################################
   
   USE global_variables
   
   IMPLICIT NONE
   ! Inputs
   INTEGER, INTENT(IN) :: i,j
   
   ColDensdisk(i,j,1) = NH 
   ColDensdisk(i,j,2) = NH2
   ColDensdisk(i,j,3) = NH2v
   ColDensdisk(i,j,4) = NCO 
   ColDensdisk(i,j,5) = NN2
   
   ColDensupdisk(i,j,1) = NH_up
   ColDensupdisk(i,j,2) = NH2_up
   ColDensupdisk(i,j,3) = NH2v_up
   ColDensupdisk(i,j,4) = NCO_up
   ColDensupdisk(i,j,5) = NN2_up
   
   ColDensdndisk(i,j,1) = NH_dn
   ColDensdndisk(i,j,2) = NH2_dn
   ColDensdndisk(i,j,3) = NH2v_dn
   ColDensdndisk(i,j,4) = NCO_dn
   ColDensdndisk(i,j,5) = NN2_dn
   
   END SUBROUTINE save_disk_columns
   !########################################################

END MODULE init_struc
