!########################################################
PROGRAM disk_main
!########################################################
   
   USE global_variables
   USE init_chem
   USE init_struc
   USE chemistry
   USE ode_solve
   USE disk_output
   
   IMPLICIT NONE
   ! --- For loops
   INTEGER :: i,j,k,l
   ! --- Abundances
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: abin
   ! --- Time
   REAL(KIND=dp) :: tbeg = 1.0_dp
   REAL(KIND=dp) :: tend = 1.0e6_dp
   REAL(KIND=dp) :: time, time_step, output_step,time_tmp,time_out
   INTEGER :: nbp = 10
   ! --- Others
   REAL(KIND=dp) :: cpu_start, cpu_end
   CHARACTER(len=11) :: tmp_name
   
   CHARACTER(len=1),DIMENSION(50) :: ar
   CHARACTER(len=1),DIMENSION(50) :: az
   
   INTEGER :: niter
   
   CHARACTER(len=80) :: filename_output
   
   !================================================================
   ! 1 - INITIALISATION
   !================================================================

   WRITE(6,*) "================================================================"
   WRITE(6,*) " 1 - INITIALIZATION"
   WRITE(6,*) "================================================================"

   CALL CPU_TIME(cpu_start)

   is_grain_reactions = 1
   is_growth = 1
   is_alloutput = 0
   
   ! Initialise the chemistry   
   CALL initial_chem 
   ALLOCATE(abin(nspe))
   abin(1:nspe) = spe(1:nspe)%ab      
               
   ! Initialise the structure
   CALL initial_struc   
   DO i = 1,nr
      nz = nztmp(i)
      DO j = 1,nz
         abdisk(1:nspe,i,j) = abin(1:nspe)
      ENDDO
   ENDDO
      
   ! Initialise variables needed by DLSODES
   CALL XSETF(0) ! No messages from DLSODES
   CALL initial_ode

   ! Set time and output time
   tbeg = tbeg*year
   tend = tend*year
   output_step = (tend/tbeg)**(1.0_dp/(nbp-1.0_dp))
   
   !================================================================
   ! 2 - SOLVE CHEMISTRY
   !================================================================
   WRITE(6,*) "================================================================"
   WRITE(6,*) " 2 - START MAIN CALCULATION                                     "
   WRITE(6,*) "================================================================"
   
   
   time = 0.0_dp

   ! Cycle over the number of outputs
   DO i = 1, nbp
      
      time_out = tbeg * output_step**(i-1.0_dp)
      time_step = time_out - time
      
      !WRITE(6,'(a,i3,a,ES13.6,a)') '--- Idx: ',i,'  Time = ',time_out/year,' yr'
      ar = ''
      
      !print*, time_out/year 
      !cycle
      
      ! Cycle mostly done to have the column densities right
      DO WHILE (time.lt.time_out)
                  
         ! Cycle over the radius
         DO j = 1,nr
            CALL system("clear && printf '\e[3J'")
            WRITE(6,*) "================================================================"
            WRITE(6,*) " ***           COMPUTING TIME-DEPENDENT CHEMISTRY           *** "
            WRITE(6,*) "================================================================"
            WRITE(6,*) " "
            WRITE(6,'(a,i3,a)') '  Iteration: ',i,' / 10'
            WRITE(6,'(a,ES10.3,a)') '  Total integration time = ',tend/year,' yr' 
            WRITE(6,'(a,ES10.3,a)') '  Current time           = ',time_out/year,' yr' 
            WRITE(6,*) " "
            l = 50*j/nr
            ar(l) = '#'
            WRITE(6,'(4x,a,50a1,a1,1i3,1a1,2x)') '0% |', ar(:),'|',100*j/nr,'%'
            WRITE(6,*) " "
            WRITE(6,*) "================================================================"
            nz = nztmp(j)
            ! Cycle over the altitude
            DO k = nz,1,-1
               abin(1:nspe) = abdisk(1:nspe,j,k)
               density = ndisk(j,k)
               Tgas = Tgdisk(j,k)
               Tdust = Tddisk(j,k)
               Av = Avdisk(j,k)
               Avup = Avupdisk(j,k)
               Avdn = Avdndisk(j,k)
               ndsigma = pi*in_ndusta2
               ndust = in_ndust * density
               ndusta = in_ndusta * density
               ndusta2 = in_ndusta2 * density
               if (i.eq.1) then
                   abin(indgr0) =ndust/density
                   abin(indgrm) = 1.0e-99_dp
                end if
               nsites = 4.0*ndsigma*surface_site_density
               G0Hab   = G0disk*(Rstar/au)**2/((rdisk(j)/au)**2 + (zdisk(j,k)/au)**2)
               G0upHab = 1.0_dp * scale_G0
               G0dnHab = 1.0_dp * scale_G0
                                             
               CALL get_taus(j,k)
               CALL get_disk_columns(j,k)             
               CALL GET_CR_RATE()
               if (CRrate.gt.2.0e-16_dp) CRrate = 2.0e-16_dp

               time_tmp = 0.0d0
               niter = 1
               DO WHILE(time_tmp.lt.time_step)        
                  istate = 1
                  DO l=1,nspe
                     satol(l) = max(atol, 1.0e-16_dp * abin(l))
                  ENDDO
                  CALL set_work_arrays(abin)
                  CALL DLSODES(get_temporal_derivatives,nspe,abin,time_tmp,time_step,itol,relative_tolerance,&
                    satol,itask,istate,iopt,rwork,lrw,iwork,liw,get_jacobian,mf)
                  where(abin.lt.1.0e-99_dp) abin=1.0e-99_dp
                  IF (istate.eq.-3) THEN
                     PRINT*, '**ERROR** Solver failed.. ISTATE = ', istate
                     STOP
                  END IF
                  if (istate.ne.2) niter = niter + 1
                  if(niter == 10) then
                     print*,'Number max of iteration reached... Using abundances from adjacent cell'
                     if (k < nz) then
                        abin(1:nspe) = abdisk(1:nspe,j,k+1)
                     else
                        abin(1:nspe) = abdisk(1:nspe,j,k-1)
                     end if
                     exit
                  end if
               END DO               
               CALL check_conservation(j,k,abin)
               abdisk(1:nspe,j,k) = abin(1:nspe)
               IF ((is_alloutput.eq.-1).or.(i.eq.nbp)) THEN
                  ratesdisk(1:nreac,j,k)=reac(1:nreac)%rate
               ENDIF
               IF(i.eq.nbp) THEN
                 CALL save_disk_columns(j,k)
               END IF
            END DO
         END DO
         time = time + time_step
      ENDDO
      first_step_done = .true.
      IF ((is_alloutput.eq.-1).or.(i.eq.nbp)) THEN
         CALL write_output1(i,time_out)
      ENDIF
      IF (i.eq.nbp) CALL write_output2(i,time_out)   
   ENDDO
   
   WRITE(6,*) "================================================================"
   WRITE(6,*) " END                                                            "
   WRITE(6,*) "================================================================"
   
   CALL CPU_TIME(cpu_end)
   
   PRINT*, 'CPU time = ', cpu_end-cpu_start
      
END PROGRAM disk_main
!########################################################
