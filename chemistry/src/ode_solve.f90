MODULE ode_solve
   
   USE global_variables
   USE chemistry
   
   PRIVATE
   PUBLIC :: initial_ode, set_work_arrays, get_jacobian, get_temporal_derivatives
   
   CONTAINS
   
   !########################################################
   SUBROUTINE initial_ode
   !########################################################   
   
   IMPLICIT NONE
   
   ! Count the number of non-zero element of the Jacobian matrix
   CALL count_nonzeros  
   
   ! Allocate work arrays
   CALL allocate_work_array
   
   ! Allocate satol, ian and jan
   ALLOCATE(satol(nspe))
   ALLOCATE(ian(nspe))
   ALLOCATE(jan(nspe))
    
   END SUBROUTINE initial_ode
   !########################################################   
   
   !########################################################
   SUBROUTINE allocate_work_array
   !########################################################

   IMPLICIT NONE

   lrw = 20 + 3 * nb_nonzeros_values*nspe + 21 * nspe
   liw = 31 + 3 * nb_nonzeros_values*nspe + 21 * nspe

   ALLOCATE(iwork(liw))
   ALLOCATE(rwork(lrw))

   iwork(1:liw) = 0
   rwork(1:lrw) = 0.d0

   RETURN
   END SUBROUTINE allocate_work_array
   !########################################################

   !########################################################
   SUBROUTINE count_nonzeros
   !########################################################
   
   IMPLICIT NONE
   
   ! Dummy parameters for restricted call of get_jacobian
   REAL(KIND=dp), DIMENSION(nspe) :: dummypdj, dummyabin
   INTEGER :: idummy
   INTEGER, PARAMETER :: dummy_n = 3
   REAL(KIND=dp), PARAMETER :: dummy_t = 0.0_dp
   REAL(KIND=dp), DIMENSION(dummy_n) :: dummy_ian, dummy_jan
   INTEGER :: max_nonzeros, numberjac
   
   INTEGER :: i
   
   ! Forced initialisation of global variables that will be needed, 
   ! especially for the 'set_constant_rates' part. We donc care about specific values,
   ! all that counts is that we can retrieve the number of non-zeros elements.
   
   max_nonzeros = 0
   dummyabin(1:nspe) = 1.e-5_dp
   
   CALL set_constant_rates
   CALL set_dependant_rates(dummyabin)
   
   DO idummy = 1,nspe
   
      CALL get_jacobian(dummy_n, dummy_t, dummyabin,idummy,dummy_ian, dummy_jan, dummypdj)
      numberjac = 0
      DO i = 1,nspe
         IF (dummypdj(i).gt.1.0e-99_dp) THEN
            numberjac = numberjac + 1
         ENDIF
      ENDDO
   
      IF (numberjac.GT.max_nonzeros) THEN
         max_nonzeros = numberjac
      ENDIF
   
   ENDDO
   
   nb_nonzeros_values = max_nonzeros
   
   RETURN
   END SUBROUTINE count_nonzeros
   !########################################################
   
   
   !########################################################
   SUBROUTINE set_work_arrays(abin)
   !########################################################
   
   IMPLICIT NONE
   
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe) :: abin
   
   INTEGER :: i,j,k
   REAL(KIND=dp), DIMENSION(nspe) :: PDJ
   INTEGER :: NNZ
   
   ! Dummy parameters for restricted call of get_jacobian
   INTEGER, PARAMETER :: dummy_n = 3
   REAL(KIND=dp), PARAMETER :: dummy_t = 0.d0
   REAL(KIND=dp), DIMENSION(dummy_n) :: dummy_ian, dummy_jan
   
   
   ! For IA and JA
   INTEGER, DIMENSION(nspe+1) :: IA
   INTEGER, DIMENSION(liw) :: JA
   
   CALL set_constant_rates
   CALL set_dependant_rates(abin)
   
   ! Initialize work arrays
   iwork(1:liw) = 0
   rwork(1:lrw) = 0.d0
   IWORK(5) = 5
   RWORK(6) = 3.154D14
   IWORK(6) = 10000
   IWORK(7) = 2
   
   IF (.NOT.(first_step_done)) THEN
      IWORK(6)=2000
   ENDIF
   
   k=1
   
   DO j = 1,nspe
   
      CALL get_jacobian(dummy_n, dummy_t, abin,j,dummy_ian, dummy_jan, pdj)
   
      ia(j)=k
   
      DO i = 1,nspe
         IF (abs(pdj(i)).gt.1.d-99) THEN
            ja(k) = i
            k=k+1
         ENDIF
      ENDDO
   
   ENDDO
   
   ia(nspe+1)=k
   
   nnz=ia(nspe+1)-1
   iwork(30+1:30+nspe+1)=IA(1:nspe+1)
   iwork(31+nspe+1:31+nspe+NNZ)=JA(1:NNZ)
   
   RETURN
   END SUBROUTINE set_work_arrays
   !########################################################
   
   !########################################################
   SUBROUTINE get_jacobian(n, t, abin, j, ian, jan, pdj)
   !########################################################
   
   IMPLICIT NONE
   
   ! Inputs
   INTEGER, INTENT(IN)                     :: n        ! number of first order ODEs.
   INTEGER, INTENT(IN)                     :: j        ! Index representing the J-th column of the jacobian
   REAL(KIND=dp), INTENT(IN)               :: t        ! the initial value of the independent variable t.
   REAL(KIND=dp), INTENT(IN), DIMENSION(n) :: ian      ! structure descriptor array of size N + 1.
   REAL(KIND=dp), INTENT(IN), DIMENSION(n) :: jan      ! structure descriptor array of size NNZ.
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe)  :: abin ! abundances
   
   ! Outputs
   REAL(KIND=dp), INTENT(OUT), DIMENSION(nspe) :: PDJ ! J-th column of df/dy
   
   ! Locals
   INTEGER                                :: no_spe             ! Index corresponding to no species
   REAL(KIND=dp), DIMENSION(nspe+1)       :: PDJ2
   INTEGER                                :: i, l, m
   INTEGER                                :: r1, r2, r3         ! Index for reactants
   INTEGER                                :: p1, p2, p3, p4, p5 ! Index for products
   INTEGER                                :: idx                ! The index of a given reaction
   
   ! Temp values to increase speed
   REAL(KIND=dp)                          :: square_density
   REAL(KIND=dp)                          :: tmp_value
   
   no_spe=nspe+1
   
   square_density = density * density
   
   PDJ2(1:nspe+1) = 0.0_dp
   
   DO i=1,nb_reactions_using_species(j)
   
      idx = relevant_reactions(i,j) ! j being the species index, given as a parameter
   
      r1 = reac(idx)%r(1)
      r2 = reac(idx)%r(2)
      r3 = reac(idx)%r(3)
   
      p1 = reac(idx)%p(1)
      p2 = reac(idx)%p(2)
      p3 = reac(idx)%p(3)
      p4 = reac(idx)%p(4)
      p5 = reac(idx)%p(5)
   
      ! One reactant only
      IF (r2.eq.no_spe) THEN
         IF (r1.eq.J) THEN
            tmp_value = reac(idx)%rate
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
      ENDIF
   
      ! Two bodies reaction
      ELSEIF (r3.eq.no_spe) THEN
         IF (r1.eq.J) THEN
            tmp_value = reac(idx)%rate * abin(r2) * density
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
            PDJ2(r2) = PDJ2(r2) - tmp_value
         ENDIF
   
         IF (r2.eq.J) THEN
            tmp_value = reac(idx)%rate * abin(r1) * density
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
            PDJ2(r2) = PDJ2(r2) - tmp_value
         ENDIF
   
      ! Three bodies reaction
      ELSE
         IF (r1.eq.J) THEN
            tmp_value = reac(idx)%rate * abin(r2) * abin(r3) * square_density
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
            PDJ2(r2) = PDJ2(r2) - tmp_value
            PDJ2(r3) = PDJ2(r3) - tmp_value
         ENDIF
   
         IF (r2.eq.J) THEN
            tmp_value = reac(idx)%rate * abin(r1) * abin(r3) * square_density
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
            PDJ2(r2) = PDJ2(r2) - tmp_value
            PDJ2(r3) = PDJ2(r3) - tmp_value
         ENDIF
   
         IF (r3.eq.J) THEN
            tmp_value = reac(idx)%rate * abin(r1) * abin(r2) * square_density
            PDJ2(p1) = PDJ2(p1) + tmp_value
            PDJ2(p2) = PDJ2(p2) + tmp_value
            PDJ2(p3) = PDJ2(p3) + tmp_value
            PDJ2(p4) = PDJ2(p4) + tmp_value
            PDJ2(p5) = PDJ2(p5) + tmp_value
            PDJ2(r1) = PDJ2(r1) - tmp_value
            PDJ2(r2) = PDJ2(r2) - tmp_value
            PDJ2(r3) = PDJ2(r3) - tmp_value
         ENDIF
   
      ENDIF
   ENDDO
      
   PDJ(1:nspe)=PDJ2(1:nspe)
   
   RETURN
   END SUBROUTINE get_jacobian
   !########################################################
   
   !########################################################
   SUBROUTINE get_temporal_derivatives(N,T,abin,abdot)
   !########################################################
   
   IMPLICIT NONE
   
   ! Inputs
   integer, intent(in) :: N !<[in] number of first order ODEs.
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe) :: abin !< [in] abundances
   REAL(KIND=dp), INTENT(IN) :: T !<[in] Not used by the code, but needed for the ODEPACK call format expected for FCHEM in dlsodes
   
   ! Outputs
   REAL(KIND=dp), intent(out), dimension(nspe) :: abdot
   ! Locals
   integer :: no_spe
   REAL(KIND=dp), dimension(nspe+1) :: YD2
   REAL(KIND=dp), dimension(nspe+1) :: YDTMP1,YDTMP2
   REAL(KIND=dp), DIMENSION(nspe)   :: abtmp !< [in] abundances
   INTEGER                                :: i
   INTEGER                                :: r1, r2, r3         ! Index for reactants
   INTEGER                                :: p1, p2, p3, p4, p5 ! Index for products
   REAL(KIND=dp)                          :: rate
   REAL(KIND=dp)                          :: square_density   
   
   no_spe = nspe + 1
   square_density = density * density
   abtmp = abin
   where(abtmp.le.1.0e-99_dp) abtmp = 1.0e-99_dp   
   call set_dependant_rates(abtmp)
   abdot(1:nspe) = 0.0_dp
   yd2(1:nspe) = 0.0_dp
   YDTMP1(1:nspe) = 0.0_dp
   YDTMP2(1:nspe) = 0.0_dp
   
   ! The differential equations are calculated in a loop here
   do i = 1,nreac
   
      r1 = reac(i)%r(1)
      r2 = reac(i)%r(2)
      r3 = reac(i)%r(3)
   
      p1 = reac(i)%p(1)
      p2 = reac(i)%p(2)
      p3 = reac(i)%p(3)
      p4 = reac(i)%p(4)
      p5 = reac(i)%p(5)
   
      ! One reactant only
      if (r2.eq.no_spe) then
         RATE = reac(i)%rate * abtmp(r1)
      else
         if (r3.eq.no_spe) then
      ! Two bodies reactions
            RATE = reac(i)%rate * abtmp(r1) * abtmp(r2) * density
         else
      ! Three bodies reactions
            RATE = reac(i)%rate*abtmp(r1) * abtmp(r2) * abtmp(r3) * square_density
         endif
      endif
   
      ! Used to compute the net rate of change in the total surface material
      IF ((reac(i)%itype.ne.40).and.(reac(i)%itype.ne.41)) THEN
         YD2(p1) = YD2(p1) + RATE
         YD2(p2) = YD2(p2) + RATE
         YD2(p3) = YD2(p3) + RATE
         YD2(p4) = YD2(p4) + RATE
         YD2(p5) = YD2(p5) + RATE
   
         YD2(r1) = YD2(r1) - RATE
         YD2(r2) = YD2(r2) - RATE
         YD2(r3) = YD2(r3) - RATE
   
         YDTMP1(p1) = YDTMP1(p1) + RATE
         YDTMP1(p2) = YDTMP1(p2) + RATE
         YDTMP1(p3) = YDTMP1(p3) + RATE
         YDTMP1(p4) = YDTMP1(p4) + RATE
         YDTMP1(p5) = YDTMP1(p5) + RATE
   
         YDTMP2(r1) = YDTMP2(r1) - RATE
         YDTMP2(r2) = YDTMP2(r2) - RATE
         YDTMP2(r3) = YDTMP2(r3) - RATE
   
      ENDIF
   
   ENDDO
      
   abdot(1:nspe) = YD2(1:nspe)
      
   RETURN
   
   END SUBROUTINE get_temporal_derivatives
   !########################################################

END MODULE ode_solve