MODULE init_chem
   
   USE global_variables
   USE shielding
   
   PRIVATE
   PUBLIC :: initial_chem, check_conservation, get_elem_ab, get_major_bearing_species

   CONTAINS

   !########################################################
   SUBROUTINE initial_chem
   !########################################################

   IMPLICIT NONE
   INTEGER :: i
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tmp_ab
   INTEGER :: dum1, dum2

   ! Read chemistry files
   CALL read_input
   
   ! Index datas
   CALL index_data
   
   ! Initialize the reaction rates
   CALL init_reaction_rates
   CALL init_relevant_reactions

   ! Get the elemental abundances
   CALL get_elem_ab(spe%ab,elem%ab)
   
   ! Recompute initial_dtg_mass_ratio to remove He
   ! In the following, initial_dtg_mass_ratio is used as a H/dust mass ratio
   DO i=1,nelem
      IF (elem(i)%name.EQ.spehe) THEN
         dtg_mass_ratio = dtg_mass_ratio * (1.0_dp + 4.0_dp * elem(i)%ab)
      ENDIF
   ENDDO

   ! **** This part may need to be recomputed based on the disk structure
   gtodn = (4.0_dp * pi * grain_density * grain_radius**3.0_dp) / (3.0_dp * dtg_mass_ratio * amu)
   spe(indgr0)%ab = 1.0_dp / gtodn
   ndsigma = pi*grain_radius**2/gtodn
   nsites = 4.0*ndsigma*surface_site_density
   
   ! ****

   ! Compute the e- abundance by checking conservation
   ALLOCATE(tmp_ab(nspe))
   tmp_ab = spe%ab
   CALL check_conservation(1,1,tmp_ab)
   spe%ab = tmp_ab
   DEALLOCATE(tmp_ab)
   
   ! Check the network
   CALL CheckNetwork
   
   ! Read shielding factors
   CALL read_visser()
   
   dum1 = 0
   dum2 = 0
   DO i = 1, nreac
      ! Constant rates
      if (reac(i)%itype == 0 .or. &
          reac(i)%itype == 4 .or. &
          reac(i)%itype == 5 .or. &
          reac(i)%itype == 6 .or. &
          reac(i)%itype == 7 .or. &
          reac(i)%itype == 8 .or. &
          reac(i)%itype == 9 .or. &
          reac(i)%itype == 11 .or. &
          reac(i)%itype == 12 .or. &
          reac(i)%itype == 13 .or. &
          reac(i)%itype == 78) then
          dum1 = dum1 + 1
       ! Variables rates
       else
          dum2 = dum2 + 1
       end if
   END DO
   ALLOCATE(idx_const_rates(dum1))
   ALLOCATE(idx_var_rates(dum2))
   dum1 = 0
   dum2 = 0
   DO i = 1, nreac
      ! Constant rates
      if (reac(i)%itype == 0 .or. &
          reac(i)%itype == 4 .or. &
          reac(i)%itype == 5 .or. &
          reac(i)%itype == 6 .or. &
          reac(i)%itype == 7 .or. &
          reac(i)%itype == 8 .or. &
          reac(i)%itype == 9 .or. &
          reac(i)%itype == 11 .or. &
          reac(i)%itype == 12 .or. &
          reac(i)%itype == 13 .or. &
          reac(i)%itype == 78) then
          dum1 = dum1 + 1
          idx_const_rates(dum1) = i
       ! Variables rates
       else
          dum2 = dum2 + 1
          idx_var_rates(dum2) = i
       end if
   END DO

   END SUBROUTINE initial_chem
   !########################################################
   
   !########################################################
   SUBROUTINE READ_INPUT
   !########################################################

   USE auxilary

   IMPLICIT NONE

   INTEGER            :: i,j,k,l
   INTEGER            :: tmp
   INTEGER            :: ntmp
   INTEGER            :: r1, r2, r3, p1, p2, p3, p4, p5
   INTEGER            :: error
   
   INTEGER            :: nline
   CHARACTER(LEN=300) :: input_file
   CHARACTER(len=300) :: dummy
   CHARACTER(len=300) :: fmt
   CHARACTER(len=300), DIMENSION(:), ALLOCATABLE   :: line

   CHARACTER(LEN=11), DIMENSION(:,:), ALLOCATABLE  :: compounds
   CHARACTER(LEN=11), DIMENSION(:,:), ALLOCATABLE  :: tmp_compounds
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_alpha
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_beta
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_gamma
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_tmp
   INTEGER, DIMENSION(:), ALLOCATABLE              :: tmp_itype
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_tmin
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: tmp_tmax
   INTEGER, DIMENSION(:), ALLOCATABLE              :: tmp_formula
   INTEGER, DIMENSION(:), ALLOCATABLE              :: tmp_id

   CHARACTER(LEN=11), DIMENSION(:), ALLOCATABLE    :: tmp_name
   INTEGER, DIMENSION(:), ALLOCATABLE              :: tmp_int
   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE      :: tmp_rea

   !--------------------------------------------------------
   ! Read chemical element
   !--------------------------------------------------------

   input_file = "data/chem/element.in"
   CALL GetNline(input_file, nline)
   IF (nelem .NE. nline) THEN
     PRINT*, '** ERROR ** The actual number of elements in element.in is different than nelem'
     STOP 
   END IF
   ALLOCATE(elem(nelem))
   ALLOCATE(line(nelem))
   CALL ReadFile(input_file, nelem, line)
   DO i = 1, nelem
      dummy = line(i)
      READ(dummy,'(a, f8.3)') elem(i)%name, elem(i)%mass
   ENDDO
   DEALLOCATE(line)

   !--------------------------------------------------------
   ! Read chemical species (gas and grains)
   !--------------------------------------------------------

   input_file = "data/chem/gas_species.in"
   CALL GetNline(input_file, nspegas)
      
   input_file = "data/chem/grain_species.in"
   CALL GetNline(input_file, nspegr)

   nspe = nspegas + nspegr
   
   ! Allocation of species characteristics and initialisation
   ALLOCATE(spe(nspe+1))
   DO i = 1, nspe+1
      spe(i)%name = ""
      spe(i)%ab = 1.0e-99_dp
      spe(i)%index = 0
      spe(i)%mass = 0.0_dp
      spe(i)%compo(1:nelem) = 0
      spe(i)%form_enthalpy = 0.0_dp
      spe(i)%bind = 0.0_dp
      spe(i)%diff = 0.0_dp
      spe(i)%charge = 0
      spe(i)%stick = 0.0_dp
   ENDDO
   
   input_file = "data/chem/gas_species.in"
   ALLOCATE(line(nspegas))
   CALL ReadFile(input_file, nspegas, line)
   write(fmt, '(a,i0,a)') '(A11,i3,', nelem, '(I3),2x)'
   DO i = 1, nspegas
      dummy = line(i)
      READ(dummy,fmt) spe(i)%name, spe(i)%charge, spe(i)%compo
   ENDDO
   DEALLOCATE(line)
      
   input_file = "data/chem/grain_species.in"
   ALLOCATE(line(nspegr))
   CALL ReadFile(input_file, nspegr, line)
   write(fmt, '(a,i0,a)') '(A11,i3,', nelem, '(I3),2x)'
   j = 1 + nspegas
   DO i = 1, nspegr
      dummy = line(i)
      READ(dummy,fmt) spe(j)%name, spe(j)%charge, spe(j)%compo
      j = j + 1
   ENDDO
   DEALLOCATE(line)
   
   
   !--------------------------------------------------------
   ! Read initial abundances
   !--------------------------------------------------------
   
   input_file = "data/chem/abundances.in"
   CALL GetNline(input_file, nline)
   ALLOCATE(line(nline))
   ALLOCATE(tmp_name(nline))
   ALLOCATE(tmp_rea(nline,1))
   CALL ReadFile(input_file, nline, line)
   DO i = 1, nline
      dummy = line(i)
      READ(dummy,'(A11,3x,E12.6)') tmp_name(i), tmp_rea(i,1)
   ENDDO
   DO i = 1, nspe
      DO j = 1, nline
         IF(spe(i)%name.EQ.tmp_name(j)) spe(i)%ab = tmp_rea(j,1)
      ENDDO
   ENDDO
   DEALLOCATE(line)
   DEALLOCATE(tmp_name)
   DEALLOCATE(tmp_rea)
   

   !--------------------------------------------------------
   ! Read surface parameters
   !--------------------------------------------------------

   input_file = "data/chem/surface_parameters.in"
   CALL GetNline(input_file, nline)
   ALLOCATE(line(nline))
   ALLOCATE(tmp_name(nline))
   ALLOCATE(tmp_int(nline))
   ALLOCATE(tmp_rea(nline,4))
   CALL ReadFile(input_file, nline, line)
   DO i = 1, nline
      dummy = line(i)
      READ(dummy,'(A11,I4,F7.0,F6.0,D8.1,27X,F8.2)') tmp_name(i),tmp_int(i),(tmp_rea(i,j),j=1,4)
   ENDDO
   DO i = 1, nspe
      DO j = 1, nline
         IF(spe(i)%name.eq.tmp_name(j)) THEN
            spe(i)%bind = tmp_rea(j,1)
            spe(i)%diff = tmp_rea(j,2)
            spe(i)%form_enthalpy = tmp_rea(j,4)
         ENDIF
      ENDDO
      if ((spe(i)%name(1:1).eq.'J'.or.spe(i)%name(1:1).eq.'K').and.&
           spe(i)%bind <= 0.0_dp) then
           PRINT*, '** WARNING ** ', trim(spe(i)%name), ' does not have a binding energy'
      end if
      IF((spe(i)%name.NE.spejh).and.(spe(i)%name(1:1).EQ.'J')) THEN
         spe(i)%diff = diff_binding_ratio_surf * spe(i)%bind
      ENDIF
      IF((spe(i)%name.NE.spekh).and.(spe(i)%name(1:1).EQ.'K')) THEN
         spe(i)%diff = diff_binding_ratio_mant * spe(i)%bind
      ENDIF
      !print*, spe(i)%name, spe(i)%bind, spe(i)%diff
   ENDDO
   DEALLOCATE(line)
   DEALLOCATE(tmp_name)
   DEALLOCATE(tmp_int)
   DEALLOCATE(tmp_rea)

   ! We redetermine nspegas and nspegr from the label of each species
   nspegas = 0
   nspegr = 0
   DO i = 1,nspe
      IF((spe(i)%name(1:1).EQ.'J').OR.(spe(i)%name(1:1).EQ.'K')) THEN
         nspegr = nspegr + 1
      ELSE
         nspegas = nspegas + 1
      ENDIF
   ENDDO
   
   !--------------------------------------------------------
   ! Read chemical reactions (gas and grains)
   !--------------------------------------------------------
   
   input_file = "data/chem/gas_reactions.in"
   CALL GetNline(input_file, nreacgas)
   
   input_file = "data/chem/grain_reactions.in"
   CALL GetNline(input_file, nreacgr)
   
   nreac = nreacgas + nreacgr
   
   ! Allocation of reactions characteristics and initialisation
   ALLOCATE(reac(nreac))
   DO i = 1,nreac
      reac(i)%comp(1:ncomp) = '           '
      reac(i)%r(1) = nspe + 1
      reac(i)%r(2) = nspe + 1
      reac(i)%r(3) = nspe + 1
      reac(i)%p(1) = nspe + 1
      reac(i)%p(2) = nspe + 1
      reac(i)%p(3) = nspe + 1
      reac(i)%p(4) = nspe + 1
      reac(i)%p(5) = nspe + 1
      reac(i)%alpha = 0.0_dp
      reac(i)%beta = 0.0_dp
      reac(i)%gamma = 0.0_dp
      reac(i)%itype = 0
      reac(i)%tmin = 0.0_dp
      reac(i)%tmax = 0.0_dp
      reac(i)%formula = 0
      reac(i)%id = 0
      reac(i)%rate = 0.0_dp
      reac(i)%br = 0.0_dp
      reac(i)%br1 = 0.0_dp
      reac(i)%br2 = 0.0_dp
   ENDDO
   
   ALLOCATE(compounds(ncomp,nreac))
   ALLOCATE(tmp_compounds(ncomp,nreac))
   ALLOCATE(tmp_alpha(nreac))
   ALLOCATE(tmp_beta(nreac))
   ALLOCATE(tmp_gamma(nreac))
   ALLOCATE(tmp_tmp(nreac))
   tmp_tmp = 0.0
   ALLOCATE(tmp_itype(nreac))
   ALLOCATE(tmp_tmin(nreac))
   ALLOCATE(tmp_tmax(nreac))
   ALLOCATE(tmp_formula(nreac))
   ALLOCATE(tmp_id(nreac))
   
   
   input_file = "data/chem/gas_reactions.in"
   ALLOCATE(line(nreacgas))
   CALL ReadFile(input_file, nreacgas, line)
   DO i = 1, nreacgas
      dummy = line(i)
      READ(dummy,'(3A11,1x,5A11,3D11.3,D9.3,14x,I3,2f7.0,i3,i6)') (tmp_compounds(j,i), j=1,ncomp), tmp_alpha(i), tmp_beta(i), &
                  tmp_gamma(i),tmp_tmp(i), tmp_itype(i), tmp_tmin(i), tmp_tmax(i), tmp_formula(i), tmp_id(i)             
   ENDDO
   !pause
   DEALLOCATE(line)
   
   input_file = "data/chem/grain_reactions.in"
   ALLOCATE(line(nreacgr))
   CALL ReadFile(input_file, nreacgr, line)
   j = nreacgas + 1
   DO i = 1, nreacgr
      dummy = line(i)
      READ(dummy,'(3A11,1x,5A11,3D11.3,D9.3,14x,I3,2f7.0,i3,i6)') (tmp_compounds(k,j), k=1,ncomp), tmp_alpha(j), tmp_beta(j), &
                  tmp_gamma(j),tmp_tmp(j), tmp_itype(j), tmp_tmin(j), tmp_tmax(j), tmp_formula(j), tmp_id(j)
      j = j + 1
   ENDDO
   DEALLOCATE(line)
   
   ! Reorder reaction file entries with the type
   k=1
   DO i = 0,99
      DO j = 1,nreac
      IF (tmp_itype(j).EQ.i) THEN
         compounds(:,k) = tmp_compounds(:,j)
         reac(k)%alpha = tmp_alpha(j)
         reac(k)%beta = tmp_beta(j)
         reac(k)%gamma = tmp_gamma(j)
         reac(k)%tmp = tmp_tmp(j)
         reac(k)%itype = tmp_itype(j)
         reac(k)%tmin = tmp_tmin(j)
         reac(k)%tmax = tmp_tmax(j)
         reac(k)%formula = tmp_formula(j)
         reac(k)%id = tmp_id(j)
         k = k+1
      ENDIF
      ENDDO
   ENDDO
   
   ! Replace the species names by blanks for non chemical species
   DO j=1,nreac
      DO i=1,ncomp
         SELECT CASE(compounds(i,j))
            CASE ('CR', 'CRP', 'Photon')
            compounds(i,j) = '           '
         END SELECT
      ENDDO
   ENDDO
   !stop
   
   ! We assign the reacants and products id based on the name of each compound
   !  - By default, non existing reactants (dummy species) 
   !    will be assigned (nspe+1)
   
   DO i = 1,nreac
      DO j = 1,nspe
         IF (compounds(1,i).EQ.spe(j)%name) reac(i)%r(1) = j
         IF (compounds(2,i).EQ.spe(j)%name) reac(i)%r(2) = j
         IF (compounds(3,i).EQ.spe(j)%name) reac(i)%r(3) = j
   
         IF (compounds(4,i).EQ.spe(j)%name) reac(i)%p(1) = j
         IF (compounds(5,i).EQ.spe(j)%name) reac(i)%p(2) = j
         IF (compounds(6,i).EQ.spe(j)%name) reac(i)%p(3) = j
         IF (compounds(7,i).EQ.spe(j)%name) reac(i)%p(4) = j
         IF (compounds(8,i).EQ.spe(j)%name) reac(i)%p(5) = j
      ENDDO
      reac(i)%comp(1:ncomp) = compounds(:,i)
   ENDDO
      
   DEALLOCATE(compounds)
   DEALLOCATE(tmp_compounds)
   DEALLOCATE(tmp_alpha)
   DEALLOCATE(tmp_beta)
   DEALLOCATE(tmp_gamma)
   DEALLOCATE(tmp_itype)
   DEALLOCATE(tmp_tmin)
   DEALLOCATE(tmp_tmax)
   DEALLOCATE(tmp_formula)
   DEALLOCATE(tmp_id)
   
   !--------------------------------------------------------
   ! Shielding functions
   !--------------------------------------------------------
   
   ! Shielding function CO and N2
   input_file = "data/chem/self_shielding_12CO.dat"
   OPEN(10,file=input_file,status="old")
   READ(10,*)
   DO i = 1, 251
      DO j = 1, 226
         READ(10,'(3es10.3)') NH2CO_file(i,j), NCO_file(i,j), thetaCO_file(i,j)
      END DO
   END DO
   CLOSE(10)
   
   input_file = "data/chem/self_shielding_N2.dat"
   OPEN(10,file=input_file,status="old")
   READ(10,*)
   DO i = 1, 226
      DO j = 1, 226
         READ(10,'(3es10.3)') NH2N2_file(i,j), NN2_file(i,j), thetaN2_file(i,j)
      END DO
   END DO
   CLOSE(10)
   
   RETURN

   END SUBROUTINE READ_INPUT
   !########################################################
   
   !########################################################
   SUBROUTINE INDEX_DATA
   !########################################################

   IMPLICIT NONE

   REAL(KIND=dp)     :: msum
   INTEGER           :: ilab, j, k, i, isptemp
   INTEGER           :: tmp
   INTEGER           :: ksum ! sum of number of primary element composing the species. If equal to 1, the current species is elemental
   REAL(KIND=dp)     :: mass_tmp !  temporary value to exchange two index in the mass array
   CHARACTER(len=11) :: name_tmp !  temporary value to exchange two index in the name array

   ! Set elements' characteristics
   ! --- Find the atomic species associated with a given element

   ilab = 1
   DO j = 1,nspe
      ksum = 0
      ! ------ Calculate species elemental_mass
      DO k = 1,nelem
         ksum = ksum + spe(j)%compo(k)
      ENDDO

      ! ------ Check for atomic species
      IF ((ksum.EQ.1).AND.(spe(j)%charge.EQ.0).AND.&
          (spe(j)%name(1:1).NE.'J').AND.(spe(j)%name(1:1).NE.'K').AND.&
          (spe(j)%name(1:1).NE.'X').and.(spe(j)%name(1:3).NE.'PAH')) THEN
         
         IF (ilab.GT.nelem) THEN
            WRITE(6, *) '*** More fundamental elements than nelem ***'
            STOP
         ENDIF

         ! --------- Save species number
         elem(ilab)%index = j
         ilab = ilab + 1
      ENDIF

      ! ------ Check for electron species number
      IF (spe(j)%name.EQ.speelec) THEN
         indel = j
      ENDIF
   ENDDO



   ! --- Re-arrange order of elements to match species_composition columns (reactions file)
   DO j = 1, nelem-1
      tmp = elem(j)%index
      IF (spe(tmp)%compo(j).NE.1) THEN
         DO k = j+1,nelem
            IF (spe(tmp)%compo(j).EQ.1) THEN
               isptemp=elem(k)%index
               mass_tmp = elem(k)%mass
               name_tmp = elem(k)%name

               elem(k)%mass = elem(j)%mass
               elem(k)%name  = elem(j)%name
               elem(k)%index=elem(j)%index

               elem(k)%mass = mass_tmp
               elem(k)%name = name_tmp
               elem(j)%index=isptemp
            ENDIF
         ENDDO
      ENDIF
   ENDDO
   
   ! Set species characteristics

   ! --- Set reference species
   DO i = 1,nspe
      ! ------ Calculate elemental_masses
      msum = 0.0_dp
      DO k = 1,nelem
         msum=msum + elem(k)%mass * spe(i)%compo(k)
      ENDDO
      spe(i)%mass = msum
      IF (spe(i)%name.EQ.speelec) spe(i)%mass = electron_mass ! electron mass in amu
      IF (spe(i)%name.EQ.spegr0 .OR. spe(i)%name.EQ.spegrm) &
          spe(i)%mass = 4.0_dp * pi * grain_radius**3.0_dp * grain_density / 3.0_dp / amu
   ENDDO

   ! Find ITYPE first and last reactions
   DO i = 0,99
      type_id_start(i)=0
      type_id_stop(i)=0
      DO j=1,nreac
         IF ((reac(j)%itype.EQ.i).AND.(type_id_start(i).EQ.0)) type_id_start(i)=j
         IF (reac(j)%itype.EQ.i) type_id_stop(i)=j
      ENDDO
   ENDDO

   ! Find the index of useful species
   DO i=1,nspe
      IF (spe(i)%name.EQ.speh) indh=i
      IF (spe(i)%name.EQ.spehp) indhp=i
      IF (spe(i)%name.EQ.speh2) indh2=i
      IF (spe(i)%name.EQ.speh2v) indh2v=i
      IF (spe(i)%name.EQ.speco) indco=i
      IF (spe(i)%name.EQ.'C@O        ') ind13co=i
      IF (spe(i)%name.EQ.spehe) indhe=i
      IF (spe(i)%name.EQ.spehep) indhep=i
      IF (spe(i)%name.EQ.spen2) indn2=i
      IF (spe(i)%name.EQ.spegr0) indgr0=i
      IF (spe(i)%name.EQ.spegrm) indgrm=i
      IF (spe(i)%name.EQ.speh2o) indh2o=i
      IF (spe(i)%name.EQ.speo) indo=i
      IF (spe(i)%name.EQ.'C+         ') indcp=i
      IF (spe(i)%name.EQ.'PAH        ') indpah=i
      IF (spe(i)%name.EQ.'PAH-       ') indpahm=i
      IF (spe(i)%name.EQ.'PAH+       ') indpahp=i
      IF (spe(i)%name.EQ.'PAHH+      ') indpahhp=i
      IF (spe(i)%name.EQ.'PAHC+      ') indpahcp=i
      
   ENDDO
   
   RETURN

   END SUBROUTINE INDEX_DATA
   !########################################################

   !########################################################
   SUBROUTINE init_reaction_rates
   !########################################################

   IMPLICIT NONE

   INTEGER :: i,j,k
   INTEGER :: n1, n2, n3
   INTEGER :: r1j,r2j,p1j,p2j,p3j,p4j,p5j
   INTEGER :: r1k,r2k,p1k,p2k,p3k,p4k,p5k
   INTEGER :: npath, nevap, badflag, atoms
   REAL(KIND=dp) :: evfrac, dhfsum, evfrac1, evfrac2
   REAL(KIND=dp) :: sum1, sum2
   REAL(KIND=dp), PARAMETER :: vib_to_dissip_freq_ratio = 1.0E-02_dp ! From Garrod et al. 2007

   ! Set the sticking coefficents
   spe%stick=0.0_dp
   DO i = 1, nspe
      IF (spe(i)%charge.EQ.0) THEN
          spe(i)%stick = 1.0_dp
      ELSE IF (spe(i)%charge.GT.0.0_dp) THEN
          spe(i)%stick = 0.0_dp
      ELSE
          spe(i)%stick = 0.0_dp
      ENDIF
      IF(spe(i)%name.eq.spegr0) spe(i)%stick = 0.0_dp
      IF(i.gt.nspegas) spe(i)%stick = 0.0_dp
   ENDDO
   
   ! Compute the branching ratios
   DO j=1,nreac
      ! ------ Initialise all branching_ratio rate factors, and get species 1 & 2
      reac(j)%br = 1.0_dp
      reac(j)%br1 = 1.0_dp
      reac(j)%br2 = 1.0_dp

      r1j = reac(j)%r(1)
      r2j = reac(j)%r(2)

      p1j = reac(j)%p(1)
      p2j = reac(j)%p(2)
      p3j = reac(j)%p(3)
      p4j = reac(j)%p(4)
      p5j = reac(j)%p(5)

      ! Brancging for surface reactions
      IF (reac(j)%itype.EQ.14) THEN
         npath=0
         ! ------ Check for branching
         DO k=1,nreac
            p1k = reac(k)%p(1)
            IF(reac(k)%itype.EQ.reac(j)%itype) THEN
               IF (((reac(j)%r(1).EQ.reac(k)%r(1)).AND.&
                    (reac(j)%r(2).EQ.reac(k)%r(2))).OR.&
                   ((reac(j)%r(2).EQ.reac(k)%r(1)).AND.&
                    (reac(j)%r(1).EQ.reac(k)%r(2)))) THEN
                   IF ((spe(p1k)%name(1:1).EQ.'J').OR.&
                       (spe(p1k)%name(1:1).EQ.'K')) npath = npath + 1
               ENDIF
            ENDIF
         ENDDO

         ! ------ Compute branching ratios of surface reactions
         IF (NPATH.EQ.0) THEN
           reac(j)%br=0.d0
         ELSE
           reac(j)%br=reac(j)%br/dble(NPATH)
         ENDIF
      

         ! ------ Factor of 2 for same species reactions
         IF (r1j.EQ.r2j) reac(j)%br=reac(j)%br/2.0_dp

         ! ------ Calculate evaporation fraction
         nevap = 0
         DO k=1,nreac

            r1k = reac(k)%r(1)
            r2k = reac(k)%r(2)

            p1k = reac(k)%p(1)
            p2k = reac(k)%p(2)
            p3k = reac(k)%p(3)
            p4k = reac(k)%p(4)
            p5k = reac(k)%p(5)

            IF ((spe(p1j)%name(1:1).EQ.'J').AND.(reac(k)%alpha.GE.1.0e-99_dp)) THEN
               IF ((spe(r1j)%name.EQ.spe(r1k)%name).AND.&
                   (spe(r2j)%name.EQ.spe(r2k)%name).AND.&
                   (spe(p1j)%name(2:).EQ.spe(p1k)%name).AND.&
                   (spe(p2j)%name(2:).EQ.spe(p2k)%name).AND.&
                   (spe(p3j)%name(2:).EQ.spe(p3k)%name).AND.&
                   (spe(p1k)%name(1:1).NE.'J')) nevap = nevap + 1
            ENDIF
            IF ((spe(p1j)%name(1:1).NE.'J').AND.(reac(j)%alpha.GE.1.0e-99_dp)) THEN
               IF ((spe(r1j)%name.EQ.spe(r1k)%name).AND.&
                   (spe(r2j)%name.EQ.spe(r2k)%name).AND.&
                   (spe(p1j)%name.EQ.spe(p1k)%name(2:)).AND.&
                   (spe(p2j)%name.EQ.spe(p2k)%name(2:)).AND.&
                   (spe(p3j)%name.EQ.spe(p3k)%name(2:)).AND.&
                   (spe(p1k)%name(1:1).EQ.'J')) nevap = nevap + 1
            ENDIF
         ENDDO

         n1=0
         n2=0
         n3=0
         DO i=nspegas+1, nspe
            IF(spe(i)%name(1:1).NE.'K') THEN
               IF (spe(p1j)%name(1:1).EQ.'J') THEN
                  IF (spe(p1j)%name.EQ.spe(i)%name) n1 = i
                  IF (spe(p2j)%name.EQ.spe(i)%name) n2 = i
                  IF (spe(p3j)%name.EQ.spe(i)%name) n3 = i
               ENDIF
               IF ((spe(p1j)%name(1:1).NE.'J').AND.&
                   (spe(p1j)%name(1:1).NE.'X')) THEN
                  IF (spe(p1j)%name.EQ.spe(i)%name(2:)) n1 = i
                  IF (spe(p2j)%name.EQ.spe(i)%name(2:)) n2 = i
                  IF (spe(p3j)%name.EQ.spe(i)%name(2:)) n3 = i
               ENDIF
            ENDIF
         ENDDO
         
         dhfsum=spe(r1j)%form_enthalpy+spe(r2j)%form_enthalpy
         IF (n1.NE.0) dhfsum=dhfsum-spe(n1)%form_enthalpy
         IF (n2.NE.0) dhfsum=dhfsum-spe(n2)%form_enthalpy
         IF (n3.NE.0) dhfsum=dhfsum-spe(n3)%form_enthalpy
         ! ------ Convert from kcal to J, from J to K
         dhfsum = dhfsum * 4.184e3_dp / 1.38054E-23_dp
         ! ------ Convert from #moles-1 to #reactions-1
         dhfsum = dhfsum / avogadro

         dhfsum = dhfsum - reac(j)%gamma

         IF (n1.NE.0) sum1 = spe(n1)%bind
         IF (n2.NE.0) sum1 = MAX(spe(n1)%bind,spe(n2)%bind)
         IF (n3.NE.0) sum1 = MAX(spe(n1)%bind,spe(n2)%bind,spe(n3)%bind)
      
         atoms = 0
         IF (n1.NE.0) THEN
            DO k=1,nelem
               atoms = atoms + spe(n1)%compo(k)
            ENDDO
         ENDIF
         
         sum2 = 1.0_dp - (sum1 / (dhfsum + 1.0e-99_dp))
         IF (atoms.EQ.2) sum2=sum2**(3.0_dp*atoms-5.0_dp)
         IF (atoms.GT.2) sum2=sum2**(3.0_dp*atoms-6.0_dp)
         sum2 = vib_to_dissip_freq_ratio * sum2
         evfrac = sum2 / (1.0_dp + sum2)

         badflag = 0
         IF (spe(r1j)%form_enthalpy.LE.-999.0) THEN
            evfrac = 0.0_dp
            evfrac1 = 0.0_dp
            evfrac2 = 0.0_dp
            badflag = badflag + 1
         ENDIF
         IF (spe(r2j)%form_enthalpy.LE.-999.0) THEN
            evfrac = 0.0_dp
            evfrac1 = 0.0_dp
            evfrac2 = 0.0_dp
            badflag = badflag + 1
         ENDIF
         IF(n1.NE.0) THEN
            IF (spe(n1)%form_enthalpy.LE.-999.0) THEN
               evfrac = 0.0_dp
               evfrac1 = 0.0_dp
               evfrac2 = 0.0_dp
               badflag = badflag + 1
            ENDIF
         ENDIF
         IF (n2.NE.0) THEN
            evfrac = 0.0_dp
            evfrac1 = 0.0_dp
            evfrac2 = 0.0_dp
            badflag = badflag + 1
         ENDIF
         IF (n3.NE.0) THEN
            evfrac = 0.0_dp
            evfrac1 = 0.0_dp
            evfrac2 = 0.0_dp
            badflag = badflag + 1
         ENDIF

         IF (evfrac.GE.1.0_dp) evfrac = 1.0_dp
         IF (evfrac.LE.0.0_dp) evfrac = 0.0_dp
         IF (nevap.EQ.0) evfrac = 0.0_dp
         IF (dhfsum.LE.0.0_dp) evfrac = 0.0_dp
         
         IF ((spe(p1j)%name(1:1).EQ.'J').OR.&
             (spe(p1j)%name(1:1).EQ.'K')) THEN
            evfrac = 1.0_dp-evfrac
            evfrac1 = 1.0_dp-evfrac1
            evfrac2 = 1.0_dp-evfrac2
         ENDIF

         reac(j)%br1 = reac(j)%br * evfrac1
         reac(j)%br2 = reac(j)%br * evfrac2
         reac(j)%br = reac(j)%br * evfrac
         

      ENDIF
   ENDDO

   RETURN

   END SUBROUTINE init_reaction_rates
   !########################################################

   !########################################################
   SUBROUTINE get_elem_ab(abin,el_ab)
   !########################################################

   IMPLICIT NONE
   INTEGER :: i,j
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe) :: abin
   REAL(KIND=dp), INTENT(OUT), DIMENSION(nelem) :: el_ab

   el_ab(1:nelem) = 0.0_dp

   DO i = 1,nspe
      DO j = 1,nelem
         if (abin(i)>1.0e-99_dp) el_ab(j) = el_ab(j) + spe(i)%compo(j) * abin(i)
      ENDDO
   ENDDO
   
   END SUBROUTINE get_elem_ab
   !########################################################

   !########################################################
   SUBROUTINE check_conservation(ir,iz,tmp_ab)
   !########################################################

   IMPLICIT NONE

   INTEGER, INTENT(IN) ::ir, iz
   REAL(KIND=dp),INTENT(INOUT) , DIMENSION(nspe) :: tmp_ab
   REAL(KIND=dp), DIMENSION(nelem) :: el_ab
   REAL(KIND=dp) :: chasum
   REAL(KIND=dp) :: dum
   INTEGER :: i, k

   CALL GET_ELEM_AB(tmp_ab,el_ab)

   ! Prevent too low abundances
   DO i=1,nspe
      IF(tmp_ab(i).LT.1.0e-99_dp) THEN
         tmp_ab(i) = 1.0e-99_dp
      ENDIF
   ENDDO

   !--- Conserve electrons
   chasum = 0.0_dp
   DO i=1,nspe
      IF (i.NE.indel) chasum = chasum + spe(i)%charge * tmp_ab(i)
   ENDDO
   IF (chasum.LE.0.0_dp) chasum = 0.0_dp
   tmp_ab(indel) = chasum

   !--- Check conservation
   DO k=1,nelem
      IF(elem(k)%ab.gt.1.0e-35_dp) THEN
      IF (abs(elem(k)%ab - el_ab(k)) / elem(k)%ab .ge. 1.0e-3_dp) THEN
         WRITE(*,'(a22,a6,a17,6es10.3)') '** WARNING ** Element ', elem(k)%name, ' is not conserved', &
         elem(k)%ab,el_ab(k), abs(elem(k)%ab - el_ab(k)) / elem(k)%ab, rdisk(ir)/au,zdisk(ir,iz)/au,&
         zdisk(ir,iz)/rdisk(ir)
      ENDIF
      ENDIF
   ENDDO

   RETURN
   END SUBROUTINE check_conservation
   !#############################################################

   !########################################################
   SUBROUTINE init_relevant_reactions
   !########################################################


   IMPLICIT NONE

   ! Locals
   INTEGER, DIMENSION(nreac, nspe+1) :: is_species_used ! For each species, tell which reactions use it or not
                                                        ! (0 if not used, 1 if used at least once)

   INTEGER :: reaction, species, idx
   INTEGER :: r1, r2, r3

   is_species_used = 0

   DO reaction=1, nreac

      r1 = reac(reaction)%r(1)
      r2 = reac(reaction)%r(2)
      r3 = reac(reaction)%r(3)

      is_species_used(reaction, r1) = 1
      is_species_used(reaction, r2) = 1
      is_species_used(reaction, r3) = 1

   ENDDO

   ! We get the total number of reactions in which each species can be involved
   ! We skip the 'nb_species+1' species that is only a fake species for "no reactant"
   ALLOCATE(nb_reactions_using_species(nspe))
   nb_reactions_using_species(1:nspe) = sum(is_species_used(1:nreac, 1:nspe), 1)

   ! What is the maximum number of reactions involving one particular species?
   max_reactions_same_species = maxval(nb_reactions_using_species(1:nspe))

   ALLOCATE(relevant_reactions(max_reactions_same_species, nspe))

   relevant_reactions(1:max_reactions_same_species, 1:nspe) = 0 ! For the extra elements (because not all species
   ! will have 'max_reactions' reactions involving it).

   ! For each species, we get the references of reactions that have it as a reactant. The number of reactions is different for each species
   ! Thus, at least one species will have a full line of meaningfull indexes. The other will have the rest of their line completed by zeros.
   DO species=1,nspe
      idx = 1
         DO reaction=1, nreac
            IF (is_species_used(reaction, species).eq.1) THEN
               relevant_reactions(idx, species) = reaction
               idx = idx + 1
            ENDIF
         ENDDO
   ENDDO

   RETURN

   END SUBROUTINE init_relevant_reactions
   !########################################################
   
   !########################################################
   SUBROUTINE CheckNetwork
   !########################################################
   
   IMPLICIT NONE
   INTEGER :: i,j,k,l
   INTEGER :: check
   INTEGER :: nprod, ndest
   REAL(kind=dp) :: left_sum, right_sum
   INTEGER, DIMENSION(3) :: rfound
   INTEGER, DIMENSION(3) :: idxr
   INTEGER, DIMENSION(5) :: pfound
   INTEGER, DIMENSION(5) :: idxp
   INTEGER :: nrt1, npt1
   INTEGER :: nrt2, npt2
   
   OPEN(10,file='check_network.log')
   
   !--- Check that: 
   !     - All gas-phase species have a grain equivalent
   !     - Grain species are well initialised
   DO i = 1, nspe
      IF (spe(i)%charge==0 .and. &
          spe(i)%name.ne.'GRAIN0     ' .and. &
          spe(i)%name.ne.'H2v        ' .and. &
          spe(i)%name(1:1) .ne.'J' .and. &
          spe(i)%name(1:1) .ne.'K') THEN
      check = 0
      DO j = 1, nspe
         IF ((spe(j)%name(1:1) .eq.'J' .or. spe(j)%name(1:1) .eq.'K').and. &
             spe(i)%name == spe(j)%name(2:)) THEN
             check = 1
          END IF
      END DO
      IF (check==0) write(10,*) trim(spe(i)%name), ' does not have a J equivalent'
      END IF
   END DO
      
   DO i = 1,nspe
      IF (spe(i)%name(1:1) .ne.'J'.and.spe(i)%name(1:1) .ne.'K'.and.&
          spe(i)%charge==0 .and. spe(i)%name.ne.'GRAIN0     ') THEN
         check = 0
         IF (spe(i)%stick > 0.0) check = 1
         IF (check==0) write(10,*) trim(spe(i)%name), ' have a sticking coef = 0'
      END IF
      IF (spe(i)%name(1:1) .eq.'J'.or.spe(i)%name(1:1) .eq.'K') THEN
         check = 0
         IF (spe(i)%bind > 0.0.and.&
             spe(i)%diff > 0.0.and.&
             spe(i)%mass > 0.0) check = 1
         IF (check==0) THEN
            write(10,*) trim(spe(i)%name), ' have important quantities missing:'
            write(10,*) '   - Ebind = ',spe(i)%bind
            write(10,*) '   - Ediff = ',spe(i)%diff
            write(10,*) '   - Mass = ',spe(i)%mass
         END IF
      END IF
   END DO
   
   ! --- Check that all species have at least one production and 
   !     and a destruction channel
   DO i = 1, nspe
      nprod = 0
      ndest = 0
      DO j = 1, nreac
         IF (any(reac(j)%comp(1:3).eq.spe(i)%name)) ndest = ndest + 1
         IF (any(reac(j)%comp(4:8).eq.spe(i)%name)) nprod = nprod + 1
      END DO
      IF (ndest==0) write(10,*) spe(i)%name, ' does not have any destruction channel'
      IF (nprod==0) write(10,*) spe(i)%name, ' does not have any production channel'
   END DO
   
   ! --- Check that all reactions are equilibrated
   DO i = 1, nreac
      !-- Loop over prime element
      DO k = 1, nelem
         !-- Reactant
         left_sum = 0.0_dp
         DO j = 1, 3
            IF (reac(i)%r(j).ne.nspe+1) THEN
               left_sum = left_sum + spe(reac(i)%r(j))%compo(k)
            END IF
         END DO
         !--
         
         !-- Products
         right_sum = 0.0_dp
         DO j = 1, 5
            IF (reac(i)%p(j).ne.nspe+1) THEN
               right_sum = right_sum + spe(reac(i)%p(j))%compo(k)
            END IF
         END DO
         !--
      
         IF (left_sum.ne.right_sum) THEN
            PRINT*, '** ERROR ** This reaction is not equilibrated in ', trim(elem(k)%name),':'
            PRINT*, reac(i)%itype, reac(i)%comp(1:3),'->  ',reac(i)%comp(4:8)
            STOP
         END IF
      END DO
   END DO
   
   
   ! --- Check for duplicated reactions
   DO i = 1, nreac
      DO j = i+1, nreac
         IF (reac(i)%itype.eq.reac(j)%itype) THEN
         
            ! --- Cycle over the reactants
            nrt1 = 0
            rfound = 0
            idxr = 0
            DO k = 1, 3
               if(reac(i)%r(k).ne.nspe+1) then
                  nrt1 = nrt1 + 1
                  nrt2 = 0
                  DO l = 1,3
                     if(reac(j)%r(l).ne.nspe+1) then
                        nrt2 = nrt2 + 1
                        IF(reac(i)%r(k)==reac(j)%r(l).and.idxr(l).ne.reac(j)%r(l)) THEN
                           rfound(k) = 1
                           idxr(l) = reac(j)%r(l)
                        END IF 
                     end if
                  END DO
               end if
            END DO
            
            ! --- Cycle over the products
            npt1 = 0
            pfound = 0
            idxp = 0
            DO k = 1, 5
               if(reac(i)%p(k).ne.nspe+1) then
                  npt1 = npt1 + 1
                  npt2 = 0
                  DO l = 1,5
                     if(reac(j)%p(l).ne.nspe+1) THEN
                        npt2 = npt2 + 1
                        IF(reac(i)%p(k)==reac(j)%p(l).and.idxp(l).ne.reac(j)%p(l)) THEN
                           pfound(k) = 1
                           idxp(l) = reac(j)%p(l)
                        END IF 
                     END IF
                  END DO
               end if
            END DO
                                    
            IF (sum(rfound)==nrt1 .and. sum(pfound)==npt1 .and. &
               nrt1 == nrt2 .and. npt1 == npt2 .and. &
               reac(i)%id .ne. reac(j)%id) THEN
               PRINT*, '** WARNING ** Reaction present twice'
               PRINT*, i,reac(i)%id, reac(i)%comp, reac(i)%itype,&
                       reac(i)%alpha,reac(i)%beta, reac(i)%gamma
               PRINT*, j,reac(j)%id, reac(j)%comp, reac(j)%itype,&
                       reac(j)%alpha,reac(j)%beta, reac(j)%gamma
               STOP
            END IF 
              
         END IF
      END DO
   END DO
   
   close(10)
   
   END SUBROUTINE CheckNetwork
   !########################################################
   
   !########################################################
   SUBROUTINE get_major_bearing_species(abin)
   !########################################################

   USE auxilary

   IMPLICIT NONE
   INTEGER :: i,j
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe) :: abin
   
   REAL(KIND=dp), DIMENSION(nspe) :: ab_main
   INTEGER, DIMENSION(nspe) :: index
   REAL(KIND=dp) :: dummy

   ab_main(:) = 0.0_dp

   DO j = 1,nelem
      DO i = 1,nspe
         ab_main(i) = spe(i)%compo(j) * abin(i)
      ENDDO
      call get_sorted_index(ab_main,index)
      dummy = 0.0_dp
      DO i = nspe,1,-1
         dummy = dummy + ab_main(index(i))
         print*, elem(j)%name, spe(index(i))%name, abin(index(i)), ab_main(index(i))*100.0/elem(j)%ab, &
         dummy*100.0/elem(j)%ab, elem(j)%ab
      END DO
      pause
   ENDDO

   END SUBROUTINE get_major_bearing_species
   !########################################################
   
END MODULE init_chem
