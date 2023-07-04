MODULE CHEMISTRY

   USE GLOBAL_VARIABLES
   
   PRIVATE
   PUBLIC :: set_constant_rates, set_dependant_rates, set_dependant_rates_3phase, shielding_co, GET_CR_RATE
   
   CONTAINS   
   
   !########################################################
   SUBROUTINE SET_CONSTANT_RATES
   !########################################################
   
   IMPLICIT NONE
   
   REAL(KIND=dp) :: ratetmp
   REAL(KIND=dp)                         :: T300, TSQ, T_EV
   INTEGER                               :: nsta, nfin
   INTEGER                               :: j, w, m, n
   INTEGER                               :: r1,r2
   INTEGER, DIMENSION(10)                :: indice
   REAL(KIND=dp), DIMENSION(10)          :: distmin, distmax
   REAL(KIND=dp)                         :: evaporation_rates
   REAL(KIND=dp)                         :: nu0
   REAL(KIND=dp)                         :: barr
   REAL(KIND=dp)                         :: tt
   REAL(KIND=dp),PARAMETER :: ampe = 1836.1_dp
   INTEGER :: ii, nn
   
   T300 = Tgas / 300.0_dp
   TSQ = SQRT(Tgas)
   tt = Tgas * 1.0e-3_dp + 1.0_dp
   
   nn = size(idx_const_rates)
   
   DO ii = 1, nn
      
   j = idx_const_rates(ii)
   reac(j)%rate = 0.0_dp
      
   SELECT CASE(reac(j)%itype)
   
   !------------------------------------------------------------
   !
   !      --- GAS PHASE REACTIONS
   !
   !------------------------------------------------------------
   
   !-------------------------------------------------------------
   !  (0) Gas phase reactions with grains (s-1)
   !      - Needs to be computed in the model because rates are 
   !        for 0.1 micro meter grains (see Wakelam 2012 - KIDA)
   !-------------------------------------------------------------
   
   CASE(0)
      r2 = reac(j)%r(2)
      reac(j)%rate = pi * grain_radius**2 * sqrt(8.0_dp*k_b*Tgas / pi / amu / spe(r2)%mass)
   
   !-------------------------------------------------------------
   !  (4-8) Bimolecular gas phase reactions
   !        Several possible formula:
   !         - Kooij (modIFied Arhenius) - formula = 3
   !         - Ionpol1 - formula = 4
   !         - Ionpol2 - formula = 5
   !-------------------------------------------------------------
   
   w = 1
   distmin(:) = 9999.0_dp
   distmax(:) = 9999.0_dp   
   CASE(4:8)  
   
      ! -- Kooij (formula = 3)--
      IF (reac(j)%formula.EQ.3) THEN
         reac(j)%rate = reac(j)%alpha * (T300**reac(j)%beta) * exp(-reac(j)%gamma / Tgas)
   
         ! Check for temperature bounderies
         IF (Tgas.LT.reac(j)%tmin) THEN
            reac(j)%rate = reac(j)%alpha * ((reac(j)%tmin / 300.0_dp)**reac(j)%beta) * exp(-reac(j)%gamma / reac(j)%tmin)
         ENDIF
   
         ! Check for the presence of several rate coefficients present in the network for the same reaction
         IF (reac(j+1)%id.EQ.reac(j)%id) THEN
            indice(w) = j
            distmin(w) = reac(j)%tmin - Tgas
            distmax(w) = Tgas - reac(j)%tmax
            w = w + 1
         ENDIF
   
         IF ((reac(j+1)%id.NE.reac(j)%id).AND.(w.NE.1)) THEN
            indice(w)  = J
            distmin(w) = reac(j)%tmin - Tgas
            distmax(w) = Tgas - reac(j)%tmax
            DO m=1,w
               n = indice(m)
               IF (Tgas.LT.reac(n)%tmin .and. n/=0) reac(n)%rate = 0.0_dp
               IF (Tgas.GT.reac(n)%tmax .and. n/=0) reac(n)%rate = 0.0_dp
            ENDDO
            IF (maxval(reac(indice(1:w))%rate).LT.1.0e-99_dp) THEN
               IF (minval(abs(distmin)).LT.minval(abs(distmax))) THEN
                  n=indice(minloc(abs(distmin),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * ((reac(n)%tmin/300.0_dp)**reac(n)%beta) * & 
                                          exp(-reac(n)%gamma / reac(n)%tmin)
               ELSE
                  n=indice(minloc(abs(distmax),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * ((reac(n)%tmax/300.0_dp)**reac(n)%beta) * &
                                          exp(-reac(n)%gamma / reac(n)%tmax)
               ENDIF
            ENDIF
            w=1
            indice(:)=0
            distmin(:)=9999.0_dp
            distmax(:)=9999.0_dp
         ENDIF
   
      ENDIF
   
      ! -- Ionpol1 (formula = 4)--
      IF (reac(j)%formula.EQ.4) THEN
         reac(j)%rate = reac(j)%alpha * reac(j)%beta * (0.62_dp+0.4767_dp * &
                        reac(j)%gamma * ((300.0_dp/Tgas)**0.5_dp))
   
         ! Check for temperature bounderies
         IF (Tgas.LT.reac(j)%tmin) THEN
            reac(j)%rate = reac(j)%alpha * reac(j)%beta * (0.62_dp+0.4767_dp * &
                           reac(j)%gamma * ((300.0_dp/reac(j)%tmin)**0.5_dp))
         ENDIF
   
         ! Check for the presence of several rate coefficients present in the network for the same reaction
         IF (reac(j+1)%id.EQ.reac(j)%id) THEN
            indice(w)=j
            distmin(w) = reac(j)%tmin - Tgas
            distmax(w) = Tgas - reac(j)%tmax
            w = w + 1
         ENDIF
   
         IF ((reac(j+1)%id.NE.reac(j)%id).AND.(W.NE.1)) THEN
            indice(w)=j
            distmin(w) = reac(j)%tmin - Tgas
            distmax(w) = Tgas - reac(j)%tmax
            DO m=1,w
               n=indice(m)
               IF (Tgas.LT.reac(n)%tmin .and. n/=0) reac(n)%rate = 0.0_dp
               IF (Tgas.GT.reac(n)%tmax .and. n/=0) reac(n)%rate = 0.0_dp
            ENDDO
   
            IF (maxval(reac(indice(1:w))%rate).LT.1.0e-99_dp) THEN
               IF (minval(abs(distmin)).LT.minval(abs(distmax))) THEN
                  n=indice(minloc(abs(distmin),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * reac(n)%beta * (0.62_dp+0.4767_dp * &
                                          reac(n)%gamma * ((300.0_dp/reac(n)%tmin)**0.5_dp))
               ELSE
                  n=indice(minloc(abs(distmax),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * reac(n)%beta * (0.62_dp+0.4767_dp * &
                                          reac(n)%gamma * ((300.0_dp/reac(n)%tmax)**0.5_dp))
               ENDIF
            ENDIF
            w=1
            indice(:)=0
            distmin(:)=9999.0_dp
            distmax(:)=9999.0_dp
         ENDIF
      ENDIF
   
      ! -- Ionpol2 (formula = 5)--
      IF (reac(j)%formula.EQ.5) THEN
   
         ! Check for temperature boundaries and apply the formula in consequence
         IF (Tgas.LT.reac(j)%tmin) THEN
            reac(j)%rate = reac(j)%alpha * reac(j)%beta * &
                          ((1.0_dp+0.0967_dp * reac(j)%gamma * (300.0_dp/reac(j)%tmin)**0.5_dp) + &
                          (reac(j)%gamma**2.0_dp * 300.0_dp / (10.526_dp * reac(j)%tmin)))
         !ELSEIF (Tgas.GT.reac(j)%tmax) THEN
         !   reac(j)%rate = reac(j)%alpha * reac(j)%beta * &
         !                 ((1.0_dp+0.0967_dp * reac(j)%gamma * (300.0_dp/reac(j)%tmax)**0.5_dp) + &
         !                 (reac(j)%gamma**2.0_dp * 300.0_dp / (10.526_dp * reac(j)%tmax)))
         ELSE
            reac(j)%rate = reac(j)%alpha * reac(j)%beta * &
                          ((1.0_dp+0.0967_dp * reac(j)%gamma * (300.0_dp/Tgas)**0.5_dp) + &
                          (reac(j)%gamma**2.0_dp * 300.0_dp / (10.526_dp * Tgas)))
         ENDIF
   
         ! Check for the presence of several rate coefficients present in the network for the same reaction
         IF (reac(j+1)%id.EQ.reac(j)%id) THEN
            indice(w)=j
            distmin(w)=reac(j)%tmin-Tgas
            distmax(w)=Tgas-reac(j)%tmax
            w = w + 1
         ENDIF
   
         IF ((reac(j+1)%id.NE.reac(j)%id).AND.(W.NE.1)) THEN
            indice(w)=j
            distmin(w)=reac(j)%tmin-Tgas
            distmax(w)=Tgas-reac(j)%tmax
            DO m=1,w
               n=indice(m)
               IF (Tgas.LT.reac(n)%tmin .and. n/=0) reac(n)%rate = 0.0_dp
               IF (Tgas.GT.reac(n)%tmax .and. n/=0) reac(n)%rate = 0.0_dp
            ENDDO
            IF (maxval(reac(indice(1:w))%rate).LT.1.0d-99) THEN
               IF (minval(abs(distmin)).LT.minval(abs(distmax))) THEN
                  n=indice(minloc(abs(distmin),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * reac(n)%beta * &
                                          ((1.0_dp+0.0967_dp * reac(n)%gamma * (300.0_dp/reac(n)%tmin)**0.5_dp) + &
                                          (reac(n)%gamma**2.0_dp * 300.0_dp / (10.526_dp * reac(n)%tmin)))
               ELSE
                  n=indice(minloc(abs(distmax),dim=1))
                  IF(n/=0) reac(n)%rate = reac(n)%alpha * reac(n)%beta * &
                                         ((1.0_dp+0.0967_dp * reac(n)%gamma * (300.0_dp/reac(n)%tmax)**0.5_dp) + &
                                         (reac(n)%gamma**2.0_dp * 300.0_dp / (10.526_dp * reac(n)%tmax)))
               ENDIF
            ENDIF
            w=1
            indice(:)=0
            distmin(:)=9999.0_dp
            distmax(:)=9999.0_dp
         ENDIF
      ENDIF
   
   !-------------------------------------------------------------
   !  (11) H2* radiative de-excitation
   !-------------------------------------------------------------
   CASE(11)
      r1 = reac(j)%r(1)
      reac(j)%rate = reac(j)%alpha

   !-------------------------------------------------------------
   !  (12) H2* collisional de-excitation
   !-------------------------------------------------------------
   CASE(12)
      r1 = reac(j)%r(1)
      r2 = reac(j)%r(2)
      IF(spe(r1)%name.eq.speh2v.and.spe(r2)%name.eq.speh) THEN
         reac(j)%rate = 10.0_dp**(-11.06_dp + 0.0555_dp/tt - 2.390_dp/tt**2)
      ELSE IF(spe(r1)%name.eq.speh2v.and.spe(r2)%name.eq.speh2) THEN
         reac(j)%rate = 10.0_dp**(-11.08_dp - 3.671_dp/tt - 2.023_dp/tt**2)
      END IF
         
   END SELECT    
   
   if (reac(j)%rate.lt.1.0e-99_dp) reac(j)%rate = 0.0_dp   
   
   END DO
   
   RETURN
   
   END SUBROUTINE SET_CONSTANT_RATES
   !########################################################
   
   !########################################################
   SUBROUTINE SET_DEPENDANT_RATES(abin)
   !########################################################
   
   USE auxilary
   
   IMPLICIT NONE
   
   ! Inputs
   REAL(KIND=dp), INTENT(IN), DIMENSION(nspe) :: abin
   
   ! Locals
   REAL(KIND=dp) :: probh2h2 ! probability for the encounter desorption process.
   REAL(KIND=dp) :: tetah2, tetah2up, tetah2dn
   REAL(KIND=dp) :: tetah2v, tetah2vup, tetah2vdn
   REAL(KIND=dp) :: tetaco, tetacoup, tetacodn
   REAL(KIND=dp) :: tetan2, tetan2up, tetan2dn
   REAL(KIND=dp) :: T300, TSQ
   INTEGER :: j, l, m, k
   INTEGER :: r1,r2
   INTEGER :: p1,p2
   INTEGER :: imod1, imod2   
   REAL(KIND=dp) :: ab_tot  !< Sum of all abundances on grain (surface+mantle)
   REAL(KIND=dp) :: ab_surf !< Sum of all abundances on grain surface
   REAL(KIND=dp) :: ab_mant !< Sum of all abundances on grain mantle
   REAL(KIND=dp) :: ab_lay  !< Abundance in 1 layer
   REAL(KIND=dp) :: sumlaysurf
   REAL(KIND=dp) :: sumlaymant
   REAL(KIND=dp) :: sumlay !< Total number of layers on the grain surface
   REAL(KIND=dp) :: abH2O   !< H2O abundance on grains surface
   REAL(KIND=dp) :: abNH3   !< NH3 abundance on grains surface
   REAL(KIND=dp) :: abCO2   !< CO2 abundance on grains surface
   REAL(KIND=dp) :: abCH4   !< CH4 abundance on grains surface
   REAL(KIND=dp) :: abCH3OH !< CH3OH abundance on grains surface
   REAL(KIND=dp) :: abCO    !< CO abundance on grains surface
   REAL(KIND=dp) :: prob_reac, prob_diff  
   REAL(KIND=dp) :: activ, barr
   REAL(KIND=dp) :: acc_prefac
   REAL(KIND=dp) :: evaporation_rates
   REAL(KIND=dp) :: diff, diff1, diff2
   REAL(KIND=dp) :: nu0, nu01, nu02
   REAL(KIND=dp) :: qtun1, qtun2
   REAL(KIND=dp) :: surf_react_proba
   REAL(KIND=dp) :: redmas
   REAL(KIND=dp) :: uvcr
   INTEGER :: idx, idsave
   REAL(KIND=dp) :: sigma
   REAL(KIND=dp) :: zee_s, a1c
   REAL(KIND=dp) :: attenu, Pabs
   REAL(KIND=dp) :: a, b, c , x, dummy, dens_tot
   REAL(KIND=dp) :: new_sigma, da
   REAL(kind=dp) :: nfac1, nfac2, nfac3
   REAL(KIND=dp) :: tmp, xx, scale_rot, scale_vib, tmp_rate
   INTEGER :: ii, nn
   
   
   T300=Tgas/300.0_dp
   TSQ=SQRT(Tgas)
   
   ab_tot   = 0.0_dp
   ab_surf  = 0.0_dp
   ab_mant  = 0.0_dp
   abCO     = 0.0_dp
   abH2O    = 0.0_dp
   abNH3    = 0.0_dp
   abCH4    = 0.0_dp
   abCH3OH  = 0.0_dp
   dens_tot = 0.0_dp
   DO j = nspegas+1,nspe
      IF (abin(j).ge.1.0e-99_dp) THEN
      IF(spe(j)%name(1:1).ne."J" .and.spe(j)%name(1:1).ne."K" ) print*, 'PROBLEM !!!!'
      ab_tot = ab_tot + abin(j)
      dens_tot = dens_tot + abin(j)*spe(j)%mass*amu*density
      IF(spe(j)%name(1:1).eq."J") ab_surf = ab_surf + abin(j)
      IF(spe(j)%name(1:1).eq."K") ab_mant = ab_mant + abin(j)
      IF(spe(j)%name.EQ.'JCO        ') abCO  = abin(j)
      IF(spe(j)%name.EQ.'JH2O       ') abH2O = abin(j)
      IF(spe(j)%name.EQ.'JNH3       ') abNH3 = abin(j)
      IF(spe(j)%name.EQ.'JCO2       ') abCO2 = abin(j)
      IF(spe(j)%name.EQ.'JCH4       ') abCH4 = abin(j)
      IF(spe(j)%name.EQ.'JCH3OH     ') abCH3OH = abin(j)
   ENDIF
   ENDDO
      
   nsites = 4.0*ndsigma*surface_site_density
   sumlay = ab_tot / nsites
   dummy = (3.0_dp*ab_tot*density/(4.0_dp*pi*ndust))**(1.0_dp/3.0_dp)
   if (sumlay.gt.dummy .and. dummy.ge.1.0_dp) then
      nsites = nsites * sumlay/dummy
   endif
   
   ab_lay = nsites
   sumlay = ab_tot / nsites
   sumlaysurf = ab_surf / nsites
   sumlaymant = ab_mant / nsites
   Nlaysurfsave = sumlaysurf
   Nlaymantsave = sumlaymant
   
   Pabs = 7.0e-3_dp
   attenu = (Pabs-1.0_dp)*((1.0_dp - Pabs)**sumlay-1.0_dp)/Pabs/sumlay
   if (attenu.gt.1.0_dp) attenu = 1.0_dp
   acc_prefac = ndsigma * spe(indH2)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(indH2)%mass)
   uvcr = CRrate / 1.3e-17_dp
   
   nn = size(idx_var_rates)
   DO ii = 1, nn
      
      j = idx_var_rates(ii)
      
      r1 = reac(j)%r(1)
      r2 = reac(j)%r(2)
      
      p1 = reac(j)%p(1)
      p2 = reac(j)%p(2)
      
      reac(j)%rate = 0.0_dp
      IF (is_grain_reactions.eq.0 .and. reac(j)%itype>=14) CYCLE
      
      SELECT CASE(reac(j)%itype)
   
   !------------------------------------------------------------
   !
   !      --- GAS PHASE REACTIONS
   !
   !------------------------------------------------------------
   
   !-------------------------------------------------------------
   !  (1) Photodissoc/ionisation with cosmic and X rays (s-1)
   !-------------------------------------------------------------
   CASE(1)
        zee_s=0.0_dp
        ! Secondary correction (Eq.A1 of Wolfire+1995)
        ! distribute according to fractional h and h2 abundances
        if(reac(j)%r(1).eq.indh) then ! H ionization by CR
            a1c=0.3908_dp*(1.0_dp-abin(indel)**0.4092_dp)**1.7592_dp
            zee_s=(35.0_dp/13.6_dp-1.0_dp)*a1c*(35.0_dp/1.0e3_dp)**(2.313_dp*abin(indel)) & 
                * (abin(indh)/(abin(indh)+abin(indh2)))
        elseif(reac(j)%r(1).eq.indh2.and.reac(j)%p(2).eq.indel) then ! H2
            a1c=0.3908_dp*(1.0_dp-abin(indel)**0.4092_dp)**1.7592_dp
            zee_s=(35.0_dp/13.6_dp-1.0_dp)*a1c*(35.0_dp/1.0e3_dp)**(2.313_dp*abin(indel)) & 
                * (abin(indh2)/(abin(indh)+abin(indh2)))
          elseif(reac(j)%r(1).eq.indhe) then ! Helium
            a1c=0.0554_dp*(1.0_dp-abin(indel)**0.4614_dp)**1.6660_dp
            zee_s=(35.0_dp/24.6_dp-1.0_dp)*a1c*(35.0_dp/1.0e3_dp)**(2.313_dp*abin(indel))
          endif 
          reac(j)%rate = reac(j)%alpha * CRrate * (1.0_dp + zee_s)
   
   
   !-------------------------------------------------------------
   !  (2) Gas phase photodissociations/ionisations by UV induced
   !      by cosmic rays.
   !      Now it takes into account the abundance of H2 but it
   !      still assumes that the albedo of the grains is 0.5,
   !      hence the 0.5 in the equation. Check Wakelam et al. (2012)
   !      to understand the formula. This is now time dependent.
   !-------------------------------------------------------------
   CASE(2)
      reac(j)%rate = reac(j)%alpha * CRrate * abin(indh2)/0.5_dp
   
   
   !-------------------------------------------------------------
   !  (3) Gas phase photodissociations/ionisations by UV of H2, CO
   !      and N2
   !-------------------------------------------------------------
   CASE(3)
      r1 = reac(j)%r(1)
      p1 = reac(j)%p(1)
      p2 = reac(j)%p(2)
   
      reac(j)%rate=reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                   reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                   reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)
   
      ! *** H2 self-shielding ***
      IF (spe(r1)%name.EQ.speh2) THEN
            call shielding_h2(tauh2,tetah2)
            call shielding_h2(tauuph2,tetah2up)
            call shielding_h2(taudnh2,tetah2dn)
            reac(j)%rate= tetah2   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                          tetah2up * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                          tetah2dn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
            if (spe(p1)%name.EQ.speh2v) then   
               reac(j)%rate = 9.0_dp * reac(j)%rate   
            endif         
   
      ENDIF
      ! ***
      
      ! *** H2* self-shielding ***
      IF (spe(r1)%name.EQ.speh2v) THEN
         call shielding_h2(tauh2v,tetah2v)
         call shielding_h2(tauuph2v,tetah2vup)
         call shielding_h2(taudnh2v,tetah2vdn)                       
         reac(j)%rate= tetah2v   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetah2vup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetah2vdn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
         ! --- Add direct continuum rate for H2*
         reac(j)%rate = reac(j)%rate + 1.0e-11_dp * (  G0Hab * exp(-reac(j)%gamma *   Av) + &
                                                     G0upHab * exp(-reac(j)%gamma * Avup) + &
                                                     G0dnHab * exp(-reac(j)%gamma * Avdn))
      ENDIF
      ! ***
      
      ! *** CO self-shielding and isotopologues ***
      IF (spe(r1)%name.EQ.speco) THEN
         call shielding_co(nh2,nco+n13co,1,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,1,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,1,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
      ENDIF
      IF (spe(r1)%name.EQ.'C@O        ') THEN
         call shielding_co(nh2,nco+n13co,2,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,2,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,2,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
      ENDIF
      IF (spe(r1)%name.EQ.'CO@        ') THEN         
         call shielding_co(nh2,nco+n13co,3,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,3,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,3,tetacodn) 
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)  
         
      ENDIF
      IF (spe(r1)%name.EQ.'C@O@       ') THEN
         call shielding_co(nh2,nco+n13co,4,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,4,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,4,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)         
      ENDIF
   
      ! *** N2 self-shielding ***
      IF (spe(r1)%name.EQ.spen2) THEN
         call shielding_n2(nh2,nn2,tetan2)
         call shielding_n2(nh2_up,nn2_up,tetan2up)
         call shielding_n2(nh2_dn,nn2_dn,tetan2dn)
         reac(j)%rate= tetan2   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetan2up * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetan2dn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)
      ENDIF
      ! ****
      

   !-------------------------------------------------------------
   !  (10) Ad-hoc formation of h2 on grains when
   !       is_grain_reactions = 0 is_h2_adhoc_form = 1
   !-------------------------------------------------------------
   CASE(10)
      IF ((IS_GRAIN_REACTIONS.eq.0).OR.(IS_H2_ADHOC_FORM.eq.1)) then
            reac(j)%rate = reac(j)%alpha * (ndsigma/5.0e-22_dp) / abin(r1)
      ELSE
            reac(j)%rate = 0.0_dp
      ENDIF
   
   !------------------------------------------------------------
   !
   !      --- GRAIN REACTIONS
   !
   !------------------------------------------------------------
      
   !-------------------------------------------------------------
   !  (99) Accretion rates on grain surfaces
   !-------------------------------------------------------------
   
   CASE(99)
      r1 = reac(j)%r(1)
      if (is_growth == 0) then
         acc_prefac = ndsigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
      else
         if (first_step_done) then
            da = sumlay*ds
            dummy = 1.0_dp + 2.0_dp*ndusta*da/ndusta2 + ndust*da**2/ndusta2
            new_sigma = ndsigma * dummy
            acc_prefac = new_sigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
         else
            acc_prefac = ndsigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
         endif
      end if
      IF((spe(r1)%name.eq.speh).or.(spe(r1)%name.eq.speh2)) THEN
         CALL sticking_special_cases(j,sumlay,acc_prefac)
      ENDIF
      reac(j)%rate = reac(j)%alpha * reac(j)%br * acc_prefac * TSQ * density
   
   !-------------------------------------------------------------
   !  (15) Thermal evaporation (s-1)
   !-------------------------------------------------------------  
   CASE(15)
      r1 = reac(j)%r(1)
      nu0 = SQRT(2.0_dp * k_b / pi**2.0_dp / amu * surface_site_density * &
                            spe(r1)%bind / spe(r1)%mass)
      evaporation_rates = nu0 * exp(-spe(r1)%bind / Tdust)
      reac(j)%rate = reac(j)%alpha * reac(j)%br * evaporation_rates
      IF (sumlay.gt.nb_active_lay)  reac(j)%rate = reac(j)%rate*nb_active_lay/sumlay
   
   
   !-------------------------------------------------------------
   !  (16) Cosmic-ray evaporation (s-1)
   !-------------------------------------------------------------
   CASE(16)
      r1 = reac(j)%r(1)
      nu0 = SQRT(2.0_dp * k_b / pi**2.0_dp / amu * surface_site_density * &
                            spe(r1)%bind / spe(r1)%mass)
      reac(j)%rate = reac(j)%alpha * reac(j)%br * (CRrate / 1.3e-17_dp) &
                     * nu0 * fe_ionisation_rate * cr_peak_duration &
                     * exp(-spe(r1)%bind / cr_peak_grain_temp)
      IF (sumlay.gt.nb_active_lay)  reac(j)%rate = reac(j)%rate*nb_active_lay/sumlay
      reac(j)%rate = 0.0_dp
   
   !-------------------------------------------------------------
   !  (66) Photodesorption by external UV
   !       1.d8 is I_ISRF-FUV from
   !       Oberg et al. 2007, ApJ, 662, 23
   !-------------------------------------------------------------
   CASE(66)
      reac(j)%rate = Yphdes / surface_site_density * 1.0e8_dp * &
                                        (  G0Hab * exp(-2.0_dp *   Av) &
                                       + G0upHab * exp(-2.0_dp * Avup) &
                                       + G0dnHab * exp(-2.0_dp * Avdn))
      IF (sumlay.gt.nb_active_lay)  reac(j)%rate = reac(j)%rate*nb_active_lay/sumlay
      IF (is_photodesorb.Eq.0) reac(j)%rate = 0.0_dp   
   
   !-------------------------------------------------------------
   !  (67) Photodesorption by CR generated UV
   !-------------------------------------------------------------   
   CASE(67)
      reac(j)%rate = 1.0e4_dp * uvcr * Yphdes / surface_site_density / 4.0_dp
      IF (sumlay.gt.nb_active_lay)  reac(j)%rate = reac(j)%rate*nb_active_lay/sumlay
      IF (is_photodesorb.Eq.0) reac(j)%rate = 0.0_dp
   
   !-------------------------------------------------------------
   !  (17 and 18) Photodiss by Cosmic rays on grain surfaces (s-1)
   !-------------------------------------------------------------
   CASE(17:18)
      r1 = reac(j)%r(1)
      reac(j)%rate = reac(j)%alpha * CRrate * abin(indh2) / 0.5_dp      
      ! Photon attenuation in the mantle with an probability per monolayer Pabs 
      r1 = reac(j)%r(1)
      if(spe(r1)%name(1:1).eq.'K') then
         reac(j)%rate = reac(j)%rate * attenu
      end if
   
   !-------------------------------------------------------------
   !  (14) Grain surface reactions
   !-------------------------------------------------------------
   CASE(14)
   
      r1 = reac(j)%r(1)
      r2 = reac(j)%r(2)
   
      imod1 = 0
      imod2 = 0
      barr = 1.0d0
   
      ! --------- Thermal hopping diffusion method
      nu01 = SQRT(2.0_dp * k_b / pi**2.0_dp / amu * surface_site_density &
                            * spe(r1)%bind / spe(r1)%mass)
      diff1 = nu01 * exp(-spe(r1)%diff / Tdust) / nsites
      qtun1 = nu01 / nsites * exp(-2.0_dp*diffusion_barrier_thickness/h_barre*&
                              SQRT(2.0_dp*amu*spe(r1)%mass*K_B*spe(r1)%diff))
      nu02 = SQRT(2.0_dp * k_b / pi**2.0_dp / amu * surface_site_density &
                            * spe(r2)%bind / spe(r2)%mass)
      diff2 = nu02 * exp(-spe(r2)%diff / Tdust) / nsites
      qtun2 = nu02 / nsites * exp(-2.0_dp*diffusion_barrier_thickness/h_barre*&
                              SQRT(2.0_dp*amu*spe(r2)%bind*K_B*spe(r2)%diff))
   
      ! --------- Check for JH,JH2
      IF (spe(r1)%name.EQ.spejh)  imod1=1
      IF (spe(r1)%name.EQ.spejh2) imod1=2
      IF (spe(r2)%name.EQ.spejh)  imod2=1
      IF (spe(r2)%name.EQ.spejh2) imod2=2
   
      ! --------- QM for JH,JH2 only - others are too heavy
      IF (imod1+imod2.NE.0) THEN
         IF (grain_tunneling_diffusion.EQ.1) THEN
            IF ((imod1.NE.0).AND.(qtun1.GT.diff1)) THEN
               diff1=qtun1
            ENDIF
            IF ((imod2.NE.0).AND.(qtun2.GT.diff2)) THEN
               diff2=qtun2
            ENDIF
         ENDIF
      ENDIF
      
      ! ------ For reactions with barrier
      IF (reac(j)%gamma.GE.1.0e-99_dp) THEN
         ! ------ Calculate quantum activation energy
         redmas = spe(r1)%mass * spe(r2)%mass / (spe(r1)%mass + spe(r2)%mass)
         surf_react_proba = 2.0_dp * chemical_barrier_thickness / h_barre * sqrt(2.0_dp * amu * redmas * k_b * reac(j)%gamma)
   
         ! --------- Calculate activation energy barrier multiplier
         activ = reac(j)%gamma / Tdust
         ! ------------ Choose fastest of classical or tunnelling
         IF (activ.GT.surf_react_proba .and. (any(reac(j)%comp(1:3).eq.spejh) .or. &
            any(reac(j)%comp(1:3).eq.spekh).or. any(reac(j)%comp(1:3).eq.spejh2) .or. &
            any(reac(j)%comp(1:3).eq.spekh2))) then
            activ = surf_react_proba
         ENDIF
         barr=exp(-activ)
         IF(is_reac_diff==1) then
            prob_reac = max(nu01,nu01) * exp(-activ)
            prob_diff = (diff1 + diff2) * nsites
   
            barr = prob_reac+prob_diff
            barr = prob_reac / barr
   
         ENDIF
      ENDIF
   
      diff = diff1 + diff2
      reac(j)%rate = reac(j)%alpha * reac(j)%br * barr * diff / density
      
      ! IF the number of mantle layer is > 1, we consider that t(diff) (the time required by a species to
      ! scan the entire grain sites) is given by the Number of sites on a layer time the number of layer times t(hop)
      IF((spe(r1)%name(1:1).eq.'K').and.(spe(r2)%name(1:1).eq.'K').and.sumlaymant.gt.1.0_dp) THEN
         reac(j)%rate = reac(j)%rate / sumlaymant
      ENDIF
   
      ! H2 formation by LH mechanism is turned off when the ad hoc formation of H2 is activated
      IF ((spe(r1)%name.EQ.spejh).AND.(spe(r2)%name.EQ.spejh)) THEN
         IF(is_h2_adhoc_form.eq.1) reac(j)%rate = 0.0_dp
      ENDIF
   
      ! "Encounter desorption" process for JH2 (Hincelin et al. 2014,A&A)
      ! The reaction JH2+JH2->JH2+H2 must be in the grain_reactions.in file to be accounted
      ! in practice the DOminant processes are really the thermal hoping and thermal desorption of H2.
      IF ((spe(r1)%name.EQ.spejh2).AND.(spe(r2)%name.EQ.spejh2)) THEN
         reac(j)%rate = reac(j)%alpha * barr * diff / density
         nu0 = SQRT(2.0_dp * k_b / pi**2.0_dp / amu * surface_site_density * &
                               spe(r1)%bind / spe(r1)%mass)
         evaporation_rates = nu0 * exp(-23.0_dp / Tdust)   
         probh2h2 = evaporation_rates / (nu0 * exp(-diff_binding_ratio_surf * 23.0_dp / &
                           Tdust) + evaporation_rates)
         reac(j)%rate=reac(j)%rate * probh2h2
      ENDIF
      
      IF (sumlay.gt.nb_active_lay)  reac(j)%rate = reac(j)%rate*(nb_active_lay/sumlay)**2
   
   
   !-------------------------------------------------------------
   !  (19-20) Photodissociations by UV photons on grain surfaces
   !   - 19 : Photodissociations by UV photons on grain surfaces
   !   - 20 : Photodissociations by UV photons on grain surfaces
   !          (when the gas-phase equivalent of the product is an
   !          ion)
   !-------------------------------------------------------------
   CASE(19:20)

      r1 = reac(j)%r(1)
      
      p1 = reac(j)%p(1)
      p2 = reac(j)%p(2)
   
      reac(j)%rate=reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                   reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                   reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)
                
      ! *** CO self-shielding and isotopologues ***
      IF (spe(r1)%name.EQ.speco) THEN
         call shielding_co(nh2,nco+n13co,1,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,1,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,1,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
      ENDIF
      IF (spe(r1)%name.EQ.'C@O        ') THEN
         call shielding_co(nh2,nco+n13co,2,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,2,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,2,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn) 
      ENDIF
      IF (spe(r1)%name.EQ.'CO@        ') THEN         
         call shielding_co(nh2,nco+n13co,3,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,3,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,3,tetacodn) 
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)  
      
      ENDIF
      IF (spe(r1)%name.EQ.'C@O@       ') THEN
         call shielding_co(nh2,nco+n13co,4,tetaco)
         call shielding_co(nh2_up,nco_up+n13co_up,4,tetacoup)
         call shielding_co(nh2_dn,nco_dn+n13co_dn,4,tetacodn)
         reac(j)%rate= tetaco   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetacoup * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetacodn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)         
      ENDIF
      ! ***
      
      ! *** N2 self-shielding ***
      IF (spe(r1)%name(2:).EQ.spen2) THEN
         reac(j)%rate= tetan2   * reac(j)%alpha *   G0Hab * exp(-reac(j)%gamma * Av) + &
                       tetan2up * reac(j)%alpha * G0upHab * exp(-reac(j)%gamma * Avup) + &
                       tetan2dn * reac(j)%alpha * G0dnHab * exp(-reac(j)%gamma * Avdn)      
      ENDIF
      ! ****
      
      ! Photon attenuation in the mantle with an probability per monolayer Pabs 
      r1 = reac(j)%r(1)
      if(spe(r1)%name(1:1).eq.'K') then
         reac(j)%rate = reac(j)%rate * attenu
      end if
   
      CASE(30:31)
         reac(j)%rate=0.0_dp
            
      END SELECT
      if (reac(j)%rate.lt.1.0e-99_dp) reac(j)%rate = 0.0_dp
   END DO
   
   RETURN
   
   END SUBROUTINE SET_DEPENDANT_RATES
   !########################################################
   
   !########################################################
   SUBROUTINE STICKING_SPECIAL_CASES(J,SUMLAY,acc_prefac)
   !########################################################
   
   IMPLICIT NONE
   
   ! Inputs
   INTEGER, INTENT(IN) :: J ! index of a given reaction
   REAL(KIND=dp), INTENT(IN) :: SUMLAY
   REAL(KIND=dp), INTENT(OUT) :: acc_prefac
   
   ! Local
   INTEGER       :: r1
   REAL(KIND=dp) :: stick, stick_ice, stick_bare
   REAL(KIND=dp) :: new_sigma, da, dummy
   
   r1 = reac(j)%r(1)
   
   IF (spe(r1)%name.eq.speh)  THEN
     stick_bare = (1.0_dp + 2.5_dp*Tgas/25.0_dp)/&
                  (1.0_dp + Tgas/25.0_dp)**2.5_dp
     stick_ice  = (1.0_dp + 2.5_dp*Tgas/52.0_dp)/&
                  (1.0_dp + Tgas/52.0_dp)**2.5_dp
   ENDIF
   
   IF (spe(r1)%name.eq.speh2) THEN
      stick_bare = 0.95_dp * (1.0_dp + 2.5_dp*Tgas/56.0_dp)/&
                   (1.0_dp + Tgas/56.0_dp)**2.5_dp
      stick_ice  = 0.76_dp * (1.0_dp + 2.5_dp*Tgas/87.0_dp)/&
                   (1.0_dp + Tgas/87.0_dp)**2.5_dp
   ENDIF
   
   IF(sumlay.le.1.0_dp) then
     stick = (1.0_dp-sumlay)*stick_bare + sumlay*stick_ice
   ELSE
     stick = stick_ice
   ENDIF
   
   if (is_growth == 0) then
      acc_prefac = ndsigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
   else
      if (first_step_done) then
         da = sumlay*ds
         dummy = 1.0_dp + 2.0_dp*ndusta*da/ndusta2 + ndust*da**2/ndusta2
         new_sigma = ndsigma * dummy
         acc_prefac = new_sigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
      else
         acc_prefac = ndsigma * spe(r1)%stick * sqrt(8.0_dp*k_b / pi / amu / spe(r1)%mass)
      endif
   end if
      
   END SUBROUTINE STICKING_SPECIAL_CASES
   !########################################################
   
   !########################################################
   subroutine shielding_h2(tau,beta)
   !########################################################

      implicit none
      real(kind=dp), intent(in)  :: tau
      real(kind=dp), intent(out) :: beta
      real(kind=dp) :: tmp, eps
      integer :: n
      real(kind=dp) :: v1, dv, b
   
      dv  = 0.13_dp * sqrt(Tgas/2.02_dp)
      dv = 1.0_dp/dv
      b = 9.2e-3_dp * dv
      v1 = 5e2_dp * dv

      if(tau.le.10) then
         beta = 0.0_dp
         eps = 1.0_dp
         n = 0
         do while (eps.gt.1e-6)
            tmp = beta
            beta = beta + ((-1.0_dp)**n * tau**n) / (fact(n)*(n+1.0_dp)**0.5  * pi**(n/2.0_dp))
            eps = abs((beta-tmp)/beta)
            n = n + 1
         end do
      else 
         beta = (log(tau/sqrt(pi))**(-0.5)/tau + (b/tau)**0.5) * &
         erfc((tau*b/pi/v1**2.0)**0.5)
      endif

   end subroutine shielding_h2
   !########################################################

   !########################################################
   real(kind=dp) function fact(n)
   !########################################################
      implicit none
      integer :: i,n
      real(kind=dp) :: p
   
      p = 1
      do i = 1, n
         p = p * i
      end do
      fact = p
   
   end function fact
   !########################################################
   
   !########################################################
   subroutine shielding_co(colH2,colCO,idx,beta)
   !########################################################
      
      USE shielding
      
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: colH2, colCO
      REAL(KIND=dp), INTENT(OUT) :: beta
      INTEGER, intent(in) :: idx
      
      INTEGER :: l,m,i,j
      REAL(KIND=dp) :: logcolH2, logcolCO
      REAL(KIND=dp) :: denom
      REAL(KIND=dp) :: x1, x2
      REAL(KIND=dp) :: xp, yp
      REAL(KIND=dp) :: y1, y2
      REAL(KIND=dp) :: rat1, rat2
      REAL(KIND=dp) :: f11, f21, f12, f22
      
      logcolH2 = log10(colH2)
      logcolCO = log10(colCO)
      
      IF (logcolCO .lt. log10(NCO_visser(1))) logcolCO = log10(NCO_visser(1))
      IF (logcolH2 .lt. log10(NH2_visser(1))) logcolH2 = log10(NH2_visser(1))
      
      IF (logcolCO .gt. log10(NCO_visser(47))) logcolCO = log10(NCO_visser(47))
      IF (logcolH2 .gt. log10(NH2_visser(42))) logcolH2 = log10(NH2_visser(42))
      
      DO l=1,46
         IF((log10(NCO_visser(l)).le.logcolCO).AND.(log10(NCO_visser(l+1)).ge.logcolCO)) THEN
            DO m=1,41
               IF((log10(NH2_visser(m)).le.logcolH2).AND.(log10(NH2_visser(m+1)).ge.logcolH2)) THEN
      
                  ! Bilinear interpolation of the CO shielding function
                  xp = logcolH2
                  yp = logcolCO
      
                  x1 = log10(NH2_visser(m))
                  x2 = log10(NH2_visser(m+1))
      
                  y1 = log10(NCO_visser(l))
                  y2 = log10(NCO_visser(l+1))
      
                  if (idx==1) then   ! 12C16O
                     f11 = shield_12c16o_visser(l,m)
                     f21 = shield_12c16o_visser(l,m+1)
                     f12 = shield_12c16o_visser(l+1,m)
                     f22 = shield_12c16o_visser(l+1,m+1)
                  else if (idx==2) then ! 13C16O
                     f11 = shield_13c16o_visser(l,m)
                     f21 = shield_13c16o_visser(l,m+1)
                     f12 = shield_13c16o_visser(l+1,m)
                     f22 = shield_13c16o_visser(l+1,m+1)
                  else if (idx==3) then ! 12C18O
                     f11 = shield_12c18o_visser(l,m)
                     f21 = shield_12c18o_visser(l,m+1)
                     f12 = shield_12c18o_visser(l+1,m)
                     f22 = shield_12c18o_visser(l+1,m+1)
                  else if (idx==4) then ! 13C18O
                     f11 = shield_13c18o_visser(l,m)
                     f21 = shield_13c18o_visser(l,m+1)
                     f12 = shield_13c18o_visser(l+1,m)
                     f22 = shield_13c18o_visser(l+1,m+1)
                  else
                     print*, 'Error in shielding CO...'
                     stop
                  end if
      
                  denom = x2 - x1
                  rat1 = (x2 - xp) * f11 + (xp - x1) * f21
                  rat1 = rat1 / denom
                  rat2 = (x2 - xp) * f12 + (xp - x1) * f22
                  rat2 = rat2 / denom
      
                  denom = y2 - y1
                  beta = (y2 - yp) * rat1 + (yp - y1) * rat2
                  beta = beta / denom
      
                  return
      
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      
   end subroutine shielding_co
   !########################################################
   
   !########################################################
   subroutine shielding_n2(colH2,colN2,beta)
   !########################################################
   
   ! *** Not very efficient : to be re-written properly
   !     - Linear interpolation in case the shielding factor is not defined
   !     - Bilinear interpolation otherwise
      
      USE shielding
      
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: colH2, colN2
      REAL(KIND=dp), INTENT(OUT) :: beta
      
      INTEGER :: l,m, i,j
      REAL(KIND=dp) :: logcolH2, logcolN2
      REAL(KIND=dp) :: denom
      REAL(KIND=dp) :: x1, x2
      REAL(KIND=dp) :: xp, yp
      REAL(KIND=dp) :: y1, y2
      REAL(KIND=dp) :: rat1, rat2
      REAL(KIND=dp) :: f11, f21, f12, f22
      
      logcolH2 = log10(colH2)
      logcolN2 = log10(colN2)
      
      beta = 1.0_dp
      
      if(logcolH2.lt.14.0) then
        i=1
      elseif(logcolH2.gt.23.0) then
        i=226
      else
        i= int(1 + (logcolH2-14.0)/0.04)
      endif
      
      if(logcolN2.lt.10.0) then
        j=1
      elseif(logcolN2.gt.19.0) then
        j=226
      else
        j= int(1 + (logcolN2-10.0)/0.04)
      endif
      
      beta = thetaN2_file(i,j)
      
      !IF(logcolN2 .lt. log10(nn2_n2(1)).and.logcolH2 .gt. log10(nH2_n2(1))) THEN
      !   IF (logcolH2 .gt. log10(nh2_n2(46))) logcolH2 = log10(nh2_n2(46))
      !   DO m=1,45
      !      IF((log10(NH2_n2(m)).le.logcolH2).AND.(log10(NH2_n2(m+1)).ge.logcolH2)) THEN
      !
      !         xp = logcolH2
      !
      !         x1 = log10(nh2_n2(m))
      !         x2 = log10(nh2_n2(m+1))
      !
      !         f11 = shielding_function_n2_nh1e20(m,1)
      !         f21 = shielding_function_n2_nh1e20(m+1,1)
      !
      !         beta = f11 + (xp - x1) * (f21-f11)/(x2-x1)
      !         
      !         return
      !         
      !      ENDIF
      !   ENDDO
      !ELSE IF(logcolN2 .gt. log10(nn2_n2(1)).and.logcolH2 .lt. log10(nH2_n2(1))) THEN
      !   IF (logcolN2 .gt. log10(nn2_n2(46))) logcolN2 = log10(nn2_n2(46))
      !   DO m=1,45
      !      IF((log10(NN2_n2(m)).le.logcolN2).AND.(log10(NN2_n2(m+1)).ge.logcolN2)) THEN
      !
      !         xp = logcolN2
      !
      !         x1 = log10(nn2_n2(m))
      !         x2 = log10(nn2_n2(m+1))
      !
      !         f11 = shielding_function_n2_nh1e20(1,m)
      !         f21 = shielding_function_n2_nh1e20(1,m+1)
      !
      !         beta = f11 + (xp - x1) * (f21-f11)/(x2-x1)
      !
      !         return
      !
      !      ENDIF
      !   ENDDO
      !ELSE
      !   IF (logcolN2 .gt. log10(nn2_n2(46))) logcolN2 = log10(nn2_n2(46))
      !   IF (logcolH2 .gt. log10(nh2_n2(46))) logcolH2 = log10(nh2_n2(46))
      !   DO l=1,45
      !      IF((log10(nn2_n2(l)).le.logcolN2).AND.(log10(nn2_n2(l+1)).ge.logcolN2)) THEN
      !         DO m=1,45
      !            IF((log10(NH2_n2(m)).le.logcolH2).AND.(log10(NH2_n2(m+1)).ge.logcolH2)) THEN
      !   
      !               ! Bilinear interpolation of the N2 shielding function
      !               xp = logcolH2
      !               yp = logcolN2
      !   
      !               x1 = log10(nh2_n2(m))
      !               x2 = log10(nh2_n2(m+1))
      !   
      !               y1 = log10(nn2_n2(l))
      !               y2 = log10(nn2_n2(l+1))
      !   
      !               f11 = shielding_function_n2_nh1e20(m,l)
      !               f21 = shielding_function_n2_nh1e20(m+1,l)
      !               f12 = shielding_function_n2_nh1e20(m,l+1)
      !               f22 = shielding_function_n2_nh1e20(m+1,l+1)
      !   
      !               denom = x2 - x1
      !               rat1 = (x2 - xp) * f11 + (xp - x1) * f21
      !               rat1 = rat1 / denom
      !               rat2 = (x2 - xp) * f12 + (xp - x1) * f22
      !               rat2 = rat2 / denom
      !   
      !               denom = y2 - y1
      !               beta = (y2 - yp) * rat1 + (yp - y1) * rat2
      !               beta = beta / denom
      !   
      !               return
      !   
      !            ENDIF
      !         ENDDO
      !      ENDIF
      !   ENDDO
      !ENDIF
      
   end subroutine shielding_n2
   !########################################################
   
   !########################################################
   SUBROUTINE GET_CR_RATE()
   !########################################################
      
      IMPLICIT NONE
      integer :: i
      real(kind=dp), dimension(10) :: ck
      real(kind=dp) :: Ntmp
      
      Ntmp = NH2_up
      if (Ntmp.lt.1.0e19) Ntmp = 1.0e19
      if (Ntmp.gt.1.0e27) Ntmp = 1.0e27
      
      ck( 1) = -3.331056497233e6_dp
      ck( 2) =  1.207744586503e6_dp
      ck( 3) = -1.913914106234e5_dp
      ck( 4) =  1.731822350618e4_dp
      ck( 5) = -9.790557206178e2_dp
      ck( 6) =  3.543830893824e1_dp
      ck( 7) = -8.034869454520e-1_dp
      ck( 8) =  1.048808593086e-2_dp
      ck( 9) = -6.188760100997e-5_dp
      ck(10) =  3.122820990797e-8_dp
      
      CRrate = 0.0_dp
      do i = 1, 10
         CRrate = CRrate + ck(i)*log10(Ntmp)**(i-1)
      end do
      CRrate = 10**CRrate
      
   END SUBROUTINE GET_CR_RATE
   !########################################################

END MODULE CHEMISTRY
