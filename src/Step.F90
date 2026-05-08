! -*-f90-*-
subroutine Step(dts)

  use GR1D_module
  use ye_of_rho
  use nulibtable
#if HAVE_LEAK_ROS
  use leakage_rosswog
#endif
  implicit none
  
  real*8 dts !time step that is passed in

  real*8 beta_rk,alpha_rk,cooling_rk_1,cooling_rk_2

  integer rkindex,i,m,gi,itc,j,k
  integer(kind=4) :: eosflag,keyerr,keytemp
  real*8 eosdummy(15)
  real*8 tempeps1(n1), tempeps2(n1)
  real*8 m1_eps_before_source(n1)
  real*8 epsin0
  
  logical nan,inf
  logical m1_rk_coupled,m1_update_diagnostics

  ! Is it time to turn on turbulence?
  if (do_turbulence) then
     if (tpb_for_turbulence .lt. 0.0d0) then
        activate_turbulence = .true.
     else if (bounce) then
        if (time .gt. t_bounce + tpb_for_turbulence) then
           activate_turbulence = .true.
        else 
           activate_turbulence = .false.  
        endif
     else 
        activate_turbulence = .false.
     endif
  endif

  !If the shock reaches the outer boundary you can safely assume that the SN exploded
  if ( (shock_radius .ge. x1(n1-ghosts1-1)) .or. &
       (shock_radius/length_gf .ge. 15000.0d5) .and. .not. explosion_reached) then
     write(*,*) "Explosion! :-)"
     open(unit=666,file=trim(adjustl(outdir))//"/explosion",status="unknown")
     write(666,*) 1
     close(666)
     explosion_reached = .true.
  endif

  !calculate v_turb for this time step
  if (activate_turbulence) then
     call Brunt_Vaisala(dts)
  endif

  !set up conserved variables
  call prim2con

  !GR, do not need sqrt_gamma
  if (GR) then
     q_hat(:,:) = q(:,:)
     q_hat_old(:,:) = q(:,:)
  else
     qold(:,:) = q(:,:)
     do m=1,n_cons
        q_hat(:,m) = sqrt_gamma(:) * q(:,m)
     enddo
     q_hat_old(:,:) = q_hat(:,:)
  endif

  if(iorder_hydro.eq.2) then
     alpha_rk = 2.0d0
     beta_rk = 1.0d0
     cooling_rk_1 = 0.5d0
     cooling_rk_2 = 0.5d0
  else if(iorder_hydro.eq.3) then
     alpha_rk = 4.0d0
     beta_rk = 3.0d0
     cooling_rk_1 = 0.25d0
     cooling_rk_2 = 1.0d0/12.0d0
  else
     alpha_rk = 1.0d0
     beta_rk = 1.0d0
     cooling_rk_1 = 1.0d0
     cooling_rk_2 = 0.0d0 !not used
  endif

  denergyloss(:) = 0.0d0
  m1_rk_coupled = do_M1.and.M1_integrate_with_hydro_rk
  m1_update_diagnostics = .false.

  if (M1_integrate_with_hydro_rk) then
     if (.not.do_M1) stop "M1_integrate_with_hydro_rk requires do_M1"
     if (.not.do_hydro) stop "M1_integrate_with_hydro_rk requires do_hydro"
     if (.not.GR) stop "M1_integrate_with_hydro_rk currently requires GR"
  endif

  if (m1_rk_coupled) then
     q_M1_prev = q_M1
     M1_hydro_source(:,:) = 0.0d0
     total_net_heating = 0.0d0
     total_net_deintdt = 0.0d0
     total_mass_gain = 0.0d0
     igain(1) = -1
     gain_radius = 0.0d0
  endif

  if(.not.do_hydro) goto 123

  do rkindex=1,iorder_hydro

     rkstep = rkindex

     gravsource(:,:) = 0.0d0
     flux_diff(:,:) = 0.0d0
     coolingsource(:,:) = 0.0d0
     presssource(:,:) = 0.0d0

     call reconstruct 
     call boundaries(0,0)

     if (activate_turbulence) then
        turb_source(:,:) = 0.0d0
        call turb_diff_terms
        call turbulence_sources
     endif
     
     if(flux_type .eq. "HLLE") then
        call flux_differences_hlle
     else
        stop "Sorry. Don't have the Riemann solver you want..."
     endif

     if(gravity_active.and.(geometry.eq.2)) then
       call gravity
     endif

     if(hydro_formulation.eq."conservative_p_source".and.(geometry.eq.2)) then
        call press_sources 
     endif
  
     if (GR.and.gravity_active) then
        if (do_nupress) then
           call nu_press_sources
        endif
     endif

#if HAVE_LEAK_ROS
     if (do_leak_ros.and.bounce) then
        call leak_rosswog
     endif
#endif

     if (m1_rk_coupled) then
        m1_update_diagnostics = rkindex.eq.iorder_hydro
        if (m1_update_diagnostics) m1_eps_before_source(:) = eps(:)
        call M1_euler_advance(dts)
        call M1_fill_hydro_source(dts,m1_update_diagnostics)
     else
        M1_hydro_source(:,:) = 0.0d0
     endif

     if (rkindex .eq. 1 ) then
        do i=ghosts1,n1-1
           
           ! rho,D
           q_hat(i,1) = q_hat_old(i,1) + dts * ( - flux_diff(i,1) )
           
           ! rho*v, S
           q_hat(i,2) = q_hat_old(i,2) + dts * ( - flux_diff(i,2) &
                + gravsource(i,2) + presssource(i,2) + &
                coolingsource(i,2) + M1_hydro_source(i,2))
           
           ! energy, tau
           q_hat(i,3) = q_hat_old(i,3) + dts * ( - flux_diff(i,3) &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3) + M1_hydro_source(i,3))
           denergyloss(i) = cooling_rk_1*coolingsource(i,3)*dts
           
           ! ye
           q_hat(i,4) = q_hat_old(i,4) + dts * ( - flux_diff(i,4) &
                + coolingsource(i,4) + M1_hydro_source(i,4))
           
        enddo
        
        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = q_hat_old(i,5) + dts * ( - flux_diff(i,5) &
                   + presssource(i,5) )
           enddo
        endif

        if(activate_turbulence) then
           do i=ghosts1,n1-1
              ! add source of turbulent energy to eps
              q_hat(i,3) = q_hat(i,3) + dts * turb_source(i,3)
                  q_hat(i,6) = q_hat_old(i,6) + dts * ( - flux_diff(i,6) &
                        + turb_source(i,6))
				  ! Sometimes v_turb**2 can become negative, which is not physical	
                  if (q_hat(i,6) .lt. 0.d0) then
                      q_hat(i,6) = 0.0d0
                  endif
           enddo
        endif
        
     elseif (rkindex .eq. 2 ) then
        do i=ghosts1,n1-1
           q_hat(i,1) = ( beta_rk * q_hat_old(i,1) + q_hat(i,1)  &
                + dts * ( - flux_diff(i,1) ) ) / alpha_rk
           
           q_hat(i,2) = ( beta_rk * q_hat_old(i,2) + q_hat(i,2)  &
                + dts * ( - flux_diff(i,2)    &
                + gravsource(i,2) & 
                + presssource(i,2) + coolingsource(i,2) &
                + M1_hydro_source(i,2)) ) / alpha_rk
           
           q_hat(i,3) = ( beta_rk * q_hat_old(i,3) + q_hat(i,3)  &
                + dts * ( - flux_diff(i,3)    &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3) + M1_hydro_source(i,3) ) ) / alpha_rk
           denergyloss(i) = denergyloss(i) + cooling_rk_2 * &
                coolingsource(i,3)*dts
           
           q_hat(i,4) = ( beta_rk * q_hat_old(i,4)           &
                + q_hat(i,4)                            &
                + dts * ( - flux_diff(i,4)    &
                + coolingsource(i,4) + M1_hydro_source(i,4) ) ) / alpha_rk
           
        enddo

        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = ( beta_rk * q_hat_old(i,5)           &
                   + q_hat(i,5)                            &
                   + dts * ( - flux_diff(i,5)    &
                   + presssource(i,5) ) ) / alpha_rk
           enddo
        endif

        if(activate_turbulence) then
           do i=ghosts1,n1-1
              q_hat(i,3) = q_hat(i,3) + dts            &
                   * turb_source(i,3) / alpha_rk              
                  q_hat(i,6) = (beta_rk * q_hat_old(i,6) + q_hat(i,6) &
                       + dts * (- flux_diff(i,6) + turb_source(i,6)))/alpha_rk
				  ! Sometimes v_turb**2 can become negative, which is not physical	
                  if (q_hat(i,6) .lt. 0.d0) then
                      q_hat(i,6) = 0.0d0
                  endif
           enddo
        endif 
       
     elseif (rkindex .eq. 3) then
        do i=ghosts1,n1-1
           q_hat(i,1) = ( q_hat_old(i,1) + 2.0d0*q_hat(i,1)  &
                + 2.0d0*dts * ( - flux_diff(i,1)) ) / 3.0d0
           
           q_hat(i,2) = ( q_hat_old(i,2) + 2.0d0*q_hat(i,2)  &
                + 2.0d0*dts * ( - flux_diff(i,2)    &
                + gravsource(i,2) & 
                + presssource(i,2) + coolingsource(i,2) &
                + M1_hydro_source(i,2)) ) / 3.0d0
           
           q_hat(i,3) = ( q_hat_old(i,3) + 2.0d0*q_hat(i,3)  &
                + 2.0d0*dts * ( - flux_diff(i,3)    &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3) + M1_hydro_source(i,3) ) ) / 3.0d0
           denergyloss(i) = denergyloss(i) + 2.0d0/3.0d0 * &
                coolingsource(i,3)*dts
           
           q_hat(i,4) = ( q_hat_old(i,4) + 2.0d0*q_hat(i,4)  &
                + 2.0d0*dts * ( - flux_diff(i,4)    &
                + coolingsource(i,4) + M1_hydro_source(i,4) ) ) / 3.0d0
        enddo
        
        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = ( q_hat_old(i,5) + 2.0d0*q_hat(i,5)  &
                   + 2.0d0*dts * ( - flux_diff(i,5)    &
                   + presssource(i,5) ) ) / 3.0d0
           enddo
        endif
        
        if(activate_turbulence) then
           do i=ghosts1,n1-1
              q_hat(i,3) = q_hat(i,3) + 2.0d0 * dts  &
                   * turb_source(i,3) / 3.0d0
                  q_hat(i,6) = (q_hat_old(i,6) + 2.0d0*q_hat(i,6) &
                       + 2.0d0*dts*( - flux_diff(i,6) &
                       + turb_source(i,6)))/3.0d0
				  ! Sometimes v_turb**2 can become negative, which is not physical	
                  if (q_hat(i,6) .lt. 0.d0) then
                      q_hat(i,6) = 0.0d0
                  endif
           enddo
        endif 
     else 
        stop 'Only iorder_hydro = 1, 2, and 3 implemented!'
     endif

     if (m1_rk_coupled) call M1_blend_rk_stage(rkindex)
     
     do m=1,n_cons
        if (GR) then
           q(:,m) = q_hat(:,m)
	else 
           q(:,m) = q_hat(:,m) / sqrt_gamma(:)
        endif
     enddo

     if (GR.and.gravity_active) then
        !find mgrav & X
        call con2GR
     endif
     
     !reconstruct primatives
     call con2prim

     ! eos update, eps fixed, find temp,entropy,cs2 etc.
     do i=ghosts1+1,n1-ghosts1
        keyerr = 0
        keytemp = 0
        tempeps1(i) = eps(i)
        call eos_full(i,rho(i),temp(i),ye(i),eps(i),press(i),pressth(i), & 
             entropy(i), &
             cs2(i), & 
             eosdummy(2),&
             eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
             eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
             eosdummy(11),eosdummy(12),eosdummy(13),nuchem(i), &
             keytemp,keyerr,eoskey,eos_rf_prec)
        tempeps2(i) = eps(i)
        if(keyerr.ne.0) then
           ! -> Issues with the EOS, this can happen around bounce
           !    and is due to very large temperature gradients
           !    in the bouncing inner core in adiabatic collapse
           !    for which the EOS was not really designed. The
           !    problems seen here should not show up for leakage/ye_of_rho
           !    runs.
           write(6,*) "############################################"
           write(6,*) "EOS PROBLEM in Step:"
           write(6,*) "timestep number: ",nt
           write(6,"(i4,1P10E15.6)") i,x1(i),rho(i)/rho_gf,temp(i),eps(i)/eps_gf,ye(i)
           write(6,*) "keyerr: ",keyerr
           call flush(6)
           if(.not.fake_neutrinos.and..not.do_leak_ros) then
              ! let's pump in a tiny bit of energy < 0.01*epsin_orig
              itc = 0
              epsin0 = eps(i)
              do while(keyerr.ne.0.and.itc.lt.10) 
                 itc = itc + 1
                 eps(i) = eps(i) + epsin0 * 1.0001d0
                 call eos_full(i,rho(i),temp(i),ye(i),eps(i),press(i),pressth(i), & 
                      entropy(i), &
                      cs2(i), & 
                      eosdummy(2),&
                      eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                      eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                      eosdummy(11),eosdummy(12),eosdummy(13),nuchem(i), &
                      keytemp,keyerr,eoskey,eos_rf_prec)
              enddo
              write(6,*) itc,keyerr
              write(6,*) "############################################"
              if(keyerr.ne.0) then
                 stop "problem in reconstruct: Step, could not be fixed."
              endif
           else
              stop "problem in reconstruct: Step"
           endif
        endif
     enddo
     
     !GR gravity updates that rely on primitive variables
     if (GR.and.gravity_active) then
        call GR_alp
        call GR_boundaries
     elseif (GR) then
        if (geometry.eq.2) then
           !GR Sedov, reflective boundaries on inside
           call GR_boundaries
        else 
           !planer, shocktube
           gi = 0
           do i=ghosts1,1,-1
              gi=gi+1
              v(i) = 0.0d0
              vp(i) = 0.0d0
              vm(i) = 0.0d0   
              W(i) = 1.0d0
           enddo
           do i=n1-ghosts1,n1
              gi=n1-ghosts1-1
              v(i) = v(gi)
              W(i) = W(gi)
           enddo
        endif
     endif

     call mass_interior

     if (do_nupress) then
        call neutrino_pressure
     endif

    call boundaries(0,0)

    if (m1_rk_coupled) then
       call M1_refresh_radiation_diagnostics
       if (m1_update_diagnostics) then
          call M1_update_deint_diagnostic(m1_eps_before_source,dts)
          dyedt_hydro(:) = (ye(:) - ye_prev(:))/dts - dyedt_neutrino(:)
       endif
    endif
    
 enddo

123 continue

 !do operator split here
 !M1
 if (do_M1.and.(.not.M1_integrate_with_hydro_rk)) then

    dyedt_hydro(:) = (ye(:) - ye_prev(:))/dts
    ye_prev(:) = ye(:)

    qold = q
    call M1_euler_advance(dts)

    if (do_hydro.or.(M1_testcase_number.eq.1.and.time.gt.0.0012d0)) then
       call M1_conservativeupdate(dts)
    endif

    dyedt_neutrino(:) = (ye(:) - ye_prev(:))/dts

 endif

 !ye of rho prescription
 if(bounce) then
    if (do_ye_of_rho.and.(time.lt.t_bounce+0.005d0)) then
       call adjust_ye
       !things slightly change so redo all the variables
       call prim2con
       call boundaries(0,0)
       if(GR) then
          call con2GR
          call GR_alp
          call GR_boundaries
       endif
    endif
 else
    if (do_ye_of_rho) then
       call adjust_ye
       !things slightly change so redo all the variables
       call prim2con
       call boundaries(0,0)
       if(GR) then
          call con2GR
          call GR_alp
          call GR_boundaries
       endif
    endif
 endif

contains

  subroutine M1_euler_advance(stage_dts)

    implicit none
    real*8, intent(in) :: stage_dts

    real*8 :: stage_implicit_factor
    integer :: ii,jj,kk

    !we need to find the new plus/minus states, GR (alp,X) boundaries are done
    call reconstruct
    call boundaries(0,0)
    !If Newtonian, need to set v to v1 for velocity terms.
    if (.not.GR) then
       v = v1
       vp = v1p
       vm = v1m
    endif

    stage_implicit_factor = 1.0d0
    q_M1_old = q_M1

    !reset source term
    M1_matter_source = 0.0d0

    !update interaction rates
    call M1_updateeas

    !reconstruct energy and flux in space and energy
    call M1_reconstruct

    !update closure variables
    call M1_closure

    !get explicit fluxes at the current stage state
    call M1_explicitterms(stage_dts,stage_implicit_factor)

    if (M1_do_backwardfix.eq.1) then
       do jj=1,number_groups
          do ii=1,number_species
             do kk=ghosts1+1,M1_imaxradii
                if ((eas(kk,ii,jj,2)+eas(kk,ii,jj,3))* &
                     (x1i(kk+1)-x1i(kk)).lt.0.01d0) then
                   B_M1(kk,ii,jj,1:2) = -flux_M1(kk,ii,jj,1:2)*0.5d0
                else
                   B_M1(kk,ii,jj,1:2) = -flux_M1(kk,ii,jj,1:2)
                endif
             enddo
          enddo
       enddo
    else
       B_M1 = -flux_M1
    endif
    C_M1 = -flux_M1_energy
    D_M1 = flux_M1_scatter

    !do implicit source solve and calculate matter source terms
    call M1_implicitstep(stage_dts,stage_implicit_factor)

    !code for backward euler explicit flux fix
    if (M1_do_backwardfix.eq.1) then
       call M1_reconstruct
       call M1_closure
       call M1_explicitterms(stage_dts,stage_implicit_factor)

       do jj=1,number_groups
          do ii=1,number_species
             do kk=ghosts1+1,M1_imaxradii
                if ((eas(kk,ii,jj,2)+eas(kk,ii,jj,3))* &
                     (x1i(kk+1)-x1i(kk)).lt.0.01d0) then
                   q_M1(kk,ii,jj,1:2) = q_M1(kk,ii,jj,1:2) - &
                        B_M1(kk,ii,jj,1:2) - flux_M1(kk,ii,jj,1:2)
                   if (q_M1(kk,ii,jj,1).lt.0.0d0) then
                      write(*,*) kk,ii,jj,q_M1(kk,ii,jj,1)
                      stop "negative en after flux correct"
                   endif

                   if (abs(q_M1(kk,ii,jj,2)/X(kk)).gt.q_M1(kk,ii,jj,1)) then
                      q_M1(kk,ii,jj,2) = X(kk)*q_M1(kk,ii,jj,2) / &
                           abs((1.0d0+1.0d-10)*q_M1(kk,ii,jj,2)/ &
                           q_M1(kk,ii,jj,1))
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine M1_euler_advance

  subroutine M1_fill_hydro_source(stage_dts,update_diagnostics)

    implicit none
    real*8, intent(in) :: stage_dts
    logical, intent(in) :: update_diagnostics

    real*8 :: dDye,dtau,oneX,maxye
    integer :: kk,maxyeloc
    logical :: passfluxtest

    M1_hydro_source(:,:) = 0.0d0
    depsdt(:) = 0.0d0
    dyedt(:) = 0.0d0
    dyedt_neutrino(:) = 0.0d0
    maxye = 0.0d0
    maxyeloc = 0

    if (update_diagnostics) then
       total_net_heating = 0.0d0
       total_net_deintdt = 0.0d0
       total_mass_gain = 0.0d0
       igain(1) = -1
       gain_radius = 0.0d0
    endif

    do kk=ghosts1+1,M1_imaxradii
       oneX = X(kk)

       passfluxtest = rho(kk)/rho_gf.lt.3.0d10
       if ((M1_matter_source(kk,3).gt.0.0d0).and. &
            (entropy(kk).gt.6.0d0).and.passfluxtest) then
          M1_matter_source(kk,3) = M1_matter_source(kk,3)*M1_heat_fac
       endif

       M1_hydro_source(kk,2) = 4.0d0*pi*M1_matter_source(kk,2)
       M1_hydro_source(kk,3) = 4.0d0*pi*M1_matter_source(kk,3)
       M1_hydro_source(kk,4) = 4.0d0*pi*M1_matter_source(kk,4)*oneX* &
            (amu_cgs*mass_gf)

       depsdt(kk) = M1_matter_source(kk,3)/rho(kk)*4.0d0*pi/eps_gf*time_gf
       dyedt(kk) = M1_hydro_source(kk,4)/q(kk,1)*time_gf
       dyedt_neutrino(kk) = M1_hydro_source(kk,4)/q(kk,1)

       if (update_diagnostics) then
          dDye = stage_dts*M1_hydro_source(kk,4)
          dtau = stage_dts*M1_hydro_source(kk,3)

          if ((dtau.gt.0.0d0).and.(entropy(kk).gt.6.0d0).and. &
               passfluxtest) then
             total_net_heating = total_net_heating + &
                  dtau*X(kk)*volume(kk)/(energy_gf*stage_dts/time_gf)
             total_mass_gain = total_mass_gain + volume(kk)*rho(kk)
             if (igain(1).lt.0) igain(1) = kk
          endif

          total_energy_absorped = total_energy_absorped + &
               dtau*volume(kk)/energy_gf/(stage_dts/time_gf)

          if (abs(dDye/q(kk,1)).gt.abs(maxye)) then
             maxye = dDye/q(kk,1)
             maxyeloc = kk
          endif
       endif
    enddo

    if (update_diagnostics.and.abs(maxye).gt.0.02d0) then
       dt_reduction_factor = dt_reduction_factor*0.9d0
       write(*,*) "Warning, ye seems unstable, reducing time step to compensate", &
            stage_dts,dt_reduction_factor,maxyeloc
    endif

    M1_matter_source(:,:) = 0.0d0

  end subroutine M1_fill_hydro_source

  subroutine M1_blend_rk_stage(stage_index)

    implicit none
    integer, intent(in) :: stage_index

    if (stage_index.eq.1) then
       return
    else if (stage_index.eq.2) then
       q_M1 = (beta_rk*q_M1_prev + q_M1)/alpha_rk
    else if (stage_index.eq.3) then
       q_M1 = (q_M1_prev + 2.0d0*q_M1)/3.0d0
    endif

  end subroutine M1_blend_rk_stage

  subroutine M1_refresh_radiation_diagnostics

    implicit none

    real*8 :: alp2,invalp2,invX,invX2,X2,W2,v2,oneW,onev,oneX
    real*8 :: sign_one,oneM1en,oneM1flux,oneeddy
    integer :: ii,jj,kk

    call M1_reconstruct
    call M1_closure

    press_nu(:) = 0.0d0
    energy_nu(:) = 0.0d0
    mom_nu(:) = 0.0d0
    ynu(:) = 0.0d0

    do kk=ghosts1+1,M1_imaxradii
       if (GR) then
          alp2 = alp(kk)*alp(kk)
          invalp2 = 1.0d0/alp2
          X2 = X(kk)*X(kk)
          oneX = X(kk)
          invX = 1.0d0/X(kk)
          invX2 = 1.0d0/X2
       else
          alp2 = 1.0d0
          invalp2 = 1.0d0
          X2 = 1.0d0
          oneX = 1.0d0
          invX = 1.0d0
          invX2 = 1.0d0
       endif

       if (v_order.eq.-1) then
          if (GR) then
             W2 = W(kk)**2
             oneW = W(kk)
             v2 = v(kk)**2
             onev = v(kk)
          else
             W2 = 1.0d0/(1.0d0-v1(kk)**2)
             oneW = sqrt(W2)
             v2 = v1(kk)**2
             onev = v1(kk)
          endif
       else if (v_order.eq.0) then
          W2 = 1.0d0
          oneW = 1.0d0
          v2 = 0.0d0
          onev = 0.0d0
       else
          stop "add in vorder"
       endif

       do ii=1,number_species_to_evolve
          do jj=1,number_groups
             oneM1en = q_M1(kk,ii,jj,1)
             oneM1flux = q_M1(kk,ii,jj,2)
             oneeddy = q_M1(kk,ii,jj,3)

             q_M1_fluid(kk,ii,jj,1) = oneM1en*W2 - &
                  2.0d0*oneM1flux*W2*onev*invX + &
                  oneeddy*oneM1en*W2*v2*invX2

             q_M1_fluid(kk,ii,jj,2) = -(oneM1en*oneW - &
                  oneM1flux*oneW*onev/oneX)*W2*onev*invX + &
                  W2*oneW*oneM1flux*invX2 - &
                  oneeddy*oneM1en*invX2**2*W2*oneW*onev*oneX

             if (ii.eq.1) sign_one = 1.0d0
             if (ii.eq.2) sign_one = -1.0d0
             if (ii.gt.2) sign_one = 0.0d0

             ynu(kk) = ynu(kk) + sign_one*q_M1_fluid(kk,ii,jj,1)* &
                  4.0d0*pi/rho(kk)*nulibtable_inv_energies(jj)* &
                  (amu_cgs*mass_gf)
             press_nu(kk) = press_nu(kk) + oneeddy*oneM1en*4.0d0*pi* &
                  invX2**2
             energy_nu(kk) = energy_nu(kk) + oneM1en*4.0d0*pi
             mom_nu(kk) = mom_nu(kk) + oneM1flux*4.0d0*pi
          enddo
       enddo
    enddo

  end subroutine M1_refresh_radiation_diagnostics

  subroutine M1_update_deint_diagnostic(eps_before,stage_dts)

    implicit none
    real*8, intent(in) :: eps_before(n1)
    real*8, intent(in) :: stage_dts

    integer :: kk

    total_net_deintdt = 0.0d0
    do kk=ghosts1+1,n1-ghosts1
       if (eps(kk).gt.eps_before(kk)) then
          if ((rho(kk)/rho_gf.lt.3.0d10).and.(entropy(kk).gt.6.0d0)) then
             total_net_deintdt = total_net_deintdt + &
                  volume(kk)*rho(kk)*(eps(kk)-eps_before(kk))/ &
                  energy_gf/(stage_dts/time_gf)
          endif
       endif
    enddo

  end subroutine M1_update_deint_diagnostic

end subroutine Step
