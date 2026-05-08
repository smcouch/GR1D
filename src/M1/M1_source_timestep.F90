!-*-f90-*-
subroutine M1_source_timestep_limit(dt_source,dt_ies,dt_energycoupling, &
     dt_realizable,limiter_kind,limiter_zone,limiter_species,limiter_group)

  use GR1D_module
  implicit none

  real*8, intent(out) :: dt_source,dt_ies,dt_energycoupling
  real*8, intent(out) :: dt_realizable
  integer, intent(out) :: limiter_kind,limiter_zone,limiter_species
  integer, intent(out) :: limiter_group

  integer, parameter :: source_none = 0
  integer, parameter :: source_ies = 1
  integer, parameter :: source_energycoupling = 2
  integer, parameter :: source_combined = 3
  real*8, parameter :: huge_dt = 1.0d99

  real*8 :: U(2*number_groups)
  real*8 :: S_ies(2*number_groups),S_energy(2*number_groups)
  real*8 :: S_total(2*number_groups)
  real*8 :: dt_local,dt_realizable_local,oneX_local
  real*8 :: lambda_ies,lambda_energy,lambda_total
  integer :: i,k,g,local_group

  dt_source = huge_dt
  dt_ies = huge_dt
  dt_energycoupling = huge_dt
  dt_realizable = huge_dt
  limiter_kind = source_none
  limiter_zone = 0
  limiter_species = 0
  limiter_group = 0

  if (.not.do_M1) return
  if (.not.M1_source_dt_limiter) return
  if (.not.include_Ielectron_exp.and..not.include_energycoupling_exp) return
  if (number_species_to_evolve.le.0.or.number_groups.le.0) return

  ! M1_explicitterms returns before these source terms at very early times.
  if (time.lt.0.000001d0) return

  if (include_Ielectron_exp) call M1_updateeas
  call M1_reconstruct
  call M1_closure

  do k=ghosts1+1,M1_imaxradii
     do i=1,number_species_to_evolve
        if (GR) then
           oneX_local = X(k)
        else
           oneX_local = 1.0d0
        endif

        do g=1,number_groups
           U(g) = q_M1(k,i,g,1)
           U(g+number_groups) = q_M1(k,i,g,2)
        enddo

        S_total = 0.0d0
        lambda_ies = 0.0d0
        lambda_energy = 0.0d0
        lambda_total = 0.0d0

        if (include_Ielectron_exp) then
           call M1_ies_source_rate_bound(k,i,U,S_ies,lambda_ies,.true.)
           call M1_source_dt_bound_from_lambda(S_ies,lambda_ies,U, &
                oneX_local,dt_local,dt_realizable_local,local_group)
           if (dt_local.lt.dt_ies) dt_ies = dt_local
           S_total = S_total + S_ies
           lambda_total = lambda_total + lambda_ies
        endif

        if (include_energycoupling_exp) then
           call M1_energycoupling_source_rate_bound(k,i,U,S_energy, &
                lambda_energy)
           call M1_source_dt_bound_from_lambda(S_energy,lambda_energy,U, &
                oneX_local,dt_local,dt_realizable_local,local_group)
           if (dt_local.lt.dt_energycoupling) dt_energycoupling = dt_local
           S_total = S_total + S_energy
           lambda_total = lambda_total + lambda_energy
        endif

        call M1_source_dt_bound_from_lambda(S_total,lambda_total,U, &
             oneX_local,dt_local,dt_realizable_local,local_group)
        if (dt_realizable_local.lt.dt_realizable) dt_realizable = &
             dt_realizable_local

        if (dt_local.lt.dt_source) then
           dt_source = dt_local
           if (include_Ielectron_exp.and.include_energycoupling_exp) then
              limiter_kind = source_combined
           else if (include_Ielectron_exp) then
              limiter_kind = source_ies
           else
              limiter_kind = source_energycoupling
           endif
           limiter_zone = k
           limiter_species = i
           limiter_group = local_group
        endif
     enddo
  enddo

end subroutine M1_source_timestep_limit

subroutine M1_source_dt_bound_from_lambda(S,lambda,U,oneX_local, &
     dt_bound,dt_realizable_bound,bound_group)

  use GR1D_module
  implicit none

  real*8, intent(in) :: S(2*number_groups)
  real*8, intent(in) :: lambda
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(in) :: oneX_local
  real*8, intent(out) :: dt_bound
  real*8, intent(out) :: dt_realizable_bound
  integer, intent(out) :: bound_group

  integer :: lambda_group

  lambda_group = 1
  if (lambda.gt.0.0d0) then
     call M1_source_dt_bound_core(S,lambda,lambda_group,U,oneX_local, &
          dt_bound,dt_realizable_bound,bound_group)
  else
     call M1_source_dt_bound_core(S,0.0d0,lambda_group,U,oneX_local, &
          dt_bound,dt_realizable_bound,bound_group)
  endif

end subroutine M1_source_dt_bound_from_lambda

subroutine M1_source_dt_bound_core(S,lambda,lambda_group,U,oneX_local, &
     dt_bound,dt_realizable_bound,bound_group)

  use GR1D_module
  implicit none

  real*8, intent(in) :: S(2*number_groups)
  real*8, intent(in) :: lambda
  integer, intent(in) :: lambda_group
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(in) :: oneX_local
  real*8, intent(out) :: dt_bound
  real*8, intent(out) :: dt_realizable_bound
  integer, intent(out) :: bound_group

  real*8, parameter :: huge_dt = 1.0d99
  real*8 :: theta,pi_pos,scale,floor_value
  real*8 :: rplus,rminus,drplusdt,drminusdt,realizable_rate
  real*8 :: species_energy,active_floor,margin_floor
  real*8 :: candidate
  integer :: row,g,theta_group,pi_group,realizable_group

  floor_value = max(M1_source_dt_floor,tiny)
  theta = 0.0d0
  pi_pos = 0.0d0
  realizable_rate = 0.0d0
  bound_group = 1
  theta_group = 1
  pi_group = 1
  realizable_group = 1
  species_energy = 0.0d0

  do g=1,number_groups
     species_energy = species_energy + max(U(g),0.0d0)
  enddo
  active_floor = max(0.0d0,M1_source_realizable_floor_abs, &
       M1_source_realizable_floor_rel*species_energy)

  do row=1,2*number_groups
     if (row.le.number_groups) then
        scale = max(abs(U(row)),floor_value)
     else
        g = row-number_groups
        scale = max(abs(U(row)),1.0d-3*abs(U(g)),floor_value)
     endif
     if (abs(S(row))/scale.gt.theta) then
        theta = abs(S(row))/scale
        if (row.le.number_groups) then
           theta_group = row
        else
           theta_group = g
        endif
     endif
  enddo

  do g=1,number_groups
     if (S(g).lt.0.0d0) then
        if (-S(g)/max(U(g),floor_value).gt.pi_pos) then
           pi_pos = -S(g)/max(U(g),floor_value)
           pi_group = g
        endif
     endif

     if (U(g).gt.active_floor) then
        margin_floor = max(floor_value, &
             M1_source_realizable_margin*max(U(g),active_floor))
        rplus = U(g) - U(g+number_groups)/oneX_local
        rminus = U(g) + U(g+number_groups)/oneX_local
        drplusdt = S(g) - S(g+number_groups)/oneX_local
        drminusdt = S(g) + S(g+number_groups)/oneX_local

        if (drplusdt.lt.0.0d0) then
           if (-drplusdt/max(rplus,margin_floor).gt.realizable_rate) then
              realizable_rate = -drplusdt/max(rplus,margin_floor)
              realizable_group = g
           endif
        endif
        if (drminusdt.lt.0.0d0) then
           if (-drminusdt/max(rminus,margin_floor).gt.realizable_rate) then
              realizable_rate = -drminusdt/max(rminus,margin_floor)
              realizable_group = g
           endif
        endif
     endif
  enddo

  dt_bound = huge_dt
  dt_realizable_bound = huge_dt
  if (lambda.gt.0.0d0) then
     candidate = M1_source_cfl_linear/lambda
     if (candidate.lt.dt_bound) then
        dt_bound = candidate
        bound_group = lambda_group
     endif
  endif
  if (theta.gt.0.0d0) then
     candidate = M1_source_cfl_fraction/theta
     if (candidate.lt.dt_bound) then
        dt_bound = candidate
        bound_group = theta_group
     endif
  endif
  if (pi_pos.gt.0.0d0) then
     candidate = M1_source_cfl_positive/pi_pos
     if (candidate.lt.dt_bound) then
        dt_bound = candidate
        bound_group = pi_group
     endif
  endif
  if (realizable_rate.gt.0.0d0) then
     candidate = M1_source_cfl_realizable/realizable_rate
     dt_realizable_bound = candidate
     if (candidate.lt.dt_bound) then
        dt_bound = candidate
        bound_group = realizable_group
     endif
  endif

end subroutine M1_source_dt_bound_core

subroutine M1_source_timestep_cache_reset

  use GR1D_module
  implicit none

  M1_source_dt_cache_valid = .false.
  M1_source_dt_cache = 1.0d99
  M1_source_dt_ies_cache = 1.0d99
  M1_source_dt_energycoupling_cache = 1.0d99
  M1_source_dt_realizable_cache = 1.0d99
  M1_source_limiter_kind_cache = 0
  M1_source_limiter_zone_cache = 0
  M1_source_limiter_species_cache = 0
  M1_source_limiter_group_cache = 0

end subroutine M1_source_timestep_cache_reset

subroutine M1_source_timestep_cache_get(dt_source,dt_ies,dt_energycoupling, &
     dt_realizable,limiter_kind,limiter_zone,limiter_species,limiter_group, &
     cache_hit)

  use GR1D_module
  implicit none

  real*8, intent(out) :: dt_source,dt_ies,dt_energycoupling
  real*8, intent(out) :: dt_realizable
  integer, intent(out) :: limiter_kind,limiter_zone,limiter_species
  integer, intent(out) :: limiter_group
  logical, intent(out) :: cache_hit

  cache_hit = M1_source_dt_cache_valid
  if (.not.cache_hit) then
     dt_source = 1.0d99
     dt_ies = 1.0d99
     dt_energycoupling = 1.0d99
     dt_realizable = 1.0d99
     limiter_kind = 0
     limiter_zone = 0
     limiter_species = 0
     limiter_group = 0
     return
  endif

  dt_source = M1_source_dt_cache
  dt_ies = M1_source_dt_ies_cache
  dt_energycoupling = M1_source_dt_energycoupling_cache
  dt_realizable = M1_source_dt_realizable_cache
  limiter_kind = M1_source_limiter_kind_cache
  limiter_zone = M1_source_limiter_zone_cache
  limiter_species = M1_source_limiter_species_cache
  limiter_group = M1_source_limiter_group_cache

end subroutine M1_source_timestep_cache_get

subroutine M1_source_timestep_cache_update(dts,implicit_factor)

  use GR1D_module
  implicit none

  real*8, intent(in) :: dts,implicit_factor

  integer, parameter :: source_none = 0
  integer, parameter :: source_ies = 1
  integer, parameter :: source_energycoupling = 2
  integer, parameter :: source_combined = 3
  real*8, parameter :: huge_dt = 1.0d99

  real*8 :: U(2*number_groups)
  real*8 :: S_ies(2*number_groups),S_energy(2*number_groups)
  real*8 :: S_total(2*number_groups)
  real*8 :: dt_source,dt_ies,dt_energycoupling,dt_realizable
  real*8 :: dt_local,dt_realizable_local,oneX_local,dt_factor
  real*8 :: lambda_ies,lambda_energy,lambda_total,safety
  integer :: i,k,g,local_group
  integer :: limiter_kind,limiter_zone,limiter_species,limiter_group

  call M1_source_timestep_cache_reset

  if (.not.do_M1) return
  if (.not.M1_source_dt_limiter) return
  if (.not.include_Ielectron_exp.and..not.include_energycoupling_exp) return
  if (number_species_to_evolve.le.0.or.number_groups.le.0) return
  if (time.lt.0.000001d0) return
  dt_factor = dts*implicit_factor
  if (dt_factor.le.0.0d0) return

  dt_source = huge_dt
  dt_ies = huge_dt
  dt_energycoupling = huge_dt
  dt_realizable = huge_dt
  limiter_kind = source_none
  limiter_zone = 0
  limiter_species = 0
  limiter_group = 0

  do k=ghosts1+1,M1_imaxradii
     do i=1,number_species_to_evolve
        if (GR) then
           oneX_local = X(k)
        else
           oneX_local = 1.0d0
        endif

        do g=1,number_groups
           U(g) = q_M1(k,i,g,1)
           U(g+number_groups) = q_M1(k,i,g,2)
        enddo

        S_total = 0.0d0
        lambda_total = 0.0d0

        if (include_Ielectron_exp) then
           do g=1,number_groups
              S_ies(g) = flux_M1_scatter(k,i,g,1)/dt_factor
              S_ies(g+number_groups) = flux_M1_scatter(k,i,g,2)/dt_factor
           enddo
           lambda_ies = M1_source_lambda_ies(k,i)
           call M1_source_dt_bound_from_lambda(S_ies,lambda_ies,U, &
                oneX_local,dt_local,dt_realizable_local,local_group)
           if (dt_local.lt.dt_ies) dt_ies = dt_local
           S_total = S_total + S_ies
           lambda_total = lambda_total + lambda_ies
        endif

        if (include_energycoupling_exp) then
           do g=1,number_groups
              S_energy(g) = -flux_M1_energy(k,i,g,1)/dt_factor
              S_energy(g+number_groups) = -flux_M1_energy(k,i,g,2)/dt_factor
           enddo
           lambda_energy = M1_source_lambda_energycoupling(k,i)
           call M1_source_dt_bound_from_lambda(S_energy,lambda_energy,U, &
                oneX_local,dt_local,dt_realizable_local,local_group)
           if (dt_local.lt.dt_energycoupling) dt_energycoupling = dt_local
           S_total = S_total + S_energy
           lambda_total = lambda_total + lambda_energy
        endif

        call M1_source_dt_bound_from_lambda(S_total,lambda_total,U, &
             oneX_local,dt_local,dt_realizable_local,local_group)
        if (dt_realizable_local.lt.dt_realizable) dt_realizable = &
             dt_realizable_local

        if (dt_local.lt.dt_source) then
           dt_source = dt_local
           if (include_Ielectron_exp.and.include_energycoupling_exp) then
              limiter_kind = source_combined
           else if (include_Ielectron_exp) then
              limiter_kind = source_ies
           else
              limiter_kind = source_energycoupling
           endif
           limiter_zone = k
           limiter_species = i
           limiter_group = local_group
        endif
     enddo
  enddo

  safety = max(tiny,min(1.0d0,M1_source_dt_cache_safety))
  M1_source_dt_cache = safety*dt_source
  M1_source_dt_ies_cache = safety*dt_ies
  M1_source_dt_energycoupling_cache = safety*dt_energycoupling
  M1_source_dt_realizable_cache = safety*dt_realizable
  M1_source_limiter_kind_cache = limiter_kind
  M1_source_limiter_zone_cache = limiter_zone
  M1_source_limiter_species_cache = limiter_species
  M1_source_limiter_group_cache = limiter_group
  M1_source_dt_cache_valid = .true.

end subroutine M1_source_timestep_cache_update
