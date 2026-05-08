!-*-f90-*-
subroutine M1_source_timestep_limit(dt_source,dt_ies,dt_energycoupling, &
     dt_realizable,limiter_kind,limiter_zone,limiter_species,limiter_group)

  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies,nulibtable_energies, &
       nulibtable_etop,nulibtable_ebottom,nulibtable_logenergies, &
       nulibtable_logetop
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
  real*8 :: A_ies(2*number_groups,2*number_groups)
  real*8 :: A_energy(2*number_groups,2*number_groups)
  real*8 :: A_total(2*number_groups,2*number_groups)
  real*8 :: dt_local,dt_realizable_local,oneX_local
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
        A_total = 0.0d0

        if (include_Ielectron_exp) then
           call finite_difference_source(source_ies,k,i,U,S_ies,A_ies)
           call source_dt_bound(S_ies,A_ies,U,oneX_local,dt_local, &
                dt_realizable_local,local_group)
           if (dt_local.lt.dt_ies) dt_ies = dt_local
           S_total = S_total + S_ies
           A_total = A_total + A_ies
        endif

        if (include_energycoupling_exp) then
           call finite_difference_source(source_energycoupling,k,i,U, &
                S_energy,A_energy)
           call source_dt_bound(S_energy,A_energy,U,oneX_local,dt_local, &
                dt_realizable_local,local_group)
           if (dt_local.lt.dt_energycoupling) dt_energycoupling = dt_local
           S_total = S_total + S_energy
           A_total = A_total + A_energy
        endif

        if (include_Ielectron_exp.and.include_energycoupling_exp) then
           call source_dt_bound(S_total,A_total,U,oneX_local,dt_local, &
                dt_realizable_local,local_group)
        else if (include_Ielectron_exp) then
           call source_dt_bound(S_total,A_total,U,oneX_local,dt_local, &
                dt_realizable_local,local_group)
        else
           call source_dt_bound(S_total,A_total,U,oneX_local,dt_local, &
                dt_realizable_local,local_group)
        endif
        if (dt_realizable_local.lt.dt_realizable) dt_realizable = dt_realizable_local

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

contains

  subroutine finite_difference_source(source_kind,zone_index,species_index,U,S0,A)

    implicit none
    integer, intent(in) :: source_kind,zone_index,species_index
    real*8, intent(in) :: U(2*number_groups)
    real*8, intent(out) :: S0(2*number_groups)
    real*8, intent(out) :: A(2*number_groups,2*number_groups)

    real*8 :: Uplus(2*number_groups),Uminus(2*number_groups)
    real*8 :: Splus(2*number_groups),Sminus(2*number_groups)
    real*8 :: h,scale,floor_value
    integer :: col

    call evaluate_source(source_kind,zone_index,species_index,U,S0)
    A = 0.0d0
    floor_value = max(M1_source_dt_floor,tiny)

    do col=1,2*number_groups
       scale = max(abs(U(col)),floor_value)
       h = 1.0d-6*scale
       Uplus = U
       Uminus = U
       Uplus(col) = Uplus(col) + h

       if (col.le.number_groups.and.Uminus(col)-h.le.floor_value) then
          call evaluate_source(source_kind,zone_index,species_index,Uplus,Splus)
          A(:,col) = (Splus(:)-S0(:))/h
       else
          Uminus(col) = Uminus(col) - h
          call evaluate_source(source_kind,zone_index,species_index,Uplus,Splus)
          call evaluate_source(source_kind,zone_index,species_index,Uminus,Sminus)
          A(:,col) = (Splus(:)-Sminus(:))/(2.0d0*h)
       endif
    enddo

  end subroutine finite_difference_source

  subroutine evaluate_source(source_kind,zone_index,species_index,U,S)

    implicit none
    integer, intent(in) :: source_kind,zone_index,species_index
    real*8, intent(in) :: U(2*number_groups)
    real*8, intent(out) :: S(2*number_groups)

    if (source_kind.eq.source_ies) then
       call evaluate_ies_source(zone_index,species_index,U,S)
    else if (source_kind.eq.source_energycoupling) then
       call evaluate_energycoupling_source(zone_index,species_index,U,S)
    else
       S = 0.0d0
    endif

  end subroutine evaluate_source

  subroutine source_dt_bound(S,A,U,oneX_local,dt_bound,dt_realizable_bound, &
       bound_group)

    implicit none
    real*8, intent(in) :: S(2*number_groups)
    real*8, intent(in) :: A(2*number_groups,2*number_groups)
    real*8, intent(in) :: U(2*number_groups)
    real*8, intent(in) :: oneX_local
    real*8, intent(out) :: dt_bound
    real*8, intent(out) :: dt_realizable_bound
    integer, intent(out) :: bound_group

    real*8 :: lambda,theta,pi_pos,row_sum,scale,floor_value
    real*8 :: rplus,rminus,drplusdt,drminusdt,realizable_rate
    real*8 :: species_energy,active_floor,margin_floor
    real*8 :: candidate
    integer :: row,col,g,lambda_group,theta_group,pi_group,realizable_group

    floor_value = max(M1_source_dt_floor,tiny)
    lambda = 0.0d0
    theta = 0.0d0
    pi_pos = 0.0d0
    realizable_rate = 0.0d0
    bound_group = 1
    lambda_group = 1
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
       row_sum = 0.0d0
       do col=1,2*number_groups
          row_sum = row_sum + abs(A(row,col))
       enddo
       if (row_sum.gt.lambda) then
          lambda = row_sum
          if (row.le.number_groups) then
             lambda_group = row
          else
             lambda_group = row-number_groups
          endif
       endif
    enddo

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

  end subroutine source_dt_bound

  subroutine evaluate_energycoupling_source(zone_index,species_index,U,S)

    implicit none
    integer, intent(in) :: zone_index,species_index
    real*8, intent(in) :: U(2*number_groups)
    real*8, intent(out) :: S(2*number_groups)

    real*8 :: div_v,dvdt_local
    real*8 :: h,dmdr,dmdt,dXdr,dXdt,Kdownrr,dWdt,dWdr,dWvuprdt,dWvuprdr
    real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6)
    real*8 :: heattermff_NL(6),velocity_coeffs(6,2)
    real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
    real*8 :: M1en_energy_fluid(number_groups),M1eddy_energy(number_groups)
    real*8 :: M1pff_energy(number_groups),M1qrrr_energy(number_groups)
    real*8 :: M1qffr_energy(number_groups),littlefactors(number_groups,6)
    real*8 :: velocity(number_groups,2),M1flux_energy_interface(number_groups,2)
    real*8 :: logdistro(number_groups),loginterface_distroj(number_groups)
    real*8 :: interface_distroj(number_groups),xi(number_groups)
    real*8 :: FL(number_groups,2),FR(number_groups,2)
    real*8 :: M1flux_diff_energy(number_groups,2)
    real*8 :: em,ep,enext,temp_term,width
    real*8 :: W2,v2,oneW,onev
    integer :: j

    S = 0.0d0

    if (GR) then
       if (zone_index.eq.ghosts1+1) then
          div_v = v(zone_index)/x1(zone_index)
       else
          div_v = (v(zone_index+1)-v(zone_index-1))/ &
               (x1(zone_index+1)-x1(zone_index-1))
       endif
       if (nt.eq.0) then
          dvdt_local = 0.0d0
       else
          dvdt_local = M1_source_dvdt(zone_index)
       endif
    else
       if (zone_index.eq.ghosts1+1) then
          div_v = v1(zone_index)/x1(zone_index)
       else
          div_v = (v1(zone_index+1)-v1(zone_index-1))/ &
               (x1(zone_index+1)-x1(zone_index-1))
       endif
       if (nt.eq.0) then
          dvdt_local = 0.0d0
       else
          dvdt_local = M1_source_dvdt(zone_index)
       endif
    endif

    if (v_order.eq.-1) then
       if (GR) then
          W2 = W(zone_index)**2
          v2 = v(zone_index)**2
          oneW = W(zone_index)
          onev = v(zone_index)
       else
          W2 = 1.0d0/(1.0d0-v1(zone_index)**2)
          v2 = v1(zone_index)**2
          oneW = sqrt(W2)
          onev = v1(zone_index)
       endif
    else if (v_order.eq.0) then
       W2 = 1.0d0
       oneW = 1.0d0
       v2 = 0.0d0
       onev = 0.0d0
       div_v = 0.0d0
       dvdt_local = 0.0d0
    else
       stop "M1_source_timestep: unsupported v_order"
    endif

    if (GR) then
       h = 1.0d0+eps(zone_index)+press(zone_index)/rho(zone_index)
       dmdr = 4.0d0*pi*x1(zone_index)**2* &
            (rho(zone_index)*h*W(zone_index)**2-press(zone_index))
       dmdt = -4.0d0*pi*x1(zone_index)**2*alp(zone_index)* &
            rho(zone_index)*h*W(zone_index)**2*v(zone_index)/X(zone_index)
       dXdr = X(zone_index)**3*(dmdr/x1(zone_index)- &
            mgrav(zone_index)/x1(zone_index)**2)
       dXdt = X(zone_index)**3*dmdt/x1(zone_index)
       Kdownrr = -X(zone_index)*dXdt/alp(zone_index)
       dWdt = W2*oneW*onev*dvdt_local
       dWdr = W2*oneW*onev*div_v
       dWvuprdt = onev/X(zone_index)*dWdt - &
            oneW*onev/X(zone_index)**2*dXdt + oneW/X(zone_index)*dvdt_local
       dWvuprdr = onev/X(zone_index)*dWdr - &
            oneW*onev/X(zone_index)**2*dXdr + oneW/X(zone_index)*div_v
    else
       dmdr = 0.0d0
       dmdt = 0.0d0
       dXdr = 0.0d0
       dXdt = 0.0d0
       Kdownrr = 0.0d0
       dWdt = W2*oneW*onev*dvdt_local
       dWdr = W2*oneW*onev*div_v
       dWvuprdt = onev*dWdt + oneW*dvdt_local
       dWvuprdr = onev*dWdr + oneW*div_v
    endif

    Z = 0.0d0
    Yupr = 0.0d0
    Xuprr = 0.0d0
    Xupff = 0.0d0
    heatterm_NL = 0.0d0
    heattermff_NL = 0.0d0

    if (GR) then
       Z(1) = 1.0d0/oneW
       Z(2) = onev/(X(zone_index)*oneW)
       Z(3) = v2/(X(zone_index)**2*oneW)
       Z(5) = v2*onev/X(zone_index)**3
       Yupr(2) = 1.0d0/(X(zone_index)**2*oneW)
       Yupr(3) = onev/(X(zone_index)**3*oneW)
       Yupr(5) = v2/X(zone_index)**4
       Xuprr(3) = 1.0d0/(X(zone_index)**4*oneW)
       Xuprr(5) = onev/X(zone_index)**5
       Xupff(4) = 1.0d0/(x1(zone_index)**2*oneW)
       Xupff(6) = onev/(X(zone_index)*x1(zone_index)**2)
       heatterm_NL(5) = 1.0d0
       heattermff_NL(6) = 1.0d0
    else
       Z(1) = 1.0d0/oneW
       Z(2) = onev/oneW
       Z(3) = v2/oneW
       Z(5) = v2*onev
       Yupr(2) = 1.0d0/oneW
       Yupr(3) = onev/oneW
       Yupr(5) = v2
       Xuprr(3) = 1.0d0/oneW
       Xuprr(5) = onev
       Xupff(4) = 1.0d0/(x1(zone_index)**2*oneW)
       Xupff(6) = onev/x1(zone_index)**2
       heatterm_NL(5) = 1.0d0
       heattermff_NL(6) = 1.0d0
    endif

    if (GR) then
       velocity_coeffs(:,1) = alp(zone_index)*(oneW*((Z(:)*onev/ &
            X(zone_index)-Yupr(:))*dphidr(zone_index) - Xuprr(:)*onev*dXdr - &
            2.0d0*Xupff(:)*onev/X(zone_index)*x1(zone_index) + &
            Xuprr(:)*Kdownrr) + Z(:)/alp(zone_index)*dWdt + Yupr(:)*dWdr - &
            X(zone_index)**2*Yupr(:)/alp(zone_index)*dWvuprdt - &
            X(zone_index)**2*Xuprr(:)*dWvuprdr)

       velocity_coeffs(:,2) = alp(zone_index)*(oneW*((Yupr(:)*onev* &
            X(zone_index)-Xuprr(:)*X(zone_index)**2)*dphidr(zone_index) - &
            heatterm_NL(:)/X(zone_index)**4*onev*dXdr - &
            2.0d0*heattermff_NL(:)/x1(zone_index)**2*onev/ &
            X(zone_index)*x1(zone_index) + heatterm_NL(:)/ &
            X(zone_index)**4*Kdownrr) + Yupr(:)*X(zone_index)**2/ &
            alp(zone_index)*dWdt + Xuprr(:)*X(zone_index)**2*dWdr - &
            Xuprr(:)*X(zone_index)**4/alp(zone_index)*dWvuprdt - &
            heatterm_NL(:)/X(zone_index)**2*dWvuprdr)
    else
       if (do_effectivepotential) then
          velocity_coeffs(:,1) = alp(zone_index)*(oneW*(-2.0d0* &
               Xupff(:)*onev*x1(zone_index)) + Z(:)/alp(zone_index)* &
               dWdt + Yupr(:)*dWdr - Yupr(:)*dWvuprdt/alp(zone_index) - &
               Xuprr(:)*dWvuprdr)
          velocity_coeffs(2,1) = velocity_coeffs(2,1) - alp(zone_index)* &
               dphidr(zone_index)

          velocity_coeffs(:,2) = alp(zone_index)*(oneW*(-2.0d0* &
               heattermff_NL(:)/x1(zone_index)**2*onev*x1(zone_index)) + &
               Yupr(:)*dWdt/alp(zone_index) + Xuprr(:)*dWdr - &
               Xuprr(:)*dWvuprdt/alp(zone_index) - heatterm_NL(:)*dWvuprdr)
          velocity_coeffs(3,2) = velocity_coeffs(3,2) - alp(zone_index)* &
               dphidr(zone_index)
       else
          velocity_coeffs(:,1) = oneW*(-2.0d0*Xupff(:)*onev*x1(zone_index)) + &
               Z(:)/alp(zone_index)*dWdt + Yupr(:)*dWdr - Yupr(:)*dWvuprdt - &
               Xuprr(:)*dWvuprdr

          velocity_coeffs(:,2) = oneW*(-2.0d0*heattermff_NL(:)/ &
               x1(zone_index)**2*onev*x1(zone_index)) + Yupr(:)*dWdt + &
               Xuprr(:)*dWdr - Xuprr(:)*dWvuprdt - heatterm_NL(:)*dWvuprdr
       endif
    endif

    do j=1,number_groups
       width = nulibtable_etop(j)-nulibtable_ebottom(j)
       M1en_energy(j) = U(j)/width
       M1flux_energy(j) = U(j+number_groups)/width
       M1eddy_energy(j) = q_M1(zone_index,species_index,j,3)
       M1pff_energy(j) = q_M1_extra(zone_index,species_index,j,1)
       M1qrrr_energy(j) = q_M1_extra(zone_index,species_index,j,2)/width
       M1qffr_energy(j) = q_M1_extra(zone_index,species_index,j,3)/width
       M1en_energy_fluid(j) = q_M1_fluid(zone_index,species_index,j,1)

       littlefactors(j,1) = M1en_energy(j)
       littlefactors(j,2) = M1flux_energy(j)
       littlefactors(j,3) = M1eddy_energy(j)*M1en_energy(j)
       littlefactors(j,4) = M1pff_energy(j)*M1en_energy(j)
       littlefactors(j,5) = M1qrrr_energy(j)
       littlefactors(j,6) = M1qffr_energy(j)
       velocity(j,1) = sum(littlefactors(j,:)*velocity_coeffs(:,1))
       velocity(j,2) = sum(littlefactors(j,:)*velocity_coeffs(:,2))
    enddo

    M1flux_energy_interface = 0.0d0
    logdistro = log(M1en_energy_fluid(:)*M1_moment_to_distro(:))
    do j=1,number_groups-1
       loginterface_distroj(j) = logdistro(j) + &
            (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
            (logdistro(j+1)-logdistro(j))/ &
            (nulibtable_logenergies(j+1)-nulibtable_logenergies(j))
       interface_distroj(j) = exp(loginterface_distroj(j))
    enddo
    j=number_groups
    loginterface_distroj(j) = logdistro(j) + &
         (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
         (logdistro(j)-logdistro(j-1))/ &
         (nulibtable_logenergies(j)-nulibtable_logenergies(j-1))
    interface_distroj(j) = exp(loginterface_distroj(j))

    xi(1) = 1.0d0
    do j=2,number_groups
       xi(j) = interface_distroj(j)/(interface_distroj(j)+ &
            interface_distroj(j-1))
    enddo

    temp_term = (nulibtable_etop(1)-nulibtable_ebottom(1))/ &
         (1.0d0-nulibtable_energies(1)/nulibtable_energies(2))*xi(1)
    FL(1,1) = velocity(1,1)*temp_term
    FL(1,2) = velocity(1,2)*temp_term
    M1flux_energy_interface(1,1) = M1flux_energy_interface(1,1) + &
         FL(1,1)/nulibtable_etop(1)
    M1flux_energy_interface(1,2) = M1flux_energy_interface(1,2) + &
         FL(1,2)/nulibtable_etop(1)

    enext = nulibtable_energies(number_groups)**2/ &
         nulibtable_energies(number_groups-1)
    do j=2,number_groups
       if (j.eq.number_groups) then
          temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
               (1.0d0-nulibtable_energies(j)/enext)*xi(j)
       else
          temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
               (1.0d0-nulibtable_energies(j)/nulibtable_energies(j+1))*xi(j)
       endif
       FL(j,1) = velocity(j,1)*temp_term
       FL(j,2) = velocity(j,2)*temp_term

       temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
            (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0)* &
            (1.0d0-xi(j))
       FR(j,1) = velocity(j,1)*temp_term
       FR(j,2) = velocity(j,2)*temp_term

       M1flux_energy_interface(j-1,1) = M1flux_energy_interface(j-1,1) + &
            FR(j,1)/nulibtable_ebottom(j)
       M1flux_energy_interface(j-1,2) = M1flux_energy_interface(j-1,2) + &
            FR(j,2)/nulibtable_ebottom(j)
       M1flux_energy_interface(j,1) = M1flux_energy_interface(j,1) + &
            FL(j,1)/nulibtable_etop(j)
       M1flux_energy_interface(j,2) = M1flux_energy_interface(j,2) + &
            FL(j,2)/nulibtable_etop(j)
    enddo

    j=1
    ep = nulibtable_etop(j)
    M1flux_diff_energy(j,1) = ep*M1flux_energy_interface(j,1)
    M1flux_diff_energy(j,2) = ep*M1flux_energy_interface(j,2)
    do j=2,number_groups
       em = nulibtable_ebottom(j)
       ep = nulibtable_etop(j)
       M1flux_diff_energy(j,1) = ep*M1flux_energy_interface(j,1) - &
            em*M1flux_energy_interface(j-1,1)
       M1flux_diff_energy(j,2) = ep*M1flux_energy_interface(j,2) - &
            em*M1flux_energy_interface(j-1,2)
    enddo

    S(1:number_groups) = -M1flux_diff_energy(:,1)
    S(number_groups+1:2*number_groups) = -M1flux_diff_energy(:,2)

  end subroutine evaluate_energycoupling_source

  subroutine evaluate_ies_source(zone_index,species_index,U,S)

    implicit none
    integer, intent(in) :: zone_index,species_index
    real*8, intent(in) :: U(2*number_groups)
    real*8, intent(out) :: S(2*number_groups)

    real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
    real*8 :: M1eddy_energy(number_groups),M1pff_energy(number_groups)
    real*8 :: M1qrrr_energy(number_groups),M1qffr_energy(number_groups)
    real*8 :: M1chi_energy(number_groups)
    real*8 :: local_M(2,2),local_J(number_groups),local_H(number_groups,2)
    real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
    real*8 :: local_Ltilde(number_groups,2,2)
    real*8 :: JoverE(number_groups),JoverF(number_groups)
    real*8 :: HoverE(number_groups,2),HoverF(number_groups,2)
    real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)
    real*8 :: ies_sourceterms(2*number_groups)
    real*8 :: h,invalp,invalp2,invX,invX2,X2,alp2,W2,v2,invr
    real*8 :: invr2,onev,oneW,onealp,oneX
    real*8 :: local_u(2),local_uup(2),local_littleh(2,2)
    real*8 :: local_littlehupup(2,2)
    real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in,ies_temp
    real*8 :: species_factor
    integer :: j,j_prime,ii,jj

    S = 0.0d0

    if (species_index.eq.3.and.number_species.eq.3) then
       species_factor = 4.0d0
    else if (species_index.eq.3) then
       stop "M1_source_timestep: unsupported species factor"
    else
       species_factor = 1.0d0
    endif

    M1en_energy = U(1:number_groups)/species_factor
    M1flux_energy = U(number_groups+1:2*number_groups)/species_factor
    M1eddy_energy = q_M1(zone_index,species_index,:,3)
    M1pff_energy = q_M1_extra(zone_index,species_index,:,1)
    M1qrrr_energy = q_M1_extra(zone_index,species_index,:,2)/species_factor
    M1qffr_energy = q_M1_extra(zone_index,species_index,:,3)/species_factor
    M1chi_energy = q_M1_extra(zone_index,species_index,:,4)

    local_M = 0.0d0
    local_J = 0.0d0
    local_H = 0.0d0
    local_L = 0.0d0
    local_Hdown = 0.0d0
    local_Ltilde = 0.0d0
    ies_sourceterms = 0.0d0

    h = 1.0d0+eps(zone_index)+press(zone_index)/rho(zone_index)

    if (GR) then
       invalp = 1.0d0/alp(zone_index)
       invalp2 = invalp**2
       alp2 = alp(zone_index)**2
       onealp = alp(zone_index)
       invX = 1.0d0/X(zone_index)
       invX2 = invX**2
       X2 = X(zone_index)**2
       oneX = X(zone_index)
    else
       if (do_effectivepotential) then
          invalp = 1.0d0/alp(zone_index)
          invalp2 = invalp**2
          alp2 = alp(zone_index)**2
          onealp = alp(zone_index)
       else
          invalp = 1.0d0
          invalp2 = 1.0d0
          alp2 = 1.0d0
          onealp = 1.0d0
       endif
       invX = 1.0d0
       invX2 = 1.0d0
       X2 = 1.0d0
       oneX = 1.0d0
    endif

    if (v_order.eq.-1) then
       if (GR) then
          W2 = W(zone_index)**2
          oneW = W(zone_index)
          v2 = v(zone_index)**2
          onev = v(zone_index)
       else
          W2 = 1.0d0/(1.0d0-v1(zone_index)**2)
          oneW = sqrt(W2)
          v2 = v1(zone_index)**2
          onev = v1(zone_index)
       endif
    else if (v_order.eq.0) then
       W2 = 0.0d0
       oneW = 1.0d0
       v2 = 0.0d0
       onev = 0.0d0
    else
       stop "M1_source_timestep: unsupported v_order"
    endif

    invr = 1.0d0/x1(zone_index)
    invr2 = invr*invr

    local_u(1) = -oneW*onealp
    local_u(2) = oneW*onev*oneX
    local_uup(1) = -local_u(1)*invalp2
    local_uup(2) = local_u(2)*invX2

    local_littleh(1,1) = -v2*W2
    local_littleh(2,1) = local_u(2)*invX2*local_u(1)
    local_littleh(1,2) = -local_u(1)*invalp2*local_u(2)
    local_littleh(2,2) = W2

    local_littlehupup(1,1) = -invalp2 + local_uup(1)*local_uup(1)
    local_littlehupup(1,2) = local_uup(1)*local_uup(2)
    local_littlehupup(2,1) = local_uup(2)*local_uup(1)
    local_littlehupup(2,2) = invX2 + local_uup(2)*local_uup(2)

    do j=1,number_groups
       local_M(1,1) = M1en_energy(j)*invalp2
       local_M(1,2) = M1flux_energy(j)*invX2*invalp
       local_M(2,1) = local_M(1,2)
       local_M(2,2) = M1eddy_energy(j)*M1en_energy(j)*invX2**2

       do ii=1,2
          do jj=1,2
             local_J(j) = local_J(j) + local_M(ii,jj)*local_u(ii)*local_u(jj)
             local_H(j,1) = local_H(j,1) - &
                  local_M(ii,jj)*local_u(ii)*local_littleh(1,jj)
             local_H(j,2) = local_H(j,2) - &
                  local_M(ii,jj)*local_u(ii)*local_littleh(2,jj)
          enddo
       enddo

       do ii=1,2
          do jj=1,2
             local_L(j,1,1) = local_L(j,1,1) + &
                  local_M(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj)* &
                  (1.5d0*M1chi_energy(j)-0.5d0) + &
                  local_littlehupup(1,1)*local_J(j)/3.0d0* &
                  (1.5d0-1.5d0*M1chi_energy(j))
             local_L(j,1,2) = local_L(j,1,2) + &
                  local_M(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj)* &
                  (1.5d0*M1chi_energy(j)-0.5d0) + &
                  local_littlehupup(1,2)*local_J(j)/3.0d0* &
                  (1.5d0-1.5d0*M1chi_energy(j))
             local_L(j,2,1) = local_L(j,2,1) + &
                  local_M(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj)* &
                  (1.5d0*M1chi_energy(j)-0.5d0) + &
                  local_littlehupup(2,1)*local_J(j)/3.0d0* &
                  (1.5d0-1.5d0*M1chi_energy(j))
             local_L(j,2,2) = local_L(j,2,2) + &
                  local_M(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj)* &
                  (1.5d0*M1chi_energy(j)-0.5d0) + &
                  local_littlehupup(2,2)*local_J(j)/3.0d0* &
                  (1.5d0-1.5d0*M1chi_energy(j))
          enddo
       enddo

       local_Hdown(j,1) = local_H(j,1)*local_littleh(1,1)*(-alp2) + &
            local_H(j,2)*local_littleh(1,2)*(-alp2)
       local_Hdown(j,2) = local_H(j,1)*local_littleh(2,1)*X2 + &
            local_H(j,2)*local_littleh(2,2)*X2

       local_Ltilde(j,1,1) = -local_L(j,1,1)*alp2 - &
            local_J(j)*local_littleh(1,1)*onethird
       local_Ltilde(j,1,2) = local_L(j,1,2)*X2 - &
            local_J(j)*local_littleh(1,2)*onethird
       local_Ltilde(j,2,1) = -local_L(j,2,1)*alp2 - &
            local_J(j)*local_littleh(2,1)*onethird
       local_Ltilde(j,2,2) = local_L(j,2,2)*X2 - &
            local_J(j)*local_littleh(2,2)*onethird

       JoverE(j) = W2*(1.0d0+M1eddy_energy(j)*v2*invX2)
       JoverF(j) = -2.0d0*W2*onev*invX

       HoverE(j,1) = oneW*local_littleh(1,1)*invalp - &
            local_littleh(1,2)*local_u(2)*M1eddy_energy(j)*invX2**2
       HoverE(j,2) = oneW*local_littleh(2,1)*invalp - &
            local_littleh(2,2)*local_u(2)*M1eddy_energy(j)*invX2**2
       HoverF(j,1) = oneW*local_littleh(1,2)*invX2 - &
            local_littleh(1,1)*local_u(2)*invX2*invalp
       HoverF(j,2) = oneW*local_littleh(2,2)*invX2 - &
            local_littleh(2,1)*local_u(2)*invX2*invalp

       LoverE(j,1,1) = local_littleh(1,1)**2*invalp2 + &
            M1eddy_energy(j)*invX2**2*local_littleh(1,2)**2
       LoverE(j,1,2) = invalp2*local_littleh(1,1)*local_littleh(2,1) + &
            M1eddy_energy(j)*invX2**2*local_littleh(1,2)*local_littleh(2,2)
       LoverE(j,2,1) = LoverE(j,1,2)
       LoverE(j,2,2) = invalp2*local_littleh(2,1)**2 + &
            M1eddy_energy(j)*invX2**2*local_littleh(2,2)**2
       LoverF(j,1,1) = 2.0d0*invX2*invalp* &
            local_littleh(1,2)*local_littleh(1,1)
       LoverF(j,1,2) = invX2*invalp*(local_littleh(1,1)* &
            local_littleh(2,2)+local_littleh(1,2)*local_littleh(2,1))
       LoverF(j,2,1) = LoverF(j,1,2)
       LoverF(j,2,2) = 2.0d0*invX2*invalp* &
            local_littleh(2,1)*local_littleh(2,2)
    enddo

    do j=1,number_groups
       do j_prime=1,number_groups
          nucubed = M1_moment_to_distro_inverse(j)
          nucubedprime = M1_moment_to_distro_inverse(j_prime)

          R0out = 0.5d0*ies(zone_index,species_index,j,j_prime,1)
          R1out = 1.5d0*ies(zone_index,species_index,j,j_prime,2)
          if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
          R0in = 0.5d0*ies(zone_index,species_index,j_prime,j,1)
          R1in = 1.5d0*ies(zone_index,species_index,j_prime,j,2)

          if (rho(zone_index).gt.5.0d12*rho_gf) then
             R0out = R0out*(5.0d12*rho_gf/rho(zone_index))**1.5d0
             R1out = R1out*(5.0d12*rho_gf/rho(zone_index))**1.5d0
             if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
             R0in = R0in*(5.0d12*rho_gf/rho(zone_index))**1.5d0
             R1in = R1in*(5.0d12*rho_gf/rho(zone_index))**1.5d0
          endif

          ies_temp = species_factor*alp2*4.0d0*pi*( &
               ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
               local_H(j,1))*R0in*JoverE(j_prime) + &
               R1in*(HoverE(j_prime,1)*((nucubed-local_J(j))* &
               local_littleh(1,1)*onethird - local_u(1)*(-invalp2)* &
               local_Hdown(j,1) - local_Ltilde(j,1,1)) + &
               HoverE(j_prime,2)*((nucubed-local_J(j))* &
               local_littleh(1,2)*onethird - local_u(1)*(-invalp2)* &
               local_Hdown(j,2) - local_Ltilde(j,1,2))))* &
               nulibtable_inv_energies(j_prime)
          ies_sourceterms(j) = ies_sourceterms(j) + &
               ies_temp*M1en_energy(j_prime)

          ies_temp = species_factor*alp2*4.0d0*pi*( &
               ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
               local_H(j,1))*R0in*JoverF(j_prime) + &
               R1in*(HoverF(j_prime,1)*((nucubed-local_J(j))* &
               local_littleh(1,1)*onethird - local_u(1)*(-invalp2)* &
               local_Hdown(j,1) - local_Ltilde(j,1,1)) + &
               HoverF(j_prime,2)*((nucubed-local_J(j))* &
               local_littleh(1,2)*onethird - local_u(1)*(-invalp2)* &
               local_Hdown(j,2) - local_Ltilde(j,1,2))))* &
               nulibtable_inv_energies(j_prime)
          ies_sourceterms(j) = ies_sourceterms(j) + &
               ies_temp*M1flux_energy(j_prime)

          ies_temp = species_factor*alp2*4.0d0*pi*( &
               -R0out*(nucubedprime-local_J(j_prime))* &
               (JoverE(j)*local_u(1)*(-invalp2)+HoverE(j,1)) + &
               R1out*(local_Hdown(j_prime,1)*(HoverE(j,1)* &
               local_u(1)*(-invalp2)+LoverE(j,1,1)) + &
               local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(1)* &
               (-invalp2)+LoverE(j,1,2))))*nulibtable_inv_energies(j_prime)
          ies_sourceterms(j) = ies_sourceterms(j) + ies_temp*M1en_energy(j)

          ies_temp = species_factor*alp2*4.0d0*pi*( &
               -R0out*(nucubedprime-local_J(j_prime))* &
               (JoverF(j)*local_u(1)*(-invalp2)+HoverF(j,1)) + &
               R1out*(local_Hdown(j_prime,1)*(HoverF(j,1)* &
               local_u(1)*(-invalp2)+LoverF(j,1,1)) + &
               local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(1)* &
               (-invalp2)+LoverF(j,1,2))))*nulibtable_inv_energies(j_prime)
          ies_sourceterms(j) = ies_sourceterms(j) + ies_temp*M1flux_energy(j)

          ies_temp = species_factor*onealp*X2*4.0d0*pi*( &
               ((nucubed-local_J(j))*(local_u(2)*invX2) - &
               local_H(j,2))*R0in*JoverE(j_prime) + &
               R1in*(HoverE(j_prime,1)*((nucubed-local_J(j))* &
               local_littleh(2,1)*onethird - local_u(2)*X2* &
               local_Hdown(j,1) - local_Ltilde(j,2,1)) + &
               HoverE(j_prime,2)*((nucubed-local_J(j))* &
               local_littleh(2,2)*onethird - local_u(2)*X2* &
               local_Hdown(j,2) - local_Ltilde(j,2,2))))* &
               nulibtable_inv_energies(j_prime)
          ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
               ies_temp*M1en_energy(j_prime)

          ies_temp = species_factor*onealp*X2*4.0d0*pi*( &
               ((nucubed-local_J(j))*(local_u(2)*invX2) - &
               local_H(j,2))*R0in*JoverF(j_prime) + &
               R1in*(HoverF(j_prime,1)*((nucubed-local_J(j))* &
               local_littleh(2,1)*onethird - local_u(2)*X2* &
               local_Hdown(j,1) - local_Ltilde(j,2,1)) + &
               HoverF(j_prime,2)*((nucubed-local_J(j))* &
               local_littleh(2,2)*onethird - local_u(2)*X2* &
               local_Hdown(j,2) - local_Ltilde(j,2,2))))* &
               nulibtable_inv_energies(j_prime)
          ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
               ies_temp*M1flux_energy(j_prime)

          ies_temp = species_factor*onealp*X2*4.0d0*pi*( &
               -R0out*(nucubedprime-local_J(j_prime))* &
               (JoverE(j)*local_u(2)*invX2+HoverE(j,2)) + &
               R1out*(local_Hdown(j_prime,1)*(HoverE(j,1)* &
               local_u(2)*invX2+LoverE(j,2,1)) + &
               local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(2)* &
               invX2+LoverE(j,2,2))))*nulibtable_inv_energies(j_prime)
          ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
               ies_temp*M1en_energy(j)

          ies_temp = species_factor*onealp*X2*4.0d0*pi*( &
               -R0out*(nucubedprime-local_J(j_prime))* &
               (JoverF(j)*local_u(2)*invX2+HoverF(j,2)) + &
               R1out*(local_Hdown(j_prime,1)*(HoverF(j,1)* &
               local_u(2)*invX2+LoverF(j,2,1)) + &
               local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(2)* &
               invX2+LoverF(j,2,2))))*nulibtable_inv_energies(j_prime)
          ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
               ies_temp*M1flux_energy(j)
       enddo
    enddo

    S = ies_sourceterms

  end subroutine evaluate_ies_source

end subroutine M1_source_timestep_limit
