!-*-f90-*-

subroutine M1_energycoupling_source_jacobian(zone_index,species_index,U,S,A, &
     need_jacobian)

  use GR1D_module
  use nulibtable, only : nulibtable_energies,nulibtable_etop, &
       nulibtable_ebottom,nulibtable_logenergies,nulibtable_logetop
  implicit none

  integer, intent(in) :: zone_index,species_index
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(out) :: S(2*number_groups)
  real*8, intent(out) :: A(2*number_groups,2*number_groups)
  logical, intent(in) :: need_jacobian

  real*8 :: div_v,dvdt_local
  real*8 :: h,dmdr,dmdt,dXdr,dXdt,Kdownrr,dWdt,dWdr,dWvuprdt,dWvuprdr
  real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6)
  real*8 :: heattermff_NL(6),velocity_coeffs(6,2)
  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1en_energy_fluid(number_groups),M1eddy_energy(number_groups)
  real*8 :: M1pff_energy(number_groups),M1qrrr_energy(number_groups)
  real*8 :: M1qffr_energy(number_groups),littlefactors(number_groups,6)
  real*8 :: velocity(number_groups,2),dvelocity(number_groups,2,2*number_groups)
  real*8 :: energy_interface(number_groups,2)
  real*8 :: denergy_interface(number_groups,2,2*number_groups)
  real*8 :: logdistro(number_groups),loginterface_distroj(number_groups)
  real*8 :: interface_distroj(number_groups),xi(number_groups)
  real*8 :: FL(number_groups,2),FR(number_groups,2)
  real*8 :: dFL(number_groups,2,2*number_groups)
  real*8 :: dFR(number_groups,2,2*number_groups)
  real*8 :: M1flux_diff_energy(number_groups,2)
  real*8 :: em,ep,enext,temp_term,width
  real*8 :: W2,v2,oneW,onev
  integer :: j,col

  S = 0.0d0
  A = 0.0d0
  dvelocity = 0.0d0
  dFL = 0.0d0
  dFR = 0.0d0
  energy_interface = 0.0d0
  denergy_interface = 0.0d0

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
     stop "M1_energycoupling_source_jacobian: unsupported v_order"
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

     if (need_jacobian) then
        dvelocity(j,1,j) = (velocity_coeffs(1,1) + &
             M1eddy_energy(j)*velocity_coeffs(3,1) + &
             M1pff_energy(j)*velocity_coeffs(4,1))/width
        dvelocity(j,2,j) = (velocity_coeffs(1,2) + &
             M1eddy_energy(j)*velocity_coeffs(3,2) + &
             M1pff_energy(j)*velocity_coeffs(4,2))/width
        dvelocity(j,1,j+number_groups) = velocity_coeffs(2,1)/width
        dvelocity(j,2,j+number_groups) = velocity_coeffs(2,2)/width
     endif
  enddo

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
  if (need_jacobian) then
     dFL(1,1,:) = dvelocity(1,1,:)*temp_term
     dFL(1,2,:) = dvelocity(1,2,:)*temp_term
  endif

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
     if (need_jacobian) then
        dFL(j,1,:) = dvelocity(j,1,:)*temp_term
        dFL(j,2,:) = dvelocity(j,2,:)*temp_term
     endif

     temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
          (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0)* &
          (1.0d0-xi(j))
     FR(j,1) = velocity(j,1)*temp_term
     FR(j,2) = velocity(j,2)*temp_term
     if (need_jacobian) then
        dFR(j,1,:) = dvelocity(j,1,:)*temp_term
        dFR(j,2,:) = dvelocity(j,2,:)*temp_term
     endif
  enddo

  do j=1,number_groups
     energy_interface(j,1) = energy_interface(j,1) + &
          FL(j,1)/nulibtable_etop(j)
     energy_interface(j,2) = energy_interface(j,2) + &
          FL(j,2)/nulibtable_etop(j)
     if (need_jacobian) then
        denergy_interface(j,1,:) = denergy_interface(j,1,:) + &
             dFL(j,1,:)/nulibtable_etop(j)
        denergy_interface(j,2,:) = denergy_interface(j,2,:) + &
             dFL(j,2,:)/nulibtable_etop(j)
     endif

     if (j.gt.1) then
        energy_interface(j-1,1) = energy_interface(j-1,1) + &
             FR(j,1)/nulibtable_ebottom(j)
        energy_interface(j-1,2) = energy_interface(j-1,2) + &
             FR(j,2)/nulibtable_ebottom(j)
        if (need_jacobian) then
           denergy_interface(j-1,1,:) = denergy_interface(j-1,1,:) + &
                dFR(j,1,:)/nulibtable_ebottom(j)
           denergy_interface(j-1,2,:) = denergy_interface(j-1,2,:) + &
                dFR(j,2,:)/nulibtable_ebottom(j)
        endif
     endif
  enddo

  j=1
  ep = nulibtable_etop(j)
  M1flux_diff_energy(j,1) = ep*energy_interface(j,1)
  M1flux_diff_energy(j,2) = ep*energy_interface(j,2)
  S(j) = -M1flux_diff_energy(j,1)
  S(j+number_groups) = -M1flux_diff_energy(j,2)
  if (need_jacobian) then
     A(j,:) = -ep*denergy_interface(j,1,:)
     A(j+number_groups,:) = -ep*denergy_interface(j,2,:)
  endif

  do j=2,number_groups
     em = nulibtable_ebottom(j)
     ep = nulibtable_etop(j)
     M1flux_diff_energy(j,1) = ep*energy_interface(j,1) - &
          em*energy_interface(j-1,1)
     M1flux_diff_energy(j,2) = ep*energy_interface(j,2) - &
          em*energy_interface(j-1,2)
     S(j) = -M1flux_diff_energy(j,1)
     S(j+number_groups) = -M1flux_diff_energy(j,2)
     if (need_jacobian) then
        do col=1,2*number_groups
           A(j,col) = -(ep*denergy_interface(j,1,col) - &
                em*denergy_interface(j-1,1,col))
           A(j+number_groups,col) = -(ep*denergy_interface(j,2,col) - &
                em*denergy_interface(j-1,2,col))
        enddo
     endif
  enddo

end subroutine M1_energycoupling_source_jacobian

subroutine M1_energycoupling_source_rate_bound(zone_index,species_index,U,S, &
     lambda_energy)

  use GR1D_module
  use nulibtable, only : nulibtable_energies,nulibtable_etop, &
       nulibtable_ebottom,nulibtable_logenergies,nulibtable_logetop
  implicit none

  integer, intent(in) :: zone_index,species_index
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(out) :: S(2*number_groups)
  real*8, intent(out) :: lambda_energy

  real*8 :: div_v,dvdt_local
  real*8 :: h,dmdr,dmdt,dXdr,dXdt,Kdownrr,dWdt,dWdr,dWvuprdt,dWvuprdr
  real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6)
  real*8 :: heattermff_NL(6),velocity_coeffs(6,2)
  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1en_energy_fluid(number_groups),M1eddy_energy(number_groups)
  real*8 :: M1pff_energy(number_groups),M1qrrr_energy(number_groups)
  real*8 :: M1qffr_energy(number_groups),littlefactors(number_groups,6)
  real*8 :: velocity(number_groups,2),velocity_bound(number_groups,2)
  real*8 :: energy_interface(number_groups,2)
  real*8 :: interface_bound(number_groups,2)
  real*8 :: logdistro(number_groups),loginterface_distroj(number_groups)
  real*8 :: interface_distroj(number_groups),xi(number_groups)
  real*8 :: FL(number_groups,2),FR(number_groups,2)
  real*8 :: FL_bound(number_groups,2),FR_bound(number_groups,2)
  real*8 :: M1flux_diff_energy(number_groups,2)
  real*8 :: em,ep,enext,temp_term,width
  real*8 :: W2,v2,oneW,onev
  integer :: j

  S = 0.0d0
  lambda_energy = 0.0d0
  energy_interface = 0.0d0
  interface_bound = 0.0d0
  velocity_bound = 0.0d0
  FL_bound = 0.0d0
  FR_bound = 0.0d0

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
     stop "M1_energycoupling_source_rate_bound: unsupported v_order"
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

     velocity_bound(j,1) = (abs(velocity_coeffs(1,1) + &
          M1eddy_energy(j)*velocity_coeffs(3,1) + &
          M1pff_energy(j)*velocity_coeffs(4,1)) + &
          abs(velocity_coeffs(2,1)))/width
     velocity_bound(j,2) = (abs(velocity_coeffs(1,2) + &
          M1eddy_energy(j)*velocity_coeffs(3,2) + &
          M1pff_energy(j)*velocity_coeffs(4,2)) + &
          abs(velocity_coeffs(2,2)))/width
  enddo

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
  FL_bound(1,1) = velocity_bound(1,1)*abs(temp_term)
  FL_bound(1,2) = velocity_bound(1,2)*abs(temp_term)

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
     FL_bound(j,1) = velocity_bound(j,1)*abs(temp_term)
     FL_bound(j,2) = velocity_bound(j,2)*abs(temp_term)

     temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
          (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0)* &
          (1.0d0-xi(j))
     FR(j,1) = velocity(j,1)*temp_term
     FR(j,2) = velocity(j,2)*temp_term
     FR_bound(j,1) = velocity_bound(j,1)*abs(temp_term)
     FR_bound(j,2) = velocity_bound(j,2)*abs(temp_term)
  enddo

  do j=1,number_groups
     energy_interface(j,1) = energy_interface(j,1) + &
          FL(j,1)/nulibtable_etop(j)
     energy_interface(j,2) = energy_interface(j,2) + &
          FL(j,2)/nulibtable_etop(j)
     interface_bound(j,1) = interface_bound(j,1) + &
          FL_bound(j,1)/nulibtable_etop(j)
     interface_bound(j,2) = interface_bound(j,2) + &
          FL_bound(j,2)/nulibtable_etop(j)

     if (j.gt.1) then
        energy_interface(j-1,1) = energy_interface(j-1,1) + &
             FR(j,1)/nulibtable_ebottom(j)
        energy_interface(j-1,2) = energy_interface(j-1,2) + &
             FR(j,2)/nulibtable_ebottom(j)
        interface_bound(j-1,1) = interface_bound(j-1,1) + &
             FR_bound(j,1)/nulibtable_ebottom(j)
        interface_bound(j-1,2) = interface_bound(j-1,2) + &
             FR_bound(j,2)/nulibtable_ebottom(j)
     endif
  enddo

  j=1
  ep = nulibtable_etop(j)
  M1flux_diff_energy(j,1) = ep*energy_interface(j,1)
  M1flux_diff_energy(j,2) = ep*energy_interface(j,2)
  S(j) = -M1flux_diff_energy(j,1)
  S(j+number_groups) = -M1flux_diff_energy(j,2)
  lambda_energy = max(lambda_energy,abs(ep)*interface_bound(j,1), &
       abs(ep)*interface_bound(j,2))

  do j=2,number_groups
     em = nulibtable_ebottom(j)
     ep = nulibtable_etop(j)
     M1flux_diff_energy(j,1) = ep*energy_interface(j,1) - &
          em*energy_interface(j-1,1)
     M1flux_diff_energy(j,2) = ep*energy_interface(j,2) - &
          em*energy_interface(j-1,2)
     S(j) = -M1flux_diff_energy(j,1)
     S(j+number_groups) = -M1flux_diff_energy(j,2)
     lambda_energy = max(lambda_energy, &
          abs(ep)*interface_bound(j,1)+abs(em)*interface_bound(j-1,1), &
          abs(ep)*interface_bound(j,2)+abs(em)*interface_bound(j-1,2))
  enddo

end subroutine M1_energycoupling_source_rate_bound

subroutine M1_ies_source_rate_bound(zone_index,species_index,U,S,lambda_ies, &
     suppress_high_density)

  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies
  implicit none

  integer, intent(in) :: zone_index,species_index
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(out) :: S(2*number_groups)
  real*8, intent(out) :: lambda_ies
  logical, intent(in) :: suppress_high_density

  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1eddy_energy(number_groups),M1chi_energy(number_groups)
  real*8 :: local_M(2,2),local_J(number_groups),local_H(number_groups,2)
  real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
  real*8 :: local_Ltilde(number_groups,2,2)
  real*8 :: JoverE(number_groups),JoverF(number_groups)
  real*8 :: HoverE(number_groups,2),HoverF(number_groups,2)
  real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)
  real*8 :: row_bound(2*number_groups)
  real*8 :: invalp,invalp2,invX,invX2,X2,alp2,W2,v2
  real*8 :: onev,oneW,onealp,oneX
  real*8 :: local_u(2),local_uup(2),local_littleh(2,2)
  real*8 :: local_littlehupup(2,2)
  real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in
  real*8 :: species_factor,ispecies_factor,kernel_fac
  real*8 :: coef,term,base1,base2,factor1,factor2
  integer :: j,j_prime,ii,jj

  S = 0.0d0
  row_bound = 0.0d0
  lambda_ies = 0.0d0

  if (species_index.eq.3.and.number_species.eq.3) then
     species_factor = 4.0d0
  else if (species_index.eq.3) then
     stop "M1_ies_source_rate_bound: unsupported species factor"
  else
     species_factor = 1.0d0
  endif
  ispecies_factor = 1.0d0/species_factor

  M1en_energy = U(1:number_groups)*ispecies_factor
  M1flux_energy = U(number_groups+1:2*number_groups)*ispecies_factor
  M1eddy_energy = q_M1(zone_index,species_index,:,3)
  M1chi_energy = q_M1_extra(zone_index,species_index,:,4)

  local_J = 0.0d0
  local_H = 0.0d0
  local_L = 0.0d0
  local_Hdown = 0.0d0
  local_Ltilde = 0.0d0

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
     stop "M1_ies_source_rate_bound: unsupported v_order"
  endif

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
     local_M = 0.0d0
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
                local_littlehupup(1,1)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,1,2) = local_L(j,1,2) + &
                local_M(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(1,2)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,2,1) = local_L(j,2,1) + &
                local_M(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(2,1)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,2,2) = local_L(j,2,2) + &
                local_M(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(2,2)*local_J(j)*onethird* &
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

        if (suppress_high_density.and.rho(zone_index).gt.5.0d12*rho_gf) then
           kernel_fac = (5.0d12*rho_gf/rho(zone_index))**1.5d0
           R0out = R0out*kernel_fac
           R1out = R1out*kernel_fac
           if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
           R0in = R0in*kernel_fac
           R1in = R1in*kernel_fac
        endif

        coef = species_factor*alp2*4.0d0*pi*nulibtable_inv_energies(j_prime)
        base1 = ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
             local_H(j,1))*R0in
        base2 = (nucubed-local_J(j))*local_littleh(1,1)*onethird - &
             local_u(1)*(-invalp2)*local_Hdown(j,1) - local_Ltilde(j,1,1)
        factor1 = (nucubed-local_J(j))*local_littleh(1,2)*onethird - &
             local_u(1)*(-invalp2)*local_Hdown(j,2) - local_Ltilde(j,1,2)

        term = base1*JoverE(j_prime) + &
             R1in*(HoverE(j_prime,1)*base2 + HoverE(j_prime,2)*factor1)
        S(j) = S(j) + coef*term*M1en_energy(j_prime)
        row_bound(j) = row_bound(j) + abs(coef*term)*ispecies_factor

        term = base1*JoverF(j_prime) + &
             R1in*(HoverF(j_prime,1)*base2 + HoverF(j_prime,2)*factor1)
        S(j) = S(j) + coef*term*M1flux_energy(j_prime)
        row_bound(j) = row_bound(j) + abs(coef*term)*ispecies_factor

        factor1 = JoverE(j)*local_u(1)*(-invalp2)+HoverE(j,1)
        factor2 = HoverE(j,1)*local_u(1)*(-invalp2)+LoverE(j,1,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(1)* &
             (-invalp2)+LoverE(j,1,2)))
        S(j) = S(j) + coef*term*M1en_energy(j)
        row_bound(j) = row_bound(j) + abs(coef*term)*ispecies_factor

        factor1 = JoverF(j)*local_u(1)*(-invalp2)+HoverF(j,1)
        factor2 = HoverF(j,1)*local_u(1)*(-invalp2)+LoverF(j,1,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(1)* &
             (-invalp2)+LoverF(j,1,2)))
        S(j) = S(j) + coef*term*M1flux_energy(j)
        row_bound(j) = row_bound(j) + abs(coef*term)*ispecies_factor

        coef = species_factor*onealp*X2*4.0d0*pi* &
             nulibtable_inv_energies(j_prime)
        base1 = ((nucubed-local_J(j))*(local_u(2)*invX2) - &
             local_H(j,2))*R0in
        base2 = (nucubed-local_J(j))*local_littleh(2,1)*onethird - &
             local_u(2)*X2*local_Hdown(j,1) - local_Ltilde(j,2,1)
        factor1 = (nucubed-local_J(j))*local_littleh(2,2)*onethird - &
             local_u(2)*X2*local_Hdown(j,2) - local_Ltilde(j,2,2)

        term = base1*JoverE(j_prime) + &
             R1in*(HoverE(j_prime,1)*base2 + HoverE(j_prime,2)*factor1)
        S(j+number_groups) = S(j+number_groups) + &
             coef*term*M1en_energy(j_prime)
        row_bound(j+number_groups) = row_bound(j+number_groups) + &
             abs(coef*term)*ispecies_factor

        term = base1*JoverF(j_prime) + &
             R1in*(HoverF(j_prime,1)*base2 + HoverF(j_prime,2)*factor1)
        S(j+number_groups) = S(j+number_groups) + &
             coef*term*M1flux_energy(j_prime)
        row_bound(j+number_groups) = row_bound(j+number_groups) + &
             abs(coef*term)*ispecies_factor

        factor1 = JoverE(j)*local_u(2)*invX2+HoverE(j,2)
        factor2 = HoverE(j,1)*local_u(2)*invX2+LoverE(j,2,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(2)* &
             invX2+LoverE(j,2,2)))
        S(j+number_groups) = S(j+number_groups) + coef*term*M1en_energy(j)
        row_bound(j+number_groups) = row_bound(j+number_groups) + &
             abs(coef*term)*ispecies_factor

        factor1 = JoverF(j)*local_u(2)*invX2+HoverF(j,2)
        factor2 = HoverF(j,1)*local_u(2)*invX2+LoverF(j,2,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(2)* &
             invX2+LoverF(j,2,2)))
        S(j+number_groups) = S(j+number_groups) + coef*term*M1flux_energy(j)
        row_bound(j+number_groups) = row_bound(j+number_groups) + &
             abs(coef*term)*ispecies_factor
     enddo
  enddo

  lambda_ies = maxval(row_bound)

end subroutine M1_ies_source_rate_bound

subroutine M1_ies_source_jacobian(zone_index,species_index,U,S,A, &
     need_jacobian,suppress_high_density)

  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies
  implicit none

  integer, intent(in) :: zone_index,species_index
  real*8, intent(in) :: U(2*number_groups)
  real*8, intent(out) :: S(2*number_groups)
  real*8, intent(out) :: A(2*number_groups,2*number_groups)
  logical, intent(in) :: need_jacobian,suppress_high_density

  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1eddy_energy(number_groups),M1chi_energy(number_groups)
  real*8 :: local_M(2,2),dM(2,2)
  real*8 :: local_J(number_groups),local_H(number_groups,2)
  real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
  real*8 :: local_Ltilde(number_groups,2,2)
  real*8 :: dJ(number_groups,2*number_groups)
  real*8 :: dH(number_groups,2,2*number_groups)
  real*8 :: dL(number_groups,2,2,2*number_groups)
  real*8 :: dHdown(number_groups,2,2*number_groups)
  real*8 :: dLtilde(number_groups,2,2,2*number_groups)
  real*8 :: JoverE(number_groups),JoverF(number_groups)
  real*8 :: HoverE(number_groups,2),HoverF(number_groups,2)
  real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)
  real*8 :: invalp,invalp2,invX,invX2,X2,alp2,W2,v2
  real*8 :: onev,oneW,onealp,oneX
  real*8 :: local_u(2),local_uup(2),local_littleh(2,2)
  real*8 :: local_littlehupup(2,2)
  real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in
  real*8 :: species_factor,ispecies_factor,kernel_fac
  real*8 :: coef,term,dterm,dmult
  real*8 :: base1,base2,dbase1,dbase2,factor1,factor2
  real*8 :: dlocal_J
  integer :: j,j_prime,ii,jj,col

  S = 0.0d0
  A = 0.0d0

  if (species_index.eq.3.and.number_species.eq.3) then
     species_factor = 4.0d0
  else if (species_index.eq.3) then
     stop "M1_ies_source_jacobian: unsupported species factor"
  else
     species_factor = 1.0d0
  endif
  ispecies_factor = 1.0d0/species_factor

  M1en_energy = U(1:number_groups)*ispecies_factor
  M1flux_energy = U(number_groups+1:2*number_groups)*ispecies_factor
  M1eddy_energy = q_M1(zone_index,species_index,:,3)
  M1chi_energy = q_M1_extra(zone_index,species_index,:,4)

  local_J = 0.0d0
  local_H = 0.0d0
  local_L = 0.0d0
  local_Hdown = 0.0d0
  local_Ltilde = 0.0d0
  dJ = 0.0d0
  dH = 0.0d0
  dL = 0.0d0
  dHdown = 0.0d0
  dLtilde = 0.0d0

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
     stop "M1_ies_source_jacobian: unsupported v_order"
  endif

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
     local_M = 0.0d0
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
                local_littlehupup(1,1)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,1,2) = local_L(j,1,2) + &
                local_M(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(1,2)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,2,1) = local_L(j,2,1) + &
                local_M(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(2,1)*local_J(j)*onethird* &
                (1.5d0-1.5d0*M1chi_energy(j))
           local_L(j,2,2) = local_L(j,2,2) + &
                local_M(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj)* &
                (1.5d0*M1chi_energy(j)-0.5d0) + &
                local_littlehupup(2,2)*local_J(j)*onethird* &
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

     if (need_jacobian) then
        do col=1,2*number_groups
           dM = 0.0d0
           if (col.eq.j) then
              dM(1,1) = invalp2*ispecies_factor
              dM(2,2) = M1eddy_energy(j)*invX2**2*ispecies_factor
           else if (col.eq.j+number_groups) then
              dM(1,2) = invX2*invalp*ispecies_factor
              dM(2,1) = dM(1,2)
           endif

           do ii=1,2
              do jj=1,2
                 dJ(j,col) = dJ(j,col) + dM(ii,jj)*local_u(ii)*local_u(jj)
                 dH(j,1,col) = dH(j,1,col) - &
                      dM(ii,jj)*local_u(ii)*local_littleh(1,jj)
                 dH(j,2,col) = dH(j,2,col) - &
                      dM(ii,jj)*local_u(ii)*local_littleh(2,jj)
              enddo
           enddo

           do ii=1,2
              do jj=1,2
                 dL(j,1,1,col) = dL(j,1,1,col) + &
                      dM(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj)* &
                      (1.5d0*M1chi_energy(j)-0.5d0) + &
                      local_littlehupup(1,1)*dJ(j,col)*onethird* &
                      (1.5d0-1.5d0*M1chi_energy(j))
                 dL(j,1,2,col) = dL(j,1,2,col) + &
                      dM(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj)* &
                      (1.5d0*M1chi_energy(j)-0.5d0) + &
                      local_littlehupup(1,2)*dJ(j,col)*onethird* &
                      (1.5d0-1.5d0*M1chi_energy(j))
                 dL(j,2,1,col) = dL(j,2,1,col) + &
                      dM(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj)* &
                      (1.5d0*M1chi_energy(j)-0.5d0) + &
                      local_littlehupup(2,1)*dJ(j,col)*onethird* &
                      (1.5d0-1.5d0*M1chi_energy(j))
                 dL(j,2,2,col) = dL(j,2,2,col) + &
                      dM(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj)* &
                      (1.5d0*M1chi_energy(j)-0.5d0) + &
                      local_littlehupup(2,2)*dJ(j,col)*onethird* &
                      (1.5d0-1.5d0*M1chi_energy(j))
              enddo
           enddo

           dHdown(j,1,col) = dH(j,1,col)*local_littleh(1,1)*(-alp2) + &
                dH(j,2,col)*local_littleh(1,2)*(-alp2)
           dHdown(j,2,col) = dH(j,1,col)*local_littleh(2,1)*X2 + &
                dH(j,2,col)*local_littleh(2,2)*X2

           dLtilde(j,1,1,col) = -dL(j,1,1,col)*alp2 - &
                dJ(j,col)*local_littleh(1,1)*onethird
           dLtilde(j,1,2,col) = dL(j,1,2,col)*X2 - &
                dJ(j,col)*local_littleh(1,2)*onethird
           dLtilde(j,2,1,col) = -dL(j,2,1,col)*alp2 - &
                dJ(j,col)*local_littleh(2,1)*onethird
           dLtilde(j,2,2,col) = dL(j,2,2,col)*X2 - &
                dJ(j,col)*local_littleh(2,2)*onethird
        enddo
     endif

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

        if (suppress_high_density.and.rho(zone_index).gt.5.0d12*rho_gf) then
           kernel_fac = (5.0d12*rho_gf/rho(zone_index))**1.5d0
           R0out = R0out*kernel_fac
           R1out = R1out*kernel_fac
           if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
           R0in = R0in*kernel_fac
           R1in = R1in*kernel_fac
        endif

        coef = species_factor*alp2*4.0d0*pi*nulibtable_inv_energies(j_prime)
        base1 = ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
             local_H(j,1))*R0in
        base2 = (nucubed-local_J(j))*local_littleh(1,1)*onethird - &
             local_u(1)*(-invalp2)*local_Hdown(j,1) - local_Ltilde(j,1,1)
        factor1 = (nucubed-local_J(j))*local_littleh(1,2)*onethird - &
             local_u(1)*(-invalp2)*local_Hdown(j,2) - local_Ltilde(j,1,2)
        term = base1*JoverE(j_prime) + &
             R1in*(HoverE(j_prime,1)*base2 + HoverE(j_prime,2)*factor1)
        S(j) = S(j) + coef*term*M1en_energy(j_prime)
        if (need_jacobian) then
           do col=1,2*number_groups
              dbase1 = ((-dJ(j,col))*(-local_u(1)*invalp2) - &
                   dH(j,1,col))*R0in
              dbase2 = -dJ(j,col)*local_littleh(1,1)*onethird - &
                   local_u(1)*(-invalp2)*dHdown(j,1,col) - &
                   dLtilde(j,1,1,col)
              dlocal_J = -dJ(j,col)*local_littleh(1,2)*onethird - &
                   local_u(1)*(-invalp2)*dHdown(j,2,col) - &
                   dLtilde(j,1,2,col)
              dterm = dbase1*JoverE(j_prime) + &
                   R1in*(HoverE(j_prime,1)*dbase2 + HoverE(j_prime,2)*dlocal_J)
              dmult = 0.0d0
              if (col.eq.j_prime) dmult = ispecies_factor
              A(j,col) = A(j,col) + coef*(dterm*M1en_energy(j_prime) + &
                   term*dmult)
           enddo
        endif

        term = base1*JoverF(j_prime) + &
             R1in*(HoverF(j_prime,1)*base2 + HoverF(j_prime,2)*factor1)
        S(j) = S(j) + coef*term*M1flux_energy(j_prime)
        if (need_jacobian) then
           do col=1,2*number_groups
              dbase1 = ((-dJ(j,col))*(-local_u(1)*invalp2) - &
                   dH(j,1,col))*R0in
              dbase2 = -dJ(j,col)*local_littleh(1,1)*onethird - &
                   local_u(1)*(-invalp2)*dHdown(j,1,col) - &
                   dLtilde(j,1,1,col)
              dlocal_J = -dJ(j,col)*local_littleh(1,2)*onethird - &
                   local_u(1)*(-invalp2)*dHdown(j,2,col) - &
                   dLtilde(j,1,2,col)
              dterm = dbase1*JoverF(j_prime) + &
                   R1in*(HoverF(j_prime,1)*dbase2 + HoverF(j_prime,2)*dlocal_J)
              dmult = 0.0d0
              if (col.eq.j_prime+number_groups) dmult = ispecies_factor
              A(j,col) = A(j,col) + coef*(dterm*M1flux_energy(j_prime) + &
                   term*dmult)
           enddo
        endif

        factor1 = JoverE(j)*local_u(1)*(-invalp2)+HoverE(j,1)
        factor2 = HoverE(j,1)*local_u(1)*(-invalp2)+LoverE(j,1,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(1)* &
             (-invalp2)+LoverE(j,1,2)))
        S(j) = S(j) + coef*term*M1en_energy(j)
        if (need_jacobian) then
           do col=1,2*number_groups
              dterm = R0out*dJ(j_prime,col)*factor1 + &
                   R1out*(dHdown(j_prime,1,col)*factor2 + &
                   dHdown(j_prime,2,col)*(HoverE(j,2)*local_u(1)* &
                   (-invalp2)+LoverE(j,1,2)))
              dmult = 0.0d0
              if (col.eq.j) dmult = ispecies_factor
              A(j,col) = A(j,col) + coef*(dterm*M1en_energy(j) + term*dmult)
           enddo
        endif

        factor1 = JoverF(j)*local_u(1)*(-invalp2)+HoverF(j,1)
        factor2 = HoverF(j,1)*local_u(1)*(-invalp2)+LoverF(j,1,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(1)* &
             (-invalp2)+LoverF(j,1,2)))
        S(j) = S(j) + coef*term*M1flux_energy(j)
        if (need_jacobian) then
           do col=1,2*number_groups
              dterm = R0out*dJ(j_prime,col)*factor1 + &
                   R1out*(dHdown(j_prime,1,col)*factor2 + &
                   dHdown(j_prime,2,col)*(HoverF(j,2)*local_u(1)* &
                   (-invalp2)+LoverF(j,1,2)))
              dmult = 0.0d0
              if (col.eq.j+number_groups) dmult = ispecies_factor
              A(j,col) = A(j,col) + coef*(dterm*M1flux_energy(j) + term*dmult)
           enddo
        endif

        coef = species_factor*onealp*X2*4.0d0*pi* &
             nulibtable_inv_energies(j_prime)
        base1 = ((nucubed-local_J(j))*(local_u(2)*invX2) - &
             local_H(j,2))*R0in
        base2 = (nucubed-local_J(j))*local_littleh(2,1)*onethird - &
             local_u(2)*X2*local_Hdown(j,1) - local_Ltilde(j,2,1)
        factor1 = (nucubed-local_J(j))*local_littleh(2,2)*onethird - &
             local_u(2)*X2*local_Hdown(j,2) - local_Ltilde(j,2,2)
        term = base1*JoverE(j_prime) + &
             R1in*(HoverE(j_prime,1)*base2 + HoverE(j_prime,2)*factor1)
        S(j+number_groups) = S(j+number_groups) + &
             coef*term*M1en_energy(j_prime)
        if (need_jacobian) then
           do col=1,2*number_groups
              dbase1 = ((-dJ(j,col))*(local_u(2)*invX2) - &
                   dH(j,2,col))*R0in
              dbase2 = -dJ(j,col)*local_littleh(2,1)*onethird - &
                   local_u(2)*X2*dHdown(j,1,col) - dLtilde(j,2,1,col)
              dlocal_J = -dJ(j,col)*local_littleh(2,2)*onethird - &
                   local_u(2)*X2*dHdown(j,2,col) - dLtilde(j,2,2,col)
              dterm = dbase1*JoverE(j_prime) + &
                   R1in*(HoverE(j_prime,1)*dbase2 + HoverE(j_prime,2)*dlocal_J)
              dmult = 0.0d0
              if (col.eq.j_prime) dmult = ispecies_factor
              A(j+number_groups,col) = A(j+number_groups,col) + &
                   coef*(dterm*M1en_energy(j_prime) + term*dmult)
           enddo
        endif

        term = base1*JoverF(j_prime) + &
             R1in*(HoverF(j_prime,1)*base2 + HoverF(j_prime,2)*factor1)
        S(j+number_groups) = S(j+number_groups) + &
             coef*term*M1flux_energy(j_prime)
        if (need_jacobian) then
           do col=1,2*number_groups
              dbase1 = ((-dJ(j,col))*(local_u(2)*invX2) - &
                   dH(j,2,col))*R0in
              dbase2 = -dJ(j,col)*local_littleh(2,1)*onethird - &
                   local_u(2)*X2*dHdown(j,1,col) - dLtilde(j,2,1,col)
              dlocal_J = -dJ(j,col)*local_littleh(2,2)*onethird - &
                   local_u(2)*X2*dHdown(j,2,col) - dLtilde(j,2,2,col)
              dterm = dbase1*JoverF(j_prime) + &
                   R1in*(HoverF(j_prime,1)*dbase2 + HoverF(j_prime,2)*dlocal_J)
              dmult = 0.0d0
              if (col.eq.j_prime+number_groups) dmult = ispecies_factor
              A(j+number_groups,col) = A(j+number_groups,col) + &
                   coef*(dterm*M1flux_energy(j_prime) + term*dmult)
           enddo
        endif

        factor1 = JoverE(j)*local_u(2)*invX2+HoverE(j,2)
        factor2 = HoverE(j,1)*local_u(2)*invX2+LoverE(j,2,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverE(j,2)*local_u(2)* &
             invX2+LoverE(j,2,2)))
        S(j+number_groups) = S(j+number_groups) + coef*term*M1en_energy(j)
        if (need_jacobian) then
           do col=1,2*number_groups
              dterm = R0out*dJ(j_prime,col)*factor1 + &
                   R1out*(dHdown(j_prime,1,col)*factor2 + &
                   dHdown(j_prime,2,col)*(HoverE(j,2)*local_u(2)* &
                   invX2+LoverE(j,2,2)))
              dmult = 0.0d0
              if (col.eq.j) dmult = ispecies_factor
              A(j+number_groups,col) = A(j+number_groups,col) + &
                   coef*(dterm*M1en_energy(j) + term*dmult)
           enddo
        endif

        factor1 = JoverF(j)*local_u(2)*invX2+HoverF(j,2)
        factor2 = HoverF(j,1)*local_u(2)*invX2+LoverF(j,2,1)
        term = -R0out*(nucubedprime-local_J(j_prime))*factor1 + &
             R1out*(local_Hdown(j_prime,1)*factor2 + &
             local_Hdown(j_prime,2)*(HoverF(j,2)*local_u(2)* &
             invX2+LoverF(j,2,2)))
        S(j+number_groups) = S(j+number_groups) + coef*term*M1flux_energy(j)
        if (need_jacobian) then
           do col=1,2*number_groups
              dterm = R0out*dJ(j_prime,col)*factor1 + &
                   R1out*(dHdown(j_prime,1,col)*factor2 + &
                   dHdown(j_prime,2,col)*(HoverF(j,2)*local_u(2)* &
                   invX2+LoverF(j,2,2)))
              dmult = 0.0d0
              if (col.eq.j+number_groups) dmult = ispecies_factor
              A(j+number_groups,col) = A(j+number_groups,col) + &
                   coef*(dterm*M1flux_energy(j) + term*dmult)
           enddo
        endif
     enddo
  enddo

end subroutine M1_ies_source_jacobian
