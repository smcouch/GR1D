!-*-f90-*-
!this routine calculates all of the explicit terms (spatial
!flux,energy flux, and scattering terms) It returns the terms in the
!flux_M1 type variables flux_M1 (for spatial), flux_M1_energy (for
!energy), flux_M1_scatter (for scattering)
subroutine M1_explicitterms(dts,implicit_factor)
  
  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies,nulibtable_ewidths, &
       nulibtable_energies,nulibtable_etop,nulibtable_ebottom, &
       nulibtable_logenergies,nulibtable_logetop
  implicit none

  !inputs:
  real*8 :: dts !time step
  real*8 :: implicit_factor !time step

  !local, spatial
  real*8 :: M1en_space(n1),M1flux_space(n1)
  real*8 :: M1chi_space(n1),M1eddy_space(n1)
  real*8 :: M1en_space_plus(n1),M1en_space_minus(n1)
  real*8 :: M1flux_space_plus(n1),M1flux_space_minus(n1)
  real*8 :: M1chi_space_plus(n1),M1chi_space_minus(n1)
  real*8 :: M1eddy_space_plus(n1),M1eddy_space_minus(n1)
  real*8 :: l_min_thin(2),l_max_thin(2)
  real*8 :: l_min_thick(2),l_max_thick(2) 
  real*8 :: l_min(3),l_max(3)
  real*8 :: p,discrim,sqrtdiscrim
  real*8 :: M1flux_interface(n1,2),M1flux_diff(n1,2)
  real*8 :: rm,rp,dx

  !local, energy
  real*8 :: div_v(n1),dvdt(n1)
  real*8 :: h,dmdr,dmdt,dXdr,dXdt,Kdownrr,dWdt,dWdr,dWvuprdt,dWvuprdr
  real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6),heattermff_NL(6)
  real*8 :: velocity_coeffs(6,2)
  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1eddy_energy(number_groups),M1pff_energy(number_groups)
  real*8 :: M1en_energy_fluid(number_groups),M1chi_energy(number_groups)
  real*8 :: M1qrrr_energy(number_groups),M1qffr_energy(number_groups)
  real*8 :: littlefactors_int(number_groups,6)
  real*8 :: littlefactors(number_groups,6)
  real*8 :: velocity(number_groups,2)
  real*8 :: velocity_top(number_groups,2)
  real*8 :: log_distro(number_groups)
  real*8 :: log_energy_interface_distroj(number_groups)
  real*8 :: energy_interface_distroj(number_groups)
  real*8 :: energy_interface_M1en(number_groups)
  real*8 :: M1flux_energy_interface(number_groups,2)
  real*8 :: em,ep,de,enext
  real*8 :: M1flux_diff_energy(number_groups,2)
  real*8 :: M1_moment_to_distro_nobinwidth(number_groups) 
  real*8 :: M1_moment_to_distro_top_nobinwidth(number_groups) 
  real*8 :: nulibtable_etop_logged(number_groups)
  real*8 :: nulibtable_energies_logged(number_groups)

  real*8 :: a_asym,kappa_inter,ipeclet_mean
  real*8 :: Jkplus1,Jk,diffusive_flux,advected_energy

  real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in,ies_temp,species_factor
  real*8 :: ies_sourceterms(2*number_groups)
  real*8 :: local_M(2,2),local_J(number_groups),local_H(number_groups,2)
  real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
  real*8 :: local_Ltilde(number_groups,2,2),JoverE(number_groups)
  real*8 :: JoverF(number_groups),HoverE(number_groups,2),HoverF(number_groups,2)
  real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)
  real*8 :: invalp,invalp2,invX,invX2,X2,alp2,W2,v2,invr,invr2,onev,oneW,onealp,oneX,oneWm,oneWp
  real*8 :: local_u(2),local_littleh(2,2),local_littlehupup(2,2),local_uup(2)

  !mueller
  real*8 :: logdistro(number_groups)
  real*8 :: loginterface_distroj(number_groups)
  real*8 :: interface_distroj(number_groups)
  real*8 :: xi(number_groups)
  real*8 :: FL(number_groups,2),FR(number_groups,2)
  real*8 :: temp_term
  real*8 :: U_source(2*number_groups),S_source(2*number_groups)
  real*8 :: lambda_source
  
  real*8 :: fluxtemp_1,fluxtemp_2
  real*8 :: limitingflux

  ! turbulence
  real*8 :: diffusive_turb_flux, grad_Enu, D_nu_turb
  real*8 :: Lambda_mixp
  
  !counters
  integer i,j,k,j_prime,jj,ii

  !reset flux to zero
  flux_M1 = 0.0d0
  flux_M1_energy = 0.0d0
  flux_M1_scatter = 0.0d0
  M1_source_lambda_ies = 0.0d0
  M1_source_lambda_energycoupling = 0.0d0
  call M1_source_timestep_cache_reset

!#################################################################
!#################################################################
!########################Spatial##################################
!#################################################################
!#################################################################

  !first spatial. This is essentially a reimann solve + corrections
  !over all the radial zones. Variables are already interpolated

  !$OMP PARALLEL DO PRIVATE(i,M1en_space,M1flux_space,M1chi_space,M1eddy_space, &
  !$OMP M1en_space_plus,M1flux_space_plus,M1chi_space_plus,M1eddy_space_plus, &
  !$OMP M1en_space_minus,M1flux_space_minus,M1chi_space_minus,M1eddy_space_minus, &
  !$OMP k,oneWm,oneWp,kappa_inter,ipeclet_mean,a_asym,l_min_thin,l_max_thin, &
  !$OMP p,discrim,sqrtdiscrim,l_min_thick,l_max_thick,l_min,l_max,Jkplus1,Jk,diffusive_flux, &
  !$OMP advected_energy,M1flux_interface,limitingflux,rm,rp,dx,M1flux_diff, &
  !$OMP diffusive_turb_flux, grad_Enu, D_nu_turb, Lambda_mixp)
  do j=1,number_groups
     do i=1,number_species_to_evolve

        M1en_space = q_M1(:,i,j,1)
        M1flux_space = q_M1(:,i,j,2)
        M1chi_space = q_M1_extra(:,i,j,4)
        M1eddy_space = q_M1(:,i,j,3)

        M1en_space_plus = q_M1p(:,i,j,1,1)
        M1en_space_minus = q_M1m(:,i,j,1,1)
        M1flux_space_plus = q_M1p(:,i,j,2,1)
        M1flux_space_minus = q_M1m(:,i,j,2,1)

        M1eddy_space_plus = q_M1p(:,i,j,3,1)
        M1eddy_space_minus = q_M1m(:,i,j,3,1)
        M1chi_space_plus = q_M1_extrap(:,i,j,1,1)
        M1chi_space_minus = q_M1_extram(:,i,j,1,1)

        !now find speeds at each interface
        do k=ghosts1,M1_imaxradii

           !only used in when GR.eq.0 and v_order.eq.-1
           oneWm = 1.0d0/sqrt(1.0d0-v1m(k+1)**2)
           oneWp = 1.0d0/sqrt(1.0d0-v1p(k)**2)

           !determine the regime for the flux calculation
           kappa_inter = sqrt((eas(k,i,j,2)+eas(k,i,j,3))*(eas(k+1,i,j,2)+eas(k+1,i,j,3)))
           if (v_order.eq.-1) then
              if (GR) then
                 ipeclet_mean = 1.0d0/(Wm(k+1)**3*(1.0d0+vm(k+1))*Xm(k+1)**2* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              else
                 ipeclet_mean = 1.0d0/(oneWm**3*(1.0d0+v1m(k+1))* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              endif
                 
           else if (v_order.eq.0) then
              if (GR) then
                 ipeclet_mean = 1.0d0/(Xm(k+1)**2* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              else
                 ipeclet_mean = 1.0d0/(kappa_inter*(x1(k+1)-x1(k))) 
              endif
           else 
              stop "add in v order, peclet number"
           endif

           !luke's term
           a_asym = tanh(ipeclet_mean)
           if (a_asym.gt.1.0d0) then
              a_asym = 1.0d0
           endif

           !speeds
           !minus interface (k+1 zone)
           !thin limit:
           if (GR) then
              l_min_thin(1) = -Xm(k+1)*alpm(k+1)
              l_max_thin(1) = Xm(k+1)*alpm(k+1)
           else
              if (do_effectivepotential) then
                 l_min_thin(1) = -alpm(k+1)
                 l_max_thin(1) = alpm(k+1)
              else
                 l_min_thin(1) = -1.0d0
                 l_max_thin(1) = 1.0d0
              endif
           endif
           
           !thick limit:
           if (v_order.eq.-1) then
              if (GR) then
                 p = alpm(k+1)*vm(k+1)*Xm(k+1)
                 discrim = alpm(k+1)**2*Xm(k+1)**2*(2.0d0*Wm(k+1)**2+1.0d0)-2.0d0*Wm(k+1)**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(1) = min((2.0d0*Wm(k+1)**2*p - sqrtdiscrim)/(2.0d0*Wm(k+1)**2+1.0d0),p)
                 l_max_thick(1) = max((2.0d0*Wm(k+1)**2*p + sqrtdiscrim)/(2.0d0*Wm(k+1)**2+1.0d0),p)
              else
                 if (do_effectivepotential) then
                    p = alpm(k+1)*v1m(k+1)
                    discrim = alpm(k+1)**2*(2.0d0*oneWm**2+1.0d0)-2.0d0*oneWm**2*p*p
                 else
                    p = v1m(k+1)
                    discrim = (2.0d0*oneWm**2+1.0d0)-2.0d0*oneWm**2*p*p
                 endif
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(1) = min((2.0d0*oneWm**2*p - sqrtdiscrim)/(2.0d0*oneWm**2+1.0d0),p)
                 l_max_thick(1) = max((2.0d0*oneWm**2*p + sqrtdiscrim)/(2.0d0*oneWm**2+1.0d0),p)
              endif


           else if (v_order.eq.0) then
              p = 0.0d0
              if (GR) then
                 discrim = alpm(k+1)**2*Xm(k+1)**2*3.0d0
              else
                 if (do_effectivepotential) then
                    discrim = alpm(k+1)**2*3.0d0
                 else
                    discrim = 3.0d0
                 endif
              endif
              sqrtdiscrim = sqrt(discrim)
              l_min_thick(1) = -sqrtdiscrim/3.0d0
              l_max_thick(1) = sqrtdiscrim/3.0d0
           else
              stop "add me in"
           endif

           !actual speed
           l_min(1) = (3.0d0*M1chi_space(k+1)-1.0d0)*0.5d0*l_min_thin(1) + &
                (1.0d0 - M1chi_space(k+1))*1.5d0*l_min_thick(1)

           l_max(1) = (3.0d0*M1chi_space(k+1)-1.0d0)*0.5d0*l_max_thin(1) + &
                (1.0d0 - M1chi_space(k+1))*1.5d0*l_max_thick(1)

           !plus interface (k zone)
           !thin limit:
           if (GR) then
              l_min_thin(2) = -Xp(k)*alpp(k)
              l_max_thin(2) = Xp(k)*alpp(k)
           else
              if (do_effectivepotential) then
                 l_min_thin(2) = -alpp(k)
                 l_max_thin(2) = alpp(k)
              else
                 l_min_thin(2) = -1.0d0
                 l_max_thin(2) = 1.0d0
              endif
           endif
              
           
           !thick limit:
           if (v_order.eq.-1) then
              if (GR) then
                 p = alpp(k)*vp(k)*Xp(k)
                 discrim = alpp(k)**2*Xp(k)**2*(2.0d0*Wp(k)**2+1.0d0)-2.0d0*Wp(k)**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(2) = min((2.0d0*Wp(k)**2*p - sqrtdiscrim)/(2.0d0*Wp(k)**2+1.0d0),p)
                 l_max_thick(2) = max((2.0d0*Wp(k)**2*p + sqrtdiscrim)/(2.0d0*Wp(k)**2+1.0d0),p)
              else
                 if (do_effectivepotential) then
                    p = alpp(k)*v1p(k)
                    discrim = alpp(k)**2*(2.0d0*oneWp**2+1.0d0)-2.0d0*oneWp**2*p*p
                 else
                    p = v1p(k)
                    discrim = (2.0d0*oneWp**2+1.0d0)-2.0d0*oneWp**2*p*p
                 endif
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(2) = min((2.0d0*oneWp**2*p - sqrtdiscrim)/(2.0d0*oneWp**2+1.0d0),p)
                 l_max_thick(2) = max((2.0d0*oneWp**2*p + sqrtdiscrim)/(2.0d0*oneWp**2+1.0d0),p)
              endif

           else if (v_order.eq.0) then
              p = 0.0d0
              if (GR) then
                 discrim = alpp(k)**2*Xp(k)**2*3.0d0
              else
                 if (do_effectivepotential) then
                    discrim = alpp(k)**2*3.0d0
                 else
                    discrim = 3.0d0
                 endif
              endif

              sqrtdiscrim = sqrt(discrim)
              l_min_thick(2) = -sqrtdiscrim/3.0d0
              l_max_thick(2) = sqrtdiscrim/3.0d0
           else
              stop "add me in"
           endif
           
           !actualy speed
           l_min(2) = (3.0d0*M1chi_space(k)-1.0d0)*0.5d0*l_min_thin(2) + &
                (1.0d0 - M1chi_space(k))*1.5d0*l_min_thick(2)
           l_max(2) = (3.0d0*M1chi_space(k)-1.0d0)*0.5d0*l_max_thin(2) + &
                (1.0d0 - M1chi_space(k))*1.5d0*l_max_thick(2)


           !check for NaNs
           if (l_min(1).ne.l_min(1)) then
              write(*,*) "NaNs in speeds 1",M1chi_space_minus(k+1),l_min_thin(1),l_min_thick(1),k,i,j
              stop
           else if(l_min(2).ne.l_min(2)) then
              write(*,*) "NaNs in speeds 2"
              stop
           else if(l_max(1).ne.l_max(1)) then
              write(*,*) "NaNs in speeds 3"
              stop
           else if(l_max(2).ne.l_max(2)) then
              write(*,*) "NaNs in speeds 4"
              stop
           endif

           l_max(3) = max(l_max(1),l_max(2))
           l_min(3) = max(l_min(1),l_min(2))

           !Lukes idea/suggestion
           if (v_order.eq.-1) then
              if (GR) then
!                 Jkplus1 = W(k+1)**2*(M1en_space(k+1)*(1.0d0+v(k+1)**2* &
!                      M1eddy_space(k+1)/X(k+1)**2)-2.0d0*M1flux_space(k+1)*v(k+1)/X(k+1))
!                 Jk = W(k)**2*(M1en_space(k)*(1.0d0+v(k)**2* &
!                      M1eddy_space(k)/X(k)**2)-2.0d0*M1flux_space(k)*v(k)/X(k))

                 !shibata's approximation to above, probably should use above.
                 Jk = 3.0d0/(2.0d0*W(k)**2+1.0d0)*((2.0d0*W(k)**2-1.0d0)*M1en_space(k) - &
                      2.0d0*W(k)**2*v(k)/X(k)*M1flux_space(k))
                 Jkplus1 = 3.0d0/(2.0d0*W(k+1)**2+1.0d0)*((2.0d0*W(k+1)**2-1.0d0)*M1en_space(k+1) - &
                      2.0d0*W(k+1)**2*v(k+1)/X(k+1)*M1flux_space(k+1))

                 if (a_asym.eq.1.0d0) then
                    diffusive_flux = 0.0d0
                 else
                    diffusive_flux = -W(k)/(3.0d0*kappa_inter*X(k)**2)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif
                 
                 if (vp(k)*vm(k+1).gt.0.0d0) then
                    if (vp(k).lt.0.0d0) then
                       advected_energy = 4.0d0*Wm(k+1)**2*vm(k+1)*Xm(k+1)*Jkplus1*onethird
                    else
                       advected_energy = 4.0d0*Wp(k)**2*vp(k)*Xp(k)*Jk*onethird
                    endif
                 else
                    advected_energy = 0.0d0
                 endif
                 
              else

                 Jk = (1.0d0/(1.0d0-v1(k)**2))*(M1en_space(k)*(1.0d0+v1(k)**2* &
                      M1eddy_space(k))-2.0d0*M1flux_space(k)*v1(k))
                 Jkplus1 = (1.0d0/(1.0d0-v1(k+1)**2))*(M1en_space(k+1)*(1.0d0+v1(k+1)**2* &
                      M1eddy_space(k+1))-2.0d0*M1flux_space(k+1)*v1(k+1))

              
                 if (a_asym.eq.1.0d0) then
                    diffusive_flux = 0.0d0
                 else
                    diffusive_flux = -1.0d0/(sqrt(1.0d0-v1(k)**2)*3.0d0*kappa_inter)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif

                 if (v1p(k)*v1m(k+1).gt.0.0d0) then
                    if (v1p(k).lt.0.0d0) then
                       advected_energy = 4.0d0*oneWm**2*v1m(k+1)*Jkplus1*onethird
                    else
                       advected_energy = 4.0d0*oneWp**2*v1p(k)*Jk*onethird
                    endif
                 else
                    advected_energy = 0.0d0
                 endif

              endif
                 
           else if (v_order.eq.0) then
              Jkplus1 = M1en_space(k+1)
              Jk = M1en_space(k)
              
              if (a_asym.eq.1.0d0) then
                 diffusive_flux = 0.0d0
              else
                 if (GR) then
                    diffusive_flux = -1.0d0/(3.0d0*kappa_inter*X(k)**2)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 else
                    diffusive_flux = -1.0d0/(3.0d0*kappa_inter)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif
              endif
                 
              advected_energy = 0.0d0

           else
              stop "add v order"
           endif

           if (activate_turbulence) then
             if (k .lt. ghosts1+1) then
                 diffusive_turb_flux = 0.0d0
             else
                 Lambda_mixp = alpha_turb * pressp(k) / (rhop(k) * dphidr(k))
                 Lambda_mixp = min(Lambda_mixp,x1(k))

                 grad_Enu = (M1en_space(k+1) - M1en_space(k))/(x1(k+1) - x1(k))
                 D_nu_turb = alpha_turb_nu * v_turbp(k) * Lambda_mixp
                 diffusive_turb_flux = - D_nu_turb * grad_Enu
                 ! what happens with GR???? Probably nothing
              endif

              if (a_asym.eq.1.0d0) diffusive_turb_flux = 0.0d0

           else
              diffusive_turb_flux = 0.0d0
           endif
                     
           ! If you use a low number of energy groups this might help the code to not
           ! crush during bounce, since it reduces the free streaming of neutrinos due
           ! to turbulence, which we know isn't very prevalent in the first 20 ms
           !if (bounce .and. time .lt. t_bounce + 0.02d0) then
           !   if ( ABS(diffusive_turb_flux) > ABS(diffusive_flux)*0.2d0 ) &
           !     diffusive_turb_flux = sign(0.2d0,diffusive_turb_flux) * ABS(diffusive_flux)
           !endif

           !flux at interface has two componants, asympotic part, free streaming part.
           M1flux_interface(k,1) = a_asym*( &
                ((l_max(3)*M1flux_space_plus(k)-&
                l_min(3)*M1flux_space_minus(k+1)) + l_max(3)*l_min(3)* &
                (M1en_space_minus(k+1)-M1en_space_plus(k)))/(l_max(3)-l_min(3)) &
                ) + &
                (1.0d0-a_asym)*( &
                diffusive_flux + advected_energy) + &
                !(1.0d0-a_asym) * diffusive_turb_flux
                (1.0d0-a_asym)**4 * diffusive_turb_flux

           !shift interface flux to geometric mean is the difference
           !is too large.  This helps the deleptonziation burst
           !stablely evolve
           if (M1en_space_plus(k)/M1en_space_minus(k+1).gt.10.0d0) then
              limitingflux = sqrt(M1eddy_space_plus(k)*M1en_space_plus(k)* &
                   M1eddy_space_minus(k+1)*M1en_space_minus(k+1))
           else
              limitingflux = (M1eddy_space_plus(k)*M1en_space_plus(k)+ &
                   M1eddy_space_minus(k+1)*M1en_space_minus(k+1))/2.0d0
           endif

           M1flux_interface(k,2) = ( &
                a_asym*((l_max(3)*M1eddy_space_plus(k)*M1en_space_plus(k)- &
                l_min(3)*M1eddy_space_minus(k+1)*M1en_space_minus(k+1)) + &
                l_max(3)*l_min(3)*(M1flux_space_minus(k+1)- &
                M1flux_space_plus(k)))/(l_max(3)-l_min(3)) &
                ) + &
                (1.0d0-a_asym)*( &
                limitingflux &
                )
           
        enddo

        do k=ghosts1+1,M1_imaxradii
           rm = x1i(k)
           rp = x1i(k+1)
           dx = (rp-rm)

           if (GR) then
              M1flux_diff(k,1) = (alpp(k)/Xp(k)**2*x1i(k+1)**2*M1flux_interface(k,1)- &
                   alpm(k)/Xm(k)**2*x1i(k)**2*M1flux_interface(k-1,1))/(dx*x1(k)**2)
              M1flux_diff(k,2) = (alpp(k)/Xp(k)**2*x1i(k+1)**2*M1flux_interface(k,2)- &
                   alpm(k)/Xm(k)**2*x1i(k)**2*M1flux_interface(k-1,2))/(dx*x1(k)**2)
           else
              if (do_effectivepotential) then
                 M1flux_diff(k,1) = (alpp(k)*x1i(k+1)**2*M1flux_interface(k,1)- &
                      alpm(k)*x1i(k)**2*M1flux_interface(k-1,1))/(dx*x1(k)**2)
                 M1flux_diff(k,2) = (alpp(k)*x1i(k+1)**2*M1flux_interface(k,2)- &
                      alpm(k)*x1i(k)**2*M1flux_interface(k-1,2))/(dx*x1(k)**2)
              else
                 M1flux_diff(k,1) = (x1i(k+1)**2*M1flux_interface(k,1)- &
                      x1i(k)**2*M1flux_interface(k-1,1))/(dx*x1(k)**2)
                 M1flux_diff(k,2) = (x1i(k+1)**2*M1flux_interface(k,2)- &
                      x1i(k)**2*M1flux_interface(k-1,2))/(dx*x1(k)**2)
              endif
           endif

           flux_M1(k,i,j,1) = dts*implicit_factor*M1flux_diff(k,1)
           flux_M1(k,i,j,2) = dts*implicit_factor*M1flux_diff(k,2)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO! end do
  
  !let things get away from IC for a few time steps
  if (time.lt.0.000001d0) return

!################################################################
!################################################################
!######################Explicit Energy###########################
!################################################################
!################################################################

  if (include_energycoupling_exp) then
     !$OMP PARALLEL DO PRIVATE(i,j,U_source,S_source,lambda_source)
     do k=ghosts1+1,M1_imaxradii
        if (nt.eq.0) then
           dvdt(k) = 0.0d0
        else
           if (GR) then
              dvdt(k) = (v(k)-v_prev(k))/dts
           else
              dvdt(k) = (v1(k)-v_prev(k))/dts
           endif

        endif
        M1_source_dvdt(k) = dvdt(k)

        do i=1,number_species_to_evolve
           do j=1,number_groups
              U_source(j) = q_M1(k,i,j,1)
              U_source(j+number_groups) = q_M1(k,i,j,2)
           enddo

           call M1_energycoupling_source_rate_bound(k,i,U_source,S_source, &
                lambda_source)
           M1_source_lambda_energycoupling(k,i) = lambda_source

           do j=1,number_groups
              if (S_source(j).ne.S_source(j)) then
                 write(*,*) S_source(j),j,i,k
                 stop "flux NaNing...1"
              else if (S_source(j+number_groups).ne. &
                   S_source(j+number_groups)) then
                 write(*,*) S_source(j+number_groups),j,i,k
                 stop "flux NaNing... 2"
              endif
           enddo

           flux_M1_energy(k,i,:,1) = -dts*implicit_factor* &
                S_source(1:number_groups)
           flux_M1_energy(k,i,:,2) = -dts*implicit_factor* &
                S_source(number_groups+1:2*number_groups)
        enddo
     enddo
     !$OMP END PARALLEL DO! end do
  endif

!################################################################
!################################################################
!###################Explicit Scattering##########################
!################################################################
!################################################################


  if (include_Ielectron_exp) then
     !$OMP PARALLEL DO PRIVATE(i,j,U_source,S_source,lambda_source,alp2,X2)
     do k=ghosts1+1,M1_imaxradii
        if (GR) then
           alp2 = alp(k)**2
           X2 = X(k)**2
        else
           if (do_effectivepotential) then
              alp2 = alp(k)**2
           else
              alp2 = 1.0d0
           endif
           X2 = 1.0d0
        endif

        do i=1,number_species_to_evolve
           do j=1,number_groups
              U_source(j) = q_M1(k,i,j,1)
              U_source(j+number_groups) = q_M1(k,i,j,2)
           enddo

           call M1_ies_source_rate_bound(k,i,U_source,S_source, &
                lambda_source,.true.)
           M1_source_lambda_ies(k,i) = lambda_source

           do j=1,number_groups
              flux_M1_scatter(k,i,j,1) = dts*implicit_factor*S_source(j)
              flux_M1_scatter(k,i,j,2) = dts*implicit_factor* &
                   S_source(j+number_groups)
              ies_sourceterm(k,i,j,1) = S_source(j)/alp2
              ies_sourceterm(k,i,j,2) = S_source(j+number_groups)/X2
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO! end do
  endif

  call M1_source_timestep_cache_update(dts,implicit_factor)

end subroutine M1_explicitterms
