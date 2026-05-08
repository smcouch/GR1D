!-*-f90-*-
!this subroutine updates the neutrino-matter interaction coefficents
subroutine M1_updateeas
  
  use GR1D_module
  use nulibtable
  implicit none

  real*8 :: xrho,xtemp,xye,xeta,xlrho,xltemp,xleta,beta_eta
  integer :: k,i,j

  real*8 :: opacity_spectrum(number_species,number_groups,2)
  real*8 :: singlespecies_opacity_spectrum(number_groups,2)

  real*8 :: inelastic_tempspectrum(number_species,number_groups,number_groups,2)
  real*8 :: singlespecies_inelastic_tempspectrum(number_groups,number_groups,2)

  real*8 :: epannihil_tempspectrum(number_species,number_groups,number_groups,4)
  real*8 :: singlespecies_epannihil_tempspectrum(number_groups,number_groups,4)

  real*8 :: energy_x

  integer :: keytemp,keyerr
  real*8 :: eosdummy(17)

  if (M1_testcase_number.eq.0.or.M1_testcase_number.eq.1) then

     !$OMP PARALLEL DO PRIVATE(xrho,xtemp,xye,xeta,xlrho,xltemp,xleta,beta_eta, &
     !$OMP opacity_spectrum,singlespecies_opacity_spectrum, &
     !$OMP keytemp,keyerr,eosdummy,inelastic_tempspectrum,singlespecies_inelastic_tempspectrum, &
     !$OMP epannihil_tempspectrum,singlespecies_epannihil_tempspectrum,energy_x,i,j)
     do k=2,M1_imaxradii+ghosts1-1
        
        xrho = rho(k)/rho_gf
        xtemp = temp(k)
        xye = ye(k)
        xlrho = log10(xrho)
        xltemp = log10(xtemp)

        if (xlrho.lt.nulibtable_logrho_min) then
           eas(k,:,:,:) = 0.0d0
           if (include_Ielectron_imp.or.include_Ielectron_exp) ies(k,:,:,:,:) = 0.0d0
           if (include_epannihil_kernels) epannihil(k,:,:,:,:) = 0.0d0
           cycle
        endif

        if (xltemp.lt.nulibtable_logtemp_min) then
           stop "M1_update_eas: temp too low"
        endif
        if (xye.lt.nulibtable_ye_min) stop "M1_update_eas: ye too low"
        if (xltemp.gt.nulibtable_logtemp_max) stop "M1_update_eas: temp too high"
        if (xye.gt.nulibtable_ye_max) stop "M1_update_eas: ye too high"

        opacity_spectrum = 0.0d0
        if (number_species_to_evolve.eq.1) then
           call nulibtable_single_species_range_energy_abs_scat(xrho,xtemp,xye,1, &
                singlespecies_opacity_spectrum,number_groups,2)
           opacity_spectrum(1,:,:) = singlespecies_opacity_spectrum(:,:)
        else if (number_species_to_evolve.eq.3) then
           call nulibtable_range_species_range_energy_abs_scat(xrho,xtemp,xye, &
                opacity_spectrum,number_species,number_groups,2)
        else
           stop "set up eas interpolation for this number of species"
        endif

        if (include_epannihil_kernels) then
           opacity_spectrum(3,:,1) = 0.0d0
        endif

        eas(k,:,:,2) = opacity_spectrum(:,:,1)
        eas(k,:,:,3) = opacity_spectrum(:,:,2)

        !Recalculate emissivity from black body using one EOS call per zone.
        keytemp = 1 !keep temperature
        keyerr = 0
#if HAVE_NUC_EOS
        call nuc_eos_full(xrho,xtemp,xye,eosdummy(1),eosdummy(2),eosdummy(3), &
             eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),eosdummy(8), &
             eosdummy(9),eosdummy(10),eosdummy(11),eosdummy(12), &
             eosdummy(13),elechem(k),eosdummy(15),eosdummy(16),eosdummy(17), &
             keytemp,keyerr,eos_rf_prec)
        if(keyerr.ne.0) then
           write(6,*) "############################################"
           write(6,*) "EOS PROBLEM in M1_updateeas.F90:"
           write(6,*) "timestep number: ",nt
           write(6,"(i5,1P3E18.9)") k,xrho,xtemp,xye
           stop "This is bad!"
        endif
#else
        stop "Need nuclear EOS for M1 transport"
#endif

        beta_eta = (elechem(k)-eosdummy(17))/xtemp
        do j=1,number_groups
           energy_x = nulibtable_energies(j)/(nulib_energy_gf*xtemp)
           eas(k,1,j,1) = M1_blackbody_emissivity_factor(j)* &
                opacity_spectrum(1,j,1)/(1.0d0+exp(energy_x-beta_eta))
           eas(k,2,j,1) = M1_blackbody_emissivity_factor(j)* &
                opacity_spectrum(2,j,1)/(1.0d0+exp(energy_x+beta_eta))
           eas(k,3,j,1) = 4.0d0*M1_blackbody_emissivity_factor(j)* &
                opacity_spectrum(3,j,1)/(1.0d0+exp(energy_x))
        enddo

        xeta = elechem(k)/xtemp
        if (include_Ielectron_exp.or.include_Ielectron_imp.or.include_epannihil_kernels) &
             xleta = log10(xeta)

        if (include_Ielectron_imp.or.include_Ielectron_exp) then
           inelastic_tempspectrum = 0.0d0
           if (xltemp.lt.nulibtable_logItemp_min) stop "M1_update_eas: Itemp too low"
           if (xleta.lt.nulibtable_logIeta_min) then
              write(*,*) xrho,xtemp,xye,xeta
              stop "M1_update_eas: Ieta too low"
           endif
           if (xltemp.gt.nulibtable_logItemp_max) stop "M1_update_eas: temp too high"
           if (xleta.gt.nulibtable_logIeta_max) then
              write(*,*) xrho,xtemp,xye,xeta,nulibtable_logIeta_max,k
              stop "M1_update_eas: Ieta too high"
           endif

           if (number_species_to_evolve.eq.1) then
              call nulibtable_inelastic_single_species_range_energy2(xtemp,xeta,1, &
                   singlespecies_inelastic_tempspectrum,number_groups,number_groups,2)
              inelastic_tempspectrum(1,:,:,:) = singlespecies_inelastic_tempspectrum(:,:,:)
           else if (number_species_to_evolve.eq.3) then
              call nulibtable_inelastic_range_species_range_energy2(xtemp,xeta, &
                   inelastic_tempspectrum,number_species,number_groups,number_groups,2)
           else
              stop "set up eas interpolation for this number of species"
           endif
           
           ies(k,:,:,:,:) = inelastic_tempspectrum(:,:,:,:)

        endif

        if (include_epannihil_kernels) then
           epannihil_tempspectrum = 0.0d0
           if (xltemp.lt.nulibtable_logItemp_min) stop "M1_update_eas: Itemp too low"
           if (xleta.lt.nulibtable_logIeta_min) then
              write(*,*) xrho,xtemp,xye,xeta
              stop "M1_update_eas: Ieta too low"
           endif
           if (xltemp.gt.nulibtable_logItemp_max) stop "M1_update_eas: temp too high"
           if (xleta.gt.nulibtable_logIeta_max) then
              write(*,*) xrho,xtemp,xye,xeta,nulibtable_logIeta_max,k
              stop "M1_update_eas: Ieta too high"
           endif

           i = 3
           call nulibtable_epannihil_single_species_range_energy2(xtemp,xeta,i, &
                singlespecies_epannihil_tempspectrum,number_groups,number_groups,4)
           epannihil_tempspectrum(i,:,:,:) = singlespecies_epannihil_tempspectrum(:,:,:)
           
           epannihil(k,:,:,:,:) = epannihil_tempspectrum(:,:,:,:)

        endif

     enddo
     !$OMP END PARALLEL DO! end do

  else if (M1_testcase_number.ge.2.and.M1_testcase_number.le.9) then
     !this is taken care of
  else
     stop "add in eas updating code for this test case"
  endif

end subroutine M1_updateeas
