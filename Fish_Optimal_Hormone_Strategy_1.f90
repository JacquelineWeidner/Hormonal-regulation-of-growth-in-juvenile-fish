!*******************************************************************************
! SVN Version ID: $Id: Hormones_v2.f90 3114 2017-03-23 09:58:15Z jwe041 $
!*******************************************************************************

    
! In this model we combine physiological, environmental and evolutionary aspects of fish growth in a state-dependent model where
! the optimal regulation of growth and survival is achieved through hormonal regulation of behaviour. 
! The model has a two-part structure. In the first part (backward), the algorithm optimizes a strategy that can be considered the 
! evolutionary adaptation to a certain environment. In the second part (forward), we investigate this adaptation by simulating 
! individuals that live in the given environment and use the calculated optimal policy. We record their trajectory of growth, hormone 
! expression and individual states. 
    
    
program Hormones_v2

  ! Module for printing results to binary files
  use binary_write
  
  ! Module containing all parameter values
  use parameter_values

  implicit none

  ! VARIABLES
  integer t, L, R, Th, G, Ox                                          ! Variables used in loops (time, length, reserves, THF, GHF, OXF)
  real(kind=PREC_REAL)  reserves, reserves_new                        ! Variables for reserves; [J]
  real(kind=PREC_REAL)  reserves_max, reserves_max_new                ! Variables for max reserves relative to length; [J]
  real(kind=PREC_REAL)  surplus_before_growth, surplus_after_growth   ! Variables for energy surplus; [J]
  real(kind=PREC_REAL)  weight, weight_somatic, weight_somatic_new    ! Variables for weight; [kg]
  real(kind=PREC_REAL)  length, length_new                            ! Variables for length; [cm]
  real(kind=PREC_REAL)  chance_of_max_length                          ! Chance of reaching max length
  real(kind=PREC_REAL)  growth                                        ! Growth defined as weight of new soma; [kg]
  real(kind=PREC_REAL)  growth_cost                                   ! Energetic cost of growing; [J]   
  real(kind=PREC_REAL)  THF, GHF, OXF                                 ! Hormone function levels
  real(kind=PREC_REAL)  intake                                        ! Intake; [J timestep -1]
  real(kind=PREC_REAL)  target_intake                                 ! Target intake set by OXF (multiples of SMR_std)  
  real(kind=PREC_REAL)  foraging_required                             ! Foraging behaviour (multiples of SMR_std)
  real(kind=PREC_REAL)  foraging_cost                                 ! Energetic cost of foraging behaviour; [J timestep -1]
  real(kind=PREC_REAL)  conversion_cost_via_reserves                  ! Cost of converting metabolites from food to reserves; [J timestep -1]
  real(kind=PREC_REAL)  conversion_cost_to_growth                     ! Cost of converting metabolites stored as reserves to new soma; [J timestep -1]
  real(kind=PREC_REAL)  L_new, R_new, L_cat, R_cat                    ! L and R category for new length and reserves (continuous)
  real(kind=PREC_REAL)  dL, dR                                        ! Integer part of new length and reserves category...
  integer intL, intR                                                  ! ...and decimal part
  real(kind=PREC_REAL)  SMR_std, SMR_somatic_std, SMR_exp_CJ, SMR_coeff_CJ, SMR_coeff, SMR_THF, SDA  
                                                                      ! Variables for calculation of metabolic rate and cost of digestion (SDA)
  real(kind=PREC_REAL)  O2max_std, O2max_THF, O2_used                 ! Variables used to calculate oxygen budget (aerobic metabolic rate)
  real(kind=PREC_REAL)  M_size, M_foraging, survival, M_O2            ! Variables used for calculation of mortality
  
  ! MATRICES
  integer, dimension(1:3, 1:(R_categories+1),  1:(L_categories+1), 1:t_max)                           :: strategy
  real(kind=PREC_REAL),  dimension(     1:(R_categories+1),  1:(L_categories+1), 1:t_max)             :: fitness
  real(kind=PREC_REAL),  dimension(     1:(Ox_categories+1), 1:(G_categories+1), 1:(Th_categories+1)) :: fitnessconsequences
  
  ! Variables / matrix: forward simulation 
  integer :: j
  real(kind=PREC_REAL), dimension(1:27, 1:t_max+1, 1:j_max) :: ind
  
  
  ! SMR (based on Clarke & Johnston 1999 - General teleost fish)
  SMR_exp_CJ = 0.80
  SMR_coeff_CJ = exp(-5.43)                                                 ![mmol O2 g-1 h-1], originally -5.43
  SMR_coeff_CJ = SMR_coeff_CJ * 434. * 24.                                  ![J g-1 d-1], converted from paper, p 895, units conversion
  SMR_coeff_CJ = SMR_coeff_CJ * (0.001**(-SMR_exp_CJ)) * t_duration         ![J kg-1 timestep-1], units conversion
  SMR_coeff = SMR_coeff_CJ * (SMR_coeff_weight**(SMR_exp_CJ-SMR_exp))       !Correction acording to new exponent, converted for a 120g-fish
    
  
  ! BACKWARDS OPTIMIZATION
  ! Initialise fitness: 1 at L_max, otherwise 0.
  fitness(:,:,:) = 0.
  fitness(:,L_categories+1,:) = 1.
  print '(4a10)', '  Timestep', '    Orexin', '    Growth', '   Thyroid'
  
  TIME_LOOP: do t = t_max-1, 1, -1
  
    !STATES
    LENGTH_LOOP: do L = 1, L_categories                                                     
      length = length_min + (L-1) * length_step                                             ! Calculate body length; [cm]
      weight_somatic = k_somatic*(length**3)                                                ! Converts length into somatic weight, reserves not considered; [kg]                                                                          
      reserves_max = (k_max_reserves*(length**3) - weight_somatic)*reserves_energy_density  ! Max amount of energy a fish of given size can store; [J]    
      M_size = M_size_coeff*(length**M_size_exp)                                            ! Size dependent mortality rate; [timestep-1]
      O2max_std = O2max_coeff*(weight_somatic**O2max_exp)                                   ! Max aerobic scope (AMR); [J timestep-1] 
      SMR_somatic_std = SMR_coeff * (weight_somatic**SMR_exp)                               ! Baseline SMR based on somatic weight; [J timestep-1]
      
      ! RESERVES
      RESERVES_LOOP: do R = 1, R_categories + 1                                             
        reserves = reserves_max * (dble(R-1)/dble(R_categories))                            ! Calculate energy content in reserves; [J]
        weight = weight_somatic + reserves/reserves_energy_density                          ! Calculates weight including soma & reserves; [kg]
        SMR_std = SMR_coeff * (weight**SMR_exp)                                             ! Baseline SMR based on total weight; [J timestep-1]

        ! DECISIONS
        THF_LOOP: do Th = 1, Th_categories + 1
          THF = THF_min + (Th-1)*THF_step                                                   ! Converts THF category to blood plasma concentration; [ng ml-1] 
          SMR_THF = SMR_std*(1.+((THF/THF_max)-0.5)*THF_SMR_effect)                         ! Caluclates new SMR under THF influence; [J timestep-1]
          O2max_THF = O2max_std*(1.+((THF/THF_max)-0.5)*THF_O2max_effect)                   ! Calculates new AMR under THF influence; [J timestep-1]
          

          GHF_LOOP: do G = 1, G_categories + 1
            GHF = GHF_min + (G-1)*GHF_step                                                  ! Converts GHF category to blood plasma concentration; [ng ml-1]
            growth = (GHF/GHF_max) * growth_max * weight_somatic                            ! Weight of new soma; [kg] 
            growth_cost = growth*soma_energy_density                                        ! Energetic cost of new soma; [J]
            weight_somatic_new = weight_somatic + growth                                    ! New somatic weight based on growth; [kg]
            length_new = (weight_somatic_new/k_somatic)**(1./3.)                            ! New length; [cm]
            L_new = ((length_new- length_min)/length_step) + 1                              ! L category for new length, continuous
            intL = min(idint(L_new),L_categories)                                           ! Integer part of new length category...
            dL = max(0.,min(L_new - dble(intL),1.))                                         ! ...and decimal part
            reserves_max_new = (k_max_reserves*(length_new**3) - weight_somatic_new)*reserves_energy_density
                      
            OXF_LOOP: do Ox = 1, Ox_categories + 1
              OXF = OXF_min + (Ox-1)*OXF_step                                               ! Converts OXF category to blood plasma cocentration; [pg ml-1] 
              target_intake = (OXF/OXF_max)*OXF_effect_on_intake                            ! Amount of energy that individual has to take up in current time step 
              intake = target_intake*SMR_somatic_std                                        ! Energy from target intake referes to as multiples of SMR; [J timestep-1]
              foraging_required = target_intake/environment_scaling                         ! Required foraging intensity given the amount of energy needed and the food availability in the environment
              foraging_cost = foraging_cost_coeff*foraging_required*SMR_std                 ! Cost of foraging; [J timestep-1]
              SDA = SDA_coeff*intake                                                        ! Cost of digestion; [J timestep-1]
              surplus_before_growth = intake - SDA - SMR_THF - foraging_cost                ! Surplus energy avilable for growth, storing, etc.; [J timestep-1]
                        
              ! Conversion cost = conversion loss to intermediate metabolites 
              if (surplus_before_growth > 0.) then                                          ! ... if more energy is taken up than used (pos surplus)
                 conversion_cost_via_reserves = surplus_before_growth * (1.-conversion_efficiency_reserves)        !Conversion cost from intake to reserves
              else    ! ... if more energy is used than taken up (neg surplus). Energy has to be drained from reserves to cover energetic costs of metabolic expenses.
                 conversion_cost_via_reserves =  abs(surplus_before_growth) / conversion_efficiency_reserves * (1.-conversion_efficiency_reserves)    !Conversion cost from reserves to metabolism                                                          !No positive intake that needs to be converted to reserves
              endif
             
              conversion_cost_to_growth = (growth_cost/conversion_efficiency_growth)*(1.-conversion_efficiency_growth)  ! Cost of converting intermediate metabolites to new soma
              reserves_new = min(reserves + surplus_before_growth - conversion_cost_via_reserves - growth_cost - conversion_cost_to_growth, reserves_max_new)
                
              R_new = (reserves_new*R_categories/reserves_max_new) + 1                      ! R category for new reserves (continuous)
              intR = max(1,min(idint(R_new),R_categories))                                  ! Integer part of new reserves category...
              dR = max(0.,min(R_new - dble(intR),1.))                                       ! ...and decimal part
              O2_used = SMR_THF + SDA + foraging_cost + conversion_cost_via_reserves + conversion_cost_to_growth
              
              ! MORTALITY
              M_foraging = (M_foraging_coeff * (foraging_required**M_foraging_exp))  
                                                                                            ! Foraging-related mortality component
              M_O2 = M_O2_coeff * (O2_used/O2max_THF)**M_O2_exp                             ! Scope-related mortality component
              survival = exp(- M_sizeindependent - M_size - (M_size*M_O2) - (M_size*M_foraging) - M_interaction*(M_size*M_O2*M_foraging))  
                                                                                            ! Probability of surviving the current timestep
                            
              if (reserves_new < 0) then                                                    ! Starvation risk (If reserves is less than 0, survival is set to 0)
                survival = 0.       
              endif
 
              ! Interpolation
              fitnessconsequences(Ox,G,Th) = survival * ( &
                (1.-dR)*(1.-dL)*fitness(intR  ,intL  ,t+1) + (   dR)*(1.-dL)*fitness(intR+1,intL  ,t+1) + &
                (1.-dR)*(   dL)*fitness(intR  ,intL+1,t+1) + (   dR)*(   dL)*fitness(intR+1,intL+1,t+1))
              
            enddo OXF_LOOP

          enddo GHF_LOOP

        enddo THF_LOOP

        ! OPTIMIZE strategy for this state combination
        strategy(1:3,R,L,t) = maxloc(fitnessconsequences(:,:,:))                            ! Finds hormonstrategy that maxmimizes fitness
        fitness(R,L,t) = fitnessconsequences(strategy(1,R,L,t), strategy(2,R,L,t), strategy(3,R,L,t))   
                                                                                            ! Stores hormone categories giving optimal fitness for current reserve, length & time step
      
      enddo RESERVES_LOOP

    enddo LENGTH_LOOP
    
    print '(4i10)', t, sum(abs(strategy(1, :,:,t)-strategy(1, :,:,t+1))), sum(abs(strategy(2, :,:,t)-strategy(2, :,:,t+1))), sum(abs(strategy(3, :,:,t)-strategy(3, :,:,t+1)))
    
  enddo TIME_LOOP
  
  ! Save fitness, fitnessconsequences, strategy to file
  call matrix_to_binary(fitness, "fitness.bin")
  call matrix_to_binary(fitnessconsequences, "fitnessconsequences.bin")
  call matrix_to_binary(strategy, "strategy.bin")

  
  ! FORWARD SIMULATION
 
  ind(:,:,:) = -1000
  j = 1
  t = 1
  
  ! Initiation
  ind(1,1,:) = length_init                                                                  ! Initial length; [cm]
  ind(2,1,:) = R_categories_init                                                            ! Initial reserve category (~R)
  
  ind(4,1,:) = k_somatic * (ind(1,1,:)**3)                                                  ! Initial somatic weight; [kg]
  ind(3,1,:) = ((k_max_reserves * (ind(1,1,:)**3) - ind(4,1,:)) * reserves_energy_density) * ((ind(2,1,:)-1.) / dble(R_categories))     
																							! Initial reserves; [J]
  ind(5,1,:) = ind(4,1,:) + ind(3,1,:) / reserves_energy_density                            ! Initial weight; [kg]
  ind(23,1,:) = 1.                                                                          ! Survival 
  ind(25,1,:) = (k_max_reserves * (ind(1,1,:)**3) - ind(4,1,:)) * reserves_energy_density   ! Max reserves an individual with given size can have; [J]
 
  ! INDIVIDES (running over all individuals in the population)
  Individ_f_LOOP: do j = 1, j_max, 1                            
      
    length = ind(1,1,j)                                                                     
    weight_somatic = ind(4,1,j)                                                             
    reserves = ind(3,1,j)                                                                   
    weight = ind(5,1,j)                                                                     
    survival = ind(23,1,j)
    reserves_max = ind(25,1,j)                                                              
    
    ! TIME (running from 1 to last timestep)
    Time_f_LOOP: do t = 1, t_max, 1                                                       
                 
      M_size = M_size_coeff * (length**M_size_exp)                                          
      O2max_std = O2max_coeff * (weight_somatic**O2max_exp)                                 
      SMR_somatic_std = SMR_coeff * (weight_somatic**SMR_exp)                               
      SMR_std = SMR_coeff * (weight**SMR_exp)                                               
       
      L_cat = ((length - length_min) / length_step) + 1                                               
      intL = min(idint(L_cat),L_categories)                                                 
      dL = max(0.,min(L_cat - dble(intL),1.))                                               
           
      R_cat = (reserves * R_categories / reserves_max) + 1                                  
      intR = max(1,min(idint(R_cat),R_categories))                                          
      dR = max(0.,min(R_cat - dble(intR),1.))                                                
           
      ! Interpolation of hormoncategories from strategy matrix
            OXF = &
                (1. - dR) * (1. - dL) * strategy(1, intR    , intL    , 1) + &
                (     dR) * (1. - dL) * strategy(1, intR + 1, intL    , 1) + & 
                (1. - dR) * (     dL) * strategy(1, intR    , intL + 1, 1) + &
                (     dR) * (     dL) * strategy(1, intR + 1, intL + 1, 1) 
      
            GHF = &
                (1. - dR) * (1. - dL) * strategy(2, intR    , intL    , 1) + &
                (     dR) * (1. - dL) * strategy(2, intR + 1, intL    , 1) + & 
                (1. - dR) * (     dL) * strategy(2, intR    , intL + 1, 1) + &
                (     dR) * (     dL) * strategy(2, intR + 1, intL + 1, 1) 
           
            THF = &
                (1. - dR) * (1. - dL) * strategy(3, intR    , intL    , 1) + &
                (     dR) * (1. - dL) * strategy(3, intR + 1, intL    , 1) + & 
                (1. - dR) * (     dL) * strategy(3, intR    , intL + 1, 1) + &
                (     dR) * (     dL) * strategy(3, intR + 1, intL + 1, 1) 
                
      !Convert from category to hormone concentration
      OXF = OXF_min + (OXF - 1.) * OXF_step          
      GHF = GHF_min + (GHF - 1.) * GHF_step    
      THF = THF_min + (THF - 1.) * THF_step   
     
      !Calculate THF - dependent variables
      SMR_THF = SMR_std * (1. + ((THF / THF_max) - 0.5) * THF_SMR_effect)                   
      O2max_THF = O2max_std * (1. + ((THF / THF_max) - 0.5) * THF_O2max_effect)             
          
      !Calculate GHF - dependent variables
      growth = (GHF / GHF_max) * growth_max * weight_somatic                                
      growth_cost = growth * soma_energy_density                                            
      weight_somatic = weight_somatic + growth                                              
      length = (weight_somatic / k_somatic) ** (1./3.)                                      
      reserves_max = (k_max_reserves * (length**3) - weight_somatic) * reserves_energy_density 
                  
      !Calculate OXF - dependent variables
      target_intake = (OXF/OXF_max)*OXF_effect_on_intake                                    
      intake = target_intake*SMR_somatic_std                                                
      foraging_required = target_intake/environment_scaling                                 
      foraging_cost = foraging_cost_coeff*foraging_required*SMR_std                         
      SDA = SDA_coeff * intake                                                              
      surplus_before_growth = intake - SDA - SMR_THF - foraging_cost                                
      
      ! CONVERSION COSTS
      if (surplus_before_growth > 0.) then
         conversion_cost_via_reserves = surplus_before_growth * (1.-conversion_efficiency_reserves)        
																							! Conversion cost from intake to reserves
      else
         conversion_cost_via_reserves =  abs(surplus_before_growth) / conversion_efficiency_reserves * (1.-conversion_efficiency_reserves)    
                                                                                            ! Conversion cost from reserves to metabolism    
	  endif
      
      conversion_cost_to_growth = (growth_cost/conversion_efficiency_growth)*(1.-conversion_efficiency_growth)  
                                                                                            ! Conversion cost from reserves to growth
      reserves = min(reserves + surplus_before_growth - conversion_cost_via_reserves - growth_cost - conversion_cost_to_growth, reserves_max)
            
      weight = weight_somatic + reserves / reserves_energy_density                          
      O2_used = SMR_THF + SDA + foraging_cost + conversion_cost_via_reserves + conversion_cost_to_growth  
      
      !Mortality 
      M_foraging = (M_foraging_coeff * (foraging_required**M_foraging_exp))  
      M_O2 = M_O2_coeff * (O2_used/O2max_THF)**M_O2_exp                                      
      survival = exp(- M_sizeindependent - M_size - (M_size*M_O2) - (M_size*M_foraging) - M_interaction*(M_size*M_O2*M_foraging))  
            
      ! If reserves are less than -100 J the individual dies (survival = 0)
      if (reserves <= -100.) then
        survival = 0. 
      endif
      
      ! If reserves was less than 0, set them to 0,
      ! This should help the individuals that make small mistakes in forward
      reserves = max(0.,reserves)
             
      ! Calculating chance of reaching max length
      L_cat = ((length - length_min) / length_step) + 1                                     ! find length category from length          
      intL = min(idint(L_cat),L_categories)                                                 ! ... integer part
      dL = max(0.,min(L_cat - dble(intL),1.))                                               ! ... decimal part
      R_cat = (reserves * R_categories / reserves_max) + 1                                  ! find reserve category from reserves
      intR = max(1,min(idint(R_cat),R_categories))                                          ! ... integer part
      dR = max(0.,min(R_cat - dble(intR),1.))                                               ! ... decimal part 
      chance_of_max_length = &                                                              ! chance of reaching max length 
          (1.-dR)*(1.-dL)*fitness(intR  , intL  , t) + &
          (1.-dR)*(   dL)*fitness(intR  , intL+1, t) + &
          (   dR)*(1.-dL)*fitness(intR+1, intL  , t) + &
          (   dR)*(   dL)*fitness(intR+1, intL+1, t) 
     
      ! SAVE VARIABLES FROM END OF TIME STEP
      ind(1,t+1,j) = length                                                                 ! Length; [cm]
      ind(2,t+1,j) = R_cat                                                                  ! Reserves; [category]
      ind(3,t+1,j) = reserves                                                               ! Reserves; [J]
      ind(4,t+1,j) = weight_somatic                                                         ! Somatic weight; [kg]
      ind(5,t+1,j) = weight                                                                 ! Total weight (soma & reserves); [kg]
      ind(23,t+1,j) = ind(23,t,j)*survival                                                  ! Probability of being alive at the beginning of next timestep
      ind(25,t+1,j) = reserves_max                                                          ! Max reserves and individual with a given size could have; [J]

      !Information stored from this timestep            
      ind(6,t,j) = OXF                                                                      ! OXF; [category]
      ind(7,t,j) = GHF                                                                      ! GHF; [category]
      ind(8,t,j) = THF                                                                      ! THF; [category]
      ind(9,t,j) = SMR_std                                                                  ! SMR without THF influence; [J]
      ind(10,t,j) = SMR_THF                                                                 ! SMR influenced by THF; [J]
      ind(11,t,j) = growth                                                                  ! New somatic weight; [kg]
      ind(12,t,j) = target_intake                                                           ! Proxy for energetic demand in current time step; [diensionless]
      ind(13,t,j) = foraging_cost                                                           ! Energetic cost of foraging behaviour; [J]
      ind(14,t,j) = intake                                                                  ! Metabolizable energy an individual takes up by foraging; [J]
      ind(15,t,j) = SDA                                                                     ! Cost of digestion; [J]
      ind(16,t,j) = conversion_cost_via_reserves + conversion_cost_to_growth                ! Conversion costs from food to reserves and from reserves to growth; [J]
      ind(17,t,j) = O2_used                                                                 ! Sum of oxygen used for metabolic processes & activity; [J]
      ind(18,t,j) = O2max_THF                                                               ! AMR; [J]
      ind(19,t,j) = chance_of_max_length                                                    ! Chance of reaching max length
      ind(20,t,j) = M_size * (365./t_duration)                                              ! Size-dependent mortality component; [year-1]
      ind(21,t,j) = M_O2_coeff * M_size * M_O2  * (365./t_duration)                         ! Scope-related mortality component; [year-1]
      ind(22,t,j) = M_size * M_foraging * (365./t_duration)                                 ! Foraging-related mortality component; [year-1]
      ind(24,t,j) = survival                                                                ! Survival; [timestep-1]
      ind(26,t,j) = M_sizeindependent * (365./t_duration)                                   ! Size-independent mortality component; [year-1]
      ind(27,t,j) = (M_interaction * M_size * M_O2 * M_foraging) * (365./t_duration)        ! Active-while-vulnerable mortality component; [year-1]
          
      
      if (survival == 0) exit Time_f_LOOP                                                   ! Exit Time_LOOP if individual is dead. Start with new individual.
      if (length >= forward_length_max) exit Time_f_LOOP                                    ! Exit Time_LOOP if individual has reached max length.
            
    enddo Time_f_LOOP
  
  enddo Individ_f_LOOP
  
  ! Write ind-matrix to binary file 
  call matrix_to_binary(ind, "ind.bin")

end program Hormones_v2