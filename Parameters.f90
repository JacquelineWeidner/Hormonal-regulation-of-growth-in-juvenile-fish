
module parameter_values

  implicit none

  ! This is the portable declaration of the real type precision kind, in
  ! terms of the number of decimal digits and the maximum range of the values. 
  ! It should guarantee replicable results on all systems and compilers.
  ! (The code below makes sure the real values have 15 digits of precision and 
  ! can range from 1E-307 to 1E307. This is equivalent to real kind 8, 64 bit.)
  
  integer, parameter, public :: PREC_REAL  = selected_real_kind(15,  307)

  ! Model Resolution
  real(kind=PREC_REAL),  parameter, public :: length_min = 10., length_max = 35., length_step = 1.       ! Max (backward) and min body length of fish; [cm]
  integer, parameter, public               :: L_categories = idnint((length_max-length_min)/length_step) ! Body length categories
  real(kind=PREC_REAL), parameter, public  :: forward_length_max = 30.                                   ! Max body length (= size for maturation used in forward); [cm]
  integer, parameter, public               :: R_categories = 10                                          ! Energy reserve categories
  real(kind=PREC_REAL),  parameter, public :: THF_min = 0., THF_max = 5., THF_step = 0.03125             ! Max and min level of thyroid hormone function (THF) in blood plasma; [ng ml-1]
  integer, parameter, public               :: Th_categories = idnint((THF_max-THF_min)/THF_step)         ! Thyroid hormone function categories
  real(kind=PREC_REAL),  parameter, public :: GHF_min = 0., GHF_max = 200., GHF_step = 1.25              ! Max and min level of growth hormone function (GHF) in blood plasma; [ng ml-1]
  integer, parameter, public               :: G_categories = idnint((GHF_max-GHF_min)/GHF_step)          ! Growth hormone function categories
  real(kind=PREC_REAL),  parameter, public :: OXF_min = 0., OXF_max = 1500., OXF_step = 9.375            ! Max and min level of orexin function (OXF) in blood plasma; [pg ml-1]
  integer, parameter, public               :: Ox_categories = idnint((OXF_max-OXF_min)/OXF_step)         ! Orexin function categories
  integer, parameter, public               :: t_max = 200                                                ! Number of timesteps
  real(kind=PREC_REAL),  parameter, public :: t_duration = 7.                                            ! Number of days in one timestep

  ! PARAMETERS
  real(kind=PREC_REAL),  parameter, public :: k_somatic = 0.85d-5, k_max_reserves = 1.2d-5  ! Fulton's condition factor for fish with empty & full reserves (units kg weight, cm length)
  real(kind=PREC_REAL),  parameter, public :: reserves_energy_density = 5.d6                ! Energy density for reserves; [J kg-1]
  real(kind=PREC_REAL),  parameter, public :: soma_energy_density     = 4.d6                ! Energy density for somatic tissue; [J kg-1]
  real(kind=PREC_REAL),  parameter, public :: SMR_coeff_weight = 0.12                       ! Weight used for correction of SMR_coeff according to new exponent (for 120g-fish); [J timestep-1 kg-beta]
  real(kind=PREC_REAL),  parameter, public :: SMR_exp = 0.7                                 ! Exponent scaling SMR proportional to weight
  real(kind=PREC_REAL),  parameter, public :: THF_SMR_effect = 0.23                         ! Coefficient for increases in SMR between min and max THF
  real(kind=PREC_REAL),  parameter, public :: THF_O2max_effect = 0.24                       ! Coefficient for increases in aerobic scope between min and max THF
  real(kind=PREC_REAL),  parameter, public :: OXF_effect_on_intake = 5.                     ! Coefficient for conversion of OXF level to intake
  real(kind=PREC_REAL),  parameter, public :: growth_max = 0.04*t_duration                  ! Max growth at max GHF as proportion of weight per week
  real(kind=PREC_REAL),  parameter, public :: conversion_efficiency_reserves = 0.85         ! Conversion efficiency for converting metabolites from intake to reserves
  real(kind=PREC_REAL),  parameter, public :: conversion_efficiency_growth = 0.75           ! Conversion efficiency for converting metabolites from reserves to somatic tissue
  real(kind=PREC_REAL),  parameter, public :: SDA_coeff = 0.15                              ! Cost of digestion (proportional to intake)
  real(kind=PREC_REAL),  parameter, public :: foraging_cost_coeff = 0.2                     ! Energetic cost of foraging behaviour
  real(kind=PREC_REAL),  parameter, public :: environment_scaling = 1.8                     ! Coefficient for scaling of food availability 

  ! AMR (Adapted from Claireaux et al. 2000)
  real(kind=PREC_REAL),  parameter, public :: Temperature = 5.                              ! Temperature; [degrees]
  ! Measured aerobic scope at given temperature (from Claireaux et al. 2000, cod, eq. 9, for 100% O2-saturation)
  real(kind=PREC_REAL),  parameter, public :: O2max_coeff = (17.29*Temperature**(-0.015*Temperature+1.062)+30.01) * (1.-exp(-0.035*100+0.34)) * 336.3 * t_duration  
                                                                                            ! [mgO2 h-1 kg-1] in paper converted to [J timestep-1 kg-1]
  real(kind=PREC_REAL),  parameter, public :: O2max_exp = SMR_exp

  !Mortality
  real(kind=PREC_REAL),  parameter, public :: M_size_coeff = 1.3*(t_duration/365.)          ! Coefficient of size-dependent mortality component
  real(kind=PREC_REAL),  parameter, public :: M_size_exp = -0.75                            ! Exponent of size-independent mortality component
  real(kind=PREC_REAL),  parameter, public :: M_sizeindependent = 0.01*(t_duration/365.)    ! Size-independent mortality component; [timestep-1]
  real(kind=PREC_REAL),  parameter, public :: M_foraging_coeff = 0.08, M_foraging_exp = 3.  ! Foraging-related mortality component (coefficient & exponent)
  real(kind=PREC_REAL),  parameter, public :: M_O2_exp = 2.7, M_O2_coeff = 0.8              ! Scope-related mortality component (coefficient & exponent)
  real(kind=PREC_REAL),  parameter, public :: M_interaction = 0.6                           ! Active-while-vulnerable mortality component (coefficient)
  
  !Parameters for FORWARD
  real(kind=PREC_REAL),  parameter, public :: R_categories_init = real(R_categories,PREC_REAL) * 0.1     
																							 ! Initial reserve category in forward simulation
  integer, parameter, public               :: j_max = 1                                      ! Max number of inidviduals in population
  integer, parameter, public               :: length_init = 10                               ! Initial length of all individuals; [cm]

end module parameter_values

