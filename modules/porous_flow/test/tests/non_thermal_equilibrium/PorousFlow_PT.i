# the example a flow through porous media with non-thermal equilibrium between the fluid
# and the matrix

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 10
  ny   = 10
  xmin = 0.0
  xmax = 1.0
  ymax = 1.0
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator 
    porous_flow_vars = 'pp T_f T_m'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
[]

[Variables]
  [./T_m]
    initial_condition = 1
  [../]
  
  [./T_f]
    initial_condition = 10
  [../]

  [./pp]
    initial_condition = 1
  [../]
[]

[Kernels]
  [./p_mass_time_derivative]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  [../]

  [./p_advection]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable  = pp
    gravity = '0 0 0'
    use_displaced_mesh = false
  [../]
  ###### Energy kernels for matrix ####
  [./matrix_heat_time_derivative]
    type = PorousFlowFluidEnergyTimeDerivative
    variable = T_m
    use_displaced_mesh = false
    base_name = matrix
  [../]
 
  [./matrix_heat_conduction]
    type = PorousFlowMatrixHeatConduction
    variable = T_m
    base_name = matrix
  [../]
  
  [./fluid_matrix_heat_exchange]
    type = PorousFlowFractureMatrixHeatExchange
    variable = T_m
    couple_temperature = T_f
    heat_change_type = FluidToMatrix
    heat_convection_coef = 1
  [../]
  ####### Energy Kernels for fluid #####
  [./fluid_heat_time_derivative]
    type = PorousFlowEnergyTimeDerivative
    variable = T_f
    use_displaced_mesh = false
    base_name = fluid
  [../]

  [./fluid_heat_conduction]
    type = PorousFlowFluidHeatConduction
    variable = T_f
    base_name = fluid
  [../]

  [./fluid_heat_convection]
    type = PorousFlowHeatAdvection 
    variable = T_f
    base_name = fluid
  [../]
  
  [./matrix_fluid_heat_exchange]
    type = PorousFlowFractureMatrixHeatExchange
    variable = T_f
    couple_temperature = T_m
    heat_convection_coef = 1
    heat_change_type = MatrixToFluid
  [../]
[]

[BCs]
  [./p_source]
    type = NeumannBC
    variable = pp
    boundary = left
    value = 1
  [../]

  [./T_fluid_soure]
    type = NeumannBC
    boundary = left
    variable = T_f
    value = 1
  [../]
[]

[Modules]
  [./FluidProperties]
   # [./water97property]
   #   type = SimpleFluidProperties
   #   bulk_modulus = 2e9
   #   density0 = 1000
   #   thermal_expansion = 0
   #   viscosity = 1e-3
   # [../]
    [./water97property]
      type = SimpleFluidProperties #Water97FluidProperties
     # bulk_modulus = 1E12
    [../]
  [../]
[]


[Materials]
  [./temperature_matrix_nodal]
    type = PorousFlowTemperature
    at_nodes = true
    base_name = matrix
    temperature = T_m
  [../]

  [./temperature_matrix_qp]
    type = PorousFlowTemperature
    at_nodes = false
    base_name = matrix
    temperature = T_m
  [../]

  [./temperature_fluid_nodal]
    type = PorousFlowTemperature
    at_nodes = true
    base_name = fluid
    temperature = T_f
  [../]

  [./temperature_fluid_qp]
    type = PorousFlowTemperature
    base_name = fluid
    at_nodes = false
    temperature = T_f
  [../]

  [./thermal_conductivity]
    type = PorousFlowThermalConductivityTensor
    modify = false
    rock_thermal_conductivity = '2.58 0 0 0 2.58 0 0 0 0'
    fluid_thermal_conductivity = '1.0 0 0 0 1 0 0 0 0'
  [../]

  [./rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1.08E3
    base_name = matrix
    density = 2.7E3
  [../]

  [./ppss]
    type = PorousFlow1PhaseFullySaturated
    at_nodes = true
    porepressure = pp
    capillary_pressure = pc
  [../]

  [./ppss_qp]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
    capillary_pressure = pc
  [../]

  [./water97property]
    type = PorousFlowSingleComponentFluid
    fp = water97property
    phase = 0
    at_nodes = true
  [../]

  [./water97property_qp]
    type = PorousFlowSingleComponentFluid
    fp = water97property
    phase = 0
    at_nodes = false
  [../]

  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
  [../]

  [./poro]
    type = PorousFlowPorosityConst
    porosity = 0.05
    at_nodes = false
  [../]

  [./poro_nodal]
    at_nodes = true 
    type = PorousFlowPorosityConst
    porosity = 0.05
  [../]

  [./diff1]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1e-9'
    tortuosity = 1.0
  [../]


  [./permeability_fracture]
   # type = PorousFlowPermeabilityConst
   # permeability = '1.05e-8 0 0 0 1.05e-8 0 0 0 0'   # 1.8e-11 = a * kf
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    PorousFlowDictator = dictator
    m = 2
    phi0 = 0.25
    k0 = 1.1E-8
    n = 3
    outputs = csv
  [../]

  [./relp]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  [../]

  [./relp_nodal]
    type = PorousFlowRelativePermeabilityConst
    at_nodes = true
    phase = 0
  [../]

 # [./vol_strain]
 #   type = PorousFlowVolumetricStrain
 # [../]
[]

[Preconditioning]
  [./basic]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
  [../]
[]


[Outputs]
  [./csv]
    type = CSV
    file_base = perm_fracture_output
  [../]

  [./exodus]
    type = Exodus
    file_base = non_thermal_equib_PT_test
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type

                         -ksp_gmres_restart -snes_ls

                         -pc_hypre_boomeramg_strong_threshold'

  petsc_options_value = 'hypre boomeramg 201 cubic 0.7'

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1
  [../]
    end_time = 10 # 9 years
    l_tol = 1e-7
    l_max_its = 500
    nl_max_its = 300
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-4
[]

