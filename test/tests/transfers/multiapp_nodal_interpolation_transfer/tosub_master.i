[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = .21
  xmax = .79
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./disp_x]
    initial_condition = 0.4
  [../]
  [./disp_y]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  petsc_options = '-snes_mf_operator -ksp_monitor'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Output]
  output_initial = true
  exodus = true
  perf_log = true
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = MooseTestApp
    execute_on = timestep
    positions = '0.2 0 0'
    input_files = tosub_sub.i
  [../]
[]

[Transfers]
  [./tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = from_master
  [../]
  [./elemental_tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = elemental_from_master
  [../]
  [./radial_tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = radial_from_master
    interp_type = radial_basis
  [../]
  [./radial_elemental_tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = radial_elemental_from_master
    interp_type = radial_basis
  [../]
  [./displaced_target_tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = displaced_target_from_master
    displaced_target_mesh = true
  [../]
  [./displaced_source_tosub]
    type = MultiAppNodalInterpolationTransfer
    direction = to_multiapp
    execute_on = timestep
    multi_app = sub
    source_variable = u
    variable = displaced_source_from_master
    displaced_source_mesh = true
  [../]
[]

