const GAUSS_LEGENDRE_C_MINUS = 1 / 2 - sqrt(3) / 6
const GAUSS_LEGENDRE_C_PLUS  = 1 / 2 + sqrt(3) / 6
const BLANES_GAUSS_LEGENDRE_C = sqrt(3 / 20)

integration_nodes(::SimulationMethod) = (0, 1)
integration_nodes(::BlochMagnusConst1) = (0,)
integration_nodes(::BlochMagnusLin2) = (0, 1)
integration_nodes(::BlochMagnusMid2) = (1//2,)
integration_nodes(::BlochMagnusLinComm2) = (0, 1)
integration_nodes(::BlochMagnusQuad2) = (0, 1//2, 1)
integration_nodes(::BlochMagnusQuad4) = (0, 1//2, 1)
integration_nodes(::BlochMagnusGL2) = (GAUSS_LEGENDRE_C_MINUS, GAUSS_LEGENDRE_C_PLUS)
integration_nodes(::BlochMagnusGL4) = (GAUSS_LEGENDRE_C_MINUS, GAUSS_LEGENDRE_C_PLUS)
integration_nodes(::BlochMagnusBGL4) = (1 / 2 - BLANES_GAUSS_LEGENDRE_C, 1//2, 1 / 2 + BLANES_GAUSS_LEGENDRE_C)
integration_nodes(::BlochMagnusBGL6) = (1 / 2 - BLANES_GAUSS_LEGENDRE_C, 1//2, 1 / 2 + BLANES_GAUSS_LEGENDRE_C)

eval_intervals_per_step(method::SimulationMethod) = max(1, length(integration_nodes(method)) - 1)
eval_intervals_per_step(::BlochMagnusMid2) = 2
eval_intervals_per_step(::BlochMagnusGL2) = 3
eval_intervals_per_step(::BlochMagnusGL4) = 3
eval_intervals_per_step(::BlochMagnusBGL4) = 4
eval_intervals_per_step(::BlochMagnusBGL6) = 4
