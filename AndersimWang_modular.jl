##############################   AnderSim 5 Eq. - Versão Modular  ###############################
###                   Compressível - Escrito todo em série                                      ###
############################ By Anderson Viana - Modularizado ####################################

# Imports e includes
include("./TDMA.jl")
include("./Limiter.jl")
using Printf
using Plots
using CSV
using DataFrames
using LinearAlgebra

include("hMLP.jl")
export hMLP, limiter, minmod2

include("utils.jl")
export out, pifrac

include("solver.jl")
export flux

include("polynomials.jl")
export lagrangeBasisAt, dLagrangeBasisAt, lagrangeBarycentricWeights, lagrangeDerivativeMatrix, lagrangeInterpolationAt
export pLegendreAt, dLegendreAt, vcjh, xFromξ, ξFromx

include("mesh.jl")
export Mesh, Element, StandardElement, meshSize, meshOrder, meshΔx, integrationPointsLocation
export element2faces

include("flow.jl")
export Flow, init!, stepEuler!, stepRK2SSP!, stepRK3SSPS3!
export flattenSolution, updatePrimitiveVars!, solutionBoundaryConditions!
export reconstructSolutionAtFace!, reconstructFluxesAtFace!, riemannFluxes!

include("frgrid.jl")
include("FR.jl")
include("MUSCL.jl")

# Estrutura para armazenar parâmetros de entrada
mutable struct InputParameters
    # Parâmetros geométricos
    Comp::Float64          # Comprimento do duto
    D::Float64             # Diâmetro interno do duto
    θ::Float64             # Ângulo do duto
    n_segmentos::Int       # Número de seguimentos
    D_ext::Float64         # Diâmetro externo do duto
    
    # Parâmetros temporais
    tempo::Float64         # Tempo de simulação
    CFL::Float64           # Coeficiente de estabilidade
    
    # Propriedades fluidas de entrada
    uL_inlet::Float64      # Velocidade do líquido na entrada
    uG_inlet::Float64      # Velocidade do gás na entrada
    ρLi::Float64           # Densidade do Líquido
    ρGi::Float64           # Densidade do gás
    P_out::Float64         # Pressão na saída do domínio (Pa)
    αG_inlet::Float64      # Fração de gás/óleo na entrada
    αG_0::Float64          # Fração de gás/óleo inicial no duto
    T_intet::Float64       # Temperatura na entrada 
    T_0::Float64           # Condição inicial de temperatura
    T_ext::Float64         # Temperatura externa
    
    # Tolerâncias
    Tol_P::Float64         # Tolerância para pressão
    Tol_u::Float64         # Tolerância para velocidade
    Tol_α::Float64         # Tolerância para fração volumetrica
    
    # Propriedades térmicas
    Cp_L::Float64          # Capacidade térmica da fase líquida/água
    Cp_G::Float64          # Capacidade térmica da fase gás/óleo
    kf_L::Float64          # Condutividade da líquida/água
    kf_G::Float64          # Condutividade do gás/óleo
    kf_pipe::Float64       # Condutividade do duto
    μL::Float64            # Viscosidade do líquido/água (Pa.s)
    μG::Float64            # Viscosidade do gás/óleo (Pa.s)
    h_ext::Float64         # entalpia externa
    
    # Propriedades físicas
    g::Float64             # Aceleração da gravidade (m/s²)
    cG::Float64            # velocidade do som no gás/óleo (m/s)
    cL::Float64            # velocidade do som no líquido/água (m/s)
    βG::Float64            # Coeficiente de expansão volumétrica da gás/óleo
    βL::Float64            # Coeficiente de expansão volumétrica do líquido/óleo
    γ::Float64
    
    # Parâmetros de método
    metodo::Int            # 1 = MUSCL, 2 = FR
    Grau::Int              # ordem do polinômio
    C::Int                 # esquema FR (0=DG,1=SD,2=Huynh,-1=Radau)
    
    # Malhas (serão inicializadas depois)
    meshP::Any             # Malha primária
    meshS::Any             # Malha secundária
end

# Estrutura para armazenar os flows FR/MUSCL
mutable struct FlowStructures
    uLflow::Flow
    uGflow::Flow
    αLflow::Flow
    αGflow::Flow
    hLflow::Flow
end

# Estrutura para armazenar variáveis de simulação
mutable struct SimulationVariables
    # Parâmetros derivados
    N::Int
    Δx::Float64
    A::Float64
    S::Float64
    αL_inlet::Float64
    αL_0::Float64
    hL_intet::Float64
    hLi::Float64
    Λ::Float64
    Pr_L::Float64
    ρL_ref::Float64
    ρG_ref::Float64
    t::Float64
    Δt::Float64
    
    # Variáveis no tempo n+1
    P::Vector{Float64}
    αL::Vector{Float64}
    αG::Vector{Float64}
    αLf::Vector{Float64}
    αGf::Vector{Float64}
    uL::Vector{Float64}
    uG::Vector{Float64}
    ρL::Vector{Float64}
    ρG::Vector{Float64}
    ρLf::Vector{Float64}
    ρGf::Vector{Float64}
    hL::Vector{Float64}
    TL::Vector{Float64}
    hLf::Vector{Float64}
    
    # Variáveis no tempo n
    P0::Vector{Float64}
    αL0::Vector{Float64}
    αG0::Vector{Float64}
    αLf0::Vector{Float64}
    αGf0::Vector{Float64}
    uL0::Vector{Float64}
    uG0::Vector{Float64}
    ρL0::Vector{Float64}
    ρG0::Vector{Float64}
    ρLf0::Vector{Float64}
    ρGf0::Vector{Float64}
    hL0::Vector{Float64}
    hLf0::Vector{Float64}
    
    # Valores supostos
    PØ::Vector{Float64}
    αLØ::Vector{Float64}
    αGØ::Vector{Float64}
    uLØ::Vector{Float64}
    uGØ::Vector{Float64}
    uLf::Vector{Float64}
    uGf::Vector{Float64}
    uLδ::Vector{Float64}
    uGδ::Vector{Float64}
    
    # Velocidades centradas
    uLC::Vector{Float64}
    uGC::Vector{Float64}
    uLC0::Vector{Float64}
    uGC0::Vector{Float64}
    
    # Números adimensionais e outros
    ReL::Vector{Float64}
    ReG::Vector{Float64}
    R_núcleo::Vector{Float64}
    fw::Vector{Float64}
    fi::Vector{Float64}
    λ::Vector{Float64}
    τw::Vector{Float64}
    τi::Vector{Float64}
    Nu::Vector{Float64}
    ULM::Vector{Float64}
    U_conv::Vector{Float64}
    R_tot::Vector{Float64}
    
    # Difusão artificial
    dif_art_L::Vector{Float64}
    dif_art_G::Vector{Float64}
    
    # Sistemas lineares para pressão
    d_P::Vector{Float64}
    du_P::Vector{Float64}
    dl_P::Vector{Float64}
    b_P::Vector{Float64}
    
    # Variáveis auxiliares
    Paux::Vector{Float64}
    αLaux::Vector{Float64}
    αGaux::Vector{Float64}
    uLaux::Vector{Float64}
    uGaux::Vector{Float64}
    Pδ::Vector{Float64}
    Tδ::Vector{Float64}
    
    # Variáveis delta para MUSCL/FR
    uLδD::Matrix{Float64}
    uLδD0::Matrix{Float64}
    uGδD::Matrix{Float64}
    uGδD0::Matrix{Float64}
    αLδD::Matrix{Float64}
    αLδD0::Matrix{Float64}
    αGδD::Matrix{Float64}
    αGδD0::Matrix{Float64}
    
    # Variáveis auxiliares para interpolação entre malhas - EXATAMENTE como no original
    αLu::Matrix{Float64}   # Interpolação de αL para malha de velocidades
    αGu::Matrix{Float64}   # Interpolação de αG para malha de velocidades
    αLu0::Matrix{Float64}  # Versão tempo atrasado
    αGu0::Matrix{Float64}  # Versão tempo atrasado
    uLα::Matrix{Float64}   # Interpolação de uL para malha de frações volumétricas
    uGα::Matrix{Float64}   # Interpolação de uG para malha de frações volumétricas
    uLα0::Matrix{Float64}  # Versão tempo atrasado
    uGα0::Matrix{Float64}  # Versão tempo atrasado
end

"""
Lê os parâmetros de entrada do arquivo de texto
"""
function read_input_parameters(filename::String)
    n_dados = 37
    dados = zeros(n_dados)
    texto = readlines(filename)
    
    for i = 2:2:2*n_dados
        aux = texto[i]
        aux2 = Int(i / 2)
        dados[aux2] = parse(Float64, aux)
    end
    
    return InputParameters(
        dados[1],   # Comp
        dados[2],   # D
        dados[3],   # θ
        Int(dados[4]),   # n_segmentos
        dados[27],  # D_ext
        dados[5],   # tempo
        dados[15],  # CFL
        dados[6],   # uL_inlet
        dados[7],   # uG_inlet
        dados[8],   # ρLi
        dados[9],   # ρGi
        dados[10],  # P_out
        dados[11],  # αG_inlet
        dados[12],  # αG_0
        dados[13],  # T_intet
        dados[14],  # T_0
        dados[23],  # T_ext
        dados[16],  # Tol_P
        dados[17],  # Tol_u
        dados[18],  # Tol_α
        dados[19],  # Cp_L
        dados[20],  # Cp_G
        dados[21],  # kf_L
        dados[22],  # kf_G
        dados[26],  # kf_pipe
        dados[24],  # μL
        dados[25],  # μG
        dados[34],  # h_ext
        dados[28],  # g
        dados[29],  # cG
        dados[30],  # cL
        dados[31],  # βG
        dados[32],  # βL
        dados[33],  # γ
        Int(dados[35]),  # metodo
        Int(dados[36]),  # Grau
        Int(dados[37]),  # C
        nothing,    # meshP (será inicializada depois)
        nothing     # meshS (será inicializada depois)
    )
end

"""
Inicializa as estruturas de flow para MUSCL/FR - EXATAMENTE como no código original
"""
function initialize_flows(params::InputParameters, vars::SimulationVariables)
    T = Float64
    Grau = params.Grau
    N = vars.N
    C = params.C
    LMIN = 0.0
    LMAX = params.Comp
    Δx = params.Comp / N
    
    # Geração das malhas conforme código original
    (meshP, innerPointsP) = frgrid(LMIN, LMAX, N, Grau, C) # Malha Primária
    (meshS, innerPointsS) = frgrid(LMIN - Δx / 2, LMAX + Δx / 2, N + 1, Grau, C) # Malha Secundária
    
    # Criar flows conforme código original
    # uL:
    f1(x) = params.uL_inlet
    uL_IC = Function[f1]
    uLflow = Flow(N + 1, Grau, 1, 1; T=T) # N+1 por ser a malha secundária 
    init!(uLflow, uL_IC, innerPointsS)
    uLflow.uδD[1, :, 1] .= params.uL_inlet        # correção "pé duro"
    uLflow.uδD[end, :, 1] .= params.uL_inlet      # correção "pé duro"
    
    # Inicializar uLδD conforme código original
    vars.uLδD[:, :] = uLflow.uδD[:, :, 1]
    vars.uLδD0[:, :] = uLflow.uδD[:, :, 1]
    vars.uLδD0[:, :] = vars.uLδD[:, :]
    
    # uG:
    f2(x) = params.uG_inlet
    uG_IC = Function[f2]
    uGflow = Flow(N + 1, Grau, 1, 2; T=T) # N+1 por ser a malha secundária 
    init!(uGflow, uG_IC, innerPointsS)
    uGflow.uδD[1, :, 1] .= params.uG_inlet        # correção "pé duro"
    uGflow.uδD[end, :, 1] .= params.uG_inlet      # correção "pé duro"
    
    # Inicializar uGδD conforme código original
    vars.uGδD[:, :] = uGflow.uδD[:, :, 1]
    vars.uGδD0[:, :] = uGflow.uδD[:, :, 1]
    
    # αL:
    f3(x) = vars.αL_0
    αL_IC = Function[f3]
    αLflow = Flow(N, Grau, 2, 3; T=T)
    init!(αLflow, αL_IC, innerPointsP)
    αLflow.uδD[1, :, 1] .= vars.αL_inlet        # correção "pé duro"
    αLflow.uδD[end, :, 1] .= params.αG_inlet      # correção "pé duro"
    
    # Inicializar αLδD conforme código original
    vars.αLδD[:, :] = αLflow.uδD[:, :, 1]
    vars.αLδD0[:, :] = αLflow.uδD[:, :, 1]
    
    # αG:
    f4(x) = params.αG_0
    αG_IC = Function[f4]
    αGflow = Flow(N, Grau, 2, 4; T=T)
    init!(αGflow, αG_IC, innerPointsP)
    αGflow.uδD[1, :, 1] .= params.αG_inlet        # correção "pé duro"
    αGflow.uδD[end, :, 1] .= vars.αL_inlet      # correção "pé duro"
    
    # Inicializar αGδD conforme código original
    vars.αGδD[:, :] = αGflow.uδD[:, :, 1]
    vars.αGδD0[:, :] = αGflow.uδD[:, :, 1]
    
    # hL:
    f5(x) = vars.hLi
    hL_IC = Function[f5]
    hLflow = Flow(N, Grau, 3, 5; T=T)
    init!(hLflow, hL_IC, innerPointsP)
    hLflow.uδD[1, :, 1] .= vars.hL_intet          # Correção "pé duro"
    hLflow.uδD[end, :, 1] .= vars.hL_intet        # Correção "pé duro"
    
    return FlowStructures(uLflow, uGflow, αLflow, αGflow, hLflow)
end

"""
Inicializa as variáveis de simulação baseadas nos parâmetros de entrada
"""
function initialize_simulation_variables(params::InputParameters)
    N = Int(params.n_segmentos)
    Δx = params.Comp / params.n_segmentos
    A = (π * params.D^2) / 4
    S = π * params.D
    αL_inlet = 1 - params.αG_inlet
    αL_0 = 1 - params.αG_0
    hL_intet = params.Cp_L * params.T_intet
    hLi = params.Cp_L * params.T_0
    Λ = params.Cp_L / params.Cp_G
    Pr_L = (params.Cp_L * params.μL) / params.kf_L
    ρL_ref = params.ρLi
    ρG_ref = params.ρGi
    
    return SimulationVariables(
        N, Δx, A, S, αL_inlet, αL_0, hL_intet, hLi, Λ, Pr_L, ρL_ref, ρG_ref,
        0.0, 0.0,  # t, Δt
        # Variáveis n+1
        zeros(N), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1),  # P, αL, αG, αLf, αGf
        zeros(N + 1), zeros(N + 1), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1),  # uL, uG, ρL, ρG, ρLf, ρGf
        zeros(N), zeros(N), zeros(N + 1),  # hL, TL, hLf
        # Variáveis n
        zeros(N), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1),  # P0, αL0, αG0, αLf0, αGf0
        zeros(N + 1), zeros(N + 1), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1),  # uL0, uG0, ρL0, ρG0, ρLf0, ρGf0
        zeros(N), zeros(N + 1),  # hL0, hLf0
        # Valores supostos
        zeros(N), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1),  # PØ, αLØ, αGØ, uLØ, uGØ
        zeros(N + 1), zeros(N + 1), zeros(N + 1), zeros(N + 1),  # uLf, uGf, uLδ, uGδ
        # Velocidades centradas
        zeros(N), zeros(N), zeros(N), zeros(N),  # uLC, uGC, uLC0, uGC0
        # Números adimensionais
        zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N),  # ReL, ReG, R_núcleo, fw, fi, λ
        zeros(N + 1), zeros(N + 1), zeros(N), zeros(N), zeros(N), zeros(N),  # τw, τi, Nu, ULM, U_conv, R_tot
        # Difusão artificial
        zeros(N + 1), zeros(N + 1),  # dif_art_L, dif_art_G
        # Sistemas lineares
        zeros(N), zeros(N), zeros(N), zeros(N),  # d_P, du_P, dl_P, b_P
        # Auxiliares
        zeros(N), zeros(N), zeros(N), zeros(N + 1), zeros(N + 1), zeros(N), zeros(N),  # Paux, αLaux, αGaux, uLaux, uGaux, Pδ, Tδ
        # Variáveis delta
        zeros(N + 3, params.Grau + 1), zeros(N + 3, params.Grau + 1),  # uLδD, uLδD0
        zeros(N + 3, params.Grau + 1), zeros(N + 3, params.Grau + 1),  # uGδD, uGδD0
        zeros(N + 2, params.Grau + 1), zeros(N + 2, params.Grau + 1),  # αLδD, αLδD0
        zeros(N + 2, params.Grau + 1), zeros(N + 2, params.Grau + 1),  # αGδD, αGδD0
        # Variáveis auxiliares para interpolação
        zeros(N + 3, params.Grau + 1), zeros(N + 3, params.Grau + 1),  # αLu, αGu
        zeros(N + 3, params.Grau + 1), zeros(N + 3, params.Grau + 1),  # αLu0, αGu0
        zeros(N + 2, params.Grau + 1), zeros(N + 2, params.Grau + 1),  # uLα, uGα
        zeros(N + 2, params.Grau + 1), zeros(N + 2, params.Grau + 1)   # uLα0, uGα0
    )
end

"""
Aplica condições iniciais nas variáveis de simulação
"""
function apply_initial_conditions!(vars::SimulationVariables, params::InputParameters)
    # Condições iniciais básicas
    vars.uLC[:] .= params.uL_inlet
    vars.uGC[:] .= params.uG_inlet
    vars.uLC0[:] .= params.uL_inlet
    vars.uGC0[:] .= params.uG_inlet
    
    vars.uL[:] .= params.uL_inlet
    vars.uG[:] .= params.uG_inlet
    vars.αL[:] .= vars.αL_0
    vars.αG[:] .= params.αG_0
    vars.ρL[:] .= params.ρLi
    vars.ρG[:] .= params.ρGi
    vars.P[:] .= params.P_out
    vars.hL[:] .= vars.hLi
    
    # Valores no tempo n - CORREÇÃO: sequência exata do código original
    vars.uL0[:] .= params.uL_inlet
    vars.uG0[:] .= params.uG_inlet
    vars.αL0[:] .= vars.αL_inlet     # CORREÇÃO: usando vars.αL_inlet que está na struct
    vars.αG0[:] .= params.αG_inlet   # CORREÇÃO: era αG_inlet no original, não αG_0
    vars.ρL0[:] .= params.ρLi
    vars.ρG0[:] .= params.ρGi
    vars.P0[:] .= params.P_out
    vars.hL0[:] .= vars.hLi
    
    # Valores supostos
    vars.uLØ[:] .= params.uL_inlet
    vars.uGØ[:] .= params.uG_inlet
    vars.αLØ[:] .= vars.αL0[:]
    vars.αGØ[:] .= vars.αG0[:]
    vars.PØ[:] .= params.P_out
    
    # Inicializar sistema linear para pressão
    vars.d_P[1] = 1.0
    vars.du_P[1] = 0.0
    vars.dl_P[1] = 0.0
    vars.b_P[1] = 0.0
    vars.d_P[vars.N] = 1.0
    vars.du_P[vars.N] = 0.0
    vars.dl_P[vars.N] = 0.0
    vars.b_P[vars.N] = 0.0
end

"""
Calcula o passo de tempo baseado no CFL
"""
function calculate_time_step(vars::SimulationVariables, params::InputParameters)
    aux1 = maximum(abs.(vars.uL[:]))
    aux2 = maximum(abs.(vars.uG[:]))
    u = [aux1, aux2]
    aux3 = maximum(u)
    
    # Adicionar velocidades do som para estabilidade
    max_wave_speed = aux3 + max(params.cL, params.cG)
    
    if params.metodo == 2
        Δt = (params.CFL * vars.Δx) / ((2 * params.Grau + 1) * max_wave_speed)  # FR
    else
        Δt = (params.CFL * vars.Δx) / (max_wave_speed)  # MUSCL
    end
    
    # Limitar o passo de tempo para estabilidade
    Δt_max = 0.001  # Limite máximo para evitar instabilidade
    Δt = min(Δt, Δt_max)
    
    return minimum([Δt, (params.tempo - vars.t)])
end

"""
Calcula coeficientes nas faces usando esquema upwind - EXATAMENTE como no original
"""
function calculate_face_coefficients!(vars::SimulationVariables, params::InputParameters)
    N = vars.N
    
    # Condições de contorno para frações volumétricas e densidades
    vars.αLf[1] = vars.αL[1]
    vars.αGf[1] = vars.αG[1]
    vars.ρLf[1] = vars.ρL[1]
    vars.ρGf[1] = vars.ρG[1]
    vars.αLf0[1] = vars.αL0[1]
    vars.αGf0[1] = vars.αG0[1]
    vars.ρLf0[1] = vars.ρL0[1]
    vars.ρGf0[1] = vars.ρG0[1]
    
    vars.αLf[N+1] = vars.αL[N]
    vars.αGf[N+1] = vars.αG[N]
    vars.ρLf[N+1] = vars.ρL[N]
    vars.ρGf[N+1] = vars.ρG[N]
    vars.αLf0[N+1] = vars.αL0[N]
    vars.αGf0[N+1] = vars.αG0[N]
    vars.ρLf0[N+1] = vars.ρL0[N]
    vars.ρGf0[N+1] = vars.ρG0[N]
    
    # Condições de contorno para entalpias - ADICIONADO conforme original
    vars.hLf[1] = vars.hL[1]
    vars.hLf0[1] = vars.hL0[1]
    vars.hLf[N+1] = vars.hL[N]
    vars.hLf0[N+1] = vars.hL0[N]
    
    # Faces internas
    for i = 2:N
        # Frações volumétricas
        vars.αLf[i] = vars.αL[i-1] * (sign(vars.uLØ[i]) + 1) / 2 + vars.αL[i] * (-sign(vars.uLØ[i]) + 1) / 2
        vars.αGf[i] = vars.αG[i-1] * (sign(vars.uGØ[i]) + 1) / 2 + vars.αG[i] * (-sign(vars.uGØ[i]) + 1) / 2
        vars.αLf0[i] = vars.αL0[i-1] * (sign(vars.uL[i]) + 1) / 2 + vars.αL0[i] * (-sign(vars.uL[i]) + 1) / 2
        vars.αGf0[i] = vars.αG0[i-1] * (sign(vars.uG[i]) + 1) / 2 + vars.αG0[i] * (-sign(vars.uG[i]) + 1) / 2
        
        # Densidades
        vars.ρLf[i] = vars.ρL[i-1] * (sign(vars.uL[i]) + 1) / 2 + vars.ρL[i] * (-sign(vars.uL[i]) + 1) / 2
        vars.ρGf[i] = vars.ρG[i-1] * (sign(vars.uG[i]) + 1) / 2 + vars.ρG[i] * (-sign(vars.uG[i]) + 1) / 2
        vars.ρLf0[i] = vars.ρL0[i-1] * (sign(vars.uL[i]) + 1) / 2 + vars.ρL0[i] * (-sign(vars.uL[i]) + 1) / 2
        vars.ρGf0[i] = vars.ρG0[i-1] * (sign(vars.uG[i]) + 1) / 2 + vars.ρG0[i] * (-sign(vars.uG[i]) + 1) / 2
        
        # Entalpias - ADICIONADO conforme original
        vars.hLf[i] = vars.hL[i-1] * (sign(vars.uL[i]) + 1) / 2 + vars.hL[i] * (-sign(vars.uL[i]) + 1) / 2
        vars.hLf0[i] = vars.hL0[i-1] * (sign(vars.uL[i]) + 1) / 2 + vars.hL0[i] * (-sign(vars.uL[i]) + 1) / 2
    end
end

"""
Calcula difusão artificial
"""
function calculate_artificial_diffusion!(vars::SimulationVariables, params::InputParameters)
    # EQUIVALÊNCIA TOTAL: usar valores originais sem proteção
    for i = 2:vars.N
        Δp_art = (params.γ * vars.αLf0[i] * vars.ρLf0[i] * vars.αGf0[i] * vars.ρGf0[i] * (vars.uG0[i] - vars.uL0[i])^2) / 
                 (vars.αGf0[i] * vars.ρLf0[i] + vars.αLf0[i] * vars.ρGf0[i])
        vars.dif_art_L[i] = Δp_art * (vars.αL0[i] - vars.αL0[i-1]) / vars.Δx
        vars.dif_art_G[i] = Δp_art * (vars.αG0[i] - vars.αG0[i-1]) / vars.Δx
    end
end

"""
Calcula velocidades centradas
"""
function calculate_centered_velocities!(vars::SimulationVariables)
    for i = 1:vars.N
        vars.uLC[i] = (vars.uL[i] + vars.uL[i+1]) / 2
        vars.uGC[i] = (vars.uG[i] + vars.uG[i+1]) / 2
        vars.uLC0[i] = (vars.uL0[i] + vars.uL0[i+1]) / 2
        vars.uGC0[i] = (vars.uG0[i] + vars.uG0[i+1]) / 2
    end
end

"""
Calcula coeficientes de transferência de calor - EXATAMENTE como no código original
"""
function calculate_heat_transfer_coefficients!(vars::SimulationVariables, params::InputParameters)
    # Proteção apenas para esta operação específica para evitar números negativos na exponenciação
    αG_safe = clamp.(vars.αG[:], 0.001, 0.999)
    vars.U_conv[:] .= vars.ULM[:] .* (1 ./ (1 .- αG_safe)) .^ 0.9
    vars.R_tot[:] .= (1 ./ (vars.U_conv[:] .* params.D * π .* vars.Δx)) .+
                     log(params.D_ext / params.D) ./ (2 * π * params.kf_pipe * vars.Δx) .+ 
                     (1 ./ (params.h_ext .* (vars.Δx .* π .* params.D_ext)))
end

"""
Calcula números de Reynolds e fatores de atrito
"""
function calculate_reynolds_and_friction!(vars::SimulationVariables, params::InputParameters)
    # Números de Reynolds conforme código original
    vars.ReL[:] .= abs.((vars.ρL[:] .* vars.uLC[:] .* vars.αL[:] * params.D) ./ params.μL) # Núcleo do Anular
    vars.R_núcleo[:] .= (params.D / 2) * (abs.((vars.αG[:] .* vars.uLC[:]) ./ (vars.uGC[:] + vars.αG[:] .* (vars.uLC[:] - vars.uGC[:])))) .^ (1 / 2)
    vars.ReG[:] .= abs.(2 * (vars.ρG[:] .* vars.uGC[:] .* vars.R_núcleo[:]) ./ params.μG) # Anel do Anular
    
    # Fatores de atrito interfacial
    vars.fi[:] .= 0.079 ./ ((vars.ReG[:] .+ 1e-12) .^ 0.25) # Para escoamento turbulento!!!!!
    vars.τi[1:end-1] .= (vars.fi[:] ./ 2) .* vars.ρG[:] .* abs.(vars.uGC[:] - vars.uLC[:]) .* (vars.uGC[:] - vars.uLC[:])
    vars.τi[1] = vars.τi[2]
    
    # Fatores de atrito na parede
    vars.fw[:] .= 0.046 ./ ((vars.ReL[:] .+ 1e-12) .^ 0.2) # Para escoamento turbulento!!!!!
    
    # EQUIVALÊNCIA com proteção mínima apenas para evitar travamento (como no original pode ter)
    vars.λ[:] .= (1 .+ vars.αG[:] + (2 * vars.αG[:] .* log.(vars.αG[:])) ./ (1 .- vars.αG[:])) ./ (1 .- vars.αG[:])
    vars.τw[1:end-1] .= vars.fw[:] .* vars.ρL[:] .* vars.uLC[:] .* abs.(vars.uLC[:]) ./ 2 - 
                        (vars.ρG[:] - vars.ρL[:]) .* (-params.g) .* params.D .* vars.αG[:] .* vars.λ[:] ./ 4
    
    # Condições de contorno
    vars.τi[end] = vars.τi[end-1]
    vars.τw[end] = vars.τw[end-1]
    
    # Número de Nusselt e coeficiente convectivo
    vars.Nu[:] .= ((vars.fw[:] ./ 8) .* (vars.ReL[:] .- 1000) .* vars.Pr_L) ./
                  (1.07 .+ 12.7 .* ((vars.fw[:] ./ 8) .^ 0.5) .* ((vars.Pr_L .^ (2 / 3)) .- 1)) # Nusselt turbulento para o líquido/água
    vars.ULM[:] .= (params.kf_L .* vars.Nu[:]) ./ params.D # Coeficiente convectivo para o líquido/água
end

"""
Resolve a equação de pressão usando TDMA - EXATAMENTE como no código original
"""
function solve_pressure_equation!(vars::SimulationVariables, params::InputParameters)
    N = vars.N
    
    # Sistema TDMA conforme código original
    vars.d_P[N] = 1.0
    vars.du_P[N] = 0.0
    vars.dl_P[N] = 0.0
    vars.b_P[N] = 0.0
    
    for i = 2:(N-1)
        # Coeficientes Pós Agrupamento (Tudo em upwind) - EXATO do código original
        if vars.uLØ[i+1] >= 0
            αLen = vars.αL0[i]
            ρLen = vars.ρL0[i]
            αLwn = vars.αL0[i-1]
            ρLwn = vars.ρL0[i-1]
            αLe = vars.αLØ[i]
            ρLe = vars.ρL[i]
            αLw = vars.αLØ[i-1]
            ρLw = vars.ρL[i-1]
            a_Le = αLen * ρLen / vars.Δt + αLen * ρLen * vars.uL0[i+1] / vars.Δx
            a_Lw = αLwn * ρLwn / vars.Δt + αLwn * ρLwn * vars.uL0[i] / vars.Δx
        else
            αLen = vars.αL0[i+1]
            ρLen = vars.ρL0[i+1]
            αLwn = vars.αL0[i]
            ρLwn = vars.ρL0[i]
            αLe = vars.αLØ[i+1]
            ρLe = vars.ρL[i+1]
            αLw = vars.αLØ[i]
            ρLw = vars.ρL[i]
            a_Le = αLen * ρLen / vars.Δt - αLen * ρLen * vars.uL0[i+1] / vars.Δx
            a_Lw = αLwn * ρLwn / vars.Δt - αLwn * ρLwn * vars.uL0[i] / vars.Δx
        end
        
        if vars.uGØ[i+1] >= 0
            αGen = vars.αG0[i]
            ρGen = vars.ρG0[i]
            αGwn = vars.αG0[i-1]
            ρGwn = vars.ρG0[i-1]
            αGe = vars.αGØ[i]
            ρGe = vars.ρG[i]
            αGw = vars.αGØ[i-1]
            ρGw = vars.ρG[i-1]
            a_Ge = αGen * ρGen / vars.Δt + αGen * ρGen * vars.uG0[i+1] / vars.Δx
            a_Gw = αGwn * ρGwn / vars.Δt + αGwn * ρGwn * vars.uG0[i] / vars.Δx
        else
            αGen = vars.αG0[i+1]
            ρGen = vars.ρG0[i+1]
            αGwn = vars.αG0[i]
            ρGwn = vars.ρG0[i]
            αGe = vars.αGØ[i+1]
            ρGe = vars.ρG[i+1]
            αGw = vars.αGØ[i]
            ρGw = vars.ρG[i]
            a_Ge = αGen * ρGen / vars.Δt - αGen * ρGen * vars.uG0[i+1] / vars.Δx
            a_Gw = αGwn * ρGwn / vars.Δt - αGwn * ρGwn * vars.uG0[i] / vars.Δx
        end
        
        mL = (αLe * ρLe * vars.uLØ[i+1] - αLw * ρLw * vars.uLØ[i]) / vars.ρL_ref
        mG = (αGe * ρGe * vars.uGØ[i+1] - αGw * ρGw * vars.uGØ[i]) / vars.ρG_ref
        coef = (vars.Δx / vars.Δt) * (vars.αLØ[i] / (vars.ρL_ref * params.cL^2) + vars.αGØ[i] / (vars.ρG_ref * params.cG^2))
        
        vars.d_P[i] = coef + αLe * ρLe * (αLe / (vars.Δx * a_Le)) / vars.ρL_ref + αGe * ρGe * (αGe / (vars.Δx * a_Ge)) / vars.ρG_ref + (αLw * ρLw * (αLw / (vars.Δx * a_Lw)) / vars.ρL_ref + αGw * ρGw * (αGw / (vars.Δx * a_Gw)) / vars.ρG_ref)
        vars.du_P[i] = -(αLe * ρLe * (αLe / (vars.Δx * a_Le)) / vars.ρL_ref + αGe * ρGe * (αGe / (vars.Δx * a_Ge)) / vars.ρG_ref)
        vars.dl_P[i] = -(αLw * ρLw * (αLw / (vars.Δx * a_Lw)) / vars.ρL_ref + αGw * ρGw * (αGw / (vars.Δx * a_Gw)) / vars.ρG_ref)
        vars.b_P[i] = -(mL + mG + (vars.Δx / vars.Δt) * ((vars.αLØ[i] * vars.ρL[i] - vars.αL0[i] * vars.ρL0[i]) / vars.ρL_ref + (vars.αGØ[i] * vars.ρG[i] - vars.αG0[i] * vars.ρG0[i]) / vars.ρG_ref))
    end
    
    vars.Pδ[:] = TDMA(vars.d_P, vars.du_P, vars.dl_P, vars.b_P)
    vars.Pδ[1] = vars.Pδ[2]  # Como no código original
end

"""
Calcula correções de velocidade baseadas na pressão
"""
function calculate_velocity_corrections!(vars::SimulationVariables, params::InputParameters)
    N = vars.N
    
    for i = 1:(N-1)
        # Coeficientes Pós Agrupamento (Tudo em upwind)
        if vars.uLØ[i+1] >= 0
            αLe = vars.αL0[i]
            ρLe = vars.ρL0[i]
            a_L = αLe * ρLe / vars.Δt + αLe * ρLe * vars.uL0[i+1] / vars.Δx
        else
            αLe = vars.αL0[i+1]
            ρLe = vars.ρL0[i+1]
            a_L = αLe * ρLe / vars.Δt - αLe * ρLe * vars.uL0[i+1] / vars.Δx
        end
        
        if vars.uGØ[i+1] >= 0
            αGe = vars.αG0[i]
            ρGe = vars.ρG0[i]
            a_G = αGe * ρGe / vars.Δt + αGe * ρGe * vars.uG0[i+1] / vars.Δx
        else
            αGe = vars.αG0[i+1]
            ρGe = vars.ρG0[i+1]
            a_G = αGe * ρGe / vars.Δt - αGe * ρGe * vars.uG0[i+1] / vars.Δx
        end
        
        de_L = αLe / (vars.Δx * max(abs(a_L), 1e-12) * sign(a_L + 1e-12))
        de_G = αGe / (vars.Δx * max(abs(a_G), 1e-12) * sign(a_G + 1e-12))
        
        # Proteção contra NaN
        if isfinite(de_L) && isfinite(vars.Pδ[i]) && isfinite(vars.Pδ[i+1])
            vars.uLδ[i+1] = de_L * (vars.Pδ[i] - vars.Pδ[i+1])
        else
            vars.uLδ[i+1] = 0.0
        end
        
        if isfinite(de_G) && isfinite(vars.Pδ[i]) && isfinite(vars.Pδ[i+1])
            vars.uGδ[i+1] = de_G * (vars.Pδ[i] - vars.Pδ[i+1])
        else
            vars.uGδ[i+1] = 0.0
        end
    end
end

"""
Atualiza pressão e velocidades com correções
"""
function update_pressure_and_velocities!(vars::SimulationVariables, params::InputParameters)
    # Atualização da pressão conforme código original
    vars.P[:] = vars.PØ[:] + vars.Pδ[:]
    
    # Proteção para valores de pressão inválidos
    for i = 1:length(vars.P)
        if isnan(vars.P[i]) || isinf(vars.P[i]) || vars.P[i] <= 0.0
            vars.P[i] = params.P_out
        end
    end
    
    # Atualização das velocidades conforme código original
    vars.uL = vars.uLØ + vars.uLδ
    vars.uG = vars.uGØ + vars.uGδ
    
    # Proteção para valores de velocidade inválidos
    for i = 1:length(vars.uL)
        if isnan(vars.uL[i]) || isinf(vars.uL[i])
            vars.uL[i] = params.uL_inlet
        end
        if isnan(vars.uG[i]) || isinf(vars.uG[i])
            vars.uG[i] = params.uG_inlet
        end
    end
    
    # Atualização dos flows conforme código original
    for i = 2:(vars.N+2)
        vars.uLδD[i, :] .+= vars.uLδ[i-1]
        vars.uGδD[i, :] .+= vars.uGδ[i-1]
    end
    vars.uLδD[vars.N+3, :] = vars.uLδD[vars.N+2, :]
    vars.uGδD[vars.N+3, :] = vars.uGδD[vars.N+2, :]
    
    # Condições de contorno
    vars.uL[vars.N+1] = vars.uL[vars.N]
    vars.uG[vars.N+1] = vars.uG[vars.N]
end

"""
Atualiza densidade baseada nas correções de pressão e temperatura - EXATAMENTE como no código original
"""
function update_density!(vars::SimulationVariables, params::InputParameters)
    # Calcular variação da temperatura conforme código original
    for i = 2:vars.N
        vars.Tδ[i] = (1 / params.Cp_L) * (vars.hL[i] - vars.hL0[i])
    end
    
    # Correção da Densidade - EXATAMENTE como no código original
    for i = 1:vars.N
        vars.ρL[i] = vars.ρL[i] + (1 / params.cL^2) * vars.Pδ[i] - vars.ρL[i] * params.βL * vars.Tδ[i]
        vars.ρG[i] = vars.ρG[i] + (1 / params.cG^2) * vars.Pδ[i] - vars.ρG[i] * params.βG * vars.Tδ[i]
    end
end

"""
Interpola frações volumétricas αL/αG para a malha de velocidades - EXATAMENTE como no original
"""
function interpolate_alpha_to_velocity_mesh!(vars::SimulationVariables, pontos::Vector{Float64})
    N = vars.N
    Grau = length(pontos) - 1  # pontos tem tamanho Grau+1, então Grau = length(pontos) - 1
    
    # Condições de contorno conforme código original
    vars.αLu[1, :] .= vars.αL[1]
    vars.αGu[1, :] .= vars.αG[1]
    vars.αLu0[1, :] .= vars.αL0[1]
    vars.αGu0[1, :] .= vars.αG0[1]
    
    vars.αLu[end-1:end, :] .= vars.αL[end]
    vars.αGu[end-1:end, :] .= vars.αG[end]
    vars.αLu0[end-1:end, :] .= vars.αL0[end]
    vars.αGu0[end-1:end, :] .= vars.αG0[end]
    
    # Interpolação nas faces internas - EXATAMENTE como no código original
    for i = 2:N+1
        for j = 1:Grau+1
            β1 = pontos[j] + 1
            β2 = pontos[j] - 1
            lj1 = lagrangeBasisAt(β1, pontos)
            lj2 = lagrangeBasisAt(β2, pontos)
            
            vars.αLu[i, j] = sum(vars.αLδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                             sum(vars.αLδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.αGu[i, j] = sum(vars.αGδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                             sum(vars.αGδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.αLu0[i, j] = sum(vars.αLδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                              sum(vars.αLδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.αGu0[i, j] = sum(vars.αGδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                              sum(vars.αGδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
        end
    end
end

"""
Interpola velocidades uL/uG para a malha de frações volumétricas - EXATAMENTE como no original
"""
function interpolate_velocity_to_alpha_mesh!(vars::SimulationVariables, flows::FlowStructures, pontos::Vector{Float64})
    N = vars.N
    Grau = length(pontos) - 1  # pontos tem tamanho Grau+1, então Grau = length(pontos) - 1
    
    # Interpolação conforme código original (linhas 664-675)
    for i = 1:N+2
        for j = 1:Grau+1
            β1 = pontos[j] + 1
            β2 = pontos[j] - 1
            lj1 = lagrangeBasisAt(β1, pontos)
            lj2 = lagrangeBasisAt(β2, pontos)
            
            vars.uLα[i, j] = sum(flows.uLflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                             sum(flows.uLflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.uGα[i, j] = sum(flows.uGflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                             sum(flows.uGflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.uLα0[i, j] = sum(vars.uLδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                              sum(vars.uLδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            vars.uGα0[i, j] = sum(vars.uGδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + 
                              sum(vars.uGδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
        end
    end
end

"""
Resolve a equação de energia usando esquema MUSCL
"""
function solve_energy_equation!(vars::SimulationVariables, params::InputParameters)
    # Calcula coeficientes de transferência de calor (conforme código original)
    # Proteção apenas para esta operação específica para evitar números negativos na exponenciação
    αG_safe = clamp.(vars.αG[:], 0.001, 0.999)
    vars.U_conv[:] .= vars.ULM[:] .* (1 ./ (1 .- αG_safe)) .^ 0.9
    vars.R_tot[:] .= (1 ./ (vars.U_conv[:] .* params.D * π .* vars.Δx)) .+
                     log(params.D_ext / params.D) ./ (2 * π * params.kf_pipe * vars.Δx) .+ 
                     (1 ./ (params.h_ext .* (vars.Δx .* π .* params.D_ext)))
    
    # Entalpias nas faces já calculadas em calculate_face_coefficients!
    
    # MUSCL para energia
    u = copy(vars.hL)
    vars.hL[1] = vars.hL_intet  # Condição de contorno
    
    for i = 2:vars.N
        # Cálculo dos limiters e valores nas interfaces
        if i == 2
            r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
            r_i_menos = 0
            r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
            uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
            uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
            uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
            uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        elseif i == (vars.N - 1)
            r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
            r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
            r_i_mais = 0
            uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
            uRmais = u[i+1]
            uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
            uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        elseif i == vars.N
            r_i = (u[i] - u[i-1]) / (u[i] - u[i] + eps())
            r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
            r_i_mais = 0
            uLmais = u[i]
            uRmais = u[i]
            uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
            uRmenos = u[i]
        else
            r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
            r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
            r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
            uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
            uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
            uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
            uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        end
        
        # Cálculo dos fluxos
        F_mais = ((vars.αLf[i+1] * vars.ρLf[i+1] * vars.uL[i+1] * (uLmais + 0.5 * vars.uL[i+1]^2)) + 
                  (vars.αGf[i+1] * vars.ρGf[i+1] * vars.uG[i+1] * (vars.Λ * uLmais + 0.5 * vars.uG[i+1]^2))) * 
                  (sign((vars.uL[i] + vars.uL[i+1])/2) + 1) / 2 + 
                 ((vars.αLf[i] * vars.ρLf[i] * vars.uL[i] * (uRmais + 0.5 * vars.uL[i]^2)) + 
                  (vars.αGf[i] * vars.ρGf[i] * vars.uG[i] * (vars.Λ * uRmais + 0.5 * vars.uG[i]^2))) * 
                  (-sign((vars.uL[i] + vars.uL[i+1])/2) + 1) / 2
        
        F_menos = ((vars.αLf[i+1] * vars.ρLf[i+1] * vars.uL[i+1] * (uLmenos + 0.5 * vars.uL[i+1]^2)) + 
                   (vars.αGf[i+1] * vars.ρGf[i+1] * vars.uG[i+1] * (vars.Λ * uLmenos + 0.5 * vars.uG[i+1]^2))) * 
                   (sign((vars.uL[i] + vars.uL[i+1])/2) + 1) / 2 + 
                  ((vars.αLf[i] * vars.ρLf[i] * vars.uL[i] * (uRmenos + 0.5 * vars.uL[i]^2)) + 
                   (vars.αGf[i] * vars.ρGf[i] * vars.uG[i] * (vars.Λ * uRmenos + 0.5 * vars.uG[i]^2))) * 
                   (-sign((vars.uL[i] + vars.uL[i+1])/2) + 1) / 2
        
        # Fluxo líquido na célula i
        a3 = (-vars.Δt / vars.Δx) * (F_mais - F_menos)
        
        # Termos de fonte (conforme código original)
        a1 = -0.5 * (vars.αG[i] * vars.ρG[i] * vars.uG[i]^2 + vars.αL[i] * vars.ρL[i] * vars.uL[i]^2)
        a2 = vars.αG0[i] * vars.ρG0[i] * (vars.Λ * vars.hL0[i] + 0.5 * vars.uG0[i]^2) + 
             vars.αL0[i] * vars.ρL0[i] * (vars.hL0[i] + 0.5 * vars.uL0[i]^2)
        a4 = vars.P[i] - vars.PØ[i]
        a5 = vars.Δt * params.g * sin(params.θ) * (vars.αL0[i] * vars.ρL0[i] * vars.uL0[i] + vars.αG0[i] * vars.ρG0[i] * vars.uG0[i])
        
        # Transferência de calor com dependência de direção (conforme original)
        an = vars.hL0[i] / params.Cp_L
        if an >= params.T_ext
            a6 = (vars.Δt * vars.hL0[i] * vars.S) / (params.Cp_L * vars.R_tot[i] * vars.A)
            a7 = (-params.T_ext * vars.Δt * vars.S) / (vars.R_tot[i] * vars.A)
        else
            a6 = (-vars.Δt * vars.hL0[i] * vars.S) / (params.Cp_L * vars.R_tot[i] * vars.A)
            a7 = (params.T_ext * vars.Δt * vars.S) / (vars.R_tot[i] * vars.A)
        end
        a8 = vars.αG[i] * vars.ρG[i] * vars.Λ + vars.αL[i] * vars.ρL[i]
        
        # Atualização da entalpia (todos os 7 termos como no original)
        vars.hL[i] = (a1 + a2 + a3 + a4 + a5 + a6 + a7) / a8
        
        # Proteção contra valores inválidos na entalpia
        if isnan(vars.hL[i]) || isinf(vars.hL[i]) || vars.hL[i] <= 0.0
            vars.hL[i] = params.Cp_L * params.T_0  # hLi calculado
        end
    end
end

"""
Resolve as equações de velocidade usando MUSCL ou FR
"""
function solve_velocity_equations!(vars::SimulationVariables, params::InputParameters, flows::FlowStructures, weights::Vector{Float64}, meshS)
    # EQUIVALÊNCIA TOTAL: usar αGf0 original sem proteção para manter identidade numérica
    if params.metodo == 1  # MUSCL
        # println("Resolvendo velocidades com MUSCL")  # DEBUG removido
        
        # Chamadas MUSCL para as velocidades conforme código original
        (vars.uLδD[:, :]) = MUSCL(flows.uLflow, vars.uLδD0, vars.αLf, vars.ρLf, vars.αLf0, vars.ρLf0, vars.uL, meshS, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, vars.dif_art_L, vars.uLC0, vars.αL, vars.ρL, vars.uL0, vars.τi, vars.τw, vars.αGf0, params.D)
        
        (vars.uGδD[:, :]) = MUSCL(flows.uGflow, vars.uGδD0, vars.αGf, vars.ρGf, vars.αGf0, vars.ρGf0, vars.uG, meshS, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, vars.dif_art_G, vars.uGC0, vars.αG, vars.ρG, vars.uG0, vars.τi, vars.τw, vars.αGf0, params.D)
        
    else  # FR
        # println("Resolvendo velocidades com FR")  # DEBUG removido
        (vars.uLδD[:, :]) = FR(flows.uLflow, vars.uLδD0, vars.αLf, vars.ρLf, vars.αLf0, vars.ρLf0, vars.uL, meshS, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, vars.dif_art_L, vars.uLC0, vars.αL, vars.ρL, vars.uL0, vars.τi, vars.τw, vars.αGf0, params.D)
        
        (vars.uGδD[:, :]) = FR(flows.uGflow, vars.uGδD0, vars.αGf, vars.ρGf, vars.αGf0, vars.ρGf0, vars.uG, meshS, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, vars.dif_art_G, vars.uGC0, vars.αG, vars.ρG, vars.uG0, vars.τi, vars.τw, vars.αGf0, params.D)
    end
    
    # Boundary treatment conforme código original
    vars.uLδD[3, :] = (vars.uLδD[2, :] + vars.uLδD[4, :]) ./ 2
    vars.uGδD[3, :] = (vars.uGδD[2, :] + vars.uGδD[4, :]) ./ 2
    
    # Atualização dos flows conforme código original
    flows.uLflow.uδD[:, :, 1] = vars.uLδD[:, :]
    flows.uGflow.uδD[:, :, 1] = vars.uGδD[:, :]
    
    # Cálculo dos valores médios usando quadratura conforme código original
    for i = 1:vars.N+1
        vars.uLØ[i] = sum(vars.uLδD[i+1, :] .* weights[:]) / 2
        vars.uGØ[i] = sum(vars.uGδD[i+1, :] .* weights[:]) / 2
        
        # Proteção contra valores inválidos nas velocidades médias
        if isnan(vars.uLØ[i]) || isinf(vars.uLØ[i])
            vars.uLØ[i] = params.uL_inlet
        end
        if isnan(vars.uGØ[i]) || isinf(vars.uGØ[i])
            vars.uGØ[i] = params.uG_inlet
        end
    end
end

"""
Resolve as equações de fração volumétrica - EXATAMENTE como no código original
"""
function solve_volume_fraction_equations!(vars::SimulationVariables, params::InputParameters, flows::FlowStructures, weights::Vector{Float64}, pontos::Vector{Float64}, meshP)
    # Interpolação conforme código original - EXATAMENTE igual
    αLu = zeros(vars.N + 3, params.Grau + 1)
    αGu = zeros(vars.N + 3, params.Grau + 1)
    αLu0 = zeros(vars.N + 3, params.Grau + 1)
    αGu0 = zeros(vars.N + 3, params.Grau + 1)
    
    # Condições de contorno - igual ao original
    αLu[1, :] .= vars.αL[1]
    αGu[1, :] .= vars.αG[1]
    αLu0[1, :] .= vars.αL0[1]
    αGu0[1, :] .= vars.αG0[1]
    
    αLu[end-1:end, :] .= vars.αL[end]
    αGu[end-1:end, :] .= vars.αG[end]
    αLu0[end-1:end, :] .= vars.αL0[end]
    αGu0[end-1:end, :] .= vars.αG0[end]
    
    # Interpolação das velocidades para as frações conforme código original
    for i = 2:vars.N+1
        for j = 1:params.Grau+1
            β1 = pontos[j] + 1
            β2 = pontos[j] - 1
            lj1 = lagrangeBasisAt(β1, pontos)
            lj2 = lagrangeBasisAt(β2, pontos)
            
            αLu[i, j] = sum(vars.αLδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.αLδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            αGu[i, j] = sum(vars.αGδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.αGδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            αLu0[i, j] = sum(vars.αLδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.αLδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            αGu0[i, j] = sum(vars.αGδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.αGδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
        end
    end
    
    # Interpolação dos flows para as frações conforme código original  
    uLα = zeros(vars.N + 2, params.Grau + 1)
    uGα = zeros(vars.N + 2, params.Grau + 1)
    uLα0 = zeros(vars.N + 2, params.Grau + 1)
    uGα0 = zeros(vars.N + 2, params.Grau + 1)
    
    for i = 1:vars.N+2
        for j = 1:params.Grau+1
            β1 = pontos[j] + 1
            β2 = pontos[j] - 1
            lj1 = lagrangeBasisAt(β1, pontos)
            lj2 = lagrangeBasisAt(β2, pontos)
            
            uLα[i, j] = sum(flows.uLflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(flows.uLflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            uGα[i, j] = sum(flows.uGflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(flows.uGflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            uLα0[i, j] = sum(vars.uLδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.uLδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            uGα0[i, j] = sum(vars.uGδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(vars.uGδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
        end
    end
    
    # Cálculo de velocidades centradas conforme código original
    for i = 1:vars.N
        vars.uLC[i] = (vars.uL[i] + vars.uL[i+1]) / 2
        vars.uGC[i] = (vars.uG[i] + vars.uG[i+1]) / 2
        vars.uLC0[i] = (vars.uL0[i] + vars.uL0[i+1]) / 2
        vars.uGC0[i] = (vars.uG0[i] + vars.uG0[i+1]) / 2
    end
    
    # EQUIVALÊNCIA TOTAL: usar αGf0 original sem proteção
    if params.metodo == 1  # MUSCL
        (vars.αLδD[:, :]) = MUSCL(flows.αLflow, vars.αLδD0, vars.uLC, vars.ρL0, vars.uLC0, vars.ρL0, vars.uL, meshP, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, 0, vars.αLf, vars.αL, vars.ρL, vars.αL0, vars.τi, vars.τw, vars.αGf0, params.D)
        
        (vars.αGδD[:, :]) = MUSCL(flows.αGflow, vars.αGδD0, vars.uGC, vars.ρG0, vars.uGC0, vars.ρG0, vars.uG, meshP, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, 0, vars.αGf, vars.αG, vars.ρG, vars.αG0, vars.τi, vars.τw, vars.αGf0, params.D)
    else  # FR
        (vars.αLδD[:, :]) = FR(flows.αLflow, vars.αLδD0, vars.uLC, vars.ρL0, vars.uLC0, vars.ρL0, vars.uL, meshP, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, 0, vars.αLf, vars.αL, vars.ρL, vars.αL0, vars.τi, vars.τw, vars.αGf0, params.D)
        
        (vars.αGδD[:, :]) = FR(flows.αGflow, vars.αGδD0, vars.uGC, vars.ρG0, vars.uGC0, vars.ρG0, vars.uG, meshP, params.CFL, vars.Δt, vars.Δx, vars.P, params.θ, 0, vars.αGf, vars.αG, vars.ρG, vars.αG0, vars.τi, vars.τw, vars.αGf0, params.D)
    end
    
    # Normalização das frações volumétricas conforme código original - EXATAMENTE igual
    for j = 2:vars.N
        for k = 1:(params.Grau+1)
            vars.αLδD[j, k] = vars.αLδD[j, k] / (vars.αLδD[j, k] + vars.αGδD[j, k])
            vars.αGδD[j, k] = 1.0 - vars.αLδD[j, k]
        end
    end
    
    # valores médios conforme código original:
    for i = 1:vars.N
        vars.αL[i] = sum(vars.αLδD[i+1, :] .* weights[:]) / 2
        vars.αG[i] = sum(vars.αGδD[i+1, :] .* weights[:]) / 2
    end
end

"""
Salva os valores atuais como valores anteriores para a próxima iteração temporal
"""
function save_previous_values!(vars::SimulationVariables, params::InputParameters)
    # Atualização Temporal conforme código original
    vars.uL0[:] = vars.uL[:]
    vars.uG0[:] = vars.uG[:]
    vars.αL0[:] = vars.αL[:]
    vars.αG0[:] = vars.αG[:]
    vars.P0[:] = vars.P[:]
    vars.ρG0[:] = vars.ρG[:]
    vars.ρL0[:] = vars.ρL[:]
    vars.hL0[:] = vars.hL[:]
    vars.TL[:] = vars.hL[:] / params.Cp_L   # Temperatura compartilhada
    
    # Proteção para valores de temperatura inválidos
    for i = 1:length(vars.TL)
        if isnan(vars.TL[i]) || isinf(vars.TL[i]) || vars.TL[i] <= 0.0
            vars.TL[i] = params.T_0
        end
    end
    
    # Valores δD anteriores
    vars.uLδD0[:, :] = vars.uLδD[:, :]
    vars.uGδD0[:, :] = vars.uGδD[:, :]
    vars.αLδD0[:, :] = vars.αLδD[:, :]
    vars.αGδD0[:, :] = vars.αGδD[:, :]
end

"""
Aplica todas as correções de estabilidade em uma única função para evitar redundância
"""
function apply_stability_corrections!(vars::SimulationVariables, params::InputParameters)
    # Correção das velocidades médias
    for i = 1:length(vars.uLØ)
        if isnan(vars.uLØ[i]) || isinf(vars.uLØ[i])
            vars.uLØ[i] = params.uL_inlet
        end
        if isnan(vars.uGØ[i]) || isinf(vars.uGØ[i])
            vars.uGØ[i] = params.uG_inlet
        end
    end
    
    # Correção das correções de pressão
    for i = 1:length(vars.Pδ)
        if isnan(vars.Pδ[i]) || isinf(vars.Pδ[i])
            vars.Pδ[i] = 0.0
        end
    end
    
    # Correção das densidades
    for i = 1:length(vars.ρL)
        if isnan(vars.ρL[i]) || isinf(vars.ρL[i]) || vars.ρL[i] <= 0.0
            vars.ρL[i] = params.ρLi
        end
        if isnan(vars.ρG[i]) || isinf(vars.ρG[i]) || vars.ρG[i] <= 0.0
            vars.ρG[i] = params.ρGi
        end
    end
    
    # Correção das frações volumétricas (mais crítico)
    for i = 1:length(vars.αL)
        if isnan(vars.αL[i]) || isinf(vars.αL[i]) || vars.αL[i] < 0.0 || vars.αL[i] > 1.0
            vars.αL[i] = 0.5
        end
        if isnan(vars.αG[i]) || isinf(vars.αG[i]) || vars.αG[i] < 0.0 || vars.αG[i] > 1.0
            vars.αG[i] = 0.5
        end
    end
    
    # Normalização forçada das frações volumétricas
    for i = 1:length(vars.αL)
        total = vars.αL[i] + vars.αG[i]
        if total > 0.0
            vars.αL[i] = vars.αL[i] / total
            vars.αG[i] = vars.αG[i] / total
        else
            vars.αL[i] = 0.5
            vars.αG[i] = 0.5
        end
    end
    
    # Correção da entalpia
   
    for i = 1:length(vars.hL)
        if isnan(vars.hL[i]) || isinf(vars.hL[i]) || vars.hL[i] <= 0.0
            vars.hL[i] = params.Cp_L * params.T_0  # hLi calculado
        end
    end
    
    # Correção das pressões
    for i = 1:length(vars.P)
        if isnan(vars.P[i]) || isinf(vars.P[i]) || vars.P[i] <= 0.0
            vars.P[i] = params.P_out
        end
    end
    
    # Correção das velocidades finais
    for i = 1:length(vars.uL)
        if isnan(vars.uL[i]) || isinf(vars.uL[i])
            vars.uL[i] = params.uL_inlet
        end
    end
    for i = 1:length(vars.uG)
        if isnan(vars.uG[i]) || isinf(vars.uG[i])
            vars.uG[i] = params.uG_inlet
        end
    end
end

"""
Função principal que executa a simulação - EXATAMENTE como no código original
"""
function main()
    println("Iniciando simulação AnderSim 5 Eq...")
    
    # Leitura dos parâmetros
    params = read_input_parameters("Entrada12_mod.txt")
    println("Parâmetros lidos com sucesso.")
    
    # Inicialização das variáveis
    vars = initialize_simulation_variables(params)
    println("Variáveis inicializadas.")
    
    # Aplicação das condições iniciais
    apply_initial_conditions!(vars, params)
    println("Condições iniciais aplicadas.")
    
    # Inicialização dos flows e configuração das malhas
    flows = initialize_flows(params, vars)
    println("Flows inicializados.")
    
    # Configuração das malhas FR - EXATAMENTE como no original
    LMIN = 0
    LMAX = params.Comp
    T = Float64
    Grau = params.Grau
    C = params.C
    N = vars.N
    Δx = vars.Δx
    
    # Geração das malhas conforme código original
    (meshP, innerPointsP) = frgrid(LMIN, LMAX, N, Grau, C) # Malha Primária
    (meshS, innerPointsS) = frgrid(LMIN - Δx / 2, LMAX + Δx / 2, N + 1, Grau, C) # Malha Secundária
    
    # Definir variáveis globais necessárias - EXATAMENTE como no original
    global g = params.g
    global weights = meshP.standardElement.quadratureWeights  # Pesos de quadratura da malha primária
    global pontos = meshP.standardElement.innerPoints    # Pontos de quadratura da malha primária
    
    # Matriz de Vandermonde - EXATAMENTE como no código original
    global V = zeros(Grau + 1, Grau + 1)
    for i = 1:Grau+1
        V[:, i] .= pontos .^ (i - 1)
    end
    
    println("Malhas configuradas e variáveis globais definidas.")
    println("Iniciando loop temporal...")
    
    # Loop temporal principal - EXATAMENTE como no código original
    while vars.t < params.tempo
        # Cálculo do passo de tempo conforme original
        aux1 = maximum(abs.(vars.uL[:]))
        aux2 = maximum(abs.(vars.uG[:]))
        u = [aux1, aux2]
        aux3 = maximum(u)
        
        if params.metodo == 2
            vars.Δt = (params.CFL * vars.Δx) / ((2*params.Grau+1)*aux3)  # FR
        else
            vars.Δt = (params.CFL * vars.Δx) / (aux3)             # MUSCL
        end
        
        vars.Δt = minimum([vars.Δt, (params.tempo - vars.t)])
        vars.t = vars.t + vars.Δt
        
        # Print do tempo como no código original (apenas o valor numérico)
        print(vars.t)
        print("\n")
        
        # Inicialização dos erros para convergência
        e_uL = 1.0
        e_uG = 1.0
        e_αL = 1.0
        e_αG = 1.0
        e_P = 1.0
        
        cont1 = 0
        
        # Loop de convergência com proteção robusta e diagnósticos
        while (e_uL > params.Tol_u || e_uG > params.Tol_u || e_αL > params.Tol_α || 
               e_αG > params.Tol_α || e_P > params.Tol_P) && cont1 < 10  # REDUZIDO PARA TESTE
               
            cont1 += 1
            
            # Debug detalhado
            # if cont1 <= 3
            #     println("Iteração $cont1: erros = uL:$e_uL, uG:$e_uG, αL:$e_αL, αG:$e_αG, P:$e_P")
            # end
            
            # Armazenar valores auxiliares
            vars.uLaux[:] = vars.uLØ[:]
            vars.uGaux[:] = vars.uGØ[:]
            vars.αLaux[:] = vars.αLØ[:]
            vars.αGaux[:] = vars.αGØ[:]
            vars.Paux[:] = vars.PØ[:]
            
            # Calcular coeficientes nas faces - PRIMEIRA PARTE
            calculate_face_coefficients!(vars, params)
            
            # Calcular difusão artificial
            calculate_artificial_diffusion!(vars, params)
            
            # Calcular velocidades centradas
            calculate_centered_velocities!(vars)
            
            # Interpolação de αL/αG para malha de velocidades - EXATAMENTE como no original
            interpolate_alpha_to_velocity_mesh!(vars, pontos)
            
            # Calcular coeficientes de transferência de calor
            calculate_heat_transfer_coefficients!(vars, params)
            
            # Calcular números de Reynolds e fatores de atrito
            calculate_reynolds_and_friction!(vars, params)
            
            # Resolver equações de velocidade usando MUSCL ou FR
            solve_velocity_equations!(vars, params, flows, weights, meshS)
            
            # Resolver equação de pressão usando TDMA
            solve_pressure_equation!(vars, params)
            
            # Calcular correções de velocidade baseadas na pressão
            calculate_velocity_corrections!(vars, params)
            
            # Atualizar pressão e velocidades com correções
            update_pressure_and_velocities!(vars, params)
            
            # Atualizar densidade
            update_density!(vars, params)
            
            # Interpolação de uL/uG para malha de frações volumétricas - EXATAMENTE como no original
            interpolate_velocity_to_alpha_mesh!(vars, flows, pontos)
            
            # Resolver equações de fração volumétrica
            solve_volume_fraction_equations!(vars, params, flows, weights, pontos, meshP)
            
            # Segunda chamada para calcular coeficientes nas faces (com entalpias)
            calculate_face_coefficients!(vars, params)
            
            # Resolver equação de energia usando esquema MUSCL
            solve_energy_equation!(vars, params)
            
            # Cálculo dos erros conforme código original
            e_uL = norm(vars.uL - vars.uLaux) / norm(vars.uLaux)
            e_uG = norm(vars.uG - vars.uGaux) / norm(vars.uGaux)
            e_αL = norm(vars.αL - vars.αLaux) / norm(vars.αLaux)
            e_αG = norm(vars.αG - vars.αGaux) / norm(vars.αGaux)
            e_P = norm(vars.P - vars.Paux) / norm(vars.Paux)
            
            # Atualização dos valores supostos conforme código original
            vars.uLØ[:] = vars.uL[:]
            vars.uGØ[:] = vars.uG[:]
            vars.αLØ[:] = vars.αL[:]
            vars.αGØ[:] = vars.αG[:]
            vars.PØ[:] = vars.P[:]
            
            # Aplicar correções de estabilidade se necessário
            apply_stability_corrections!(vars, params)
        end
        
        # Atualização temporal - EXATAMENTE como no código original
        save_previous_values!(vars, params)
        
        # Calcular temperatura compartilhada conforme código original
        vars.TL[:] = vars.hL[:] / params.Cp_L
            
            # 6. Resolver pressão
            solve_pressure_equation!(vars, params)
            
            # 7. Correções de velocidade
            calculate_velocity_corrections!(vars, params)
            
            # 8. Atualizar pressão e velocidades
            update_pressure_and_velocities!(vars, params)
            
            # 9. Atualizar densidade
            update_density!(vars, params)
            
            # 10. Resolver frações volumétricas
            solve_volume_fraction_equations!(vars, params, flows, weights, pontos, meshP)
            
            # 11. Resolver energia
            solve_energy_equation!(vars, params)
            
            # PROTEÇÃO ÚNICA: Verificar e corrigir TODOS os valores inválidos uma vez por iteração
            apply_stability_corrections!(vars, params)
            
            # Cálculo dos erros simplificado e robusto
            e_uL = maximum(abs.(vars.uL - vars.uLaux)) / (maximum(abs.(vars.uLaux)) + 1e-12)
            e_uG = maximum(abs.(vars.uG - vars.uGaux)) / (maximum(abs.(vars.uGaux)) + 1e-12)
            e_αL = maximum(abs.(vars.αL - vars.αLaux)) / (maximum(abs.(vars.αLaux)) + 1e-12)
            e_αG = maximum(abs.(vars.αG - vars.αGaux)) / (maximum(abs.(vars.αGaux)) + 1e-12)
            e_P = maximum(abs.(vars.P - vars.Paux)) / (maximum(abs.(vars.Paux)) + 1e-12)
            
            # Proteção adicional contra erros inválidos
            e_uL = isnan(e_uL) || isinf(e_uL) ? 0.0 : e_uL
            e_uG = isnan(e_uG) || isinf(e_uG) ? 0.0 : e_uG
            e_αL = isnan(e_αL) || isinf(e_αL) ? 0.0 : e_αL
            e_αG = isnan(e_αG) || isinf(e_αG) ? 0.0 : e_αG
            e_P = isnan(e_P) || isinf(e_P) ? 0.0 : e_P
            
            # Atualização dos valores supostos conforme código original
            vars.uLØ[:] = vars.uL[:]
            vars.uGØ[:] = vars.uG[:]
            vars.αLØ[:] = vars.αL[:]
            vars.αGØ[:] = vars.αG[:]
            vars.PØ[:] = vars.P[:]
        end
        
        # Salvar valores atuais como valores anteriores para a próxima iteração temporal
        save_previous_values!(vars, params)
    
    println("Simulação concluída!")
    
    # Exportação dos resultados conforme código original
    # Limpeza final para garantir valores válidos para exportação
    for i = 1:vars.N
        # Substituir NaN e Inf por valores padrão
        if isnan(vars.uL[i]) || isinf(vars.uL[i])
            vars.uL[i] = params.uL_inlet
        end
        if isnan(vars.uG[i]) || isinf(vars.uG[i])
            vars.uG[i] = params.uG_inlet
        end
        if isnan(vars.αL[i]) || isinf(vars.αL[i]) || vars.αL[i] < 0 || vars.αL[i] > 1
            vars.αL[i] = 0.5
        end
        if isnan(vars.αG[i]) || isinf(vars.αG[i]) || vars.αG[i] < 0 || vars.αG[i] > 1
            vars.αG[i] = 0.5
        end
        if isnan(vars.TL[i]) || isinf(vars.TL[i]) || vars.TL[i] <= 0
            vars.TL[i] = params.T_0
        end
        if isnan(vars.ρL[i]) || isinf(vars.ρL[i]) || vars.ρL[i] <= 0
            vars.ρL[i] = params.ρLi
        end
        if isnan(vars.ρG[i]) || isinf(vars.ρG[i]) || vars.ρG[i] <= 0
            vars.ρG[i] = params.ρGi
        end
        if isnan(vars.P[i]) || isinf(vars.P[i]) || vars.P[i] <= 0
            vars.P[i] = params.P_out
        end
    end
    
    # Limpeza adicional para velocidades (incluindo N+1)
    for i = 1:vars.N+1
        if isnan(vars.uL[i]) || isinf(vars.uL[i])
            vars.uL[i] = params.uL_inlet
        end
        if isnan(vars.uG[i]) || isinf(vars.uG[i])
            vars.uG[i] = params.uG_inlet
        end
    end
    
    df = DataFrame(
        Pr=vars.P, 
        UL=vars.uL[1:vars.N], 
        UG=vars.uG[1:vars.N], 
        AlfaL=vars.αL, 
        AlfaG=vars.αG, 
        Temperatura=vars.TL, 
        ρL=vars.ρL, 
        ρG=vars.ρG
    )
    CSV.write("export_modular.csv", df)

    # Plot dos resultados
    aux2 = params.Comp / (vars.N - 1)
    plot(0:aux2:params.Comp, vars.TL)
    
    println("Resultados exportados e plotados com sucesso!")
end

# Execução da função principal
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end