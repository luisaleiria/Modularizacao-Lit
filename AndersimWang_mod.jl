##############################   AnderSim 5 Eq.  ###############################
###                   Compressível - Escrito todo em série                   ###
############################ By Anderson Viana #################################
#
include("./TDMA.jl")
include("./Limiter.jl")
using Printf
using Plots
using CSV
using DataFrames
using LinearAlgebra

include("hMLP.jl")
#export MLPu, LimiterMLP, trb_cell, modes
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
#
# Leitura do Arquivo de Entrada
n_dados = 37 # número de dados de entrada
dados = zeros(n_dados) # armazenador dos dados de entrada
texto = readlines("Entrada12_mod.txt")
for i = 2:2:2*n_dados
    aux = texto[i]
    aux2 = Int(i / 2)
    global dados[aux2] = parse(Float64, aux)
end
#
# Variáveis de Entrada
Comp = dados[1]           # Comprimento do duto
D = dados[2]              # Diâmetro interno do duto
θ = dados[3]              # Ângulo do duto
n_segmentos = dados[4]    # Número de seguimentos
tempo = dados[5]          # Tempo de simulação
uL_inlet = dados[6]       # Velocidade do líquido na entrada
uG_inlet = dados[7]       # Velocidade do gás na entrada
ρLi = dados[8]            # Densidade do Líquido
ρGi = dados[9]            # Densidade do gás
P_out = dados[10]         # Pressão na saída do domínio (Pa)
αG_inlet = dados[11]      # Fração de gás/óleo na entrada
αG_0 = dados[12]          # Fração de gás/óleo inicial no duto
T_intet = dados[13]       # Temperatura na entrada 
T_0 = dados[14]           # Condição inicial de temperatura   
CFL = dados[15]           # Coeficiente de estabilidade
Tol_P = dados[16]         # Tolerância para pressão
Tol_u = dados[17]         # Tolerância para velocidade
Tol_α = dados[18]         # Tolerância para fração volumetrica
#
# Variáveis de Entrada Adicionais (Energia)
Cp_L = dados[19]             # Capacidade térmica da fase líquida/água
Cp_G = dados[20]             # Capacidade térmica da fase gás/óleo
hL_intet = Cp_L * T_intet      # Entalpia de entrada fase líquida/água
hLi = Cp_L * T_0               # Entalpia inicial líquida/água
kf_L = dados[21]             # Condutividade da líquida/água
kf_G = dados[22]             # Condutividade do gás/óleo
T_ext = dados[23]            # Temperatura externa
μL = dados[24]               # Viscosidade do líquido/água (Pa.s)
μG = dados[25]               # Viscosidade do gás/óleo (Pa.s)
kf_pipe = dados[26]          # Condutividade do duto  
D_ext = dados[27]            # Diâmetro externo do duto
#
# Parâmetros
Λ = Cp_L / Cp_G                   # Razão entre os calores específicos
Pr_L = (Cp_L * μL) / kf_L     # Número de Prandtl para o líquido/água
#Pr_G = (Cp_G * μG) / kf_G    # Número de Prandtl para o gás
#
g = dados[28]                 # Aceleração da gravidade (m/s²)
cG = dados[29]                # velocidade do som no gás/óleo (m/s)
cL = dados[30]                # velocidade do som no líquido/água (m/s)
βG = dados[31]                 # Coeficiente de expansão volumétrica da gás/óleo
βL = dados[32]                 # Coeficiente de expansão volumétrica do líquido/óleo
γ = dados[33]
h_ext = dados[34]              # entalpia externa           
# Novos parâmetros vindos do TXT
metodo = Int(dados[35])  # 1 = MUSCL, 2 = FR
Grau = Int(dados[36])     # ordem do polinômio
C = Int(dados[37])        # esquema FR (0=DG,1=SD,2=Huynh,-1=Radau)

ρL_ref = ρLi                   # 2190400 # Densidade de referência
ρG_ref = ρGi                   # 10^(-10) # 120689.7
N = Int(n_segmentos)           # Número de nós da malha principal
A = (pi * D^2) / 4             # Área da seção (m²) - nesta versão é constante!
αL_inlet = 1 - αG_inlet        # Fração de líquido na entrada
αL_0 = 1 - αG_0                # fração de líquido inicial
u_inlet = [uL_inlet, uG_inlet] # Vetor das velocidades
Δx = Comp / n_segmentos        # Tamano do passo espacial
aux = maximum(u_inlet)
A = (π * D^2) / 4                  # Área da secção transvessal
S = π * D                        # Perímetro molhado
global αG_inlet
global αL_inlet
global uL_inlet
global uG_inlet
global hL_intet
global hLi
t = 0                       # Tempo inicial
LMIN = 0                    # Comprimento mínimo
LMAX = Comp                 # Comprimento máximo

T = Float64           # define precision
# Grau = 1  # (lido do TXT)
# C = -1
  # (lido do TXT)
# Variáveis no tempo n+1
P = zeros(N)         # Pressão na malha principal
αL = zeros(N)        # Fração do líquido na malha principal
αG = zeros(N)        # Fração do gás na malha principal
αLf = zeros(N + 1)     # Fração do líquido na face
αGf = zeros(N + 1)     # Fração do gás na face
uL = zeros(N + 1)
uG = zeros(N + 1)
ρL = zeros(N)
ρG = zeros(N)
ρLf = zeros(N + 1)
ρGf = zeros(N + 1)
hL = zeros(N)        # Entalapia da fase líquido/água  
TL = zeros(N)        # Temperatura na malha principal (como a temperatura é a mesma para as fase, escolhemos a temperatura do líquido)
hLf = zeros(N + 1)     # Entalapia da fase líquido/água na face
# 
# Variáveis no tempo n
P0 = zeros(N)         # Pressão no tempo atrasado 
αL0 = zeros(N)
αG0 = zeros(N)
αLf0 = zeros(N + 1)
αGf0 = zeros(N + 1)
uL0 = zeros(N + 1)
uG0 = zeros(N + 1)
ρL0 = zeros(N)
ρG0 = zeros(N)
ρLf0 = zeros(N + 1)
ρGf0 = zeros(N + 1)
hL0 = zeros(N)          # Entalapia da fase líquido/água atrasado
hLf0 = zeros(N + 1)       # Entalapia da fase líquido/água na face atrasado

## Valor suposto
PØ = zeros(N)
αLØ = zeros(N)
αGØ = zeros(N)
uLØ = zeros(N + 1)
uGØ = zeros(N + 1)
uLf = zeros(N)
uGf = zeros(N)
uLδ = zeros(N)        # Valor para correção
uGδ = zeros(N)
#
uLC = zeros(N)
uGC = zeros(N)
uLC0 = zeros(N)
uGC0 = zeros(N)

# números adimensionais E OUTROS
ReL = zeros(N)            # Número de Reynods para o Líquido
ReG = zeros(N)            # Número de Reynods para o Líquido 
R_núcleo = zeros(N)
fw = zeros(N)
fi = zeros(N)
λ = zeros(N)
τw = zeros(N + 1)
τi = zeros(N + 1)
Nu = zeros(N)             # Nusselt para a fase 2
ULM = zeros(N)            # Coeficiente convectivo para a fase 2
U_conv = zeros(N)         # Coeficiente convectivo da mistura
R_tot = zeros(N)          # Resistência a transferência de calor total

#
# Variáveis conserváveis
fi0 = zeros(N)
#
dif_art_L = zeros(N + 1)
dif_art_G = zeros(N + 1)
#

# Condições Iniciais
uLC[:] .= uL_inlet
uGC[:] .= uG_inlet
uLC0[:] .= uL_inlet
uGC0[:] .= uG_inlet

# Todos os elementos do vetores
uL[:] .= uL_inlet
uG[:] .= uG_inlet
αL[:] .= αL_0
αG[:] .= αG_0
ρL[:] .= ρLi
ρG[:] .= ρGi
P[:] .= P_out
hL[:] .= hLi             # Coloco hL0 
#
# Mapa de Valores
# Iniciais
uL0[:] .= uL_inlet
uG0[:] .= uG_inlet
αL0[:] .= αL_inlet
αG0[:] .= αG_inlet
ρL0[:] .= ρLi
ρG0[:] .= ρGi
P0[:] .= P_out
hL0[:] .= hLi            # Coloco hL0 


# Supostos
uLØ[:] .= uL_inlet
uGØ[:] .= uG_inlet
αLØ[:] .= αL0[:]
αGØ[:] .= αG0[:]
PØ[:] .= P_out
# Corretivos
uLδ = zeros(N + 1)
uGδ = zeros(N + 1)
Pδ = zeros(N)
Tδ = zeros(N)
#
# Inicialização dos erros relativos
e_uL = 1.0
e_uG = 1.0
e_αL = 1.0
e_αG = 1.0
e_P = 1.0
#
Paux = zeros(N)
αLaux = zeros(N)
αGaux = zeros(N)
uLaux = zeros(N + 1)
uGaux = zeros(N + 1)
##
##
# Para o FR:
(meshP, innerPointsP) = frgrid(LMIN, LMAX, N, Grau, C) # Malha Primária
(meshS, innerPointsS) = frgrid(LMIN - Δx / 2, LMAX + Δx / 2, N + 1, Grau, C) # Malha Secundária
# Initialize flow -----------------------
# uL:
f1(x) = uL_inlet
uL_IC = Function[f1]
uLflow = Flow(N + 1, Grau, 1, 1; T=T) # N+1 por ser a malha secundária 
init!(uLflow, uL_IC, innerPointsS)
uLflow.uδD[1, :, 1] .= uL_inlet        # correção "pé duro"
uLflow.uδD[end, :, 1] .= uL_inlet      # correção "pé duro"
uLδD = zeros(N + 3, Grau + 1)
uLδD0 = zeros(N + 3, Grau + 1)
uLδD[:, :] = uLflow.uδD[:, :, 1]
uLδD0[:, :] = uLδD[:, :]
# uG:
f2(x) = uG_inlet
uG_IC = Function[f2]
uGflow = Flow(N + 1, Grau, 1, 2; T=T) # N+1 por ser a malha secundária 
init!(uGflow, uG_IC, innerPointsS)
uGflow.uδD[1, :, 1] .= uG_inlet        # correção "pé duro"
uGflow.uδD[end, :, 1] .= uG_inlet      # correção "pé duro"
uGδD = zeros(N + 3, Grau + 1)
uGδD0 = zeros(N + 3, Grau + 1)
uGδD[:, :] = uGflow.uδD[:, :, 1]
uGδD0[:, :] = uGflow.uδD[:, :, 1]
# αL:
f3(x) = αL_0
αL_IC = Function[f3]
αLflow = Flow(N, Grau, 2, 3; T=T)
init!(αLflow, αL_IC, innerPointsP)
αLflow.uδD[1, :, 1] .= αL_inlet        # correção "pé duro"
αLflow.uδD[end, :, 1] .= αG_inlet      # correção "pé duro"
αLδD = zeros(N + 2, Grau + 1)
αLδD0 = zeros(N + 2, Grau + 1)
αLδD[:, :] = αLflow.uδD[:, :, 1]
αLδD0[:, :] = αLflow.uδD[:, :, 1]
# αG:
f4(x) = αG_0
αG_IC = Function[f4]
αGflow = Flow(N, Grau, 2, 4; T=T)
init!(αGflow, αG_IC, innerPointsP)
αGflow.uδD[1, :, 1] .= αG_inlet        # correção "pé duro"
αGflow.uδD[end, :, 1] .= αL_inlet      # correção "pé duro"
αGδD = zeros(N + 2, Grau + 1)
αGδD0 = zeros(N + 2, Grau + 1)
αGδD[:, :] = αGflow.uδD[:, :, 1]
αGδD0[:, :] = αGflow.uδD[:, :, 1]

# hL:
f5(x) = hLi                                                         #-> foi
hL_IC = Function[f5]                                                  #-> foi
hLflow = Flow(N, Grau, 3, 5; T=T)                                   #-> foi
init!(hLflow, hL_IC, innerPointsP)                                    #-> foi
hLflow.uδD[1, :, 1] .= hL_intet          # Correção "pé duro"           -> foi
hLflow.uδD[end, :, 1] .= hL_intet        # Correção "pé duro"           -> foi
hLδD = zeros(N + 2, Grau + 1)                                              #-> foi
hLδD0 = zeros(N + 2, Grau + 1)                                             #-> foi
hLδD[:, :] = hLflow.uδD[:, :, 1]                                         #-> foi
hLδD0[:, :] = hLflow.uδD[:, :, 1]                                        #-> foi
# Confirmar com Anderson se tá ok (meta de sexta- fer)



##
# Matriz de Vandermonde:
V = zeros(Grau + 1, Grau + 1)
for i = 1:Grau+1
    V[:, i] .= pontos .^ (i - 1)
end
##
# Vetores para resolução dos sistemas Lineares - Pressão
d_P = zeros(N) # diagonal da matriz A do sistema Ax = b
du_P = zeros(N) # diagonal inferior
dl_P = zeros(N) # diagonal superior
b_P = zeros(N)
# Condições de Contorno - Pressão
d_P[1] = 1.0
du_P[1] = 0
dl_P[1] = 0.0
b_P[1] = 0.0 #
d_P[N] = 1.0
du_P[N] = 0.0
dl_P[N] = 0.0
b_P[N] = 0.0       # Condição de Dirichlet - Não há variação - Lembrando que esse
# sistema resolve a variação da pressão!
#


αLu = zeros(N + 3, Grau + 1)
αGu = zeros(N + 3, Grau + 1)
αLu0 = zeros(N + 3, Grau + 1)
αGu0 = zeros(N + 3, Grau + 1)

uLα = zeros(N + 2, Grau + 1)
uGα = zeros(N + 2, Grau + 1)
uLα0 = zeros(N + 2, Grau + 1)
uGα0 = zeros(N + 2, Grau + 1)

# Processamento
while t < tempo # botar o tempo
    aux1 = maximum(abs.(uL[:]))
    aux2 = maximum(abs.(uG[:]))
    u = [aux1, aux2]
    aux3 = maximum(u) #+ cL
if metodo == 2
    global Δt = (CFL * Δx) / ((2*Grau+1)*aux3)  # FR
else
    global Δt = (CFL * Δx) / (aux3)             # MUSCL
end

    global t
    global Δt = minimum([Δt, (tempo - t)])
    global t = t + Δt
    print(t)
    print("\n")
    # Variáveis Globais
    # Colocar as variáveis implementadas como glogais
    global uL0
    global uG0
    global αL0
    global αG0
    global ρG0
    global ρL0
    global P0
    global uLØ
    global uGØ
    global αLØ
    global αGØ
    global PØ
    global uL
    global uG
    global αL
    global αG
    global ρG
    global ρL
    global P
    global hL

    #
    # Erros
    e_uL = 1
    e_uG = 1
    e_αL = 1
    e_αG = 1
    e_P = 1
    #
    #
    global cont1 = 0 
    while (e_uL > Tol_u || e_uG > Tol_u || e_αL > Tol_α || e_αG > Tol_α || e_P > Tol_P) && cont1 < 1000

        #print(cont)
        #print("\n")
        cont1 += 1
        #
        uLaux[:] = uLØ[:]
        uGaux[:] = uGØ[:]
        αLaux[:] = αLØ[:]
        αGaux[:] = αGØ[:]
        Paux[:] = PØ[:]
        #

        # Velocidades Supostas -------------------------------------------------
        ## preciso linearisar lá na equação de energia para o esquema de volumes finitos##
        # Coeficientes nas Faces:
        αLf[1] = αL[1]
        αGf[1] = αG[1]
        ρLf[1] = ρL[1]
        ρGf[1] = ρG[1]
        αLf0[1] = αL0[1]
        αGf0[1] = αG0[1]
        ρLf0[1] = ρL0[1]
        ρGf0[1] = ρG0[1]
        αLf[N+1] = αL[N]
        αGf[N+1] = αG[N]
        ρLf[N+1] = ρL[N]
        ρGf[N+1] = ρG[N]
        αLf0[N+1] = αL0[N]
        αGf0[N+1] = αG0[N]
        ρLf0[N+1] = ρL0[N]
        ρGf0[N+1] = ρG0[N]
        for i = 2:N
            αLf[i] = αL[i-1] * (sign(uLØ[i]) + 1) / 2 + αL[i] * (-sign(uLØ[i]) + 1) / 2
            αGf[i] = αG[i-1] * (sign(uGØ[i]) + 1) / 2 + αG[i] * (-sign(uGØ[i]) + 1) / 2
            αLf0[i] = αL0[i-1] * (sign(uL[i]) + 1) / 2 + αL0[i] * (-sign(uL[i]) + 1) / 2 #(αL0[i-1]+αL0[i])/2 #(αL0[i-1]+αL0[i])/2 #
            αGf0[i] = αG0[i-1] * (sign(uG[i]) + 1) / 2 + αG0[i] * (-sign(uG[i]) + 1) / 2 #(αG0[i-1]+αG0[i])/2 #(αG0[i-1]+αG0[i])/2 #
            ρLf[i] = ρL[i-1] * (sign(uL[i]) + 1) / 2 + ρL[i] * (-sign(uL[i]) + 1) / 2
            ρGf[i] = ρG[i-1] * (sign(uG[i]) + 1) / 2 + ρG[i] * (-sign(uG[i]) + 1) / 2
            ρLf0[i] = ρL0[i-1] * (sign(uL[i]) + 1) / 2 + ρL0[i] * (-sign(uL[i]) + 1) / 2
            ρGf0[i] = ρG0[i-1] * (sign(uG[i]) + 1) / 2 + ρG0[i] * (-sign(uG[i]) + 1) / 2
        end

        for i = 2:N
            Δp_art = (γ * αLf0[i] * ρLf0[i] * αGf0[i] * ρGf0[i] * (uG0[i] - uL0[i])^2) / (αGf0[i] * ρLf0[i] + αLf0[i] * ρGf0[i])
            dif_art_L[i] = Δp_art * (αL0[i] - αL0[i-1]) / Δx
            dif_art_G[i] = Δp_art * (αG0[i] - αG0[i-1]) / Δx
        end


        for i = 1:N
            uLC[i] = (uL[i] + uL[i+1]) / 2 #uL[i]*(sign((uL[i]+uL[i+1])/2)+1)/2 + uL[i+1]*(-sign((uL[i]+uL[i+1])/2)+1)/2
            uGC[i] = (uG[i] + uG[i+1]) / 2 #uG[i]*(sign((uG[i]+uG[i+1])/2)+1)/2 + uG[i+1]*(-sign((uG[i]+uG[i+1])/2)+1)/2
            uLC0[i] = (uL0[i] + uL0[i+1]) / 2 #uL0[i]*(sign((uL0[i]+uL0[i+1])/2)+1)/2 + uL0[i+1]*(-sign((uL0[i]+uL0[i+1])/2)+1)/2
            uGC0[i] = (uG0[i] + uG0[i+1]) / 2 #uG0[i]*(sign((uG0[i]+uG0[i+1])/2)+1)/2 + uG0[i+1]*(-sign((uG0[i]+uG0[i+1])/2)+1)/2
        end



        αLu[1, :] .= αL[1]
        αGu[1, :] .= αG[1]
        αLu0[1, :] .= αL0[1]
        αGu0[1, :] .= αG0[1]

        αLu[end-1:end, :] .= αL[end]
        αGu[end-1:end, :] .= αG[end]
        αLu0[end-1:end, :] .= αL0[end]
        αGu0[end-1:end, :] .= αG0[end]

        for i = 2:N+1
            for j = 1:Grau+1

                β1 = pontos[j] + 1
                β2 = pontos[j] - 1
                lj1 = lagrangeBasisAt(β1, pontos)
                lj2 = lagrangeBasisAt(β2, pontos)

                αLu[i, j] = sum(αLδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(αLδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                αGu[i, j] = sum(αGδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(αGδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                αLu0[i, j] = sum(αLδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(αLδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                αGu0[i, j] = sum(αGδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(αGδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))

            end
        end

        ReL[:] .= abs.((ρL[:] .* uLC[:] .* αL[:] * D) ./ μL) # Núcleo do Anular
        R_núcleo[:] .= (D / 2) * (abs.((αG[:] .* uLC[:]) ./ (uGC[:] + αG[:] .* (uLC[:] - uGC[:])))) .^ (1 / 2)
        ReG[:] .= abs.(2 * (ρG[:] .* uGC[:] .* R_núcleo[:]) ./ μG) # Anel do Anular

        fi[:] .= 0.079 ./ ((ReG[:] .+ 1e-12) .^ 0.25) # Para escoamento turbulento!!!!!
        τi[1:end-1] .= (fi[:] ./ 2) .* ρG[:] .* abs.(uGC[:] - uLC[:]) .* (uGC[:] - uLC[:])
        τi[1] = τi[2]

        fw[:] .= 0.046 ./ ((ReL[:] .+ 1e-12) .^ 0.2) # Para escoamento turbulento!!!!!
        λ[:] .= (1 .+ αG[:] + (2 * αG[:] .* log.(αG[:])) ./ (1 .- αG[:])) ./ (1 .- αG[:])
        τw[1:end-1] .= fw[:] .* ρL[:] .* uLC[:] .* abs.(uLC[:]) ./ 2 - (ρG[:] - ρL[:]) .* (-g) .* D .* αG[:] .* λ[:] ./ 4

        τi[end] = τi[end-1]
        τw[end] = τw[end-1]

        Nu[:] .= ((fw[:] ./ 8) .* (ReL[:] .- 1000) .* Pr_L) ./
                 (1.07 .+ 12.7 .* ((fw[:] ./ 8) .^ 0.5) .* ((Pr_L .^ (2 / 3)) .- 1))  # Nusselt tubulento para o líquido/água  
        ULM[:] .= (kf_L .* Nu[:]) ./ D                                   # Coeficiente convectivo para o líquido/água 



        if metodo == 1
            (uLδD[:, :]) = MUSCL(uLflow, uLδD0, αLf, ρLf, αLf0, ρLf0, uL, meshS, CFL, Δt, Δx, P, θ, dif_art_L, uLC0, αL, ρL, uL0, τi, τw, αGf0, D)
            
            (uGδD[:, :]) = MUSCL(uGflow, uGδD0, αGf, ρGf, αGf0, ρGf0, uG, meshS, CFL, Δt, Δx, P, θ, dif_art_G, uGC0, αG, ρG, uG0, τi, τw, αGf0, D)
        else
            (uLδD[:, :]) = FR(uLflow, uLδD0, αLf, ρLf, αLf0, ρLf0, uL, meshS, CFL, Δt, Δx, P, θ, dif_art_L, uLC0, αL, ρL, uL0, τi, τw, αGf0, D)

            (uGδD[:, :]) = FR(uGflow, uGδD0, αGf, ρGf, αGf0, ρGf0, uG, meshS, CFL, Δt, Δx, P, θ, dif_art_G, uGC0, αG, ρG, uG0, τi, τw, αGf0, D)
        end

        uLδD[3, :] = (uLδD[2, :] + uLδD[4, :]) ./ 2
        uGδD[3, :] = (uGδD[2, :] + uGδD[4, :]) ./ 2

        uLflow.uδD[:, :, 1] = uLδD[:, :]
        uGflow.uδD[:, :, 1] = uGδD[:, :]

        # valores médios:
        for i = 1:N+1
            uLØ[i] = sum(uLδD[i+1, :] .* weights[:]) / 2
            uGØ[i] = sum(uGδD[i+1, :] .* weights[:]) / 2
        end

        # Pressão --------------------------------------------------------------
        #
        d_P[N] = 1.0
        du_P[N] = 0.0
        dl_P[N] = 0.0
        b_P[N] = 0.0
        #
        for i = 2:(N-1)
            #
            # Coeficientes Pós Agrupamento (Tudo em upwind)
            if uLØ[i+1] >= 0
                αLen = αL0[i]
                ρLen = ρL0[i]
                αLwn = αL0[i-1]
                ρLwn = ρL0[i-1]
                αLe = αLØ[i]
                ρLe = ρL[i]
                αLw = αLØ[i-1]
                ρLw = ρL[i-1]
                a_Le = αLen * ρLen / Δt + αLen * ρLen * uL0[i+1] / Δx
                a_Lw = αLwn * ρLwn / Δt + αLwn * ρLwn * uL0[i] / Δx
            else
                αLen = αL0[i+1]
                ρLen = ρL0[i+1]
                αLwn = αL0[i]
                ρLwn = ρL0[i]
                αLe = αLØ[i+1]
                ρLe = ρL[i+1]
                αLw = αLØ[i]
                ρLw = ρL[i]
                a_Le = αLen * ρLen / Δt - αLen * ρLen * uL0[i+1] / Δx
                a_Lw = αLwn * ρLwn / Δt - αLwn * ρLwn * uL0[i] / Δx
            end
            if uGØ[i+1] >= 0
                αGen = αG0[i]
                ρGen = ρG0[i]
                αGwn = αG0[i-1]
                ρGwn = ρG0[i-1]
                αGe = αGØ[i]
                ρGe = ρG[i]
                αGw = αGØ[i-1]
                ρGw = ρG[i-1]
                a_Ge = αGen * ρGen / Δt + αGen * ρGen * uG0[i+1] / Δx
                a_Gw = αGwn * ρGwn / Δt + αGwn * ρGwn * uG0[i] / Δx
            else
                αGen = αG0[i+1]
                ρGen = ρG0[i+1]
                αGwn = αG0[i]
                ρGwn = ρG0[i]
                αGe = αGØ[i+1]
                ρGe = ρG[i+1]
                αGw = αGØ[i]
                ρGw = ρG[i]
                a_Ge = αGen * ρGen / Δt - αGen * ρGen * uG0[i+1] / Δx
                a_Gw = αGwn * ρGwn / Δt - αGwn * ρGwn * uG0[i] / Δx
            end
            mL = (αLe * ρLe * uLØ[i+1] - αLw * ρLw * uLØ[i]) / ρL_ref
            mG = (αGe * ρGe * uGØ[i+1] - αGw * ρGw * uGØ[i]) / ρG_ref
            coef = (Δx / Δt) * (αLØ[i] / (ρL_ref * cL^2) + αGØ[i] / (ρG_ref * cG^2))
            #
            d_P[i] = coef + αLe * ρLe * (αLe / (Δx * a_Le)) / ρL_ref + αGe * ρGe * (αGe / (Δx * a_Ge)) / ρG_ref + (αLw * ρLw * (αLw / (Δx * a_Lw)) / ρL_ref + αGw * ρGw * (αGw / (Δx * a_Gw)) / ρG_ref)
            du_P[i] = -(αLe * ρLe * (αLe / (Δx * a_Le)) / ρL_ref + αGe * ρGe * (αGe / (Δx * a_Ge)) / ρG_ref)
            dl_P[i] = -(αLw * ρLw * (αLw / (Δx * a_Lw)) / ρL_ref + αGw * ρGw * (αGw / (Δx * a_Gw)) / ρG_ref)
            b_P[i] = -(mL + mG + (Δx / Δt) * ((αLØ[i] * ρL[i] - αL0[i] * ρL0[i]) / ρL_ref + (αGØ[i] * ρG[i] - αG0[i] * ρG0[i]) / ρG_ref))
        end
        Pδ = TDMA(d_P, du_P, dl_P, b_P)
        #Pδ[1] = 2 * Pδ[2] - Pδ[3]
        Pδ[1] = Pδ[2] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAUSOU UNDERSHOOT
        #

        #
        # Correções --------------------------------------------------------
        for i = 1:(N-1)
            # Coeficientes Pós Agrupamento (Tudo em upwind)
            if uLØ[i+1] >= 0
                αLe = αL0[i]
                ρLe = ρL0[i]
                a_L = αLe * ρLe / Δt + αLe * ρLe * uL0[i+1] / Δx
            else
                αLe = αL0[i+1]
                ρLe = ρL0[i+1]
                a_L = αLe * ρLe / Δt - αLe * ρLe * uL0[i+1] / Δx
            end
            if uGØ[i+1] >= 0
                αGe = αG0[i]
                ρGe = ρG0[i]
                a_G = αGe * ρGe / Δt + αGe * ρGe * uG0[i+1] / Δx
            else
                αGe = αG0[i+1]
                ρGe = ρG0[i+1]
                a_G = αGe * ρGe / Δt - αGe * ρGe * uG0[i+1] / Δx
            end
            de_L = αLe / (Δx * a_L)
            de_G = αGe / (Δx * a_G)
            uLδ[i+1] = de_L * (Pδ[i] - Pδ[i+1])
            uGδ[i+1] = de_G * (Pδ[i] - Pδ[i+1])
        end
        #
        P[:] = PØ[:] + Pδ[:] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        uL = uLØ + uLδ
        uG = uGØ + uGδ

        #
        for i = 2:(N+2)
            uLflow.uδD[i, :, 1] .+= uLδ[i-1]
            uGflow.uδD[i, :, 1] .+= uGδ[i-1]
        end
        uLflow.uδD[N+3, :, 1] = uLflow.uδD[N+2, :, 1]
        uGflow.uδD[N+3, :, 1] = uGflow.uδD[N+2, :, 1]
        #
        # Contorno:
        uL[N+1] = uL[N]
        uG[N+1] = uG[N]
        #

        # Variação da temperatuara  


        for i = 2:N
            Tδ[i] = (1 / Cp_L) * (hL[i] - hL0[i])
        end
       # Tδ[2] = 2*Tδ[3] - Tδ[4]
        #
        #Correção da Densidade: Para o caso two-fluid five equations de Eduardo
        for i = 1:N
            global ρL[i] = ρL[i] + (1 / cL^2) * Pδ[i] - ρL[i] * βL * Tδ[i]
            global ρG[i] = ρG[i] + (1 / cG^2) * Pδ[i] - ρG[i] * βG * Tδ[i]
        end
        #

        # Alfas ------------------------------------------------------------
        #

        for i = 1:N+2
            for j = 1:Grau+1

                β1 = pontos[j] + 1
                β2 = pontos[j] - 1
                lj1 = lagrangeBasisAt(β1, pontos)
                lj2 = lagrangeBasisAt(β2, pontos)

                uLα[i, j] = sum(uLflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(uLflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                uGα[i, j] = sum(uGflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(uGflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                uLα0[i, j] = sum(uLδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(uLδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                uGα0[i, j] = sum(uGδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(uGδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
            end
        end

        # Coeficientes nos Centros:
        for i = 1:N
            uLC[i] = (uL[i] + uL[i+1]) / 2#uL[i]*(sign((uL[i]+uL[i+1])/2)+1)/2 + uL[i+1]*(-sign((uL[i]+uL[i+1])/2)+1)/2
            uGC[i] = (uG[i] + uG[i+1]) / 2#uG[i]*(sign((uG[i]+uG[i+1])/2)+1)/2 + uG[i+1]*(-sign((uG[i]+uG[i+1])/2)+1)/2
            uLC0[i] = (uL0[i] + uL0[i+1]) / 2#uL0[i]*(sign((uL0[i]+uL0[i+1])/2)+1)/2 + uL0[i+1]*(-sign((uL0[i]+uL0[i+1])/2)+1)/2
            uGC0[i] = (uG0[i] + uG0[i+1]) / 2#uG0[i]*(sign((uG0[i]+uG0[i+1])/2)+1)/2 + uG0[i+1]*(-sign((uG0[i]+uG0[i+1])/2)+1)/2
        end
        #

        if metodo == 1
            (αLδD[:, :]) = MUSCL(αLflow, αLδD0, uLC, ρL0, uLC0, ρL0, uL, meshP, CFL, Δt, Δx, P, θ, 0, αLf, αL, ρL, αL0, τi, τw, αGf0, D)

            (αGδD[:, :]) = MUSCL(αGflow, αGδD0, uGC, ρG0, uGC0, ρG0, uG, meshP, CFL, Δt, Δx, P, θ, 0, αGf, αG, ρG, αG0, τi, τw, αGf0, D)
        else
            (αLδD[:, :]) = FR(αLflow, αLδD0, uLC, ρL0, uLC0, ρL0, uL, meshP, CFL, Δt, Δx, P, θ, 0, αLf, αL, ρL, αL0, τi, τw, αGf0, D)

            (αGδD[:, :]) = FR(αGflow, αGδD0, uGC, ρG0, uGC0, ρG0, uG, meshP, CFL, Δt, Δx, P, θ, 0, αGf, αG, ρG, αG0, τi, τw, αGf0, D)
        end
        #
        for j = 2:N
            for k = 1:(Grau+1)
                αLδD[j, k] = αLδD[j, k] / (αLδD[j, k] + αGδD[j, k])
                αGδD[j, k] = 1.0 - αLδD[j, k]
            end
        end
        #
        # valores médios:
        for i = 1:N
            αL[i] = sum(αLδD[i+1, :] .* weights[:]) / 2
            αG[i] = sum(αGδD[i+1, :] .* weights[:]) / 2
        end



        # Coeficientes nas Faces:
        αLf[1] = αL[1]
        αGf[1] = αG[1]
        ρLf[1] = ρL[1]
        ρGf[1] = ρG[1]
        αLf0[1] = αL0[1]
        αGf0[1] = αG0[1]
        ρLf0[1] = ρL0[1]
        ρGf0[1] = ρG0[1]
        hLf[1] = hL[1]
        hLf0[1] = hL0[1]      # entalpias nas interfaces de entrada

        αLf[N+1] = αL[N]
        αGf[N+1] = αG[N]
        ρLf[N+1] = ρL[N]
        ρGf[N+1] = ρG[N]
        αLf0[N+1] = αL0[N]
        αGf0[N+1] = αG0[N]
        ρLf0[N+1] = ρL0[N]
        ρGf0[N+1] = ρG0[N]
        hLf[N+1] = hL[N]
        hLf0[N+1] = hL0[N]  # entalpias na última interface

        for i = 2:N
            # Frações volumétricas
            αLf[i] = αL[i-1] * (sign(uLØ[i]) + 1) / 2 + αL[i] * (-sign(uLØ[i]) + 1) / 2
            αGf[i] = αG[i-1] * (sign(uGØ[i]) + 1) / 2 + αG[i] * (-sign(uGØ[i]) + 1) / 2
            αLf0[i] = αL0[i-1] * (sign(uL[i]) + 1) / 2 + αL0[i] * (-sign(uL[i]) + 1) / 2
            αGf0[i] = αG0[i-1] * (sign(uG[i]) + 1) / 2 + αG0[i] * (-sign(uG[i]) + 1) / 2

            # Densidades
            ρLf[i] = ρL[i-1] * (sign(uL[i]) + 1) / 2 + ρL[i] * (-sign(uL[i]) + 1) / 2
            ρGf[i] = ρG[i-1] * (sign(uG[i]) + 1) / 2 + ρG[i] * (-sign(uG[i]) + 1) / 2
            ρLf0[i] = ρL0[i-1] * (sign(uL[i]) + 1) / 2 + ρL0[i] * (-sign(uL[i]) + 1) / 2
            ρGf0[i] = ρG0[i-1] * (sign(uG[i]) + 1) / 2 + ρG0[i] * (-sign(uG[i]) + 1) / 2

            # Entalpias
            hLf[i] = hL[i-1] * (sign(uL[i]) + 1) / 2 + hL[i] * (-sign(uL[i]) + 1) / 2
            hLf0[i] = hL0[i-1] * (sign(uL[i]) + 1) / 2 + hL0[i] * (-sign(uL[i]) + 1) / 2
        end


        ######### EQUAÇÃO DE ENERGIA #########

        U_conv[:] .= ULM[:] .* (1 ./ (1 .- αG[:])) .^ 0.9             # Cálculo do coeficiente de transferência conv. da mistura
        R_tot[:] .= (1 ./ (U_conv[:] .* D * π .* Δx)) .+
                    log(D_ext / D) ./ (2 * π * kf_pipe * Δx) .+ (1 ./ (h_ext .* (Δx .* π .* D_ext)))
        # Restência térmica total

        # MUSCL:
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        u = zeros(N)
        u[:] = hL[:]


        hL[1] = hL_intet # Condição de contorno

        for i = 2:N

            if i == 2
                r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
                r_i_menos = 0
                r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
                uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
                uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
            elseif i == (N - 1)
                r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
                r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
                r_i_mais = 0
                uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                uRmais = u[i+1]
                uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
            elseif i == (N)
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

            F_mais = ((αLf[i+1] * ρLf[i+1] * uL[i+1] * (uLmais + 0.5 * uL[i+1]^2))+(αGf[i+1] * ρGf[i+1] * uG[i+1] * (Λ * uLmais + 0.5 * uG[i+1]^2)))*(sign((uL[i]+uL[i+1])/2)+1)/2 
                    + ((αLf[i] * ρLf[i] * uL[i] * (uRmais + 0.5 * uL[i]^2))+(αGf[i] * ρGf[i] * uG[i] * (Λ * uRmais + 0.5 * uG[i]^2)))*(-sign((uL[i]+uL[i+1])/2)+1)/2 # ROE

            F_menos = ((αLf[i+1] * ρLf[i+1] * uL[i+1] * (uLmenos + 0.5 * uL[i+1]^2))+(αGf[i+1] * ρGf[i+1] * uG[i+1] * (Λ * uLmenos + 0.5 * uG[i+1]^2)))*(sign((uL[i]+uL[i+1])/2)+1)/2 
                    + ((αLf[i] * ρLf[i] * uL[i] * (uRmenos + 0.5 * uL[i]^2))+(αGf[i] * ρGf[i] * uG[i] * (Λ * uRmenos + 0.5 * uG[i]^2)))*(-sign((uL[i]+uL[i+1])/2)+1)/2 # ROE

            #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
               
            # Fluxo líquido na célula i
            a3 = (-Δt / Δx) * (F_mais - F_menos)


            a1 = -0.5 * (αG[i] * ρG[i] * uG[i]^2 + αL[i] * ρL[i] * uL[i]^2)
            a2 = αG0[i] * ρG0[i] * (Λ * hL0[i] + 0.5 * uG0[i]^2) +
                 αL0[i] * ρL0[i] * (hL0[i] + 0.5 * uL0[i]^2)
            a4 = P[i] - PØ[i]
            a5 = Δt * g * sin(θ) * (αL0[i] * ρL0[i] * uL0[i] + αG0[i] * ρG0[i] * uG0[i])
            an = hL0[i] / Cp_L
            if an >= T_ext
                a6 = (Δt * hL0[i] * S) / (Cp_L * R_tot[i] * A)
                a7 = (-T_ext * Δt * S) / (R_tot[i] * A)
            else
                a6 = (-Δt * hL0[i] * S) / (Cp_L * R_tot[i] * A)
                a7 = (T_ext * Δt * S) / (R_tot[i] * A)
            end
            a8 = αG[i] * ρG[i] * Λ + αL[i] * ρL[i]
            hL[i] = (a1 + a2 + a3 + a4 + a5 + a6 + a7) / a8 # Entalpia 
            #
        end
        # Erros:  ##############################################################
        e_uL = norm(uL - uLaux) / norm(uLaux)
        e_uG = norm(uG - uGaux) / norm(uGaux)
        e_αL = norm(αL - αLaux) / norm(αLaux)
        e_αG = norm(αG - αGaux) / norm(αGaux)
        e_P = norm(P - Paux) / norm(Paux)
        #
        # Atualização dos Supostos -----------------------------------------
        uLØ[:] = uL[:]
        uGØ[:] = uG[:]
        αLØ[:] = αL[:]
        αGØ[:] = αG[:]
        PØ[:] = P[:]

    end
    #
    # Atualização Temporal -------------------------------------------------

    uL0[:] = uL[:]
    uG0[:] = uG[:]
    αL0[:] = αL[:]
    αG0[:] = αG[:]
    P0[:] = P[:]
    ρG0 = ρG
    ρL0 = ρL
    hL0[:] = hL[:]
    TL[:] = hL[:] / Cp_L   # Temperatura comartilhada
    #
    #
    uLδD0[:, :] = uLδD[:, :]
    uGδD0[:, :] = uGδD[:, :]
    αLδD0[:, :] = αLδD[:, :]
    αGδD0[:, :] = αGδD[:, :]

end
#

df = DataFrame(Pr=P, UL=uL[1:N], UG=uG[1:N], AlfaL=αL, AlfaG=αG, Temperatura=TL, ρL=ρL, ρG=ρG)
CSV.write("export.csv", df)
#
aux2 = Comp / (N - 1)
plot(0:aux2:Comp, TL)
