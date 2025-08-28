# AnderSim 5 Eq. - Versão Modularizada

Esta é a versão modularizada do código AnderSim para simulação de escoamento bifásico compressível.

## Estrutura do Código

### Arquivos Principais

- `AndersimWang_modular.jl` - Código principal modularizado
- `AndersimWang_mod.jl` - Código original (versão monolítica)
- `exemplo_uso.jl` - Exemplo de como usar o código modularizado

### Dependências

O código requer os seguintes arquivos:
- `TDMA.jl` - Algoritmo tri-diagonal
- `Limiter.jl` - Limitadores para MUSCL
- `hMLP.jl` - Funções MLP
- `utils.jl` - Funções utilitárias
- `solver.jl` - Solucionadores
- `polynomials.jl` - Funções polinomiais
- `mesh.jl` - Estruturas de malha
- `flow.jl` - Estruturas de fluxo
- `frgrid.jl` - Malhas FR
- `FR.jl` - Métodos Flux Reconstruction
- `MUSCL.jl` - Métodos MUSCL
- `Entrada12_mod.txt` - Arquivo de entrada com parâmetros

## Estrutura Modular

### 1. Estruturas de Dados

#### `InputParameters`
Contém todos os parâmetros de entrada da simulação:
- Parâmetros geométricos (comprimento, diâmetro, ângulo)
- Parâmetros temporais (tempo de simulação, CFL)
- Propriedades dos fluidos
- Tolerâncias de convergência
- Propriedades térmicas e físicas
- Configuração do método numérico

#### `SimulationVariables`
Contém todas as variáveis da simulação:
- Variáveis no tempo atual e anterior
- Valores supostos para iteração
- Números adimensionais
- Sistemas lineares
- Variáveis auxiliares

### 2. Funções Principais

#### Inicialização
- `read_input_parameters()` - Lê parâmetros do arquivo de entrada
- `initialize_simulation_variables()` - Inicializa estruturas de dados
- `apply_initial_conditions!()` - Aplica condições iniciais

#### Cálculos Auxiliares
- `calculate_time_step()` - Calcula passo de tempo baseado no CFL
- `calculate_face_coefficients!()` - Calcula coeficientes nas faces
- `calculate_artificial_diffusion!()` - Calcula difusão artificial
- `calculate_centered_velocities!()` - Calcula velocidades centradas
- `calculate_reynolds_and_friction!()` - Calcula Reynolds e atrito

#### Resolução de Equações
- `solve_pressure_equation!()` - Resolve equação de pressão (TDMA)
- `calculate_velocity_corrections!()` - Calcula correções de velocidade
- `update_pressure_and_velocities!()` - Atualiza pressão e velocidades
- `update_density!()` - Atualiza densidade
- `solve_energy_equation!()` - Resolve equação de energia (MUSCL)
- `solve_velocity_equations!()` - Resolve equações de velocidade
- `solve_volume_fraction_equations!()` - Resolve frações volumétricas

#### Função Principal
- `main()` - Executa simulação completa

## Como Usar

### Execução Simples
```julia
julia AndersimWang_modular.jl
```

### Usando o Exemplo
```julia
# Simulação completa
julia exemplo_uso.jl

# Teste de leitura de parâmetros
julia exemplo_uso.jl parametros

# Teste de inicialização
julia exemplo_uso.jl init
```

### Uso Programático
```julia
include("AndersimWang_modular.jl")

# Leitura de parâmetros
params = read_input_parameters("Entrada12_mod.txt")

# Inicialização
vars = initialize_simulation_variables(params)
apply_initial_conditions!(vars, params)

# Execução da simulação
main()
```


## Notas de Implementação

- As estruturas são mutáveis para permitir modificações durante a simulação
- Funções com `!` modificam seus argumentos (convenção Julia)
- Tratamento de casos especiais (contornos, divisão por zero, etc.)
