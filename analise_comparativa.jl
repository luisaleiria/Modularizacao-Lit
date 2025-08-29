"""
Análise Comparativa AnderSim: Original vs Modular
================================================
Este script gera gráficos comparativos entre os resultados dos códigos
original e modular do AnderSim, analisando as diferenças encontradas.
"""

using CSV
using DataFrames
using Plots
using Statistics

function load_results()
    """Carrega os resultados dos dois códigos"""
    println("📊 Carregando resultados...")
    
    # Carregar dados
    original = CSV.read("export.csv", DataFrame)
    modular = CSV.read("export_modular.csv", DataFrame)
    
    println("✅ Dados carregados:")
    println("   Original: $(nrow(original)) pontos")
    println("   Modular:  $(nrow(modular)) pontos")
    
    return original, modular
end

function calculate_differences(original, modular)
    """Calcula diferenças entre os resultados"""
    println("🔍 Calculando diferenças...")
    
    differences = DataFrame()
    
    for col in names(original)
        if col in names(modular)
            diff_abs = abs.(original[!, col] - modular[!, col])
            diff_rel = diff_abs ./ abs.(original[!, col]) * 100
            
            differences[!, "$(col)_abs"] = diff_abs
            differences[!, "$(col)_rel"] = diff_rel
            
            println("   $col: máx=$(maximum(diff_rel))%, média=$(mean(diff_rel))%")
        end
    end
    
    return differences
end

function plot_comparison_results(original, modular, differences)
    """Gera gráficos comparativos"""
    println("📈 Gerando gráficos...")
    
    # Configurar posições (assumindo espaçamento uniforme)
    N = nrow(original)
    positions = 1:N
    
    # 1. Gráfico de velocidades
    p1 = plot(positions, original.UL, label="UL Original", lw=2, color=:blue)
    plot!(p1, positions, modular.UL, label="UL Modular", lw=2, color=:lightblue, linestyle=:dash)
    plot!(p1, positions, original.UG, label="UG Original", lw=2, color=:red)
    plot!(p1, positions, modular.UG, label="UG Modular", lw=2, color=:pink, linestyle=:dash)
    xlabel!(p1, "Posição")
    ylabel!(p1, "Velocidade [m/s]")
    title!(p1, "Comparação de Velocidades")
    
    # 2. Gráfico de densidades
    p2 = plot(positions, original.ρL, label="ρL Original", lw=2, color=:blue)
    plot!(p2, positions, modular.ρL, label="ρL Modular", lw=2, color=:lightblue, linestyle=:dash)
    plot!(p2, positions, original.ρG, label="ρG Original", lw=2, color=:red)
    plot!(p2, positions, modular.ρG, label="ρG Modular", lw=2, color=:pink, linestyle=:dash)
    xlabel!(p2, "Posição")
    ylabel!(p2, "Densidade [kg/m³]")
    title!(p2, "Comparação de Densidades")
    
    # 3. Gráfico de frações volumétricas
    p3 = plot(positions, original.AlfaL, label="αL Original", lw=2, color=:green)
    plot!(p3, positions, modular.AlfaL, label="αL Modular", lw=2, color=:lightgreen, linestyle=:dash)
    plot!(p3, positions, original.AlfaG, label="αG Original", lw=2, color=:orange)
    plot!(p3, positions, modular.AlfaG, label="αG Modular", lw=2, color=:yellow, linestyle=:dash)
    xlabel!(p3, "Posição")
    ylabel!(p3, "Fração Volumétrica [-]")
    title!(p3, "Comparação de Frações Volumétricas")
    
    # 4. Gráfico de temperatura e pressão
    p4 = plot(positions, original.Temperatura, label="T Original", lw=2, color=:purple)
    plot!(p4, positions, modular.Temperatura, label="T Modular", lw=2, color=:magenta, linestyle=:dash)
    xlabel!(p4, "Posição")
    ylabel!(p4, "Temperatura [°C]")
    title!(p4, "Comparação de Temperatura")
    
    # 5. Gráfico de diferenças relativas - Velocidades
    p5 = plot(positions, differences.UL_rel, label="UL", lw=2, color=:blue)
    plot!(p5, positions, differences.UG_rel, label="UG", lw=2, color=:red)
    xlabel!(p5, "Posição")
    ylabel!(p5, "Diferença Relativa [%]")
    title!(p5, "Diferenças Relativas - Velocidades")
    
    # 6. Gráfico de diferenças relativas - Densidades
    p6 = plot(positions, differences.ρL_rel, label="ρL", lw=2, color=:blue)
    plot!(p6, positions, differences.ρG_rel, label="ρG", lw=2, color=:red)
    xlabel!(p6, "Posição")
    ylabel!(p6, "Diferença Relativa [%]")
    title!(p6, "Diferenças Relativas - Densidades")
    
    # Combinar gráficos
    comparison_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1200,800))
    differences_plot = plot(p5, p6, layout=(1,2), size=(1000,400))
    
    # Salvar gráficos
    savefig(comparison_plot, "comparacao_resultados.png")
    savefig(differences_plot, "diferencas_relativas.png")
    
    println("✅ Gráficos salvos:")
    println("   📊 comparacao_resultados.png")
    println("   📈 diferencas_relativas.png")
    
    return comparison_plot, differences_plot
end

function generate_statistics_table(original, modular, differences)
    """Gera tabela de estatísticas"""
    println("📋 Gerando estatísticas...")
    
    stats = DataFrame(
        Variavel = String[],
        Original_Min = Float64[],
        Original_Max = Float64[],
        Original_Media = Float64[],
        Modular_Min = Float64[],
        Modular_Max = Float64[],
        Modular_Media = Float64[],
        Diff_Max_Abs = Float64[],
        Diff_Max_Rel = Float64[],
        Diff_Media_Rel = Float64[]
    )
    
    variables = ["Pr", "UL", "UG", "AlfaL", "AlfaG", "Temperatura", "ρL", "ρG"]
    
    for var in variables
        if var in names(original) && var in names(modular)
            push!(stats, [
                var,
                minimum(original[!, var]),
                maximum(original[!, var]),
                mean(original[!, var]),
                minimum(modular[!, var]),
                maximum(modular[!, var]),
                mean(modular[!, var]),
                maximum(differences[!, "$(var)_abs"]),
                maximum(differences[!, "$(var)_rel"]),
                mean(differences[!, "$(var)_rel"])
            ])
        end
    end
    
    # Salvar tabela
    CSV.write("estatisticas_comparacao.csv", stats)
    println("✅ Estatísticas salvas em: estatisticas_comparacao.csv")
    
    return stats
end

function main()
    """Função principal"""
    println("🎯 ANÁLISE COMPARATIVA ANDERSIM")
    println("=" ^ 50)
    
    # Carregar dados
    original, modular = load_results()
    
    # Calcular diferenças
    differences = calculate_differences(original, modular)
    
    # Gerar gráficos
    comp_plot, diff_plot = plot_comparison_results(original, modular, differences)
    
    # Gerar estatísticas
    stats = generate_statistics_table(original, modular, differences)
    
    # Resumo final
    println("\n📊 RESUMO DA ANÁLISE:")
    println("=" ^ 30)
    
    max_diffs = []
    for var in ["UL", "UG", "ρL", "ρG", "Temperatura"]
        if "$(var)_rel" in names(differences)
            max_diff = maximum(differences[!, "$(var)_rel"])
            mean_diff = mean(differences[!, "$(var)_rel"])
            push!(max_diffs, max_diff)
            println("   $var: máx=$(round(max_diff, digits=2))%, média=$(round(mean_diff, digits=2))%")
        end
    end
    
    overall_max = maximum(max_diffs)
    println("\n🎯 Maior diferença encontrada: $(round(overall_max, digits=2))%")
    
    if overall_max < 5
        println("✅ Diferenças pequenas - códigos equivalentes")
    elseif overall_max < 25
        println("⚠️  Diferenças moderadas - típicas de códigos modularizados")
    else
        println("❌ Diferenças significativas - revisar implementação")
    end
    
    println("\n✅ Análise completa!")
    println("📁 Arquivos gerados:")
    println("   📊 comparacao_resultados.png")
    println("   📈 diferencas_relativas.png") 
    println("   📋 estatisticas_comparacao.csv")
    
    return stats, comp_plot, diff_plot
end

# Executar análise
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
