#!/usr/bin/env julia
"""
Script para comparar os resultados do AndersimWang_mod.jl vs AndersimWang_modular.jl
"""

using CSV
using DataFrames
using Statistics

println("=== COMPARAÇÃO DOS RESULTADOS ANDERSIM ===")

# Carrega os dados
try
    original = CSV.read("export.csv", DataFrame)
    modular = CSV.read("export_modular.csv", DataFrame)
    
    println("✓ Arquivos carregados com sucesso!")
    println("  - Original: $(nrow(original)) linhas")
    println("  - Modular:  $(nrow(modular)) linhas")
    
    # Verifica se têm o mesmo número de linhas
    if nrow(original) != nrow(modular)
        println("❌ ERRO: Número de linhas diferente!")
        println("   Original: $(nrow(original)), Modular: $(nrow(modular))")
        exit(1)
    end
    
    # Verifica colunas
    cols_orig = names(original)
    cols_mod = names(modular)

    println("\n--- COLUNAS ---")
    println("Original: $cols_orig")
    println("Modular:  $cols_mod")
    
    if cols_orig != cols_mod
        println("⚠️  AVISO: Nomes das colunas são diferentes!")
    else
        println("✓ Colunas idênticas")
    end
    
    # Comparação numérica
    println("\n--- COMPARAÇÃO NUMÉRICA ---")
    
    tolerancia = 1e-10
    diferencas_encontradas = false
    
    for col in cols_orig
        if col in cols_mod
            diff = abs.(original[!, col] - modular[!, col])
            max_diff = maximum(diff)
            mean_diff = mean(diff)
            
            if max_diff > tolerancia
                diferencas_encontradas = true
                println("❌ $col: Diferença máxima = $max_diff, média = $mean_diff")
                
                # Mostra onde estão as maiores diferenças
                idx_max = argmax(diff)
                println("   Maior diferença na linha $idx_max:")
                println("   Original: $(original[idx_max, col])")
                println("   Modular:  $(modular[idx_max, col])")
            else
                println("✓ $col: Máx diff = $max_diff (OK)")
            end
        end
    end
    
    # Resumo final
    println("\n=== RESUMO ===")
    if !diferencas_encontradas
        println("🎉 RESULTADOS IDÊNTICOS!")
        println("   Os dois códigos produzem exatamente os mesmos resultados!")
    else
        println("⚠️  DIFERENÇAS ENCONTRADAS!")
        println("   Verifique os valores acima para identificar discrepâncias.")
    end
    
    # Estatísticas básicas
    println("\n--- ESTATÍSTICAS BÁSICAS (MODULAR) ---")
    for col in cols_mod
        if eltype(modular[!, col]) <: Number
            vals = modular[!, col]
            println("$col: min=$(minimum(vals)), max=$(maximum(vals)), média=$(mean(vals))")
        end
    end
    
catch e
    println("❌ ERRO ao carregar arquivos:")
    println("   $e")
    println("\nVerifique se os arquivos export.csv e export_modular.csv existem.")
end
