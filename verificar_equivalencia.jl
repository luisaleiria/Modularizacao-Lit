#!/usr/bin/env julia

using CSV, DataFrames

println("=== VERIFICAÇÃO DE EQUIVALÊNCIA ENTRE CÓDIGO ORIGINAL E MODULAR ===")
println()

# Executar código original
println("1. Executando código original...")
run(`julia AndersimWang_mod.jl`)
run(`mv export.csv export_original_final.csv`)

# Executar código modular
println("\n2. Executando código modular...")
run(`julia AndersimWang_modular.jl`)
run(`mv export.csv export_modular_final.csv`)

# Comparar resultados
println("\n3. Comparando resultados...")

try
    df_original = CSV.read("export_original_final.csv", DataFrame)
    df_modular = CSV.read("export_modular_final.csv", DataFrame)
    
    println("Dados originais: $(size(df_original)) linhas x colunas")
    println("Dados modulares: $(size(df_modular)) linhas x colunas")
    
    if size(df_original) == size(df_modular)
        println("✅ Dimensões idênticas!")
        
        # Comparar valores numéricos
        diferenca_maxima = 0.0
        total_diferencas = 0
        
        for col in names(df_original)
            if eltype(df_original[!, col]) <: Number && eltype(df_modular[!, col]) <: Number
                diff = abs.(df_original[!, col] .- df_modular[!, col])
                max_diff = maximum(diff)
                
                if max_diff > diferenca_maxima
                    diferenca_maxima = max_diff
                end
                
                diferencas_significativas = sum(diff .> 1e-10)
                total_diferencas += diferencas_significativas
                
                if diferencas_significativas > 0
                    println("⚠️  Coluna $col: $diferencas_significativas diferenças significativas (max: $max_diff)")
                end
            end
        end
        
        if total_diferencas == 0
            println("✅ SUCESSO: Os resultados são IDÊNTICOS!")
            println("   Diferença máxima encontrada: $diferenca_maxima")
        else
            println("⚠️  ATENÇÃO: Encontradas $total_diferencas diferenças numéricas")
            println("   Diferença máxima: $diferenca_maxima")
            if diferenca_maxima < 1e-8
                println("   (As diferenças são muito pequenas, provavelmente devido à precisão numérica)")
            end
        end
        
    else
        println("❌ Dimensões diferentes!")
    end
    
catch e
    println("❌ Erro ao comparar arquivos: $e")
end

println("\n=== VERIFICAÇÃO CONCLUÍDA ===")
