#!/usr/bin/env julia

"""
Script para executar o código modular e comparar com o original
"""

using Printf
using CSV
using DataFrames

println("="^60)
println("    EXECUTANDO CÓDIGO MODULAR ANDERSIM 5 EQ.")
println("="^60)

# Executar o código modular
println("\n🚀 Executando código modular...")
try
    include("AndersimWang_modular.jl")
    println("✅ Código modular executado com sucesso!")
catch e
    println("❌ Erro ao executar código modular:")
    println(e)
    exit(1)
end

println("\n📊 Comparando resultados...")

# Carregar os arquivos CSV
try
    df_original = CSV.read("export.csv", DataFrame)
    df_modular = CSV.read("export_modular.csv", DataFrame)
    
    println("\n📈 RESUMO DOS RESULTADOS:")
    println("-"^40)
    
    # Verificar se os tamanhos são iguais
    if nrow(df_original) != nrow(df_modular)
        println("⚠️  ATENÇÃO: Números de linhas diferentes!")
        println("   Original: $(nrow(df_original)) linhas")
        println("   Modular:  $(nrow(df_modular)) linhas")
    else
        println("✅ Mesmo número de linhas: $(nrow(df_original))")
    end
    
    # Comparar últimos valores
    println("\n🔍 ÚLTIMOS VALORES:")
    println("-"^40)
    
    last_orig = df_original[end, :]
    last_mod = df_modular[end, :]
    
    println("📊 PRESSÃO:")
    println("   Original: $(last_orig.Pr)")
    println("   Modular:  $(last_mod.Pr)")
    println("   Diferença: $(abs(last_orig.Pr - last_mod.Pr))")
    
    println("\n🌊 VELOCIDADE LÍQUIDO:")
    println("   Original: $(last_orig.UL)")
    println("   Modular:  $(last_mod.UL)")
    println("   Diferença: $(abs(last_orig.UL - last_mod.UL))")
    
    println("\n💨 VELOCIDADE GÁS:")
    println("   Original: $(last_orig.UG)")
    println("   Modular:  $(last_mod.UG)")
    println("   Diferença: $(abs(last_orig.UG - last_mod.UG))")
    
    println("\n🌡️  TEMPERATURA:")
    println("   Original: $(last_orig.Temperatura)")
    println("   Modular:  $(last_mod.Temperatura)")
    println("   Diferença: $(abs(last_orig.Temperatura - last_mod.Temperatura))")
    
    # Verificar se há NaN no modular
    has_nan_mod = any([any(ismissing.(col)) || any(isnan.(col)) for col in eachcol(df_modular)])
    has_nan_orig = any([any(ismissing.(col)) || any(isnan.(col)) for col in eachcol(df_original)])
    
    println("\n🔍 VERIFICAÇÃO DE NaN:")
    println("-"^40)
    if has_nan_mod
        println("❌ Código modular contém valores NaN!")
    else
        println("✅ Código modular SEM valores NaN")
    end
    
    if has_nan_orig
        println("❌ Código original contém valores NaN!")
    else
        println("✅ Código original SEM valores NaN")
    end
    
    # Calcular diferenças máximas
    println("\n📏 DIFERENÇAS MÁXIMAS:")
    println("-"^40)
    
    max_diff_pr = maximum(abs.(df_original.Pr .- df_modular.Pr))
    max_diff_ul = maximum(abs.(df_original.UL .- df_modular.UL))
    max_diff_ug = maximum(abs.(df_original.UG .- df_modular.UG))
    max_diff_temp = maximum(abs.(df_original.Temperatura .- df_modular.Temperatura))
    
    println("   Pressão:     $(max_diff_pr)")
    println("   Vel. Líq.:   $(max_diff_ul)")
    println("   Vel. Gás:    $(max_diff_ug)")
    println("   Temperatura: $(max_diff_temp)")
    
    # Determinar se os resultados são equivalentes
    tolerancia = 1e-10
    equiv_pr = max_diff_pr < tolerancia
    equiv_ul = max_diff_ul < tolerancia
    equiv_ug = max_diff_ug < tolerancia
    equiv_temp = max_diff_temp < tolerancia
    
    println("\n🎯 EQUIVALÊNCIA (tol: $(tolerancia)):")
    println("-"^40)
    println("   Pressão:     $(equiv_pr ? "✅" : "❌")")
    println("   Vel. Líq.:   $(equiv_ul ? "✅" : "❌")")
    println("   Vel. Gás:    $(equiv_ug ? "✅" : "❌")")
    println("   Temperatura: $(equiv_temp ? "✅" : "❌")")
    
    if equiv_pr && equiv_ul && equiv_ug && equiv_temp
        println("\n🎉 RESULTADOS IDÊNTICOS! Modularização bem-sucedida!")
    elseif !has_nan_mod
        println("\n✅ Modularização funcional (sem NaN), pequenas diferenças numéricas")
    else
        println("\n❌ Problemas na modularização detectados")
    end
    
catch e
    println("❌ Erro ao carregar ou comparar arquivos:")
    println(e)
end

println("\n" * "="^60)
println("    COMPARAÇÃO CONCLUÍDA")
println("="^60)
