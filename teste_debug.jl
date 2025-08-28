# Teste simples para comparar valores específicos
include("AndersimWang_modular.jl")

# Executa o código original
println("Executando código modular...")
main()

# Lê os resultados
using CSV, DataFrames
df_modular = CSV.read("export_modular.csv", DataFrame)
df_original = CSV.read("export.csv", DataFrame)

println("\nComparando primeiros valores:")
println("Original - Pressão[1]: ", df_original.Pr[1])
println("Modular  - Pressão[1]: ", df_modular.Pr[1])
println("Original - uL[1]: ", df_original.UL[1])
println("Modular  - uL[1]: ", df_modular.UL[1])
println("Original - αL[1]: ", df_original.AlfaL[1])
println("Modular  - αL[1]: ", df_modular.AlfaL[1])

println("\nComparando últimos valores:")
println("Original - Pressão[end]: ", df_original.Pr[end])
println("Modular  - Pressão[end]: ", df_modular.Pr[end])
println("Original - uL[end]: ", df_original.UL[end])
println("Modular  - uL[end]: ", df_modular.UL[end])
println("Original - αL[end]: ", df_original.AlfaL[end])
println("Modular  - αL[end]: ", df_modular.AlfaL[end])
