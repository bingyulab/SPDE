include("experiment.jl")

# ============================================================================
# EXPERIMENT 1: Noise Strength vs Mean ATP:ADP Ratio
# ============================================================================

function scan_noise_strength_experiment()
    println("\n" * "="^80)
    println("EXPERIMENT: Noise Strength vs Mean ATP:ADP Ratio")
    println("="^80)
    
    # Define noise strength ranges for each type
    noise_configs = Dict(
        :additive => 0.01:0.01:0.2,
        :multiplicative => 0.1:0.1:1.0,
        :state_dependent => 0.1:0.1:1.0,        
        :jump => 0.001:0.002:0.02
    )
    
    results = Dict()
    
    for (noise_type, strength_range) in noise_configs
        println("\nScanning $noise_type noise...")
        
        df_res = DataFrame(
            noise_strength = Float64[],
            mean_ratio = Float64[],
            std_ratio = Float64[]
        )
        
        for σ in strength_range
            println("  Noise strength σ = $σ : ")
            
            # Run simulation with modified noise
            if noise_type == :jump
                # For jump noise, interpret σ as jump size, keep rate constant
                sol, _ = simulate_model(noise_type; 
                                        λ_jump = σ,
                                        σ_jump = σ,
                                        tspan=(0.0, 4000.0))
            elseif noise_type == :additive
                sol, _ = simulate_model(noise_type; σ_additive = σ, tspan=(0.0, 4000.0))
            elseif noise_type == :multiplicative
                sol, _ = simulate_model(noise_type; σ_multiplicative = σ, tspan=(0.0, 4000.0))
            elseif noise_type == :state_dependent
                sol, _ = simulate_model(noise_type; σ_calcium = σ, tspan=(0.0, 4000.0))
            else
                continue  # skip unknown types
            end
            
            df_sol = DataFrame(sol)
            ratio = calculate_atp_adp_ratio(df_sol, 1000)
            
            if !isnothing(ratio) && length(ratio) > 10
                mean_r = mean(ratio)
                std_r = std(ratio)
                
                push!(df_res, (σ, mean_r, std_r))
                println("✓ mean=$(round(mean_r, digits=3)), std=$(round(std_r, digits=4))")
            else
                push!(df_res, (σ, NaN, NaN))
                println("⚠️  Invalid ratio data")
            end
        end
        
        results[noise_type] = df_res
    end
    
    # Plot results
    println("\nGenerating plots...")
    mkpath("imgs/bio/noise_analysis/")
    
    for (nt, df) in results
        valid_idx = .!isnan.(df.mean_ratio)
        if sum(valid_idx) > 0
            p = plot(xlabel="Noise Strength (σ)", ylabel="Mean ATP:ADP Ratio",
                     title="$nt: Noise Strength vs ATP:ADP Ratio",
                     legend=false, size=(700, 500))
            
            plot!(p, df.noise_strength[valid_idx], df.mean_ratio[valid_idx],
                  marker=:circle, linewidth=2, markersize=6,
                  ribbon=df.std_ratio[valid_idx], fillalpha=0.3,
                  color=get(color_map, nt, :blue))
            
            savefig(p, "imgs/bio/noise_analysis/noise_strength_vs_ratio_$(nt).png")
            println("Saved: imgs/bio/noise_analysis/noise_strength_vs_ratio_$(nt).png")
        end
    end
    
    # Also create combined plot for comparison
    p_combined = plot(xlabel="Noise Strength (σ)", ylabel="Mean ATP:ADP Ratio",
                      title="All Noise Types: Noise Strength vs ATP:ADP Ratio",
                      legend=:topright, size=(900, 600))
    
    for (nt, df) in results
        valid_idx = .!isnan.(df.mean_ratio)
        if sum(valid_idx) > 0
            plot!(p_combined, df.noise_strength[valid_idx], df.mean_ratio[valid_idx],
                  label=string(nt), marker=:circle, linewidth=2,
                  ribbon=df.std_ratio[valid_idx], fillalpha=0.2,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_combined, "imgs/bio/noise_analysis/noise_strength_vs_ratio_combined.png")
    println("Saved: imgs/bio/noise_analysis/noise_strength_vs_ratio_combined.png")
end

# ============================================================================
# EXPERIMENT 2: IP3 Scan with ALL analyses (SEPARATE PLOTS + Important Variables)
# ============================================================================
function scan_ip3_experiment()
    println("\n" * "="^80)
    println("EXPERIMENT 2-4: IP3 Scan & Kramers Escape Analysis")
    println("="^80)

    ip3_values = 0.1:0.1:2.0
    results_all = Dict()

    for noise_type in noise_list
        println("\nProcessing $noise_type...")
        
        # DataFrame with all metrics
        df_res = DataFrame(
            ip3 = Float64[],
            atp_adp_ratio_mean = Float64[],
            atp_adp_ratio_std = Float64[],
            dominant_freq = Float64[],
            freq_variance = Float64[],
            isi_entropy = Float64[],
            isi_mean = Float64[],
            escape_rate = Float64[],
            n_peaks = Int[]
        )
        
        # Add columns for important variables
        for var in important_variables
            df_res[!, var] = Float64[]
        end
        
        for val in ip3_values
            print("  IP3 = $val : ")
            
            try
                sol, _ = simulate_model(noise_type; tspan=(0.0, 4000.0), ip3_val=val)
                df_sol = DataFrame(sol)
                
                # Initialize row data
                row_data = Dict{String, Any}(
                    "ip3" => val,
                    "atp_adp_ratio_mean" => NaN,
                    "atp_adp_ratio_std" => NaN,
                    "dominant_freq" => NaN,
                    "freq_variance" => NaN,
                    "isi_entropy" => NaN,
                    "isi_mean" => NaN,
                    "escape_rate" => NaN,
                    "n_peaks" => 0
                )
                
                # Extract important variables (last 1000 points mean)
                for var in important_variables
                    var_data = get_last_n_points(df_sol, var, 1000)
                    if !isnothing(var_data)
                        valid_data = filter(isfinite, var_data)
                        row_data[var] = length(valid_data) > 0 ? mean(valid_data) : NaN
                    else
                        row_data[var] = NaN
                    end
                end
                
                # Calculate ATP:ADP ratio (last 1000 points)
                ratio = calculate_atp_adp_ratio(df_sol, 1000)
                
                if !isnothing(ratio) && length(ratio) > 50
                    row_data["atp_adp_ratio_mean"] = mean(ratio)
                    row_data["atp_adp_ratio_std"] = std(ratio)
                    
                    # Task 3: Frequency analysis
                    dom_freq, freq_var = estimate_frequency(ratio, 1.0)
                    row_data["dominant_freq"] = dom_freq
                    row_data["freq_variance"] = freq_var
                    
                    # Task 4: Kramers escape analysis
                    peaks = detect_peaks(ratio; threshold_percentile=75.0)
                    row_data["n_peaks"] = length(peaks)
                    
                    isi = calculate_isi(peaks, 1.0)
                    if length(isi) > 3
                        row_data["isi_mean"] = mean(isi)
                        row_data["isi_entropy"] = calculate_entropy(isi)
                        row_data["escape_rate"] = calculate_escape_rate(isi)
                    end
                    
                    println("✓ Ratio=$(round(row_data["atp_adp_ratio_mean"], digits=2)), " *
                            "Peaks=$(length(peaks)), Rate=$(round(row_data["escape_rate"], digits=4))")
                else
                    println("⚠️  Insufficient ratio data")
                end
                
                push!(df_res, row_data)
                
            catch e
                println("❌ Error: $e")
                row_data = Dict{String, Any}(
                    "ip3" => val,
                    "atp_adp_ratio_mean" => NaN,
                    "atp_adp_ratio_std" => NaN,
                    "dominant_freq" => NaN,
                    "freq_variance" => NaN,
                    "isi_entropy" => NaN,
                    "isi_mean" => NaN,
                    "escape_rate" => NaN,
                    "n_peaks" => 0
                )
                for var in important_variables
                    row_data[var] = NaN
                end
                push!(df_res, row_data)
            end
        end
        
        results_all[noise_type] = df_res
    end

    # =========================================================================
    # PLOTTING - CONSOLIDATED SUBPLOTS
    # =========================================================================
    mkpath("imgs/bio/scan/")
    
    # =========================================================================
    # IMPORTANT VARIABLES vs IP3 (RESTORED)
    # =========================================================================
    println("\n--- Plotting Important Variables vs IP3 ---")
    # Define legend position for each variable
    legend_positions = Dict(
        "adpc" => :bottomright,
        "cac" => :bottomright,
        "atpc" => :topright,
        "caer" => :topright,
        "psi" => :topright,
        "pyrm" => :topright
    )
    for var in important_variables
        leg_pos = get(legend_positions, var, :topright)
        p = plot(xlabel="IP3 (μM)", ylabel="$var (steady state)", 
                 title="Steady State $var vs IP3", 
                 legend=leg_pos, size=(900, 600),
                 left_margin=10Plots.mm, bottom_margin=8Plots.mm)
        
        for (nt, df) in results_all
            if nrow(df) > 0 && var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df[!, "ip3"][valid_idx], df[!, var][valid_idx], 
                          label=string(nt), 
                          marker=:circle, 
                          markersize=5,
                          linewidth=2,
                          color=get(color_map, nt, :auto))
                end
            end
        end
        
        savefig(p, "imgs/bio/scan/compare_$(var)_ip3.png")
        println("Saved: imgs/bio/scan/compare_$(var)_ip3.png")
    end
    
    # Combined overview of important variables
    println("\nCreating combined overview plot for important variables...")
    plots_combined = []
    for (i, var) in enumerate(important_variables)
        show_legend = (i == 1)

        p = plot(xlabel="IP3 (μM)", ylabel=var, title=var, 
                legend = show_legend ? :topright : false,
                 left_margin=8Plots.mm, bottom_margin=6Plots.mm)
        for (nt, df) in results_all
            if nrow(df) > 0 && var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df[!, "ip3"][valid_idx], df[!, var][valid_idx], 
                          label=string(nt),
                          linewidth=2, marker=:circle, markersize=3,
                          color=get(color_map, nt, :auto))
                end
            end
        end
        push!(plots_combined, p)
    end
    
    n_vars = length(plots_combined)
    n_cols = 3
    n_rows = ceil(Int, n_vars / n_cols)
    p_all = plot(plots_combined..., layout=(n_rows, n_cols), size=(500*n_cols, 350*n_rows))
    savefig(p_all, "imgs/bio/scan/all_variables_overview.png")
    println("Saved: imgs/bio/scan/all_variables_overview.png")

    # ----- TASK 2: IP3 vs Mean ATP:ADP Ratio -----
    println("\n--- Task 2: IP3 vs ATP:ADP Ratio ---")
    
    # Separate plots per noise type
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            p = plot(xlabel="IP3 (μM)", ylabel="Mean ATP:ADP Ratio",
                    title="$nt: IP3 vs ATP:ADP Ratio",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
            plot!(p, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                marker=:circle, linewidth=2, markersize=5,
                ribbon=df.atp_adp_ratio_std[valid], fillalpha=0.3,
                color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task2_ip3_vs_ratio_$(nt).png")
        end
    end

    # Combined plot
    p2_combined = plot(xlabel="IP3 (μM)", ylabel="Mean ATP:ADP Ratio",
                    title="IP3 vs ATP:ADP Ratio (All Noise Types)",
                    legend=:topright, size=(900, 600),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            plot!(p2_combined, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                label=string(nt), marker=:circle, linewidth=2,
                color=get(color_map, nt, :auto))
        end
    end
    savefig(p2_combined, "imgs/bio/scan/task2_ip3_vs_ratio_combined.png")
    println("Saved: imgs/bio/scan/task2_ip3_vs_ratio_combined.png")

    # =========================================================================
    # TASK 4: IP3 vs Escape Rate and IP3 vs Entropy (RESTORED)
    # =========================================================================
    println("\n--- Task 4: IP3 vs Escape Rate ---")
    
    # Escape Rate vs IP3 (separate per noise type)
    for (nt, df) in results_all
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            p = plot(xlabel="IP3 (μM)", ylabel="Escape Rate (1/s)",
                    title="$nt: Escape Rate vs IP3",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
            plot!(p, df.ip3[valid], df.escape_rate[valid],
                marker=:circle, linewidth=2, markersize=5,
                color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task4_escape_rate_vs_ip3_$(nt).png")
            println("Saved: imgs/bio/scan/task4_escape_rate_vs_ip3_$(nt).png")
        end
    end
    
    # Combined Escape Rate vs IP3
    p_rate_ip3 = plot(xlabel="IP3 (μM)", ylabel="Escape Rate (1/s)",
                      title="Escape Rate vs IP3 (All Noise Types)",
                      legend=:topright, size=(900, 600),
                      left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    for (nt, df) in results_all
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            plot!(p_rate_ip3, df.ip3[valid], df.escape_rate[valid],
                  label=string(nt), marker=:circle, linewidth=2,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_rate_ip3, "imgs/bio/scan/task4_escape_rate_vs_ip3_combined.png")
    println("Saved: imgs/bio/scan/task4_escape_rate_vs_ip3_combined.png")
    
    println("\n--- Task 4: IP3 vs ISI Entropy ---")
    
    # ISI Entropy vs IP3 (separate per noise type)
    for (nt, df) in results_all
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            p = plot(xlabel="IP3 (μM)", ylabel="ISI Entropy (bits)",
                    title="$nt: ISI Entropy vs IP3",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
            plot!(p, df.ip3[valid], df.isi_entropy[valid],
                marker=:circle, linewidth=2, markersize=5,
                color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task4_entropy_vs_ip3_$(nt).png")
            println("Saved: imgs/bio/scan/task4_entropy_vs_ip3_$(nt).png")
        end
    end
    
    # Combined ISI Entropy vs IP3
    p_entropy_ip3 = plot(xlabel="IP3 (μM)", ylabel="ISI Entropy (bits)",
                         title="ISI Entropy vs IP3 (All Noise Types)",
                         legend=:bottomright, size=(900, 600),
                         left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    for (nt, df) in results_all
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            plot!(p_entropy_ip3, df.ip3[valid], df.isi_entropy[valid],
                  label=string(nt), marker=:circle, linewidth=2,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_entropy_ip3, "imgs/bio/scan/task4_entropy_vs_ip3_combined.png")
    println("Saved: imgs/bio/scan/task4_entropy_vs_ip3_combined.png")

    # ----- Filter out :none for stochastic analysis -----
    stochastic_noises = filter(x -> x != :none, collect(keys(results_all)))

    # ----- 3a. Dominant Frequency (1 × N layout) -----
    subplots_freq = []

    for nt in stochastic_noises
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        
        p = plot(xlabel="IP3 (μM)", ylabel="Freq (Hz)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 2
            x = df.ip3[valid]
            y = df.dominant_freq[valid]
            
            plot!(p, x, y, marker=:circle, linewidth=2, markersize=4,
                color=get(color_map, nt, :blue))
            
            # Detect jumps
            if length(y) > 3
                d_freq = abs.(diff(y))
                threshold = 2 * median(d_freq)
                jumps = findall(d_freq .> threshold)
                for j in jumps
                    vline!(p, [x[j]], linestyle=:dash, color=:red, alpha=0.5)
                end
            end
            
            # Pattern annotation
            if length(y) > 5
                corr = cor(x, y)
                max_idx = argmax(y)
                if max_idx > 2 && max_idx < length(y) - 1
                    pattern = "peaked"
                elseif corr > 0.5
                    pattern = "↗"
                elseif corr < -0.5
                    pattern = "↘"
                else
                    pattern = "~"
                end
                annotate!(p, :topright, text(pattern, 10))
            end
        end
        push!(subplots_freq, p)
    end

    n_plots = length(subplots_freq)
    if n_plots > 0
        p_freq_all = plot(subplots_freq..., layout=(1, n_plots),
                        size=(350*n_plots, 400),
                        plot_title="Dominant Frequency vs IP3")
        savefig(p_freq_all, "imgs/bio/scan/task3_frequency_analysis.png")
        println("Saved: imgs/bio/scan/task3_frequency_analysis.png")
    end

    # ----- 3b. Frequency Variance (1 × N layout) -----
    subplots_fvar = []

    for nt in stochastic_noises
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        valid = .!isnan.(df.freq_variance)
        
        p = plot(xlabel="IP3 (μM)", ylabel="Freq Var",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 2
            x = df.ip3[valid]
            y = df.freq_variance[valid]
            
            plot!(p, x, y, marker=:circle, linewidth=2, markersize=4,
                color=get(color_map, nt, :blue))
            
            mean_var = mean(filter(isfinite, y))
            label = mean_var < 0.01 ? "regular" : "irregular"
            annotate!(p, :topright, text(label, 9))
        end
        push!(subplots_fvar, p)
    end

    n_plots = length(subplots_fvar)
    if n_plots > 0
        p_fvar_all = plot(subplots_fvar..., layout=(1, n_plots),
                        size=(350*n_plots, 400),
                        plot_title="Frequency Variance vs IP3 (Low=Regular, High=Irregular)")
        savefig(p_fvar_all, "imgs/bio/scan/task3_freq_variance_analysis.png")
        println("Saved: imgs/bio/scan/task3_freq_variance_analysis.png")
    end
    
    # ----- 4a. Escape Rate vs Energy (with single colorbar) -----
    subplots_rate = []
    all_ip3_vals = Float64[]

    for nt in stochastic_noises
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            append!(all_ip3_vals, df.ip3[valid])
        end
    end

    ip3_min = length(all_ip3_vals) > 0 ? minimum(all_ip3_vals) : 0.1
    ip3_max = length(all_ip3_vals) > 0 ? maximum(all_ip3_vals) : 2.0

    for (idx, nt) in enumerate(stochastic_noises)
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& (df.escape_rate .> 0)
                
        p = plot(xlabel="ATP:ADP (Energy)", ylabel="Rate (1/s)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=10Plots.mm, bottom_margin=10Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 2
            x = df.atp_adp_ratio_mean[valid]
            y = df.escape_rate[valid]
            
            scatter!(p, x, y, markersize=6,
                    zcolor=df.ip3[valid], color=:viridis,
                    clims=(ip3_min, ip3_max),
                    colorbar=false, alpha=0.8, markerstrokewidth=0)
            
            # Trend line
            valid_xy = isfinite.(x) .& isfinite.(y)
            if sum(valid_xy) > 2
                X_mat = hcat(ones(sum(valid_xy)), x[valid_xy])
                coeffs = X_mat \ y[valid_xy]
                x_line = range(minimum(x[valid_xy]), maximum(x[valid_xy]), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x[valid_xy], y[valid_xy])
                slope_sign = coeffs[2] < 0 ? "↘" : "↗"
                annotate!(p, :topright, text("r=$(round(corr, digits=2))$slope_sign", 9))
            end
        end
        push!(subplots_rate, p)
    end

    n_plots = length(subplots_rate)
    if n_plots > 0
        # Create a dummy plot just for the colorbar
        p_cbar = scatter([NaN], [NaN], zcolor=[ip3_min], clims=(ip3_min, ip3_max),
                        color=:viridis, colorbar=true, colorbar_title="IP3 (μM)",
                        framestyle=:none, label="", markersize=0,
                        left_margin=0Plots.mm, right_margin=0Plots.mm)
        
        widths = vcat(ones(n_plots), [0.4])
        widths = widths ./ sum(widths)
        all_plots = vcat(subplots_rate, [p_cbar])
        
        p_rate_combined = plot(all_plots..., 
                                layout=grid(1, n_plots+1, widths=widths),
                                size=(350*n_plots + 80, 400),
                                plot_title="Escape Rate vs Energy (E↑→Rate↓ = Kramers)")
        savefig(p_rate_combined, "imgs/bio/scan/task4_escape_rate_vs_energy.png")
        println("Saved: imgs/bio/scan/task4_escape_rate_vs_energy.png")
    end

    # ----- 4b. Entropy vs Energy (with single colorbar) -----
    subplots_entropy = []

    for (idx, nt) in enumerate(stochastic_noises)
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        
        valid = .!isnan.(df.isi_entropy) .& .!isnan.(df.atp_adp_ratio_mean)
                
        p = plot(xlabel="ATP:ADP (Energy)", ylabel="Entropy (bits)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=10Plots.mm, bottom_margin=10Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 2
            x = df.atp_adp_ratio_mean[valid]
            y = df.isi_entropy[valid]
            
            scatter!(p, x, y, markersize=6,
                    zcolor=df.ip3[valid], color=:viridis,
                    clims=(ip3_min, ip3_max),
                    colorbar=false, alpha=0.8, markerstrokewidth=0)
            
            valid_xy = isfinite.(x) .& isfinite.(y)
            if sum(valid_xy) > 2
                X_mat = hcat(ones(sum(valid_xy)), x[valid_xy])
                coeffs = X_mat \ y[valid_xy]
                x_line = range(minimum(x[valid_xy]), maximum(x[valid_xy]), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x[valid_xy], y[valid_xy])
                slope_sign = coeffs[2] < 0 ? "↘" : "↗"
                annotate!(p, :topright, text("r=$(round(corr, digits=2))$slope_sign", 9))
            end
        end
        push!(subplots_entropy, p)
    end

    n_plots = length(subplots_entropy)
    if n_plots > 0
        p_cbar = scatter([NaN], [NaN], zcolor=[ip3_min], clims=(ip3_min, ip3_max),
                        color=:viridis, colorbar=true, colorbar_title="IP3 (μM)",
                        framestyle=:none, label="", markersize=0,
                        left_margin=0Plots.mm, right_margin=0Plots.mm)
        widths = vcat(ones(n_plots), [0.4])
        widths = widths ./ sum(widths)
        all_plots = vcat(subplots_entropy, [p_cbar])
        
        p_entropy_combined = plot(all_plots..., 
                                layout=grid(1, n_plots+1, widths=widths),
                                size=(350*n_plots + 80, 400),
                                plot_title="Entropy vs Energy (IP3↑→Entropy↑ = More random)")
        savefig(p_entropy_combined, "imgs/bio/scan/task4_entropy_vs_energy.png")
        println("Saved: imgs/bio/scan/task4_entropy_vs_energy.png")
    end

    # ----- 4c. Kramers-Arrhenius (with single colorbar) -----
    subplots_arrhenius = []

    for (idx, nt) in enumerate(stochastic_noises)
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
                
        p = plot(xlabel="1/Energy", ylabel="log(Rate)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=10Plots.mm, bottom_margin=10Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 2
            inv_energy = 1.0 ./ df.atp_adp_ratio_mean[valid]
            log_rate = log.(df.escape_rate[valid])
            
            valid_xy = isfinite.(inv_energy) .& isfinite.(log_rate)
            
            if sum(valid_xy) > 2
                x = inv_energy[valid_xy]
                y = log_rate[valid_xy]
                
                scatter!(p, x, y, markersize=6,
                        zcolor=df.ip3[valid][valid_xy], color=:viridis,
                        clims=(ip3_min, ip3_max),
                        colorbar=false, alpha=0.8, markerstrokewidth=0)
                
                X_mat = hcat(ones(length(x)), x)
                coeffs = X_mat \ y
                x_line = range(minimum(x), maximum(x), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x, y)
                interpretation = coeffs[2] > 0 ? "E↑→R↓" : "E↑→R↑"
                annotate!(p, :topright, text("r=$(round(corr, digits=2))\n$interpretation", 8))
            end
        end
        push!(subplots_arrhenius, p)
    end

    n_plots = length(subplots_arrhenius)
    if n_plots > 0        
        p_cbar = scatter([NaN], [NaN], zcolor=[ip3_min], clims=(ip3_min, ip3_max),
                        color=:viridis, colorbar=true, colorbar_title="IP3 (μM)",
                        framestyle=:none, label="", markersize=0,
                        left_margin=0Plots.mm, right_margin=0Plots.mm)
        widths = vcat(ones(n_plots), [0.4])
        widths = widths ./ sum(widths)
        all_plots = vcat(subplots_arrhenius, [p_cbar])
        
        p_arr_combined = plot(all_plots..., 
                            layout=grid(1, n_plots+1, widths=widths),
                            size=(350*n_plots + 80, 400),
                            plot_title="Kramers-Arrhenius (Positive slope = E↑→Rate↓)")
        savefig(p_arr_combined, "imgs/bio/scan/task4_kramers_arrhenius.png")
        println("Saved: imgs/bio/scan/task4_kramers_arrhenius.png")
    end

    println("\n--- TASK 5: Investigating Variance at High IP3 ---")

    # ----- 5a. Coefficient of Variation -----
    subplots_cv = []

    for nt in stochastic_noises
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        
        valid = .!isnan.(df.atp_adp_ratio_mean) .& .!isnan.(df.atp_adp_ratio_std) .& 
                (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="IP3 (μM)", ylabel="CV (σ/μ)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 0
            cv = df.atp_adp_ratio_std[valid] ./ df.atp_adp_ratio_mean[valid]
            
            plot!(p, df.ip3[valid], cv, marker=:circle, linewidth=2, markersize=4,
                color=get(color_map, nt, :blue))
            
            vline!(p, [1.5], linestyle=:dash, color=:red, alpha=0.7)
            
            high_ip3_idx = df.ip3[valid] .> 1.5
            if sum(high_ip3_idx) > 0 && sum(.!high_ip3_idx) > 0
                cv_low = mean(cv[.!high_ip3_idx])
                cv_high = mean(cv[high_ip3_idx])
                if cv_high > cv_low * 1.2
                    ratio = round(cv_high / cv_low, digits=1)
                    annotate!(p, :topright, text("$(ratio)x↑", 9, :red))
                end
            end
        end
        push!(subplots_cv, p)
    end

    n_plots = length(subplots_cv)
    if n_plots > 0
        p_cv = plot(subplots_cv..., layout=(1, n_plots),
                    size=(350*n_plots, 400),
                    plot_title="Coefficient of Variation vs IP3 (Red=IP3=1.5)")
        savefig(p_cv, "imgs/bio/scan/task5_cv_analysis.png")
        println("Saved: imgs/bio/scan/task5_cv_analysis.png")
    end

    # ----- 5b. Standard Deviation -----
    subplots_std = []

    for nt in stochastic_noises
        if !haskey(results_all, nt)
            continue
        end
        df = results_all[nt]
        
        valid = .!isnan.(df.atp_adp_ratio_std)
        
        p = plot(xlabel="IP3 (μM)", ylabel="Std Dev",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                top_margin=5Plots.mm, right_margin=5Plots.mm)
        
        if sum(valid) > 0
            plot!(p, df.ip3[valid], df.atp_adp_ratio_std[valid],
                marker=:circle, linewidth=2, markersize=4,
                color=get(color_map, nt, :blue))
            
            vline!(p, [1.5], linestyle=:dash, color=:red, alpha=0.7)
        end
        push!(subplots_std, p)
    end

    n_plots = length(subplots_std)
    if n_plots > 0
        p_std = plot(subplots_std..., layout=(1, n_plots),
                    size=(350*n_plots, 400),
                    plot_title="Std Dev of ATP:ADP vs IP3")
        savefig(p_std, "imgs/bio/scan/task5_std_analysis.png")
        println("Saved: imgs/bio/scan/task5_std_analysis.png")
    end

    # ----- SIMPLIFIED SUMMARY GRID (4 panels only) -----
    println("\n--- Creating Simplified Summary Grid (4 panels) ---")

    p1 = plot(xlabel="IP3 (μM)", ylabel="Freq (Hz)", title="Oscillation Frequency", 
              legend=:topleft, size=(400, 350),
              left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    p2 = plot(xlabel="IP3 (μM)", ylabel="Freq Var", title="Frequency Variability", 
              legend=false, size=(400, 350),
              left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    p3 = plot(xlabel="IP3 (μM)", ylabel="Rate (1/s)", title="Escape Rate", 
              legend=false, size=(400, 350),
              left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    p4 = plot(xlabel="IP3 (μM)", ylabel="Entropy (bits)", title="ISI Entropy", 
              legend=false, size=(400, 350),
              left_margin=10Plots.mm, bottom_margin=8Plots.mm)

    for (nt, df) in results_all
        c = get(color_map, nt, :auto)
        
        # Oscillation Frequency
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            plot!(p1, df.ip3[valid], df.dominant_freq[valid], 
                  label=string(nt), color=c, linewidth=2, marker=:circle, markersize=3)
        end
        
        # Frequency Variability
        valid = .!isnan.(df.freq_variance)
        if sum(valid) > 0
            plot!(p2, df.ip3[valid], df.freq_variance[valid], color=c, linewidth=2, 
                  marker=:circle, markersize=3)
        end
        
        # Escape Rate
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            plot!(p3, df.ip3[valid], df.escape_rate[valid], color=c, linewidth=2,
                  marker=:circle, markersize=3)
        end
        
        # ISI Entropy
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            plot!(p4, df.ip3[valid], df.isi_entropy[valid], color=c, linewidth=2,
                  marker=:circle, markersize=3)
        end
    end

    p_summary = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
    savefig(p_summary, "imgs/bio/scan/kramers_summary.png")
    println("Saved: imgs/bio/scan/kramers_summary.png")

    println("\n" * "="^80)
    println("IP3 SCAN COMPLETE")
    println("="^80)

    return results_all
end

# ============================================================================
# EXPERIMENT 3 : ISI Distribution & Kramers Escape
# ============================================================================

function analyze_isi_distribution(; ip3_val=0.7)
    println("\n" * "="^80)
    println("DETAILED ISI DISTRIBUTION ANALYSIS (IP3=$ip3_val)")
    println("="^80)
    
    mkpath("imgs/bio/isi_analysis/")
    
    isi_results = Dict()
    
    for noise_type in noise_list
        println("\nAnalyzing $noise_type...")
        
        try
            sol, _ = simulate_model(noise_type; tspan=(0.0, 6000.0), ip3_val=ip3_val)
            df_sol = DataFrame(sol)
            
            # Use last 2000 points for statistics
            ratio = calculate_atp_adp_ratio(df_sol, 2000)
            
            if isnothing(ratio) || length(ratio) < 100
                println("  ⚠️  Insufficient ratio data")
                continue
            end
            
            peaks = detect_peaks(ratio; threshold_percentile=75.0)
            println("  Found $(length(peaks)) peaks (escape events)")
            
            if length(peaks) < 5
                println("  ⚠️  Too few escape events")
                continue
            end
            
            isi = calculate_isi(peaks, 1.0)
            
            if length(isi) < 5
                continue
            end
            
            # Kramers escape statistics
            isi_mean = mean(isi)
            isi_std = std(isi)
            isi_cv = isi_std / isi_mean  # CV ~ 1 for Poisson process
            isi_entropy = calculate_entropy(isi; n_bins=15)
            escape_rate = calculate_escape_rate(isi)
            energy = mean(ratio)
            
            isi_results[noise_type] = Dict(
                :isi => isi,
                :mean => isi_mean,
                :std => isi_std,
                :cv => isi_cv,
                :entropy => isi_entropy,
                :escape_rate => escape_rate,
                :energy => energy,
                :n_peaks => length(peaks)
            )
            
            println("  ISI: μ=$(round(isi_mean, digits=2))s, σ=$(round(isi_std, digits=2))s")
            println("  CV: $(round(isi_cv, digits=3)) (≈1 for Poisson)")
            println("  Entropy: $(round(isi_entropy, digits=3)) bits")
            println("  Escape Rate: $(round(escape_rate, digits=4)) Hz")
            println("  Energy (ATP:ADP): $(round(energy, digits=3))")
            
            # Individual histogram with exponential fit overlay
            p_hist = histogram(isi, bins=20, normalize=:pdf,
                              xlabel="Inter-Spike Interval (s)",
                              ylabel="Probability Density",
                              title="$noise_type: ISI Distribution\n(λ=$(round(escape_rate,digits=3)) Hz, H=$(round(isi_entropy,digits=2)) bits)",
                              legend=:topright, fillalpha=0.7, size=(700, 500),
                              color=get(color_map, noise_type, :blue), label="Data")
            
            # Overlay exponential fit (Kramers prediction)
            x_exp = range(minimum(isi), maximum(isi), length=100)
            y_exp = escape_rate .* exp.(-escape_rate .* x_exp)
            plot!(p_hist, x_exp, y_exp, linewidth=3, color=:red, 
                  linestyle=:dash, label="Exponential fit")
            
            vline!(p_hist, [isi_mean], color=:black, linewidth=2, 
                   linestyle=:dot, label="Mean")
            
            savefig(p_hist, "imgs/bio/isi_analysis/isi_histogram_$(noise_type).png")
            println("  Saved histogram")
            
        catch e
            println("  ❌ Error: $e")
        end
    end
    
    # ========== COMBINED ANALYSIS PLOTS ==========
    
    if length(isi_results) > 0
        # 1. Overlay histograms
        p_overlay = plot(xlabel="ISI (s)", ylabel="Density",
                        title="ISI Distribution Comparison",
                        legend=:topright, size=(900, 600))
        
        for (nt, res) in isi_results
            histogram!(p_overlay, res[:isi], bins=20, normalize=:pdf, alpha=0.4,
                      label="$nt (λ=$(round(res[:escape_rate], digits=3)))",
                      color=get(color_map, nt, :auto))
        end
        savefig(p_overlay, "imgs/bio/isi_analysis/isi_comparison.png")
        
        # 2. Kramers relationship: log(escape_rate) vs 1/energy (Arrhenius-like)
        p_arrhenius = plot(xlabel="1/Energy (1/ATP:ADP)", ylabel="log(Escape Rate)",
                          title="Kramers-Arrhenius Plot",
                          legend=:topright, size=(800, 600))
        
        for (nt, res) in isi_results
            inv_energy = 1.0 / res[:energy]
            log_rate = log(res[:escape_rate])
            scatter!(p_arrhenius, [inv_energy], [log_rate],
                    label=string(nt), markersize=12, markerstrokewidth=2,
                    color=get(color_map, nt, :auto))
        end
        savefig(p_arrhenius, "imgs/bio/isi_analysis/kramers_arrhenius.png")
        println("Saved: imgs/bio/isi_analysis/kramers_arrhenius.png")
        
        # 3. Entropy vs Energy (scatter)
        p_ee = plot(xlabel="Energy (ATP:ADP)", ylabel="ISI Entropy (bits)",
                   title="Entropy-Energy Relationship",
                   legend=:topright, size=(800, 600))
        
        for (nt, res) in isi_results
            scatter!(p_ee, [res[:energy]], [res[:entropy]],
                    label=string(nt), markersize=12, markerstrokewidth=2,
                    color=get(color_map, nt, :auto))
        end
        savefig(p_ee, "imgs/bio/isi_analysis/entropy_vs_energy.png")
        
        # 4. Summary bar charts
        noise_types_sorted = sort(collect(keys(isi_results)))
        
        p_bars = plot(layout=(2, 2), size=(1000, 800))
        
        # Escape rate
        rates = [isi_results[nt][:escape_rate] for nt in noise_types_sorted]
        bar!(p_bars[1], string.(noise_types_sorted), rates,
             ylabel="Escape Rate (Hz)", title="Escape Rate", legend=false,
             color=[get(color_map, nt, :gray) for nt in noise_types_sorted])
        
        # Entropy
        entropies = [isi_results[nt][:entropy] for nt in noise_types_sorted]
        bar!(p_bars[2], string.(noise_types_sorted), entropies,
             ylabel="Entropy (bits)", title="ISI Entropy", legend=false,
             color=[get(color_map, nt, :gray) for nt in noise_types_sorted])
        
        # CV (should be ~1 for Poisson)
        cvs = [isi_results[nt][:cv] for nt in noise_types_sorted]
        bar!(p_bars[3], string.(noise_types_sorted), cvs,
             ylabel="CV", title="Coefficient of Variation", legend=false,
             color=[get(color_map, nt, :gray) for nt in noise_types_sorted])
        hline!(p_bars[3], [1.0], linestyle=:dash, color=:red, label="Poisson")
        
        # Energy
        energies = [isi_results[nt][:energy] for nt in noise_types_sorted]
        bar!(p_bars[4], string.(noise_types_sorted), energies,
             ylabel="ATP:ADP", title="Mean Energy", legend=false,
             color=[get(color_map, nt, :gray) for nt in noise_types_sorted])
        
        savefig(p_bars, "imgs/bio/isi_analysis/kramers_summary_bars.png")
        println("Saved: imgs/bio/isi_analysis/kramers_summary_bars.png")
        
        # Print summary table
        println("\n" * "="^70)
        println("KRAMERS ESCAPE SUMMARY")
        println("="^70)
        @printf("%-15s %8s %8s %8s %8s %8s\n", 
                "Noise Type", "Rate", "Mean τ", "CV", "Entropy", "Energy")
        println("-"^70)
        
        for nt in noise_types_sorted
            res = isi_results[nt]
            @printf("%-15s %8.4f %8.2f %8.3f %8.3f %8.3f\n",
                   string(nt), res[:escape_rate], res[:mean], 
                   res[:cv], res[:entropy], res[:energy])
        end
        println("="^70)
        println("Note: CV ≈ 1 indicates Poisson-like escape (Kramers regime)")
    end
    
    return isi_results
end


function scan_bifurcation_experiment()
    println("\n" * "="^80)
    println("EXPERIMENT: Bifurcation Analysis - Extended IP3 Range")
    println("="^80)
    
    # Extended IP3 range to capture potential bifurcation
    ip3_values = 0.1:0.5:1000.0  # Extended to IP3 = 5.0 μM
    
    results_all = Dict()
    
    for noise_type in [:none]
        println("\nProcessing $noise_type...")
        
        df_res = DataFrame(
            ip3 = Float64[],
            atp_adp_ratio_mean = Float64[],
            atp_adp_ratio_std = Float64[],
            atp_adp_ratio_min = Float64[],
            atp_adp_ratio_max = Float64[],
            oscillation_amplitude = Float64[],
            dominant_freq = Float64[],
            n_peaks = Int[],
            is_oscillating = Bool[],
            cac_mean = Float64[],
            cac_std = Float64[],
            caer_mean = Float64[],
            caer_std = Float64[]
        )
        
        for val in ip3_values
            print("  IP3 = $val : ")
            
            try
                sol, _ = simulate_model(noise_type; tspan=(0.0, 4000.0), ip3_val=val)
                df_sol = DataFrame(sol)
                
                # Initialize row
                row_data = Dict{String, Any}(
                    "ip3" => val,
                    "atp_adp_ratio_mean" => NaN,
                    "atp_adp_ratio_std" => NaN,
                    "atp_adp_ratio_min" => NaN,
                    "atp_adp_ratio_max" => NaN,
                    "oscillation_amplitude" => NaN,
                    "dominant_freq" => NaN,
                    "n_peaks" => 0,
                    "is_oscillating" => false,
                    "cac_mean" => NaN,
                    "cac_std" => NaN,
                    "caer_mean" => NaN,
                    "caer_std" => NaN
                )
                
                # Calculate ATP:ADP ratio
                ratio = calculate_atp_adp_ratio(df_sol, 1000)
                
                if !isnothing(ratio) && length(ratio) > 50
                    valid_ratio = filter(isfinite, ratio)
                    if length(valid_ratio) > 10
                        row_data["atp_adp_ratio_mean"] = mean(valid_ratio)
                        row_data["atp_adp_ratio_std"] = std(valid_ratio)
                        row_data["atp_adp_ratio_min"] = minimum(valid_ratio)
                        row_data["atp_adp_ratio_max"] = maximum(valid_ratio)
                        row_data["oscillation_amplitude"] = maximum(valid_ratio) - minimum(valid_ratio)
                        
                        # Frequency analysis
                        dom_freq, _ = estimate_frequency(valid_ratio, 1.0)
                        row_data["dominant_freq"] = dom_freq
                        
                        # Peak detection for oscillation characterization
                        peaks = detect_peaks(valid_ratio; threshold_percentile=75.0)
                        row_data["n_peaks"] = length(peaks)
                        
                        # Determine if oscillating: need multiple peaks and significant amplitude
                        cv = row_data["atp_adp_ratio_std"] / row_data["atp_adp_ratio_mean"]
                        row_data["is_oscillating"] = (length(peaks) >= 3) && (cv > 0.1)
                    end
                end
                
                # Extract calcium concentrations
                cac_data = get_last_n_points(df_sol, "cac", 1000)
                if !isnothing(cac_data)
                    valid_cac = filter(isfinite, cac_data)
                    if length(valid_cac) > 10
                        row_data["cac_mean"] = mean(valid_cac)
                        row_data["cac_std"] = std(valid_cac)
                    end
                end
                
                caer_data = get_last_n_points(df_sol, "caer", 1000)
                if !isnothing(caer_data)
                    valid_caer = filter(isfinite, caer_data)
                    if length(valid_caer) > 10
                        row_data["caer_mean"] = mean(valid_caer)
                        row_data["caer_std"] = std(valid_caer)
                    end
                end
                
                # Status output
                osc_status = row_data["is_oscillating"] ? "oscillating" : "stable"
                println("✓ Ratio=$(round(row_data["atp_adp_ratio_mean"], digits=2)), " *
                        "Amp=$(round(row_data["oscillation_amplitude"], digits=2)), " *
                        "Peaks=$(row_data["n_peaks"]), $osc_status")
                
                push!(df_res, row_data)
                
            catch e
                println("❌ Error: $e")
                row_data = Dict{String, Any}(
                    "ip3" => val,
                    "atp_adp_ratio_mean" => NaN,
                    "atp_adp_ratio_std" => NaN,
                    "atp_adp_ratio_min" => NaN,
                    "atp_adp_ratio_max" => NaN,
                    "oscillation_amplitude" => NaN,
                    "dominant_freq" => NaN,
                    "n_peaks" => 0,
                    "is_oscillating" => false,
                    "cac_mean" => NaN,
                    "cac_std" => NaN,
                    "caer_mean" => NaN,
                    "caer_std" => NaN
                )
                push!(df_res, row_data)
            end
        end
        
        results_all[noise_type] = df_res
    end
    
    # =========================================================================
    # PLOTTING
    # =========================================================================
    mkpath("imgs/bio/bifurcation/")
    
    # ----- 1. Bifurcation Diagram: ATP:ADP ratio vs IP3 -----
    println("\n--- Creating Bifurcation Diagrams ---")
    
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            p = plot(xlabel="IP3 (μM)", ylabel="ATP:ADP Ratio",
                    title="$nt: Bifurcation Diagram",
                    legend=:topright, size=(900, 600),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
            
            # Plot mean with min-max envelope
            plot!(p, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                  label="Mean", linewidth=2, color=get(color_map, nt, :blue))
            
            # Add min-max envelope as ribbon
            if sum(.!isnan.(df.atp_adp_ratio_min[valid])) > 0
                plot!(p, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                      ribbon=(df.atp_adp_ratio_mean[valid] .- df.atp_adp_ratio_min[valid],
                              df.atp_adp_ratio_max[valid] .- df.atp_adp_ratio_mean[valid]),
                      fillalpha=0.3, label="", color=get(color_map, nt, :blue))
            end
            
            # Mark oscillating vs stable regions
            osc_idx = valid .& df.is_oscillating
            stable_idx = valid .& .!df.is_oscillating
            
            if sum(osc_idx) > 0
                scatter!(p, df.ip3[osc_idx], df.atp_adp_ratio_mean[osc_idx],
                        marker=:circle, markersize=6, color=:green, alpha=0.7,
                        label="Oscillating")
            end
            if sum(stable_idx) > 0
                scatter!(p, df.ip3[stable_idx], df.atp_adp_ratio_mean[stable_idx],
                        marker=:square, markersize=6, color=:red, alpha=0.7,
                        label="Stable")
            end
            
            savefig(p, "imgs/bio/bifurcation/bifurcation_$(nt).png")
            println("Saved: imgs/bio/bifurcation/bifurcation_$(nt).png")
        end
    end
    
    # ----- 2. Combined Bifurcation Diagram -----
    p_combined = plot(xlabel="IP3 (μM)", ylabel="ATP:ADP Ratio",
                      title="Bifurcation Diagram (All Noise Types)",
                      legend=:topright, size=(1000, 700),
                      left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            plot!(p_combined, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=3,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_combined, "imgs/bio/bifurcation/bifurcation_combined.png")
    println("Saved: imgs/bio/bifurcation/bifurcation_combined.png")
    
    # ----- 3. Oscillation Amplitude vs IP3 -----
    p_amp = plot(xlabel="IP3 (μM)", ylabel="Oscillation Amplitude",
                 title="Oscillation Amplitude vs IP3",
                 legend=:topright, size=(900, 600),
                 left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        valid = .!isnan.(df.oscillation_amplitude)
        if sum(valid) > 0
            plot!(p_amp, df.ip3[valid], df.oscillation_amplitude[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=4,
                  color=get(color_map, nt, :auto))
        end
    end
    
    # Add horizontal line at amplitude = 0.1 (threshold for oscillation)
    hline!(p_amp, [0.1], linestyle=:dash, color=:gray, label="Threshold", alpha=0.7)
    
    savefig(p_amp, "imgs/bio/bifurcation/oscillation_amplitude.png")
    println("Saved: imgs/bio/bifurcation/oscillation_amplitude.png")
    
    # ----- 4. Frequency vs IP3 -----
    p_freq = plot(xlabel="IP3 (μM)", ylabel="Dominant Frequency (Hz)",
                  title="Oscillation Frequency vs IP3",
                  legend=:topright, size=(900, 600),
                  left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            plot!(p_freq, df.ip3[valid], df.dominant_freq[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=4,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_freq, "imgs/bio/bifurcation/frequency_vs_ip3.png")
    println("Saved: imgs/bio/bifurcation/frequency_vs_ip3.png")
    
    # ----- 5. Number of Peaks (Oscillation Count) vs IP3 -----
    p_peaks = plot(xlabel="IP3 (μM)", ylabel="Number of Peaks",
                   title="Peak Count vs IP3 (0 = No Oscillation)",
                   legend=:topright, size=(900, 600),
                   left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        valid = df.n_peaks .>= 0  # All valid
        if sum(valid) > 0
            plot!(p_peaks, df.ip3[valid], df.n_peaks[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=4,
                  color=get(color_map, nt, :auto))
        end
    end
    
    # Add horizontal line at n_peaks = 3 (minimum for oscillation)
    hline!(p_peaks, [3], linestyle=:dash, color=:gray, label="Min for oscillation", alpha=0.7)
    
    savefig(p_peaks, "imgs/bio/bifurcation/peak_count_vs_ip3.png")
    println("Saved: imgs/bio/bifurcation/peak_count_vs_ip3.png")
    
    # ----- 6. Calcium Dynamics: CaC and CaER vs IP3 -----
    p_ca = plot(layout=(1, 2), size=(1200, 500),
                left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    # CaC
    for (nt, df) in results_all
        valid = .!isnan.(df.cac_mean)
        if sum(valid) > 0
            plot!(p_ca[1], df.ip3[valid], df.cac_mean[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=3,
                  color=get(color_map, nt, :auto))
        end
    end
    plot!(p_ca[1], xlabel="IP3 (μM)", ylabel="[Ca²⁺]c (μM)", 
          title="Cytosolic Calcium vs IP3", legend=:topleft)
    
    # CaER
    for (nt, df) in results_all
        valid = .!isnan.(df.caer_mean)
        if sum(valid) > 0
            plot!(p_ca[2], df.ip3[valid], df.caer_mean[valid],
                  label="", linewidth=2, marker=:circle, markersize=3,
                  color=get(color_map, nt, :auto))
        end
    end
    plot!(p_ca[2], xlabel="IP3 (μM)", ylabel="[Ca²⁺]ER (μM)", 
          title="ER Calcium vs IP3", legend=false)
    
    savefig(p_ca, "imgs/bio/bifurcation/calcium_vs_ip3.png")
    println("Saved: imgs/bio/bifurcation/calcium_vs_ip3.png")
    
    # ----- 7. Phase Diagram: Oscillating vs Stable Regions -----
    println("\n--- Creating Phase Diagram ---")
    
    p_phase = plot(xlabel="IP3 (μM)", ylabel="Noise Type",
                   title="Oscillation Phase Diagram\n(Green = Oscillating, Red = Stable)",
                   size=(1000, 400), legend=false,
                   left_margin=15Plots.mm, bottom_margin=10Plots.mm)
    
    noise_types_ordered = [:none, :additive, :multiplicative, :state_dependent, :jump]
    
    for (y_idx, nt) in enumerate(noise_types_ordered)
        if haskey(results_all, nt)
            df = results_all[nt]
            for (i, row) in enumerate(eachrow(df))
                color = row.is_oscillating ? :green : :red
                scatter!(p_phase, [row.ip3], [y_idx], 
                        marker=:square, markersize=12, color=color, alpha=0.8)
            end
        end
    end
    
    yticks!(p_phase, 1:length(noise_types_ordered), string.(noise_types_ordered))
    
    savefig(p_phase, "imgs/bio/bifurcation/phase_diagram.png")
    println("Saved: imgs/bio/bifurcation/phase_diagram.png")
    
    # ----- 8. Summary: 4-panel bifurcation overview -----
    println("\n--- Creating Summary Plot ---")
    
    p1 = plot(xlabel="IP3 (μM)", ylabel="ATP:ADP Mean", title="Bifurcation Diagram",
              legend=:topleft, left_margin=8Plots.mm, bottom_margin=6Plots.mm)
    p2 = plot(xlabel="IP3 (μM)", ylabel="Amplitude", title="Oscillation Amplitude",
              legend=false, left_margin=8Plots.mm, bottom_margin=6Plots.mm)
    p3 = plot(xlabel="IP3 (μM)", ylabel="Frequency (Hz)", title="Dominant Frequency",
              legend=false, left_margin=8Plots.mm, bottom_margin=6Plots.mm)
    p4 = plot(xlabel="IP3 (μM)", ylabel="# Peaks", title="Peak Count",
              legend=false, left_margin=8Plots.mm, bottom_margin=6Plots.mm)
    
    for (nt, df) in results_all
        c = get(color_map, nt, :auto)
        
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            plot!(p1, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                  label=string(nt), linewidth=2, marker=:circle, markersize=2, color=c)
        end
        
        valid = .!isnan.(df.oscillation_amplitude)
        if sum(valid) > 0
            plot!(p2, df.ip3[valid], df.oscillation_amplitude[valid],
                  linewidth=2, marker=:circle, markersize=2, color=c)
        end
        
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            plot!(p3, df.ip3[valid], df.dominant_freq[valid],
                  linewidth=2, marker=:circle, markersize=2, color=c)
        end
        
        plot!(p4, df.ip3, df.n_peaks, linewidth=2, marker=:circle, markersize=2, color=c)
    end
    
    # Add threshold lines
    hline!(p2, [0.1], linestyle=:dash, color=:gray, alpha=0.5)
    hline!(p4, [3], linestyle=:dash, color=:gray, alpha=0.5)
    
    p_summary = plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 900))
    savefig(p_summary, "imgs/bio/bifurcation/bifurcation_summary.png")
    println("Saved: imgs/bio/bifurcation/bifurcation_summary.png")
    
    println("\n" * "="^80)
    println("BIFURCATION ANALYSIS COMPLETE")
    println("="^80)
    
    # Print summary table
    println("\nBifurcation Summary:")
    println("-"^60)
    for (nt, df) in results_all
        osc_count = sum(df.is_oscillating)
        total = nrow(df)
        
        # Find approximate bifurcation point (where oscillation stops)
        osc_to_stable = findfirst(i -> df.is_oscillating[i] && !df.is_oscillating[i+1], 1:nrow(df)-1)
        stable_to_osc = findfirst(i -> !df.is_oscillating[i] && df.is_oscillating[i+1], 1:nrow(df)-1)
        
        println("$nt:")
        println("  Oscillating: $osc_count / $total IP3 values")
        if !isnothing(osc_to_stable)
            println("  Oscillation stops at IP3 ≈ $(df.ip3[osc_to_stable]) - $(df.ip3[osc_to_stable+1]) μM")
        end
        if !isnothing(stable_to_osc)
            println("  Oscillation starts at IP3 ≈ $(df.ip3[stable_to_osc]) - $(df.ip3[stable_to_osc+1]) μM")
        end
    end
    
    return results_all
end

# ============================================================================
# RUN ALL EXPERIMENTS
# ============================================================================

# Run the updated IP3 scan
println("\n" * "="^80)
println("RUNNING COMPREHENSIVE IP3 SCAN")
println("="^80)

scan_noise_strength_experiment() # Task 1
scan_ip3_experiment()            # Task 2, 3, 4
analyze_isi_distribution()       # Task 4 (Distribution)
scan_bifurcation_experiment()    # Task 5

println("\n" * "="^80)
println("ALL EXPERIMENTS COMPLETE!")
println("="^80)
println("\nResults saved in:")
println("  - imgs/bio/scan/")
println("  - imgs/bio/isi_analysis/")
println("  - imgs/bio/noise_analysis/")