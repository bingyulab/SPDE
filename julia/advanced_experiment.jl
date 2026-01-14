include("experiment.jl")
using JLD2: @save, @load

# Note: For IP3 equals to 0.1, 0.2, the system is unstable with state_dependent noise. 
# for other IP3 values, the system may also be unstable.
# Till IP3 = 0.6, the system becomes stable with jump noise.

# ============================================================================
# CONSOLIDATED EXPERIMENT: Single IP3 Scan with All Analyses
# ============================================================================

"""
Run a single comprehensive IP3 scan that collects all data needed for:
- Task 1: Noise strength analysis
- Task 2: IP3 vs Mean ATP:ADP Ratio
- Task 3: Frequency analysis
- Task 4: ISI distribution & Kramers escape
- Task 5: Bifurcation analysis
- Bifurcation point dynamics

All analyses use the SAME simulation data - no redundant runs.
"""

function run_comprehensive_analysis(; 
    ip3_range = 0.1:0.1:2.0,
    n_ensemble::Int = 10,
    tspan_total::Float64 = 8000.0,
    min_success_rate::Float64 = 0.3,
)
    println("\n" * "="^80)
    println("COMPREHENSIVE IP3 ANALYSIS (CONSOLIDATED)")
    println("  IP3 range: $(first(ip3_range)) to $(last(ip3_range))")
    println("  Ensemble size: $n_ensemble runs per condition")
    println("  Simulation time: $(tspan_total)s per run")
    println("="^80)

    # =========================================================================
    # STEP 1: RUN ALL SIMULATIONS ONCE AND STORE RAW DATA
    # =========================================================================
    println("\n--- PHASE 1: Running Simulations ---")
    
    # Store raw solution data for reuse
    raw_data = Dict{Symbol, Dict{Float64, Vector{DataFrame}}}()
    stability_report = Dict{Symbol, Dict{Float64, Tuple{Int, Int}}}()  # (success, total)
    
    for noise_type in noise_list
        println("\nProcessing $noise_type...")
        raw_data[noise_type] = Dict{Float64, Vector{DataFrame}}()
        stability_report[noise_type] = Dict{Float64, Tuple{Int, Int}}()
                
        for val in ip3_range
            print("  IP3 = $val : ")
            raw_data[noise_type][val] = DataFrame[]
            
            n_success = 0
            n_total = n_ensemble

            for run in 1:n_ensemble
                current_seeds = Dict(k => seed_dict[k] + run * 100 for k in keys(seed_dict))
                sol, _ = simulate_model(noise_type; 
                                        tspan=(0.0, tspan_total), 
                                        ip3_val=val, 
                                        seeds=current_seeds)
                # Check if solution is valid (not just converged, but reasonable)
                df_sol = DataFrame(sol)
                
                # Validity checks
                is_valid = true

                # Check for NaN/Inf
                for col in names(df_sol)
                    if col != "timestamp" && any(!isfinite, df_sol[!, col])
                        is_valid = false
                        break
                    end
                end
                
                # Check for negative concentrations (biological constraint)
                for var in ["cac", "caer", "atpc", "adpc"]
                    col_idx = df_find_column(df_sol, var)
                    if !isnothing(col_idx) && any(x -> x < -1e-6, df_sol[!, col_idx])
                        is_valid = false
                        break
                    end
                end
                
                # Check for unreasonably large values (explosion)
                caer_idx = df_find_column(df_sol, "caer")
                if !isnothing(caer_idx) && maximum(df_sol[!, caer_idx]) > 1000
                    is_valid = false
                end
                
                if is_valid
                    push!(raw_data[noise_type][val], df_sol)
                    n_success += 1
                end
            end
            stability_report[noise_type][val] = (n_success, n_total)
            
            # Status indicator
            success_rate = n_success / n_total
            if success_rate >= 0.8
                status = "✓"
            elseif success_rate >= min_success_rate
                status = "⚠"
            else
                status = "✗"
            end
            
            println("$status ($n_success/$n_total successful)")
        end
    end
        
    # =========================================================================
    # STEP 2: EXTRACT ALL METRICS FROM STORED DATA
    # =========================================================================
    println("\n--- PHASE 2: Extracting Metrics ---")
    
    results_all = Dict{Symbol, DataFrame}()
    
    for noise_type in noise_list
        println("\nExtracting metrics for $noise_type...")
        
        df_res = DataFrame(
            ip3 = Float64[],
            # Task 2: ATP:ADP ratio
            atp_adp_ratio_mean = Float64[],
            atp_adp_ratio_std = Float64[],
            atp_adp_ratio_sem = Float64[],
            atp_adp_ratio_min = Float64[],
            atp_adp_ratio_max = Float64[],
            # Task 3: Frequency
            dominant_freq = Float64[],
            freq_variance = Float64[],
            # Task 4: ISI/Kramers
            isi_entropy = Float64[],
            isi_entropy_std = Float64[],
            isi_mean = Float64[],
            isi_std = Float64[],
            isi_cv = Float64[],
            escape_rate = Float64[],
            escape_rate_std = Float64[],
            n_peaks = Int[],
            n_events = Int[],
            # Task 5: Bifurcation
            oscillation_amplitude = Float64[],
            is_oscillating = Bool[]
        )
        
        # Add columns for important variables
        for var in important_variables
            df_res[!, var] = Float64[]
            df_res[!, "$(var)_std"] = Float64[]
        end
        
        for val in ip3_range
            dfs = get(raw_data[noise_type], val, DataFrame[])
            
            if isempty(dfs)
                # Push NaN row
                row = create_nan_row(val)
                push!(df_res, row)
                continue
            end
            
            # Aggregate metrics across ensemble
            row = extract_all_metrics(dfs, val)
            push!(df_res, row)
        end
        
        results_all[noise_type] = df_res
    end
    
    @save "results/advanced_analysis.jld2" results_all raw_data stability_report
    println("Saved results to results/advanced_analysis.jld2")
    return results_all, raw_data, stability_report
end 

function generate_all_plots(results_all, raw_data, ip3_range, n_ensemble::Int=10)
    # =========================================================================
    # STEP 3: GENERATE ALL PLOTS
    # =========================================================================
    println("\n--- PHASE 3: Generating Plots ---")
    
    mkpath("imgs/bio/scan/")
    mkpath("imgs/bio/isi_analysis/")
    mkpath("imgs/bio/bifurcation/")
    
    # Task 2: IP3 vs ATP:ADP Ratio
    plot_task2_atp_adp_ratio(results_all, n_ensemble)
    
    # Task 3: Frequency Analysis
    plot_task3_frequency(results_all)
    
    # Task 4: ISI & Kramers Analysis
    plot_task4_kramers(results_all, raw_data, ip3_range)
    
    # Task 5: Bifurcation Analysis
    plot_bifurcation_analysis(raw_data, results_all, ip3_range)
    
    # Important Variables vs IP3
    plot_important_variables(results_all)
   
    # Summary plots
    plot_summary(results_all)
    
    println("\n" * "="^80)
    println("COMPREHENSIVE ANALYSIS COMPLETE")
    println("="^80)    
end

"""
Print a summary of which IP3/noise combinations are stable
"""
function print_stability_summary(stability_report, min_success_rate)
    println("\n" * "="^80)
    println("STABILITY SUMMARY")
    println("="^80)
    
    # Header
    ip3_values = sort(collect(keys(first(values(stability_report)))))
    
    println("\nSuccess rates (✓ ≥80%, ⚠ ≥$(Int(min_success_rate*100))%, ✗ <$(Int(min_success_rate*100))%):\t")
    
    # Print table header
    print(rpad("Noise Type", 18))
    for ip3 in ip3_values
        print(rpad("$(ip3)", 8))
    end
    println()
    println("-"^(18 + 8*length(ip3_values)))
    
    # Print each noise type
    unstable_combinations = []
    
    for nt in noise_list
        if !haskey(stability_report, nt)
            continue
        end
        
        print(rpad(string(nt), 18))
        
        for ip3 in ip3_values
            if haskey(stability_report[nt], ip3)
                n_success, n_total = stability_report[nt][ip3]
                rate = n_success / n_total
                
                if rate >= 0.8
                    symbol = "✓"
                elseif rate >= min_success_rate
                    symbol = "⚠"
                else
                    symbol = "✗"
                    push!(unstable_combinations, (nt, ip3, rate))
                end
                
                print(rpad("$symbol$(Int(round(rate*100)))%", 8))
            else
                print(rpad("-", 8))
            end
        end
        println()
    end
    
    # Report unstable combinations
    if !isempty(unstable_combinations)
        println("\n⚠️  UNSTABLE COMBINATIONS (excluded from analysis):")
        for (nt, ip3, rate) in unstable_combinations
            println("   - $nt at IP3=$ip3 ($(Int(round(rate*100)))% success)")
        end
    end
    
    println("="^80)
end

# ============================================================================
# HELPER: Create NaN row for failed simulations
# ============================================================================
function create_nan_row(ip3_val::Float64)
    row = Dict{String, Any}(
        "ip3" => ip3_val,
        "atp_adp_ratio_mean" => NaN,
        "atp_adp_ratio_std" => NaN,
        "atp_adp_ratio_sem" => NaN,
        "atp_adp_ratio_min" => NaN,
        "atp_adp_ratio_max" => NaN,
        "dominant_freq" => NaN,
        "freq_variance" => NaN,
        "isi_entropy" => NaN,
        "isi_entropy_std" => NaN,
        "isi_mean" => NaN,
        "isi_std" => NaN,
        "isi_cv" => NaN,
        "escape_rate" => NaN,
        "escape_rate_std" => NaN,
        "n_peaks" => 0,
        "n_events" => 0,
        "oscillation_amplitude" => NaN,
        "is_oscillating" => false
    )
    for var in important_variables
        row[var] = NaN
        row["$(var)_std"] = NaN
    end
    return row
end

# ============================================================================
# HELPER: Extract all metrics from ensemble of DataFrames
# ============================================================================
function extract_all_metrics(dfs::Vector{DataFrame}, ip3_val::Float64)
    # If no valid data, return NaN row
    if isempty(dfs)
        return create_nan_row(ip3_val)
    end

    # CRITICAL FIX: Use TIME-BASED window (last 4000 seconds), not point-based
    analysis_window = 4000.0  # seconds

    # Storage for ensemble aggregation
    ratio_means = Float64[]
    ratio_stds = Float64[]
    ratio_mins = Float64[]
    ratio_maxs = Float64[]
    freq_vals = Float64[]
    freq_var_vals = Float64[]
    entropy_vals = Float64[]
    rate_vals = Float64[]
    peak_counts = Int[]
    isi_means = Float64[]
    isi_stds = Float64[]
    all_isi = Float64[]  # Pooled ISI for detailed analysis
    var_means = Dict(var => Float64[] for var in important_variables)
    
    for df_sol in dfs
        ratio = calculate_atp_adp_ratio(df_sol, analysis_window)
        
        if !isnothing(ratio) && length(ratio) > 10
            # FIX: Store mean of the ratio vector, not ratio of means
            push!(ratio_means, mean(ratio))
            push!(ratio_stds, std(ratio))
            push!(ratio_mins, minimum(ratio))
            push!(ratio_maxs, maximum(ratio))
            
            # Frequency analysis
            freq, fvar = estimate_frequency(ratio, 1.0)
            if isfinite(freq) && freq > 0
                push!(freq_vals, freq)
                push!(freq_var_vals, fvar)
            end
            
            # Peak detection
            peaks = detect_peaks(ratio; threshold_percentile=75.0)
            push!(peak_counts, length(peaks))
            
            # ISI analysis (use full windowed data for better statistics)
            events = detect_escape_events(ratio; low_percentile=25.0, high_percentile=75.0)
            if length(events) >= 3
                isi = calculate_isi(events, 1.0)
                append!(all_isi, isi)
                push!(isi_means, mean(isi))
                push!(isi_stds, std(isi))
                push!(entropy_vals, calculate_entropy(isi; n_bins=15))
                push!(rate_vals, 1.0 / mean(isi))
            end
        end
        
        # Important variables (using windowed data)
        # Apply consistent TIME-BASED windowing
        t_max = maximum(df_sol.timestamp)
        t_start = max(0.0, t_max - analysis_window)
        mask = df_sol.timestamp .>= t_start
        
        if sum(mask) < 10
            continue
        end
        
        df_windowed = df_sol[mask, :]
        
        for var in important_variables
            col_idx = df_find_column(df_windowed, var)
            if !isnothing(col_idx)
                var_data = df_windowed[!, col_idx]
                valid_data = filter(isfinite, var_data)
                if length(valid_data) > 0
                    push!(var_means[var], mean(valid_data))
                end
            end
        end
    end
    
    # Aggregate results
    n_runs = length(ratio_means)

    # Check if we have enough data
    if length(ratio_means) == 0
        println("    ⚠️  No valid runs for IP3=$ip3_val")
        return create_nan_row(ip3_val)
    end
    
    # Oscillation detection
    mean_amplitude = n_runs > 0 ? mean(ratio_maxs) - mean(ratio_mins) : NaN
    is_oscillating = mean_amplitude > 0.1 && (n_runs > 0 ? mean(peak_counts) >= 3 : false)
    
    row = Dict{String, Any}(
        "ip3" => ip3_val,
        "atp_adp_ratio_mean" => n_runs > 0 ? mean(ratio_means) : NaN,
        "atp_adp_ratio_std" => n_runs > 0 ? mean(ratio_stds) : NaN,
        "atp_adp_ratio_sem" => n_runs > 1 ? std(ratio_means) / sqrt(n_runs) : NaN,
        "atp_adp_ratio_min" => n_runs > 0 ? mean(ratio_mins) : NaN,
        "atp_adp_ratio_max" => n_runs > 0 ? mean(ratio_maxs) : NaN,
        "dominant_freq" => length(freq_vals) > 0 ? mean(freq_vals) : NaN,
        "freq_variance" => length(freq_var_vals) > 0 ? mean(freq_var_vals) : NaN,
        "isi_entropy" => length(entropy_vals) > 0 ? mean(entropy_vals) : NaN,
        "isi_entropy_std" => length(entropy_vals) > 1 ? std(entropy_vals) : NaN,
        "isi_mean" => length(isi_means) > 0 ? mean(isi_means) : NaN,
        "isi_std" => length(isi_stds) > 0 ? mean(isi_stds) : NaN,
        "isi_cv" => length(isi_means) > 0 && length(isi_stds) > 0 ? 
                    mean(isi_stds) / mean(isi_means) : NaN,
        "escape_rate" => length(rate_vals) > 0 ? mean(rate_vals) : NaN,
        "escape_rate_std" => length(rate_vals) > 1 ? std(rate_vals) : NaN,
        "n_peaks" => length(peak_counts) > 0 ? round(Int, mean(peak_counts)) : 0,
        "n_events" => length(all_isi),
        "oscillation_amplitude" => mean_amplitude,
        "is_oscillating" => is_oscillating
    )
    
    for var in important_variables
        vals = var_means[var]
        row[var] = length(vals) > 0 ? mean(vals) : NaN
        row["$(var)_std"] = length(vals) > 1 ? std(vals) : NaN
    end
    return row
end

# ============================================================================
# PLOTTING FUNCTIONS (Modular)
# ============================================================================

function plot_task2_atp_adp_ratio(results_all, n_ensemble)
    println("\n--- Task 2: IP3 vs ATP:ADP Ratio ---")
    
    # Define IP3 range filters for each noise type
    ip3_min_filter = Dict(
        :none => 0.1,
        :additive => 0.2,           # Start from 0.2
        :multiplicative => 0.1,
        :state_dependent => 0.3,    # Start from 0.3
        :jump => 0.1
    )

    # Separate plots per noise type
    for (nt, df) in results_all
        # Apply IP3 filter for this noise type
        ip3_min = get(ip3_min_filter, nt, 0.1)
        ip3_filter = df.ip3 .>= ip3_min
        valid = .!isnan.(df.atp_adp_ratio_mean) .& ip3_filter

        if sum(valid) > 0
            std_vals = coalesce.(df.atp_adp_ratio_std[valid], 0.0)
            y_min_data = minimum(df.atp_adp_ratio_mean[valid] .- std_vals)
            y_max_data = maximum(df.atp_adp_ratio_mean[valid] .+ std_vals)
            y_min = min(0.0, y_min_data * 1.1)  # Include 0 or go lower if data is negative
            y_max = y_max_data * 1.1
            
            p = plot(xlabel="IP3 (μM)", ylabel="Mean ATP:ADP Ratio",
                    title="$nt: IP3 vs Mean ATP:ADP Ratio (n=$n_ensemble runs)",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                    ylims=(y_min, y_max))
            plot!(p, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                marker=:circle, linewidth=2, markersize=5,
                ribbon=df.atp_adp_ratio_std[valid], fillalpha=0.3,
                color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task2_ip3_vs_ratio_$(nt).png")
            println("Saved: imgs/bio/scan/task2_ip3_vs_ratio_$(nt).png")
        end
    end

    # Combined plot
    y_min = 0.0
    y_max = 0.0
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            std_vals = coalesce.(df.atp_adp_ratio_std[valid], 0.0)
            y_min = min(y_min, minimum(df.atp_adp_ratio_mean[valid] .- std_vals))
            y_max = max(y_max, maximum(df.atp_adp_ratio_mean[valid] .+ std_vals))
        end
    end
    y_min = min(0.0, y_min * 1.1)
    y_max *= 1.1

    p_combined = plot(xlabel="IP3 (μM)", ylabel="Mean ATP:ADP Ratio",
                    title="IP3 vs ATP:ADP Ratio (All Noise Types)",
                    legend=:topright, size=(900, 600),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                    ylims=(y_min, y_max))
    for (nt, df) in results_all
        valid = .!isnan.(df.atp_adp_ratio_mean)
        if sum(valid) > 0
            plot!(p_combined, df.ip3[valid], df.atp_adp_ratio_mean[valid],
                label=string(nt), marker=:circle, linewidth=2,
                color=get(color_map, nt, :auto))
        end
    end
    savefig(p_combined, "imgs/bio/scan/task2_ip3_vs_ratio_combined.png")
    println("Saved: imgs/bio/scan/task2_ip3_vs_ratio_combined.png")
end

function plot_task3_frequency(results_all)
    println("\n--- Task 3: Frequency Analysis ---")
    
    stochastic_noises = filter(x -> x != :none, collect(keys(results_all)))
    
    # ----- Separate plots per noise type -----
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            y_max = maximum(df.dominant_freq[valid]) * 1.2
            p = plot(xlabel="IP3 (μM)", ylabel="Dominant Frequency (Hz)",
                    title="$nt: Dominant Frequency vs IP3",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                    ylims=(0, y_max))
            plot!(p, df.ip3[valid], df.dominant_freq[valid],
                  marker=:circle, linewidth=2, markersize=5,
                  color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task3_frequency_$(nt).png")
            println("Saved: imgs/bio/scan/task3_frequency_$(nt).png")
        end
    end
    
    # ----- Combined: Dominant Frequency (1 × N layout) -----
    subplots_freq = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        
        y_max = sum(valid) > 0 ? maximum(df.dominant_freq[valid]) * 1.2 : 1.0
        p = plot(xlabel="IP3 (μM)", ylabel="Freq (Hz)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                ylims=(0, y_max))
        
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
    
    if length(subplots_freq) > 0
        p_freq = plot(subplots_freq..., layout=(1, length(subplots_freq)),
                     size=(350*length(subplots_freq), 400),
                     plot_title="Dominant Frequency vs IP3")
        savefig(p_freq, "imgs/bio/scan/task3_frequency_analysis.png")
        println("Saved: imgs/bio/scan/task3_frequency_analysis.png")
    end
    
    # ----- Separate plots for Frequency Variance -----
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.freq_variance)
        if sum(valid) > 0
            y_max = maximum(df.freq_variance[valid]) * 1.2
            p = plot(xlabel="IP3 (μM)", ylabel="Frequency Variance",
                    title="$nt: Frequency Variance vs IP3",
                    legend=false, size=(700, 500),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                    ylims=(0, y_max))
            plot!(p, df.ip3[valid], df.freq_variance[valid],
                  marker=:circle, linewidth=2, markersize=5,
                  color=get(color_map, nt, :blue))
            savefig(p, "imgs/bio/scan/task3_freq_variance_$(nt).png")
            println("Saved: imgs/bio/scan/task3_freq_variance_$(nt).png")
        end
    end
    
    # ----- Combined: Frequency Variance (1 × N layout) -----
    subplots_fvar = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.freq_variance)

        y_max = sum(valid) > 0 ? maximum(df.freq_variance[valid]) * 1.2 : 1.0
        p = plot(xlabel="IP3 (μM)", ylabel="Freq Var",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                ylims=(0, y_max))
        
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
    
    if length(subplots_fvar) > 0
        p_fvar = plot(subplots_fvar..., layout=(1, length(subplots_fvar)),
                     size=(350*length(subplots_fvar), 400),
                     plot_title="Frequency Variance vs IP3 (Low=Regular, High=Irregular)")
        savefig(p_fvar, "imgs/bio/scan/task3_freq_variance_analysis.png")
        println("Saved: imgs/bio/scan/task3_freq_variance_analysis.png")
    end
end

function plot_task4_kramers(results_all, raw_data, ip3_range)
    println("\n--- Task 4: ISI & Kramers Escape Analysis ---")
    
    stochastic_noises = filter(x -> x != :none, collect(keys(results_all)))
    all_noises = collect(keys(results_all))
    
    # =========================================================================
    # SEPARATE PLOTS: Escape Rate vs IP3 (per noise type)
    # =========================================================================
    for nt in stochastic_noises
        df = results_all[nt]
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
    p_rate = plot(xlabel="IP3 (μM)", ylabel="Escape Rate (1/s)",
                  title="Escape Rate vs IP3 (All Noise Types)",
                  legend=:topright, size=(900, 600),
                  left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        if nt == :none; continue; end
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            plot!(p_rate, df.ip3[valid], df.escape_rate[valid],
                  label=string(nt), marker=:circle, linewidth=2,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_rate, "imgs/bio/scan/task4_escape_rate_vs_ip3_combined.png")
    println("Saved: imgs/bio/scan/task4_escape_rate_vs_ip3_combined.png")
    
    # =========================================================================
    # SEPARATE PLOTS: ISI Entropy vs IP3 (per noise type)
    # =========================================================================
    for nt in all_noises
        df = results_all[nt]
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
    p_entropy = plot(xlabel="IP3 (μM)", ylabel="ISI Entropy (bits)",
                     title="ISI Entropy vs IP3 (All Noise Types)",
                     legend=:bottomright, size=(900, 600),
                     left_margin=10Plots.mm, bottom_margin=8Plots.mm)
    
    for (nt, df) in results_all
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            plot!(p_entropy, df.ip3[valid], df.isi_entropy[valid],
                  label=string(nt), marker=:circle, linewidth=2,
                  color=get(color_map, nt, :auto))
        end
    end
    savefig(p_entropy, "imgs/bio/scan/task4_entropy_vs_ip3_combined.png")
    println("Saved: imgs/bio/scan/task4_entropy_vs_ip3_combined.png")
    
    # =========================================================================
    # NEW: TWO-ROW LAYOUT - All noise types together
    # Row 1: log(ATP:ADP) as X-axis, Rate/Entropy as Y-axis, IP3 as colorbar
    # Row 2: IP3 as X-axis, Rate/Entropy as Y-axis, log(ATP:ADP) as colorbar
    # =========================================================================
    println("\n--- Kramers Analysis: Two-Row Layout (All Noise Types) ---")
    
    # Get ranges for consistent colorbars
    all_energy_vals = Float64[]
    all_ip3_vals = Float64[]
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.atp_adp_ratio_mean) .& (df.atp_adp_ratio_mean .> 0)
        if sum(valid) > 0
            append!(all_energy_vals, df.atp_adp_ratio_mean[valid])
            append!(all_ip3_vals, df.ip3[valid])
        end
    end
    
    energy_min = length(all_energy_vals) > 0 ? minimum(all_energy_vals) : 1.0
    energy_max = length(all_energy_vals) > 0 ? maximum(all_energy_vals) : 10.0
    ip3_min = length(all_ip3_vals) > 0 ? minimum(all_ip3_vals) : 0.1
    ip3_max = length(all_ip3_vals) > 0 ? maximum(all_ip3_vals) : 2.0
    
    # Use log scale for energy
    log_energy_min = log10(max(energy_min, 0.1))
    log_energy_max = log10(max(energy_max, 1.0))
    
    n_noise = length(stochastic_noises)
    
    if n_noise == 0
        println("  No stochastic noise types available")
        return
    end
    
    raw_widths = vcat(ones(n_noise), [0.4])
    widths = raw_widths ./ sum(raw_widths)  # Normalize!
    
    # =========================================================================
    # ESCAPE RATE: Two-row layout
    # Row 1: log(ATP:ADP) as X-axis, IP3 as colorbar
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar
    # =========================================================================
    
    # First, calculate global y_max for escape rate across all noise types
    rate_y_max = 0.0
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        if sum(valid) > 0
            rate_y_max = max(rate_y_max, maximum(df.escape_rate[valid]))
        end
    end
    rate_y_max *= 1.1

    # Row 1: log(ATP:ADP) as X-axis, IP3 as colorbar
    row1_plots = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="log₁₀(ATP:ADP)", ylabel="Escape Rate (1/s)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm, 
                ylims=(0, rate_y_max))
        
        if sum(valid) > 2
            x_log_energy = log10.(df.atp_adp_ratio_mean[valid])
            y_rate = df.escape_rate[valid]
            z_ip3 = df.ip3[valid]
            
            scatter!(p, x_log_energy, y_rate, markersize=7,
                    zcolor=z_ip3, color=:viridis,
                    clims=(ip3_min, ip3_max), colorbar=false, alpha=0.9)
            
            # Trend line (Kramers: should be negative slope)
            if length(x_log_energy) > 2
                X_mat = hcat(ones(length(x_log_energy)), x_log_energy)
                coeffs = X_mat \ y_rate
                x_line = range(minimum(x_log_energy), maximum(x_log_energy), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x_log_energy, y_rate)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row1_plots, p)
    end
    
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar
    row2_plots = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="IP3 (μM)", ylabel="Escape Rate (1/s)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm,
                ylims=(0, rate_y_max))
        
        if sum(valid) > 2
            x_ip3 = df.ip3[valid]
            y_rate = df.escape_rate[valid]
            log_energy = log10.(df.atp_adp_ratio_mean[valid])
            
            scatter!(p, x_ip3, y_rate, markersize=7,
                    zcolor=log_energy, color=:plasma,
                    clims=(log_energy_min, log_energy_max), colorbar=false, alpha=0.9)
            
            # Trend line
            if length(x_ip3) > 2
                X_mat = hcat(ones(length(x_ip3)), x_ip3)
                coeffs = X_mat \ y_rate
                x_line = range(minimum(x_ip3), maximum(x_ip3), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:black)
                
                corr = cor(x_ip3, y_rate)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row2_plots, p)
    end
    
    # Create colorbars
    p_cbar1 = scatter([NaN], [NaN], zcolor=[ip3_min], clims=(ip3_min, ip3_max),
                     color=:viridis, colorbar=true, colorbar_title="IP3 (μM)",
                     framestyle=:none, label="", markersize=0)
    
    p_cbar2 = scatter([NaN], [NaN], zcolor=[log_energy_min], clims=(log_energy_min, log_energy_max),
                     color=:plasma, colorbar=true, colorbar_title="log₁₀(ATP:ADP)",
                     framestyle=:none, label="", markersize=0)
    
    # Combine into 2-row layout
    if n_noise > 0
        all_plots = vcat(row1_plots, [p_cbar1], row2_plots, [p_cbar2])
                
        p_dual_row = plot(all_plots..., 
                         layout=grid(2, n_noise+1, widths=widths),
                         size=(350*n_noise + 80, 800),
                         plot_title="Escape Rate: Row1=f(Energy), Row2=f(IP3)",
                         plot_titlefontsize=11)
        savefig(p_dual_row, "imgs/bio/scan/task4_escape_rate_two_row.png")
        println("Saved: imgs/bio/scan/task4_escape_rate_two_row.png")
    end
    
    # =========================================================================
    # ISI ENTROPY: Two-row layout
    # Row 1: log(ATP:ADP) as X-axis, IP3 as colorbar
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar
    # =========================================================================
    # First, calculate global y_max for escape rate across all noise types
    rate_y_max = 0.0
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.isi_entropy) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.isi_entropy .> 0) .& (df.atp_adp_ratio_mean .> 0)
        if sum(valid) > 0
            rate_y_max = max(rate_y_max, maximum(df.isi_entropy[valid]))
        end
    end
    rate_y_max *= 1.1
    
    # Row 1: log(ATP:ADP) as X-axis, IP3 as colorbar
    row1_entropy = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.isi_entropy) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="log₁₀(ATP:ADP)", ylabel="ISI Entropy (bits)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm,
                ylims=(0, rate_y_max))
        
        if sum(valid) > 2
            x_log_energy = log10.(df.atp_adp_ratio_mean[valid])
            y_entropy = df.isi_entropy[valid]
            z_ip3 = df.ip3[valid]
            
            scatter!(p, x_log_energy, y_entropy, markersize=7,
                    zcolor=z_ip3, color=:viridis,
                    clims=(ip3_min, ip3_max), colorbar=false, alpha=0.9)
            
            if length(x_log_energy) > 2
                X_mat = hcat(ones(length(x_log_energy)), x_log_energy)
                coeffs = X_mat \ y_entropy
                x_line = range(minimum(x_log_energy), maximum(x_log_energy), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x_log_energy, y_entropy)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row1_entropy, p)
    end
    
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar
    row2_entropy = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.isi_entropy) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="IP3 (μM)", ylabel="ISI Entropy (bits)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm,
                ylims=(0, rate_y_max))
        
        if sum(valid) > 2
            x_ip3 = df.ip3[valid]
            y_entropy = df.isi_entropy[valid]
            log_energy = log10.(df.atp_adp_ratio_mean[valid])
            
            scatter!(p, x_ip3, y_entropy, markersize=7,
                    zcolor=log_energy, color=:plasma,
                    clims=(log_energy_min, log_energy_max), colorbar=false, alpha=0.9)
            
            if length(x_ip3) > 2
                X_mat = hcat(ones(length(x_ip3)), x_ip3)
                coeffs = X_mat \ y_entropy
                x_line = range(minimum(x_ip3), maximum(x_ip3), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:black)
                
                corr = cor(x_ip3, y_entropy)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row2_entropy, p)
    end
    
    if n_noise > 0
        all_entropy_plots = vcat(row1_entropy, [p_cbar1], row2_entropy, [p_cbar2])
        
        p_entropy_dual = plot(all_entropy_plots..., 
                             layout=grid(2, n_noise+1, widths=widths),
                             size=(350*n_noise + 80, 800),
                             plot_title="ISI Entropy: Row1=f(Energy), Row2=f(IP3)",
                             plot_titlefontsize=11)
        savefig(p_entropy_dual, "imgs/bio/scan/task4_entropy_two_row.png")
        println("Saved: imgs/bio/scan/task4_entropy_two_row.png")
    end
    
    # =========================================================================
    # KRAMERS-ARRHENIUS: Two-row layout (log(Rate))
    # Row 1: log(ATP:ADP) as X-axis, IP3 as colorbar - THIS IS THE ARRHENIUS PLOT
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar
    # =========================================================================
    
    # Row 1: log(ATP:ADP) as X-axis (Arrhenius plot: log(Rate) vs Energy)
    row1_kramers = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="log₁₀(ATP:ADP)", ylabel="log(Rate)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm)
        
        if sum(valid) > 2
            x_log_energy = log10.(df.atp_adp_ratio_mean[valid])
            y_log_rate = log.(df.escape_rate[valid])
            z_ip3 = df.ip3[valid]
            
            finite_mask = isfinite.(y_log_rate)
            if sum(finite_mask) > 2
                x_log_energy = x_log_energy[finite_mask]
                y_log_rate = y_log_rate[finite_mask]
                z_ip3 = z_ip3[finite_mask]
                
                scatter!(p, x_log_energy, y_log_rate, markersize=7,
                        zcolor=z_ip3, color=:viridis,
                        clims=(ip3_min, ip3_max), colorbar=false, alpha=0.9)
                
                # Arrhenius: log(Rate) = A - E/kT → negative slope means Kramers
                X_mat = hcat(ones(length(x_log_energy)), x_log_energy)
                coeffs = X_mat \ y_log_rate
                x_line = range(minimum(x_log_energy), maximum(x_log_energy), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:red)
                
                corr = cor(x_log_energy, y_log_rate)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row1_kramers, p)
    end
    
    # Row 2: IP3 as X-axis, log(ATP:ADP) as colorbar  
    row2_kramers = []
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        
        p = plot(xlabel="IP3 (μM)", ylabel="log(Rate)",
                title="$nt", legend=false, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=3Plots.mm)
        
        if sum(valid) > 2
            x_ip3 = df.ip3[valid]
            y_log_rate = log.(df.escape_rate[valid])
            energy = df.atp_adp_ratio_mean[valid]
            
            finite_mask = isfinite.(y_log_rate)
            if sum(finite_mask) > 2
                x_ip3 = x_ip3[finite_mask]
                y_log_rate = y_log_rate[finite_mask]
                energy = energy[finite_mask]
                log_energy = log10.(energy)
                
                scatter!(p, x_ip3, y_log_rate, markersize=7,
                        zcolor=log_energy, color=:plasma,
                        clims=(log_energy_min, log_energy_max), colorbar=false, alpha=0.9)
                
                X_mat = hcat(ones(length(x_ip3)), x_ip3)
                coeffs = X_mat \ y_log_rate
                x_line = range(minimum(x_ip3), maximum(x_ip3), length=50)
                y_line = coeffs[1] .+ coeffs[2] .* x_line
                plot!(p, x_line, y_line, linestyle=:dash, linewidth=2, color=:black)
                
                corr = cor(x_ip3, y_log_rate)
                annotate!(p, :bottomright, text("r=$(round(corr, digits=2))", 8))
            end
        end
        push!(row2_kramers, p)
    end
    
    if n_noise > 0
        all_kramers_plots = vcat(row1_kramers, [p_cbar1], row2_kramers, [p_cbar2])
        
        p_kramers_dual = plot(all_kramers_plots..., 
                             layout=grid(2, n_noise+1, widths=widths),
                             size=(350*n_noise + 80, 800),
                             plot_title="Kramers-Arrhenius: Row1=log(R) vs log(E), Row2=log(R) vs IP3",
                             plot_titlefontsize=10)
        savefig(p_kramers_dual, "imgs/bio/scan/task4_kramers_two_row.png")
        println("Saved: imgs/bio/scan/task4_kramers_two_row.png")
    end
    
    # =========================================================================
    # Summary correlation table
    # =========================================================================
    println("\n--- Kramers Correlation Summary ---")
    println("-"^85)
    println(rpad("Noise Type", 18), rpad("logE-Rate", 12), rpad("logE-logR", 12), 
            rpad("IP3-Rate", 12), rpad("Kramers?", 18))
    println("-"^85)
    
    for nt in stochastic_noises
        df = results_all[nt]
        valid = .!isnan.(df.escape_rate) .& .!isnan.(df.atp_adp_ratio_mean) .& 
                (df.escape_rate .> 0) .& (df.atp_adp_ratio_mean .> 0)
        
        if sum(valid) > 2
            log_energy = log10.(df.atp_adp_ratio_mean[valid])
            rate = df.escape_rate[valid]
            log_rate = log.(rate)
            ip3 = df.ip3[valid]
            
            finite_mask = isfinite.(log_rate)
            
            corr_logE_rate = cor(log_energy, rate)
            corr_logE_logR = sum(finite_mask) > 2 ? cor(log_energy[finite_mask], log_rate[finite_mask]) : NaN
            corr_ip3_rate = cor(ip3, rate)
            
            # Kramers criterion: higher energy barrier → lower escape rate
            # In Arrhenius form: log(Rate) = A - E/kT → negative correlation
            kramers_status = if !isnan(corr_logE_logR) && corr_logE_logR < -0.3
                "✓ Yes (E↑→R↓)"
            elseif !isnan(corr_logE_logR) && corr_logE_logR > 0.3
                "✗ Anti-Kramers"
            else
                "~ Weak/None"
            end
            
            println(rpad(string(nt), 18), 
                    rpad("$(round(corr_logE_rate, digits=3))", 12),
                    rpad(isnan(corr_logE_logR) ? "N/A" : "$(round(corr_logE_logR, digits=3))", 12),
                    rpad("$(round(corr_ip3_rate, digits=3))", 12),
                    kramers_status)
        else
            println(rpad(string(nt), 18), rpad("N/A", 12), rpad("N/A", 12), rpad("N/A", 12), "Insufficient data")
        end
    end
    println("-"^85)
    
    # =========================================================================
    # ISI Histograms - FIXED for deterministic case
    # =========================================================================
    println("\n--- ISI Histograms ---")
    mkpath("imgs/bio/isi_analysis/")
    
    test_ip3_values = [0.7, 0.8, 0.9, 1.0, 1.2, 1.4]
    
    for noise_type in noise_list
        if !haskey(raw_data, noise_type)
            println("  $noise_type: No raw data available")
            continue
        end
        
        # Find best IP3 value with most ISI events
        best_ip3 = nothing
        best_isi = Float64[]
        
        for test_ip3 in test_ip3_values
            dfs = get(raw_data[noise_type], test_ip3, DataFrame[])
            if isempty(dfs); continue; end
            
            test_isi = Float64[]
            for df_sol in dfs
                ratio = calculate_atp_adp_ratio(df_sol, 4000.0)
                if isnothing(ratio) || length(ratio) < 100; continue; end
                
                events = detect_escape_events(ratio; low_percentile=25.0, high_percentile=75.0)
                if length(events) >= 2
                    isi = calculate_isi(events, 1.0)
                    valid_isi = filter(x -> x > 5.0 && x < 1000.0, isi)
                    append!(test_isi, valid_isi)
                end
            end
            
            if length(test_isi) > length(best_isi)
                best_isi = test_isi
                best_ip3 = test_ip3
            end
        end
        
        if length(best_isi) < 3
            println("  $noise_type: Too few ISI samples ($(length(best_isi))) across all IP3 values")
            continue
        end
        
        all_isi = best_isi
        representative_ip3 = best_ip3
        
        isi_mean = mean(all_isi)
        isi_std = std(all_isi)
        isi_cv = isi_std / max(isi_mean, 1e-10)
        
        # Check if deterministic (very low CV) - use stricter threshold
        is_deterministic = noise_type == :none
        
        if is_deterministic
            # For deterministic: use FFT to find true period instead of event detection
            println("  $noise_type: Deterministic system - using FFT for period estimation")
            
            # Get a clean trajectory for FFT analysis
            dfs = get(raw_data[noise_type], representative_ip3, DataFrame[])
            if !isempty(dfs)
                df_sol = first(dfs)
                ratio = calculate_atp_adp_ratio(df_sol, 4000.0)
                
                if !isnothing(ratio) && length(ratio) > 100
                    # Use FFT to find dominant period
                    dt = 1.0  # Assuming 1s sampling
                    n = length(ratio)
                    
                    # Remove mean and apply window
                    ratio_centered = ratio .- mean(ratio)
                    
                    # FFT
                    fft_result = abs.(fft(ratio_centered))
                    freqs = fftfreq(n, 1/dt)
                    
                    # Only positive frequencies, skip DC
                    pos_mask = freqs .> 0.001  # Skip very low frequencies
                    pos_freqs = freqs[pos_mask]
                    pos_power = fft_result[pos_mask]
                    
                    if length(pos_power) > 0
                        # Find dominant frequency
                        max_idx = argmax(pos_power)
                        dominant_freq = pos_freqs[max_idx]
                        true_period = 1.0 / dominant_freq
                        
                        println("    FFT dominant frequency: $(round(dominant_freq, digits=4)) Hz")
                        println("    True period: $(round(true_period, digits=2)) s")
                        println("    Event-based mean: $(round(isi_mean, digits=2)) ± $(round(isi_std, digits=2)) s")
                        println("    CV from events: $(round(isi_cv, digits=4))")
                        
                        # Create a clean visualization
                        p_det = plot(layout=(1, 2), size=(1400, 500), 
                                    plot_title="$noise_type: Deterministic Oscillation Analysis (IP3=$representative_ip3)")
                        
                        # Left: Power spectrum
                        plot!(p_det[1], pos_freqs[1:min(100, length(pos_freqs))], 
                              pos_power[1:min(100, length(pos_power))],
                              xlabel="Frequency (Hz)", ylabel="Power",
                              title="Power Spectrum\tDominant: $(round(dominant_freq, digits=4)) Hz (T=$(round(true_period, digits=1))s)",
                              legend=false, linewidth=2, color=:blue,
                              left_margin=10Plots.mm, bottom_margin=8Plots.mm)
                        vline!(p_det[1], [dominant_freq], color=:red, linewidth=2, linestyle=:dash)
                        
                        # Right: Single bar showing the true period
                        bar!(p_det[2], ["Period"], [true_period],
                             yerror=[isi_std],  # Use event-based std as uncertainty
                             xlabel="", ylabel="Period (s)",
                             title="Oscillation Period\t$(round(true_period, digits=2)) s (CV=$(round(isi_cv, digits=4)))",
                             legend=false, color=get(color_map, noise_type, :blue),
                             alpha=0.7, bar_width=0.5,
                             left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                             ylims=(0, true_period * 1.5))
                        
                        # Add annotation
                        annotate!(p_det[2], 0.5, true_period * 1.3, 
                                 text("Regular oscillations\t(no stochastic noise)", 10, :center))
                        
                        savefig(p_det, "imgs/bio/isi_analysis/isi_fft_$(noise_type).png")
                        println("Saved: imgs/bio/isi_analysis/isi_fft_$(noise_type).png")
                    end
                end
            end
        end
        # Stochastic case - normal histogram with fits
        k_gamma = (isi_mean / max(isi_std, 1e-10))^2
        theta_gamma = isi_std^2 / max(isi_mean, 1e-10)
        lambda_exp = 1.0 / max(isi_mean, 1e-10)
        
        interpretation = if k_gamma < 1.5
            "Exponential-like (Poisson)"
        elseif k_gamma < 3
            "Gamma (mild regularity)"
        else
            "Regular oscillations"
        end
        
        p_hist = histogram(all_isi, bins=min(25, max(5, length(all_isi)÷3)), normalize=:pdf,
                            xlabel="Inter-Spike Interval (s)",
                            ylabel="Probability Density",
                            title="$noise_type: ISI Distribution (IP3=$representative_ip3, n=$(length(all_isi)))\t" *
                                "CV=$(round(isi_cv, digits=2)), k=$(round(k_gamma, digits=1)) - $interpretation",
                            legend=:topright, fillalpha=0.7, size=(800, 550),
                            color=get(color_map, noise_type, :blue), label="Data",
                            top_margin=10Plots.mm)
        
        x_fit = range(max(0.0, minimum(all_isi) * 0.9), maximum(all_isi) * 1.1, length=100)
        
        if k_gamma > 0 && theta_gamma > 0 && isfinite(k_gamma) && isfinite(theta_gamma) && k_gamma < 1000
            gamma_dist = Gamma(k_gamma, theta_gamma)
            y_gamma = pdf.(gamma_dist, x_fit)
            if all(isfinite.(y_gamma))
                plot!(p_hist, x_fit, y_gamma, linewidth=3, color=:green, 
                        linestyle=:solid, label="Gamma (k=$(round(k_gamma, digits=1)))")
            end
        end
        
        y_exp = lambda_exp .* exp.(-lambda_exp .* x_fit)
        plot!(p_hist, x_fit, y_exp, linewidth=2, color=:red, 
                linestyle=:dash, label="Exponential")
        
        vline!(p_hist, [isi_mean], color=:black, linewidth=2, 
                linestyle=:dot, label="Mean = $(round(isi_mean, digits=1))s")
        
        savefig(p_hist, "imgs/bio/isi_analysis/isi_histogram_$(noise_type).png")
        println("Saved: imgs/bio/isi_analysis/isi_histogram_$(noise_type).png (IP3=$representative_ip3)")

    end
    
    # Combined ISI Comparison (skip deterministic for clarity)
    println("\n--- Combined ISI Comparison ---")
    p_isi_compare = plot(xlabel="ISI (s)", ylabel="Probability Density",
                        title="ISI Distribution Comparison (Stochastic Noise Types)",
                        legend=:topright, size=(900, 600),
                        left_margin=10Plots.mm, bottom_margin=8Plots.mm, top_margin=8Plots.mm)
    
    has_data = false
    for noise_type in noise_list
        if noise_type == :none; continue; end  # Skip deterministic
        if !haskey(raw_data, noise_type); continue; end
        
        all_isi = Float64[]
        for ip3 in [0.7, 0.8, 0.9, 1.0]
            dfs = get(raw_data[noise_type], ip3, DataFrame[])
            for df_sol in dfs
                ratio = calculate_atp_adp_ratio(df_sol, 4000.0)
                if isnothing(ratio); continue; end
                
                events = detect_escape_events(ratio; low_percentile=25.0, high_percentile=75.0)
                if length(events) >= 2
                    isi = calculate_isi(events, 1.0)
                    valid_isi = filter(x -> x > 5.0 && x < 1000.0, isi)
                    append!(all_isi, valid_isi)
                end
            end
        end
        
        if length(all_isi) >= 5
            isi_cv = std(all_isi) / mean(all_isi)
            k_gamma = (mean(all_isi) / std(all_isi))^2
            histogram!(p_isi_compare, all_isi, bins=min(15, max(5, length(all_isi)÷3)), 
                      normalize=:pdf, alpha=0.4,
                      label="$noise_type (CV=$(round(isi_cv, digits=2)), k=$(round(k_gamma, digits=1)))",
                      color=get(color_map, noise_type, :auto))
            has_data = true
        end
    end
    
    if has_data
        savefig(p_isi_compare, "imgs/bio/isi_analysis/isi_comparison.png")
        println("Saved: imgs/bio/isi_analysis/isi_comparison.png")
    else
        println("  No ISI data available for comparison plot")
    end
end


# ============================================================================
# BIFURCATION ANALYSIS (Deterministic Only)
# ============================================================================
"""
Comprehensive bifurcation analysis for deterministic system:
1. Show bifurcation phenomenon (ATP:ADP ratio, frequency, amplitude vs IP3)
2. Compare key variables before, during, and after bifurcation
3. Analyze CV and std of ATP:ADP ratio across ensemble
"""
function plot_bifurcation_analysis(raw_data, results_all, ip3_range)
    println("\n" * "="^80)
    println("BIFURCATION ANALYSIS (Deterministic System)")
    println("="^80)
    
    mkpath("imgs/bio/bifurcation/")
    noise_type = :none
    
    if !haskey(raw_data, noise_type)
        println("⚠️ No deterministic raw data available.")
        return
    end

    # =========================================================================
    # PART 0: PRE-CALCULATE METRICS ON CONSISTENT WINDOW
    # =========================================================================
    # We calculate ALL metrics (Ratio, Amp, Freq, Peaks) on the last 2000s 
    # to ensure every plot tells the same story.
    println("Pre-calculating steady-state metrics (Windowed)...")
    
    bif_data = DataFrame(
        ip3 = Float64[],
        ratio_mean = Float64[],
        ratio_min = Float64[],
        ratio_max = Float64[],
        is_oscillating = Bool[],
        amp = Float64[],
        freq = Float64[],
        peaks = Int[]
    )

    # Use a consistent window for both parts
    analysis_window = 1000.0 
    
    sorted_keys = sort(collect(keys(raw_data[noise_type])))
    
    for ip3 in sorted_keys
        dfs = raw_data[noise_type][ip3]
        if isempty(dfs); continue; end
        
        df_sol = first(dfs) 
        t_max = maximum(df_sol.timestamp)
        
        ratio = calculate_atp_adp_ratio(df_sol, analysis_window)
        
        if isnothing(ratio) || length(ratio) < 10
            continue
        end

        r_mean = mean(ratio)
        r_min = minimum(ratio)
        r_max = maximum(ratio)
        r_amp = r_max - r_min
        
        # 2. Oscillation Detection
        is_osc = r_amp > (0.01 * r_mean) && r_amp > 0.05
        
        # 3. Frequency & Peaks (Simple Robust Counting)
        # Find peaks in the Ratio signal to count oscillations
        n_peaks = 0
        if length(ratio) > 3
            # Simple local maxima detector
            for i in 2:(length(ratio)-1)
                if ratio[i] > ratio[i-1] && ratio[i] > ratio[i+1]
                    n_peaks += 1
                end
            end
        end
        
        # Estimate Frequency (Hz) = (Peaks - 1) / Duration
        # Get the time window duration
        t_start = max(0.0, t_max - analysis_window)
        duration = t_max - t_start
        freq = (n_peaks > 1 && duration > 0) ? (n_peaks - 1) / duration : 0.0
        
        push!(bif_data, (ip3, r_mean, r_min, r_max, is_osc, r_amp, freq, n_peaks))
    end

    # =========================================================================
    # PART 1: Bifurcation Phenomenon (Using Recalculated Data)
    # =========================================================================
    println("\n--- Part 1: Bifurcation Phenomenon (Steady State) ---")
    
    if nrow(bif_data) > 0
        # ----- 1a. ATP:ADP Ratio Bifurcation Diagram -----
        y_max = maximum(bif_data.ratio_max) * 1.1

        p_bif = plot(xlabel="IP3 (μM)", ylabel="ATP:ADP Ratio",
                    title="Bifurcation Diagram (Steady State)",
                    legend=:topright, size=(1000, 600),
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                    ylims=(0, y_max))
        
        # Mean
        plot!(p_bif, bif_data.ip3, bif_data.ratio_mean,
              label="Mean", linewidth=3, color=:blue)
        
        # Min/Max Envelope
        plot!(p_bif, bif_data.ip3, bif_data.ratio_max,
              label="Range", linewidth=1.0, linestyle=:dash, color=:blue, alpha=0.4)
        plot!(p_bif, bif_data.ip3, bif_data.ratio_min,
              linewidth=1.0, linestyle=:dash, color=:blue, alpha=0.4, label="")
        
        # Fill
        plot!(p_bif, bif_data.ip3, bif_data.ratio_mean,
              ribbon=(bif_data.ratio_mean .- bif_data.ratio_min,
                      bif_data.ratio_max .- bif_data.ratio_mean),
              fillalpha=0.2, label="", color=:blue, linewidth=0)
        
        # Markers
        osc_rows = bif_data[bif_data.is_oscillating, :]
        stab_rows = bif_data[.!bif_data.is_oscillating, :]
        
        if nrow(osc_rows) > 0
            scatter!(p_bif, osc_rows.ip3, osc_rows.ratio_mean,
                    marker=:circle, markersize=6, color=:green, label="Oscillating")
        end
        if nrow(stab_rows) > 0
            scatter!(p_bif, stab_rows.ip3, stab_rows.ratio_mean,
                    marker=:square, markersize=6, color=:red, label="Stable")
        end
        
        savefig(p_bif, "imgs/bio/bifurcation/bifurcation_ratio.png")
        println("Saved: imgs/bio/bifurcation/bifurcation_ratio.png")
        
        # ----- 1b. Amplitude vs IP3 -----
        y_max_amp = maximum(bif_data.amp) * 1.2
        p_amp = plot(bif_data.ip3, bif_data.amp,
                    xlabel="IP3 (μM)", ylabel="Amplitude (Max-Min)", title="Oscillation Amplitude",
                    legend=false, size=(900, 500),
                    marker=:circle, linewidth=2, color=:purple,
                    left_margin=10Plots.mm, bottom_margin=8Plots.mm)
        hline!(p_amp, [0.05], linestyle=:dash, color=:gray)
        savefig(p_amp, "imgs/bio/bifurcation/bifurcation_amplitude.png")
        println("Saved: imgs/bio/bifurcation/bifurcation_amplitude.png")

        # ----- 1c. Frequency vs IP3 (RESTORED) -----
        # Filter out zero frequency (stable points) to show trend clearly
        valid_freq = bif_data.freq .> 0
        if sum(valid_freq) > 0
            y_max_freq = maximum(bif_data.freq[valid_freq]) * 1.2            
            p_freq = plot(bif_data.ip3[valid_freq], bif_data.freq[valid_freq],
                        xlabel="IP3 (μM)", ylabel="Frequency (Hz)", title="Oscillation Frequency",
                        legend=false, size=(900, 500),
                        marker=:circle, linewidth=2, color=:orange,
                        left_margin=10Plots.mm, bottom_margin=8Plots.mm, ylims=(0, y_max_freq))
            savefig(p_freq, "imgs/bio/bifurcation/bifurcation_frequency.png")
            println("Saved: imgs/bio/bifurcation/bifurcation_frequency.png")
        end

        # ----- 1d. Peak Count vs IP3 (RESTORED) -----
        y_max_peaks = maximum(bif_data.peaks) * 1.2
        p_peaks = plot(bif_data.ip3, bif_data.peaks,
                      xlabel="IP3 (μM)", ylabel="Peak Count (in 2000s)", title="Peak Count",
                      legend=false, size=(900, 500),
                      marker=:circle, linewidth=2, color=:teal,
                      left_margin=10Plots.mm, bottom_margin=8Plots.mm)
        savefig(p_peaks, "imgs/bio/bifurcation/bifurcation_peaks.png")
        println("Saved: imgs/bio/bifurcation/bifurcation_peaks.png")
        
        # ----- 1e. Summary Panel -----
        y_max_ratio = maximum(bif_data.ratio_max) * 1.1
        p1 = plot(bif_data.ip3, bif_data.ratio_mean, title="Ratio Range", ylabel="Ratio", legend=false, ylims=(0, y_max_ratio))
        plot!(p1, bif_data.ip3, bif_data.ratio_mean, ribbon=(bif_data.ratio_mean .- bif_data.ratio_min, bif_data.ratio_max .- bif_data.ratio_mean), fillalpha=0.3)
        
        p2 = plot(bif_data.ip3, bif_data.amp, title="Amplitude", ylabel="Amp", legend=false)
        
        p3 = plot(bif_data.ip3[valid_freq], bif_data.freq[valid_freq], title="Frequency", ylabel="Hz", legend=false, color=:orange)
        
        p4 = plot(bif_data.ip3, bif_data.peaks, title="Peaks", ylabel="Count", legend=false, color=:teal)
        
        p_summary = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800), margin=5Plots.mm)
        savefig(p_summary, "imgs/bio/bifurcation/bifurcation_summary.png")
        println("Saved: imgs/bio/bifurcation/bifurcation_summary.png")
    end

    # =========================================================================
    # PART 2: Key Variables Comparison (Consistent Window)
    # =========================================================================
    println("\n" * "="^40)
    println("--- Part 2: Key Variables Comparison ---")
    println("="^40)
    
    target_ip3s = [1.2, 1.4, 1.6] 
    
    # Map targets to closest keys
    selected_ip3s = Float64[]
    for t in target_ip3s
        if !isempty(sorted_keys)
            closest = sorted_keys[argmin(abs.(sorted_keys .- t))]
            push!(selected_ip3s, closest)
        end
    end
    unique!(selected_ip3s)
    println("Selected IP3 points: $selected_ip3s")

    # Generate Plots
    key_var_plots = Plots.Plot[]
    
    for (vi, var) in enumerate(important_variables)
        println("Processing $var...")
        
        p = plot(; title=var, xlabel="Time (s)", ylabel=var,
                 legend=(vi==1 ? :topright : false),
                 left_margin=12Plots.mm, bottom_margin=10Plots.mm,
                 yformatter=:scientific) 
        
        has_data = false
        colors = [:blue, :red, :green]
        styles = [:solid, :dash, :dot]
        
        previous_traces = [] 

        for (i, ip3) in enumerate(selected_ip3s)
            if !haskey(raw_data[noise_type], ip3); continue; end
            
            dfs = raw_data[noise_type][ip3]
            if isempty(dfs); continue; end
            df_sol = first(dfs)
            
            col_idx = df_find_column(df_sol, var)
            if isnothing(col_idx); continue; end

            # --- CONSISTENT WINDOWING MATCH ---
            t_max = maximum(df_sol.timestamp)
            t_start = max(0.0, t_max - analysis_window) # Exact same window as Part 1
            mask = df_sol.timestamp .>= t_start
            
            t_data = df_sol.timestamp[mask]
            y_data = df_sol[mask, col_idx]
            
            if length(t_data) < 2; continue; end
            
            t_shifted = t_data .- minimum(t_data)
            
            # Check overlap
            is_identical = false
            for (prev_ip3, prev_y) in previous_traces
                if length(prev_y) == length(y_data) && sum(abs.(prev_y .- y_data)) < 1e-9
                    is_identical = true
                end
            end
            push!(previous_traces, (ip3, y_data))
            
            c = colors[mod1(i, length(colors))]
            s = styles[mod1(i, length(styles))]
            w = is_identical ? 3.0 : 1.5
            
            plot!(p, t_shifted, y_data, 
                  label=(vi==1 ? "IP3=$(round(ip3, digits=2))" : ""), 
                  linewidth=w, alpha=0.8, color=c, linestyle=s)
            has_data = true
        end
        
        if has_data
            push!(key_var_plots, p)
        end
    end
    
    if !isempty(key_var_plots)
        n_plots = length(key_var_plots)
        n_cols = 3
        n_rows = ceil(Int, n_plots / n_cols)
        
        p_combined = plot(key_var_plots..., layout=(n_rows, n_cols), 
                          size=(500*n_cols, 350*n_rows),
                          plot_title="Comparison (Window: Last $(Int(analysis_window))s)",
                          margin=10Plots.mm)
                          
        savefig(p_combined, "imgs/bio/bifurcation/key_variables_comparison.png")
        println("Saved: imgs/bio/bifurcation/key_variables_comparison.png")
    end

    println("\n" * "="^80)
    println("BIFURCATION ANALYSIS COMPLETE")
    println("="^80)
end

function plot_important_variables(results_all)
    println("\n--- Important Variables vs IP3 ---")
    
    legend_positions = Dict(
        "adpc" => :bottomright, "cac" => :bottomright,
        "atpc" => :topright, "caer" => :topright,
        "psi" => :bottomright, "pyrm" => :topright
    )
    
    for var in important_variables
        leg_pos = get(legend_positions, var, :topright)
        # Calculate y_max across all noise types
        y_max = 0.0
        for (nt, df) in results_all
            if var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    y_max = max(y_max, maximum(df[!, var][valid_idx]))
                end
            end
        end
        y_max *= 1.1
        
        p = plot(xlabel="IP3 (μM)", ylabel="$var (steady state)", 
                 title="Steady State $var vs IP3", 
                 legend=leg_pos, size=(900, 600),
                 left_margin=10Plots.mm, bottom_margin=8Plots.mm, 
                 ylims=(var != "psi" ? (0, y_max) : (:auto)))
        
        for (nt, df) in results_all
            if var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df.ip3[valid_idx], df[!, var][valid_idx], 
                          label=string(nt), marker=:circle, markersize=5,
                          linewidth=2, color=get(color_map, nt, :auto))
                end
            end
        end
        
        savefig(p, "imgs/bio/scan/compare_$(var)_ip3.png")
        println("Saved: imgs/bio/scan/compare_$(var)_ip3.png")
    end
    
    # Combined overview
    plots_combined = []
    for (i, var) in enumerate(important_variables)
        # Calculate y_max for this variable
        y_max = 0.0
        for (nt, df) in results_all
            if var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    y_max = max(y_max, maximum(df[!, var][valid_idx]))
                end
            end
        end
        y_max *= 1.1
        p = plot(xlabel="IP3 (μM)", ylabel=var, title=var, 
                legend=(i == 1) ? :topright : false,
                left_margin=8Plots.mm, bottom_margin=6Plots.mm, ylims=(0, y_max))
        
        for (nt, df) in results_all
            if var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df.ip3[valid_idx], df[!, var][valid_idx], 
                          label=string(nt), linewidth=2, marker=:circle, markersize=3,
                          color=get(color_map, nt, :auto))
                end
            end
        end
        push!(plots_combined, p)
    end
    
    n_cols = 3
    n_rows = ceil(Int, length(plots_combined) / n_cols)
    p_all = plot(plots_combined..., layout=(n_rows, n_cols), size=(500*n_cols, 350*n_rows))
    savefig(p_all, "imgs/bio/scan/all_variables_overview.png")
    println("Saved: imgs/bio/scan/all_variables_overview.png")
end

function plot_summary(results_all)
    println("\n--- Summary Plots ---")
    # Calculate y_max for each subplot
    freq_max, fvar_max, rate_max, entropy_max = 0.0, 0.0, 0.0, 0.0
    for (nt, df) in results_all
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            freq_max = max(freq_max, maximum(df.dominant_freq[valid]))
        end
        
        valid = .!isnan.(df.freq_variance)
        if sum(valid) > 0
            fvar_max = max(fvar_max, maximum(df.freq_variance[valid]))
        end
        
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            rate_max = max(rate_max, maximum(df.escape_rate[valid]))
        end
        
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            entropy_max = max(entropy_max, maximum(df.isi_entropy[valid]))
        end
    end

    p1 = plot(xlabel="IP3 (μM)", ylabel="Freq (Hz)", title="Oscillation Frequency", 
              legend=:topleft, size=(400, 350), ylims=(0, freq_max * 1.1))
    p2 = plot(xlabel="IP3 (μM)", ylabel="Freq Var", title="Frequency Variability", 
              legend=false, size=(400, 350), ylims=(0, fvar_max * 1.1))
    p3 = plot(xlabel="IP3 (μM)", ylabel="Rate (1/s)", title="Escape Rate", 
              legend=false, size=(400, 350), ylims=(0, rate_max * 1.1))
    p4 = plot(xlabel="IP3 (μM)", ylabel="Entropy (bits)", title="ISI Entropy", 
              legend=false, size=(400, 350), ylims=(0, entropy_max * 1.1))

    for (nt, df) in results_all
        c = get(color_map, nt, :auto)
        
        valid = .!isnan.(df.dominant_freq) .& (df.dominant_freq .> 0)
        if sum(valid) > 0
            plot!(p1, df.ip3[valid], df.dominant_freq[valid], 
                  label=string(nt), color=c, linewidth=2, marker=:circle, markersize=3)
        end
        
        valid = .!isnan.(df.freq_variance)
        if sum(valid) > 0
            plot!(p2, df.ip3[valid], df.freq_variance[valid], color=c, linewidth=2, 
                  marker=:circle, markersize=3)
        end
        
        valid = .!isnan.(df.escape_rate) .& (df.escape_rate .> 0)
        if sum(valid) > 0
            plot!(p3, df.ip3[valid], df.escape_rate[valid], color=c, linewidth=2,
                  marker=:circle, markersize=3)
        end
        
        valid = .!isnan.(df.isi_entropy)
        if sum(valid) > 0
            plot!(p4, df.ip3[valid], df.isi_entropy[valid], color=c, linewidth=2,
                  marker=:circle, markersize=3)
        end
    end

    p_summary = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
    savefig(p_summary, "imgs/bio/scan/kramers_summary.png")
    println("Saved: imgs/bio/scan/kramers_summary.png")
end

# ============================================================================
# EXPERIMENT 1: Noise Strength vs Mean ATP:ADP Ratio (Kept separate)
# ============================================================================
function scan_noise_strength_experiment()
    println("\n" * "="^80)
    println("EXPERIMENT 1: Noise Strength vs Mean ATP:ADP Ratio")
    println("="^80)
    
    noise_configs = Dict(
        :additive => 0.01:0.01:0.07, 
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
            print("  σ = $σ : ")
            
            try
                if noise_type == :jump
                    sol, _ = simulate_model(noise_type; σ_jump=σ, tspan=(0.0, 4000.0))
                elseif noise_type == :additive
                    sol, _ = simulate_model(noise_type; σ_additive=σ, tspan=(0.0, 4000.0))
                elseif noise_type == :multiplicative
                    sol, _ = simulate_model(noise_type; σ_multiplicative=σ, tspan=(0.0, 4000.0))
                elseif noise_type == :state_dependent
                    sol, _ = simulate_model(noise_type; σ_calcium=σ, tspan=(0.0, 4000.0))
                else
                    continue
                end
                
                df_sol = DataFrame(sol)
                ratio = calculate_atp_adp_ratio(df_sol, 1000.0)
                
                if !isnothing(ratio) && length(ratio) > 10
                    push!(df_res, (σ, mean(ratio), std(ratio)))
                    println("✓")
                else
                    push!(df_res, (σ, NaN, NaN))
                    println("⚠️ Invalid")
                end
            catch e
                push!(df_res, (σ, NaN, NaN))
                println("❌ Error")
            end
        end
        
        results[noise_type] = df_res
    end
    
    # Plot
    mkpath("imgs/bio/noise_analysis/")
    
    # =========================================================================
    # SEPARATE PLOTS per noise type
    # =========================================================================
    for (nt, df) in results
        valid_idx = .!isnan.(df.mean_ratio)
        y_max = maximum(df.mean_ratio[valid_idx] .+ coalesce.(df.std_ratio[valid_idx], 0.0)) * 1.1            
        if sum(valid_idx) > 0
            p = plot(xlabel="Noise Strength (σ)", ylabel="Mean ATP:ADP Ratio",
                     title="$nt: Noise Strength vs ATP:ADP Ratio",
                     legend=false, size=(700, 500),
                     left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                     ylims=(0, y_max))
            
            plot!(p, df.noise_strength[valid_idx], df.mean_ratio[valid_idx],
                  marker=:circle, linewidth=2, markersize=6,
                  ribbon=df.std_ratio[valid_idx], fillalpha=0.3,
                  color=get(color_map, nt, :blue))
            
            savefig(p, "imgs/bio/noise_analysis/noise_strength_vs_ratio_$(nt).png")
            println("Saved: imgs/bio/noise_analysis/noise_strength_vs_ratio_$(nt).png")
        end
    end
    
    # =========================================================================
    # COMBINED PLOT
    # =========================================================================
    y_max = 0.0
    for (nt, df) in results
        valid_idx = .!isnan.(df.mean_ratio)
        if sum(valid_idx) > 0
            y_max = max(y_max, maximum(df.mean_ratio[valid_idx] .+ coalesce.(df.std_ratio[valid_idx], 0.0)))
        end
    end
    y_max *= 1.1
    p_combined = plot(xlabel="Noise Strength (σ)", ylabel="Mean ATP:ADP Ratio",
                      title="Noise Strength vs ATP:ADP Ratio",
                      legend=:topright, size=(900, 600),
                      left_margin=10Plots.mm, bottom_margin=8Plots.mm,
                      ylims=(0, y_max))
    
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
    
    return results
end


function print_data_statistics(raw_data::Dict{Symbol, Dict{Float64, Vector{DataFrame}}})
    println("\n" * "="^130)
    println("DATA STATISTICS REPORT (Timestamps & Biological Validity Check)")
    println("="^130)
    
    analysis_window = 4000.0  # Same window as metric extraction
    
    # Table Header
    header_fmt = "%-8s | %-12s | %-5s | %-10s | %-10s | %-10s | %-8s | %-8s | %-12s"
    row_fmt    = "%-8.1f | %-12s | %-5d | %-10.1f | %-10.1f | %-10.1f | %-8d | %-8d | %-12s"
    
    println(Printf.format(Printf.Format(header_fmt), 
        "IP3", "Noise Type", "Count", "T(end-2)", "T(end-1)", "T(end)", "NegRatio", "NegConc", "Status"))
    println("-"^130)

    # Track problematic combinations for summary
    problematic_cases = []
    
    # 1. Get sorted keys for consistent display
    all_ip3s = Set{Float64}()
    for (nt, subdict) in raw_data
        union!(all_ip3s, keys(subdict))
    end
    sorted_ip3s = sort(collect(all_ip3s))
    sorted_noise = sort(collect(keys(raw_data)))

    # 2. Iterate and inspect
    for ip3 in sorted_ip3s
        for nt in sorted_noise
            
            if haskey(raw_data, nt) && haskey(raw_data[nt], ip3)
                dfs = raw_data[nt][ip3]
                count = length(dfs)
                
                if count > 0
                    # INSPECT THE FIRST RUN AS A REPRESENTATIVE SAMPLE for timestamps
                    df_sample = first(dfs)
                    times = df_sample.timestamp
                    n_points = length(times)
                    
                    # Safe extraction of last 3 points
                    t_end   = n_points > 0 ? times[end] : NaN
                    t_last2 = n_points > 1 ? times[end-1] : NaN
                    t_last3 = n_points > 2 ? times[end-2] : NaN
                    
                    # =========================================================
                    # NEW: Check ALL runs for negative values in analysis window
                    # =========================================================
                    n_negative_ratio_runs = 0
                    n_negative_conc_runs = 0
                    
                    for (run_idx, df) in enumerate(dfs)
                        # Apply same window as metric extraction
                        t_max = maximum(df.timestamp)
                        t_start = max(0.0, t_max - analysis_window)
                        mask = df.timestamp .>= t_start
                        
                        if sum(mask) < 10
                            continue
                        end
                        
                        # Get ATPC and ADPC columns
                        atpc_idx = df_find_column(df, "atpc")
                        adpc_idx = df_find_column(df, "adpc")
                        
                        if !isnothing(atpc_idx) && !isnothing(adpc_idx)
                            atpc_windowed = df[mask, atpc_idx]
                            adpc_windowed = df[mask, adpc_idx]
                            
                            # Check for negative concentrations (biological violation)
                            has_neg_atpc = any(x -> x < 0, atpc_windowed)
                            has_neg_adpc = any(x -> x < 0, adpc_windowed)
                            
                            if has_neg_atpc || has_neg_adpc
                                n_negative_conc_runs += 1
                            end
                            
                            # Calculate ratio and check for negative/invalid values
                            # Negative ratio can occur if one concentration goes negative
                            valid_mask = (adpc_windowed .> 0) .& isfinite.(atpc_windowed) .& isfinite.(adpc_windowed)
                            
                            if sum(valid_mask) > 0
                                ratio = atpc_windowed[valid_mask] ./ adpc_windowed[valid_mask]
                                
                                # Check for negative ratios (happens when atpc < 0)
                                if any(x -> x < 0, ratio)
                                    n_negative_ratio_runs += 1
                                end
                            else
                                # All points invalid - count as problematic
                                n_negative_ratio_runs += 1
                            end
                        end
                    end
                    
                    # Determine status
                    status = "OK"
                    if n_points < 3
                        status = "FEW PTS"
                    elseif (t_end - t_last2) > 500.0
                        status = "HUGE JUMP"
                    elseif t_end < 100.0
                        status = "EARLY END"
                    elseif n_negative_conc_runs > 0
                        status = "NEG CONC!"
                    elseif n_negative_ratio_runs > 0
                        status = "NEG RATIO!"
                    end
                    
                    # Track problematic cases
                    if n_negative_conc_runs > 0 || n_negative_ratio_runs > 0
                        push!(problematic_cases, (ip3, nt, count, n_negative_ratio_runs, n_negative_conc_runs))
                    end
                    
                    println(Printf.format(Printf.Format(row_fmt), 
                        ip3, string(nt), count, t_last3, t_last2, t_end, 
                        n_negative_ratio_runs, n_negative_conc_runs, status))
                else
                    println(Printf.format(Printf.Format(row_fmt), 
                        ip3, string(nt), 0, NaN, NaN, NaN, 0, 0, "EMPTY"))
                end
            else
                println(Printf.format(Printf.Format(row_fmt), 
                    ip3, string(nt), 0, NaN, NaN, NaN, 0, 0, "MISSING"))
            end
        end
        # Separator between IP3 groups
        if !isempty(sorted_noise)
            println("-"^130)
        end
    end
    
    # =========================================================================
    # SUMMARY: Problematic cases
    # =========================================================================
    if !isempty(problematic_cases)
        println("\n" * "="^80)
        println("⚠️  BIOLOGICAL VALIDITY WARNINGS (Last $(Int(analysis_window))s window)")
        println("="^80)
        println("\nThe following combinations have runs with negative concentrations or ratios:")
        println("This indicates numerical instability violating biological constraints.\n")
        
        println(rpad("IP3", 8), rpad("Noise Type", 15), rpad("Total Runs", 12), 
                rpad("Neg Ratio", 12), rpad("Neg Conc", 12), "Severity")
        println("-"^70)
        
        for (ip3, nt, total, neg_ratio, neg_conc) in problematic_cases
            severity = if neg_conc > total/2
                "CRITICAL"
            elseif neg_conc > 0
                "HIGH"
            elseif neg_ratio > total/2
                "MODERATE"
            else
                "LOW"
            end
            
            println(rpad("$ip3", 8), rpad(string(nt), 15), rpad("$total", 12),
                    rpad("$neg_ratio", 12), rpad("$neg_conc", 12), severity)
        end
        
        println("\n" * "-"^70)
        println("Legend:")
        println("  NegRatio = Runs where ATP:ADP ratio contains negative values")
        println("  NegConc  = Runs where ATP or ADP concentration went negative")
        println("  CRITICAL = >50% runs have negative concentrations")
        println("  HIGH     = Some runs have negative concentrations")
        println("  MODERATE = >50% runs have negative ratios (but concentrations OK)")
        println("  LOW      = Few runs affected")
        println("\nRecommendation: Exclude CRITICAL/HIGH cases from analysis or")
        println("                increase simulation stability (reduce noise, smaller dt)")
    else
        println("\n✓ All runs passed biological validity check (no negative concentrations)")
    end
    
    println("="^130 * "\n")
    
    return problematic_cases
end

"""
Validate results_all and raw_data for biological validity.
Checks:
1. Sufficient data points in analysis window
2. No negative ATP:ADP ratios (biological impossibility)
3. No negative concentrations (ATPC, ADPC, CaC, CaER)
4. Reasonable value ranges (not exploded)

Returns a detailed report and flags problematic combinations.
"""
function validate_biological_consistency(results_all::Dict{Symbol, DataFrame}, 
                                         raw_data::Dict{Symbol, Dict{Float64, Vector{DataFrame}}};
                                         analysis_window::Float64 = 4000.0,
                                         verbose::Bool = true)
    println("\n" * "="^100)
    println("BIOLOGICAL VALIDITY CHECK")
    println("="^100)
    
    validation_report = Dict{Symbol, DataFrame}()
    problematic_combinations = []
    
    for noise_type in noise_list
        if !haskey(raw_data, noise_type)
            continue
        end
        
        if verbose
            println("\n--- Checking $noise_type ---")
        end
        
        df_report = DataFrame(
            ip3 = Float64[],
            n_runs = Int[],
            n_valid_runs = Int[],
            n_negative_ratio = Int[],
            n_negative_atpc = Int[],
            n_negative_adpc = Int[],
            n_insufficient_points = Int[],
            n_exploded = Int[],
            mean_ratio = Float64[],
            std_ratio = Float64[],
            min_ratio = Float64[],
            max_ratio = Float64[],
            status = String[]
        )
        
        sorted_ip3 = sort(collect(keys(raw_data[noise_type])))
        
        for ip3 in sorted_ip3
            dfs = raw_data[noise_type][ip3]
            n_runs = length(dfs)
            
            if n_runs == 0
                push!(df_report, (ip3, 0, 0, 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN, "NO DATA"))
                continue
            end
            
            # Counters
            n_valid = 0
            n_neg_ratio = 0
            n_neg_atpc = 0
            n_neg_adpc = 0
            n_insufficient = 0
            n_exploded = 0
            
            all_ratios = Float64[]
            
            for (run_idx, df) in enumerate(dfs)
                # Apply analysis window
                t_max = maximum(df.timestamp)
                t_start = max(0.0, t_max - analysis_window)
                mask = df.timestamp .>= t_start
                n_points = sum(mask)
                
                # Check 1: Sufficient points
                if n_points < 100
                    n_insufficient += 1
                    continue
                end
                
                # Get columns
                atpc_idx = df_find_column(df, "atpc")
                adpc_idx = df_find_column(df, "adpc")
                
                if isnothing(atpc_idx) || isnothing(adpc_idx)
                    n_insufficient += 1
                    continue
                end
                
                atpc = df[mask, atpc_idx]
                adpc = df[mask, adpc_idx]
                
                # Check 2: Negative concentrations
                has_neg_atpc = any(x -> x < -1e-10, atpc)
                has_neg_adpc = any(x -> x < -1e-10, adpc)
                
                if has_neg_atpc
                    n_neg_atpc += 1
                end
                if has_neg_adpc
                    n_neg_adpc += 1
                end
                
                # Check 3: Exploded values (unreasonably large)
                max_reasonable = 100.0  # Concentrations shouldn't exceed ~100 in this model
                has_exploded = maximum(abs.(atpc)) > max_reasonable || maximum(abs.(adpc)) > max_reasonable
                if has_exploded
                    n_exploded += 1
                end
                
                # Calculate ratio (only for positive values)
                valid_mask = (atpc .> 0) .& (adpc .> 0) .& isfinite.(atpc) .& isfinite.(adpc)
                
                if sum(valid_mask) < 10
                    n_insufficient += 1
                    continue
                end
                
                ratio = atpc[valid_mask] ./ adpc[valid_mask]
                
                # Check 4: Negative ratios (shouldn't happen if both positive, but check anyway)
                if any(x -> x < 0, ratio)
                    n_neg_ratio += 1
                end
                
                # Only count as valid if no issues
                if !has_neg_atpc && !has_neg_adpc && !has_exploded && sum(valid_mask) >= n_points * 0.9
                    n_valid += 1
                    append!(all_ratios, [mean(ratio)])  # Store mean ratio for this run
                end
            end
            
            # Calculate aggregate statistics from valid runs only
            if length(all_ratios) > 0
                mean_r = mean(all_ratios)
                std_r = length(all_ratios) > 1 ? std(all_ratios) : 0.0
                min_r = minimum(all_ratios)
                max_r = maximum(all_ratios)
            else
                mean_r, std_r, min_r, max_r = NaN, NaN, NaN, NaN
            end
            
            # Determine status
            status = if n_valid == n_runs
                "✓ OK"
            elseif n_valid >= n_runs * 0.7
                "⚠ PARTIAL"
            elseif n_valid > 0
                "⚠ UNSTABLE"
            else
                "✗ INVALID"
            end
            
            # Track problematic cases
            if n_neg_atpc > 0 || n_neg_adpc > 0 || n_neg_ratio > 0 || n_valid < n_runs * 0.5
                push!(problematic_combinations, (noise_type, ip3, n_valid, n_runs, status))
            end
            
            push!(df_report, (ip3, n_runs, n_valid, n_neg_ratio, n_neg_atpc, n_neg_adpc, 
                             n_insufficient, n_exploded, mean_r, std_r, min_r, max_r, status))
        end
        
        validation_report[noise_type] = df_report
        
        # Print summary for this noise type
        if verbose
            n_total_ip3 = nrow(df_report)
            n_ok = sum(df_report.status .== "✓ OK")
            n_partial = sum(startswith.(df_report.status, "⚠"))
            n_invalid = sum(df_report.status .== "✗ INVALID")
            
            println("  Total IP3 values: $n_total_ip3")
            println("  ✓ OK: $n_ok, ⚠ Partial/Unstable: $n_partial, ✗ Invalid: $n_invalid")
            
            # Show problematic IP3 values
            prob_rows = df_report[df_report.status .!= "✓ OK", :]
            if nrow(prob_rows) > 0
                println("\n  Problematic IP3 values:")
                for row in eachrow(prob_rows)
                    println("    IP3=$(row.ip3): $(row.n_valid_runs)/$(row.n_runs) valid, " *
                           "neg_atpc=$(row.n_negative_atpc), neg_adpc=$(row.n_negative_adpc), " *
                           "exploded=$(row.n_exploded) → $(row.status)")
                end
            end
        end
    end
    
    # =========================================================================
    # SUMMARY TABLE
    # =========================================================================
    println("\n" * "="^100)
    println("SUMMARY: RECOMMENDED IP3 RANGES FOR EACH NOISE TYPE")
    println("="^100)
    
    recommended_ranges = Dict{Symbol, Tuple{Float64, Float64}}()
    
    for noise_type in noise_list
        if !haskey(validation_report, noise_type)
            continue
        end
        
        df = validation_report[noise_type]
        valid_ip3 = df[df.status .== "✓ OK", :ip3]
        
        if length(valid_ip3) > 0
            ip3_min = minimum(valid_ip3)
            ip3_max = maximum(valid_ip3)
            recommended_ranges[noise_type] = (ip3_min, ip3_max)
            println("  $noise_type: IP3 ∈ [$ip3_min, $ip3_max] ($(length(valid_ip3)) valid points)")
        else
            # Fall back to partial validity
            partial_ip3 = df[startswith.(df.status, "⚠"), :ip3]
            if length(partial_ip3) > 0
                ip3_min = minimum(partial_ip3)
                ip3_max = maximum(partial_ip3)
                recommended_ranges[noise_type] = (ip3_min, ip3_max)
                println("  $noise_type: IP3 ∈ [$ip3_min, $ip3_max] ($(length(partial_ip3)) partial, use with caution)")
            else
                println("  $noise_type: ⚠️ NO VALID IP3 RANGE - all simulations unstable!")
            end
        end
    end
    
    # =========================================================================
    # COMPARISON WITH results_all
    # =========================================================================
    println("\n" * "="^100)
    println("CHECKING results_all FOR INVALID VALUES")
    println("="^100)
    
    for noise_type in noise_list
        if !haskey(results_all, noise_type)
            continue
        end
        
        df = results_all[noise_type]
        
        # Check for negative mean ratios
        neg_ratio_rows = df[.!isnan.(df.atp_adp_ratio_mean) .& (df.atp_adp_ratio_mean .< 0), :]
        
        # Check for huge standard deviations (sign of instability)
        huge_std_rows = df[.!isnan.(df.atp_adp_ratio_std) .& (df.atp_adp_ratio_std .> 5.0), :]
        
        # Check for unreasonable mean ratios
        unreasonable_rows = df[.!isnan.(df.atp_adp_ratio_mean) .& (df.atp_adp_ratio_mean .> 20.0), :]
        
        if nrow(neg_ratio_rows) > 0 || nrow(huge_std_rows) > 0 || nrow(unreasonable_rows) > 0
            println("\n⚠️  $noise_type has invalid values in results_all:")
            
            if nrow(neg_ratio_rows) > 0
                println("  NEGATIVE RATIOS at IP3: $(neg_ratio_rows.ip3)")
            end
            
            if nrow(huge_std_rows) > 0
                println("  HUGE STD (>5) at IP3: $(huge_std_rows.ip3)")
                for row in eachrow(huge_std_rows)
                    println("    IP3=$(row.ip3): mean=$(round(row.atp_adp_ratio_mean, digits=2)) ± $(round(row.atp_adp_ratio_std, digits=2))")
                end
            end
            
            if nrow(unreasonable_rows) > 0
                println("  UNREASONABLE MEAN (>20) at IP3: $(unreasonable_rows.ip3)")
            end
        else
            println("  $noise_type: ✓ All values in results_all appear reasonable")
        end
    end
    
    println("\n" * "="^100)
    
    return validation_report, problematic_combinations, recommended_ranges
end

"""
Recalculate all metrics from existing raw_data without re-running simulations.
Use this when you've fixed a bug in the metric calculation.
"""
function recalculate_results_from_raw_data(raw_data::Dict{Symbol, Dict{Float64, Vector{DataFrame}}}, 
                                           ip3_range)
    println("\n" * "="^80)
    println("RECALCULATING RESULTS FROM RAW DATA")
    println("="^80)
    
    results_all = Dict{Symbol, DataFrame}()
    
    for noise_type in noise_list
        if !haskey(raw_data, noise_type)
            println("⚠️ No data for $noise_type")
            continue
        end
        
        println("\nRecalculating metrics for $noise_type...")
        
        df_res = DataFrame(
            ip3 = Float64[],
            atp_adp_ratio_mean = Float64[],
            atp_adp_ratio_std = Float64[],
            atp_adp_ratio_sem = Float64[],
            atp_adp_ratio_min = Float64[],
            atp_adp_ratio_max = Float64[],
            dominant_freq = Float64[],
            freq_variance = Float64[],
            isi_entropy = Float64[],
            isi_entropy_std = Float64[],
            isi_mean = Float64[],
            isi_std = Float64[],
            isi_cv = Float64[],
            escape_rate = Float64[],
            escape_rate_std = Float64[],
            n_peaks = Int[],
            n_events = Int[],
            oscillation_amplitude = Float64[],
            is_oscillating = Bool[]
        )
        
        # Add columns for important variables
        for var in important_variables
            df_res[!, var] = Float64[]
            df_res[!, "$(var)_std"] = Float64[]
        end
        
        for val in ip3_range
            dfs = get(raw_data[noise_type], val, DataFrame[])
            
            if isempty(dfs)
                row = create_nan_row(val)
                push!(df_res, row)
                continue
            end
            
            # Recalculate metrics
            row = extract_all_metrics(dfs, val)
            push!(df_res, row)
        end
        
        results_all[noise_type] = df_res
    end
    
    println("\n✓ Recalculation complete!")
    return results_all
end

# ============================================================================
# RUN ALL EXPERIMENTS (CONSOLIDATED)
# ============================================================================

println("\n" * "="^80)
println("RUNNING CONSOLIDATED ANALYSIS")
println("="^80)

# Task 1: Noise strength scan (separate, different parameter)
if !isfile("imgs/bio/noise_analysis/noise_strength_vs_ratio_combined.png")
    scan_noise_strength_experiment()
end

# Tasks 2-5 + Bifurcation: Single comprehensive run
stable_ip3_ranges = Dict(
    :none => 0.1:0.1:2.0,           # Deterministic is usually stable
    :additive => 0.1:0.1:2.0,       # Usually stable
    :multiplicative => 0.1:0.1:2.0, # Usually stable  
    :state_dependent => 0.3:0.1:2.0, # Skip 0.1, 0.2 (known unstable)
    :jump => 0.6:0.1:2.0            # Skip 0.1-0.5 (known unstable)
)
n_ensemble = 10
ip3_range = 0.1:0.1:2.0
if isfile("results/advanced_analysis.jld2")
    println("\nLoading previous comprehensive analysis results...")
    @load "results/advanced_analysis.jld2" results_all raw_data stability_report
    
    validation_report, problematic_combinations, recommended_ranges = 
        validate_biological_consistency(results_all, raw_data)
    
    # RECALCULATE results from raw_data to fix any metric bugs
    # println("\nRecalculating metrics from raw data...")
    # results_all = recalculate_results_from_raw_data(raw_data, ip3_range)
    
    # Save updated results
    # @save "results/advanced_analysis.jld2" results_all raw_data stability_report
    # println("Saved recalculated results to results/advanced_analysis.jld2")
    
    verify_atp_adp_ratio_computation()
else
    println("\nRunning comprehensive analysis (this may take a while)...")
    results_all, raw_data, stability_report = run_comprehensive_analysis(
        ip3_range = ip3_range,  
        n_ensemble = n_ensemble,
        tspan_total = 8000.0,
    )
end

# Print stability summary
# print_stability_summary(stability_report, 0.3)
print_data_statistics(raw_data)

generate_all_plots(results_all, raw_data, ip3_range, n_ensemble)

println("\n" * "="^80)
println("ALL EXPERIMENTS COMPLETE!")
println("="^80)
println("\nResults saved in:")
println("  - imgs/bio/scan/")
println("  - imgs/bio/isi_analysis/")
println("  - imgs/bio/noise_analysis/")
println("  - imgs/bio/bifurcation/")