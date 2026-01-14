using FFTW
using StatsBase

"""
Verification function to check ATP:ADP ratio calculation
Prints diagnostic information to help debug the ratio calculation
"""
function verify_atp_adp_ratio_computation()
    println("\n" * "="^80)
    println("VERIFICATION: ATP:ADP RATIO CALCULATION")
    println("="^80)
    
    # Check a few IP3 values across different noise types
    test_cases = [
        (:none, 0.3),
        (:none, 1.0),
        (:none, 1.4),
        (:none, 1.8),
        (:additive, 1.4),
        (:multiplicative, 1.4),
        (:jump, 1.4),
        (:state_dependent, 1.4),
        (:state_dependent, 0.3)
    ]
    
    for (nt, ip3) in test_cases
        if !haskey(raw_data, nt) || !haskey(raw_data[nt], ip3)
            println("⚠️  $nt, IP3=$ip3: No data")
            continue
        end
        
        dfs = raw_data[nt][ip3]
        if isempty(dfs)
            println("⚠️  $nt, IP3=$ip3: Empty dataframe list")
            continue
        end
        
        # FIX: Calculate mean across ALL runs, not just first run
        all_run_ratios = Float64[]
        
        for (run_idx, df) in enumerate(dfs)
            atpc_idx = df_find_column(df, "atpc")
            adpc_idx = df_find_column(df, "adpc")
            
            if isnothing(atpc_idx) || isnothing(adpc_idx)
                continue
            end
            
            # Apply 4000s window
            window = 4000.0
            t_max = maximum(df.timestamp)
            t_start = max(0.0, t_max - window)
            mask = df.timestamp .>= t_start
            
            atpc_windowed = df[mask, atpc_idx]
            adpc_windowed = df[mask, adpc_idx]
            
            # Mean of ratios for this run
            valid_mask = (atpc_windowed .> 0) .& (adpc_windowed .> 0) .& 
                         isfinite.(atpc_windowed) .& isfinite.(adpc_windowed)
            if sum(valid_mask) > 10
                ratio_vector = atpc_windowed[valid_mask] ./ adpc_windowed[valid_mask]
                push!(all_run_ratios, mean(ratio_vector))
            end
        end
        
        # First run details (for display)
        df = first(dfs)
        atpc_idx = df_find_column(df, "atpc")
        adpc_idx = df_find_column(df, "adpc")
        
        t_max = maximum(df.timestamp)
        n_total = nrow(df)
        window = 4000.0
        t_start = max(0.0, t_max - window)
        mask = df.timestamp .>= t_start
        n_windowed = sum(mask)
        
        atpc_windowed = df[mask, atpc_idx]
        adpc_windowed = df[mask, adpc_idx]
        
        atpc_mean = mean(atpc_windowed)
        adpc_mean = mean(adpc_windowed)
        ratio_of_means = atpc_mean / adpc_mean
        
        valid_mask = (atpc_windowed .> 0) .& (adpc_windowed .> 0) .& 
                     isfinite.(atpc_windowed) .& isfinite.(adpc_windowed)
        ratio_vector = atpc_windowed[valid_mask] ./ adpc_windowed[valid_mask]
        mean_of_ratios_run1 = mean(ratio_vector)
        
        # Ensemble mean
        ensemble_mean = length(all_run_ratios) > 0 ? mean(all_run_ratios) : NaN
        ensemble_std = length(all_run_ratios) > 1 ? std(all_run_ratios) : NaN
        
        # Get stored value
        if haskey(results_all, nt)
            res_df = results_all[nt]
            res_row = res_df[res_df.ip3 .== ip3, :]
            if nrow(res_row) > 0
                stored_ratio = res_row[1, :atp_adp_ratio_mean]
            else
                stored_ratio = NaN
            end
        else
            stored_ratio = NaN
        end
        
        println("\n$nt, IP3=$ip3:")
        println("  Total time: $t_max s, Total points: $n_total")
        println("  Window: last $window s, Points in window: $n_windowed")
        println("  ATPC (windowed, run1): mean=$(round(atpc_mean, digits=3)), range=$(round(minimum(atpc_windowed), digits=3))-$(round(maximum(atpc_windowed), digits=3))")
        println("  ADPC (windowed, run1): mean=$(round(adpc_mean, digits=3)), range=$(round(minimum(adpc_windowed), digits=3))-$(round(maximum(adpc_windowed), digits=3))")
        println("  Ratio of means (WRONG): $(round(ratio_of_means, digits=3))")
        println("  Mean of ratios (run1): $(round(mean_of_ratios_run1, digits=3))")
        println("  Ensemble mean ($(length(all_run_ratios)) runs): $(round(ensemble_mean, digits=3)) ± $(isnan(ensemble_std) ? "N/A" : round(ensemble_std, digits=3))")
        println("  Stored ratio: $(isnan(stored_ratio) ? "NaN" : round(stored_ratio, digits=3))")
        
        # Compare ensemble mean to stored (not single run)
        if !isnan(stored_ratio) && !isnan(ensemble_mean) && abs(ensemble_mean - stored_ratio) > 0.01
            println("  ⚠️  MISMATCH! Difference: $(round(abs(ensemble_mean - stored_ratio), digits=3))")
        else
            println("  ✓ Values match (ensemble mean ≈ stored)")
        end
    end
    
    println("\n" * "="^80)
end

# ============================================================================
# ANALYSIS HELPER FUNCTIONS
# ============================================================================
"""
Helper: find dataframe column index by matching name substring (case-insensitive)
"""
function df_find_column(df, varname)
    for (i, nm) in enumerate(names(df))
        if occursin(varname, lowercase(string(nm))) 
            return i
        end
    end
    return nothing
end

"""
Extract last N SECONDS (not points!) from solution DataFrame for a given variable
Uses time-based windowing for consistency across different sampling rates.
"""
function get_last_n_seconds(df::DataFrame, var_name::String, window_seconds::Float64=4000.0)
    col_idx = df_find_column(df, var_name)
    if isnothing(col_idx)
        return nothing
    end
    
    # Time-based filtering
    t_max = maximum(df.timestamp)
    t_start = max(0.0, t_max - window_seconds)
    mask = df.timestamp .>= t_start
    
    if sum(mask) < 10
        return nothing
    end
    
    return df[mask, col_idx]
end

"""
DEPRECATED: Use get_last_n_seconds instead for consistency
"""
function get_last_n_points(df::DataFrame, var_name::String, n_points::Int=1000)
    col_idx = df_find_column(df, var_name)
    if isnothing(col_idx)
        return nothing
    end
    
    data = df[!, col_idx]
    n_available = length(data)
    start_idx = max(1, n_available - n_points + 1)
    return data[start_idx:end]
end

"""
Calculate ATP:ADP ratio from solution DataFrame using TIME-BASED windowing
Default: last 4000 seconds of an 8000s simulation
"""
function calculate_atp_adp_ratio(df::DataFrame, window_seconds::Float64=4000.0)
    atpc_data = get_last_n_seconds(df, "atpc", window_seconds)
    adpc_data = get_last_n_seconds(df, "adpc", window_seconds)
    
    if isnothing(atpc_data) || isnothing(adpc_data)
        return nothing
    end
    
    # Filter valid values
    valid_mask = (atpc_data .> 0) .& (adpc_data .> 0) .& 
                 isfinite.(atpc_data) .& isfinite.(adpc_data)
    
    if sum(valid_mask) < 10
        return nothing
    end
    
    ratio = atpc_data[valid_mask] ./ adpc_data[valid_mask]
    return ratio
end

"""
Detect complete oscillation cycles (escape events) using dual threshold crossing.
An escape event is: 
  1. Start below low_threshold (25th percentile)
  2. Cross above high_threshold (75th percentile) 
  3. Return below low_threshold (complete cycle)
Returns indices of cycle completions (when signal returns to baseline).
"""
function detect_escape_events(data::Vector{Float64}; 
                              low_percentile::Float64=25.0,
                              high_percentile::Float64=75.0)
    n = length(data)
    if n < 10
        return Int[]
    end
    
    # Dual thresholds
    low_threshold = quantile(data, low_percentile / 100.0)
    high_threshold = quantile(data, high_percentile / 100.0)
    
    escape_times = Int[]
    
    # State machine: 0 = at baseline, 1 = rising (crossed low), 2 = at peak (crossed high)
    state = 0
    
    for i in 2:n
        if state == 0
            # At baseline, waiting to cross low threshold upward
            if data[i] > low_threshold && data[i-1] <= low_threshold
                state = 1  # Started rising
            end
        elseif state == 1
            # Rising, waiting to cross high threshold
            if data[i] > high_threshold && data[i-1] <= high_threshold
                state = 2  # Reached peak region
            elseif data[i] < low_threshold && data[i-1] >= low_threshold
                # Fell back before reaching peak - reset
                state = 0
            end
        elseif state == 2
            # At peak, waiting to return to baseline
            if data[i] < low_threshold && data[i-1] >= low_threshold
                # Completed full cycle!
                push!(escape_times, i)
                state = 0  # Ready for next cycle
            end
        end
    end
    
    return escape_times
end

"""
Detect peaks/spikes in a time series using local maxima
"""
function detect_peaks(data::Vector{Float64}; threshold_percentile::Float64=75.0)
    n = length(data)
    if n < 3
        return Int[]
    end
    
    # Calculate threshold
    threshold = quantile(data, threshold_percentile/100.0)
    
    peaks = Int[]
    for i in 2:(n-1)
        if data[i] > data[i-1] && data[i] > data[i+1] && data[i] > threshold
            push!(peaks, i)
        end
    end
    
    return peaks
end


"""
Calculate Inter-Spike Intervals
"""
function calculate_isi(peaks::Vector{Int}, dt::Float64=1.0)
    if length(peaks) < 2; return Float64[]; end
    return Float64.(diff(peaks)) .* dt
end

"""
Calculate Shannon entropy of ISI distribution
For Kramers escape: 
- Low entropy = regular periodic oscillations (deterministic limit)
- High entropy = random escape times (strong noise regime)
"""
function calculate_entropy(data::Vector{Float64}; n_bins::Int=20)
    if length(data) < 2; return NaN; end
    
    data_range = maximum(data) - minimum(data)
    if data_range ≈ 0; return 0.0; end
    
    counts = zeros(Int, n_bins)
    for val in data
        bin_idx = min(n_bins, max(1, ceil(Int, (val - minimum(data)) / data_range * n_bins)))
        counts[bin_idx] += 1
    end
    
    probs = counts ./ sum(counts)
    
    entropy = 0.0
    for p in probs
        if p > 0
            entropy -= p * log2(p)
        end
    end
    return entropy
end


"""
Estimate oscillation frequency and variance using FFT
Returns (dominant_frequency, frequency_variance)
"""
function estimate_frequency(data::Vector{Float64}, dt::Float64=1.0)
    n = length(data)
    if n < 10
        return NaN, NaN
    end
    
    # Remove mean (DC component)
    data_centered = data .- mean(data)
    
    # FFT
    fft_result = abs.(fft(data_centered))
    freqs = fftfreq(n, 1/dt)
    
    # Only positive frequencies (skip DC at index 1)
    pos_idx = 2:(n÷2)
    pos_freqs = abs.(freqs[pos_idx])
    pos_power = fft_result[pos_idx]
    
    if length(pos_freqs) < 1
        return NaN, NaN
    end
    
    # Find dominant frequency
    peak_idx = argmax(pos_power)
    dominant_freq = pos_freqs[peak_idx]
    
    # Calculate frequency variance (spectral spread around dominant)
    total_power = sum(pos_power)
    if total_power > 0
        weighted_freq = sum(pos_freqs .* pos_power) / total_power
        freq_variance = sum(pos_power .* (pos_freqs .- weighted_freq).^2) / total_power
    else
        freq_variance = NaN
    end
    
    return dominant_freq, freq_variance
end

"""
Calculate Kramers escape rate from ISI distribution
Escape rate λ = 1/⟨τ⟩ where τ is waiting time
"""
function calculate_escape_rate(isi::Vector{Float64})
    return 1.0 / mean(isi)
end

"""
Calculate escape rate with error estimate from multiple runs
Returns (mean_rate, std_rate)
"""
function calculate_escape_rate_with_error(isi_samples::Vector{Vector{Float64}})
    rates = Float64[]
    for isi in isi_samples
        if length(isi) > 2
            push!(rates, 1.0 / mean(isi))
        end
    end
    if length(rates) == 0
        return NaN, NaN
    end
    return mean(rates), std(rates)
end