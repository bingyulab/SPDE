using FFTW
using StatsBase

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
Extract last N points from solution DataFrame for a given variable
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
Calculate ATP:ADP ratio from solution DataFrame
"""
function calculate_atp_adp_ratio(df::DataFrame, n_points::Int=1000)
    atpc_data = get_last_n_points(df, "atpc", n_points)
    adpc_data = get_last_n_points(df, "adpc", n_points)
    
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
    if length(isi) < 2
        return NaN
    end
    return 1.0 / mean(isi)
end