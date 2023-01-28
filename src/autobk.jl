using Turing, GLM, DataFrames

"""
    get_nknots(Rbkg::Float64, Δk::Float64)

The maximum number of knots is 2*Rbkg*Δk/π, and the conservative estimate is to 
use one less than the highest integer less than or equal to this value. 

> Newville, M., P. Līviņš, š Y. Yacoby, J. J. Rehr, and E. A. Stern. 
> "Near-edge x-ray-absorption fine structure of Pb: 
> A comparison of theory and experiment." Physical Review B 47, 
> no. 21 (1993): 14126.
"""
get_nknots(Rbkg::Float64, Δk::Float64) = floor(Int64, 2 * Rbkg * Δk / π) - 1

"""
    find_E0(x::Vector{Float64}, y::Vector{Float64}, roi_min::Float64, 
    roi_max::Float64)

Finds the inflection point (here finding the peak of the first derivative 
within a region of interest). Currently not as many checks as in IFEFFIT

> Newville, M., P. Līviņš, š Y. Yacoby, J. J. Rehr, and E. A. Stern. 
> "Near-edge x-ray-absorption fine structure of Pb: 
> A comparison of theory and experiment." Physical Review B 47, 
> no. 21 (1993): 14126.
"""
function find_E0(x::Vector{Float64}, y::Vector{Float64}, roi_min::Float64, 
    roi_max::Float64)
    xs = (x[1:(end - 1)] .+ x[2:end]) ./ 2
    ys = (y[1:(end - 1)] .- y[2:end]) ./ (x[1:(end - 1)] .- x[2:end])
    istart = findfirst(x -> x >= roi_min, xs)
    istop = findlast(x -> x <= roi_max, xs)
    return xs[istart + argmax(ys[istart:istop]) - 1]
end

"""
    fit_pre_edge(x::Vector{Float64}, y::Vector{Float64}, emin::Float64, 
    emax::Float64)

Fit a line to for the pre-edge region.

> Newville, M., P. Līviņš, š Y. Yacoby, J. J. Rehr, and E. A. Stern. 
> "Near-edge x-ray-absorption fine structure of Pb: 
> A comparison of theory and experiment." Physical Review B 47, 
> no. 21 (1993): 14126.
"""
function fit_pre_edge(x::Vector{Float64}, y::Vector{Float64}, emin::Float64, 
    emax::Float64)
    istart = findfirst(E -> E >= emin, x)
    istop = findlast(E -> E <= emax, x)
    return lm(@formula(y ~ x), 
        DataFrame(x = x[istart:istop], y = y[istart:istop]))
end

"""
    fit_post_edge(x::Vector{Float64}, y::Vector{Float64}, emin::Float64, 
    emax::Float64)

Fit a quadratic to for the post-edge region.

> Newville, M., P. Līviņš, š Y. Yacoby, J. J. Rehr, and E. A. Stern. 
> "Near-edge x-ray-absorption fine structure of Pb: 
> A comparison of theory and experiment." Physical Review B 47, 
> no. 21 (1993): 14126.
"""
function fit_post_edge(x::Vector{Float64}, y::Vector{Float64}, emin::Float64, 
    emax::Float64)
    istart = findfirst(E -> E >= emin, x)
    istop = findlast(E -> E <= emax, x)
    return lm(@formula(y ~ x + x^2), 
        DataFrame(x = x[istart:istop], y = y[istart:istop]))
end




"""
    get_norm_μ(μx, μy, pre_emin, pre_emax, post_emin, post_emax)

"""
function get_norm_μ(μx, μy, E0, pre_fit, post_fit)
    
    
end

@model function autobk()

end
