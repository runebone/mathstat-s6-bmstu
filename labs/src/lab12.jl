# module main

# export x, Mmax, Mmin, R, n, μ, S², m, fig, Δ

x = [
3.38, 1.21, 1.85, 2.24, 4.17, 2.99, 4.81, 2.71, 2.70, 4.41,
3.21, 3.15, 2.77, 4.05, 3.89, 1.56, 2.78, 2.04, 2.82, 3.28,
2.63, 1.89, 3.57, 3.15, 3.80, 5.40, 3.25, 2.04, 2.61, 5.06,
2.87, 2.66, 4.80, 3.86, 0.09, 2.45, 2.40, 2.14, 1.69, 2.36,

5.44, 2.77, 1.94, 2.55, 3.97, 1.88, 3.01, 4.21, 4.74, 2.02,
2.38, 2.46, 3.51, 2.89, 1.57, 3.53, 0.77, 3.31, 3.58, 2.77,
3.61, 3.71, 2.38, 3.06, 4.29, 4.76, 1.69, 1.59, 3.21, 2.74,
3.99, 3.53, 3.52, 2.84, 1.21, 2.82, 4.34, 3.65, 2.22, 2.87,

3.14, 3.58, 1.96, 3.41, 3.85, 1.96, 3.02, 4.22, 3.10, 2.68,
3.67, 1.70, 5.47, 5.02, 2.52, 3.09, 2.19, 4.44, 2.33, 2.27,
3.34, 3.05, 4.35, 3.58, 3.43, 4.49, 3.57, 3.20, 1.53, 3.53,
3.53, 1.27, 3.40, 4.53, 2.21, 3.28, 3.50, 2.01, 3.30, 1.86,
]

Mmax = findmax(x)[1]
Mmin = findmin(x)[1]

R = Mmax - Mmin
n = length(x)
μ = 1/n * sum(x)

S² = 1/(n-1) * sum((x -> (x - μ)^2).(x))

m = Int(floor(log2(n)) + 2)
Δ = R / m

using Printf

@printf "Mmin = %f\nMmax = %f\nR = %f\nn = %d\nμ = %f\nS² = %f\nm = %d\nΔ = %f\n" Mmin Mmax R n μ S² m Δ

icnt = zeros(m+1)
for i = 1:m+1
    prev = Mmin + Δ * (i - 1)
    next = prev + Δ
    for j = 1:n
        if prev ≤ x[j] ≤ next
            icnt[i] += 1
        end
    end
    # В интервале [ 6.810,  7.510)  7 элементов.
    @printf "В интервале [ %4.3f, %4.3f) %d элементов.\n" prev next icnt[i]
end

using GLMakie

fig = Figure()
# ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 1])

# hist!(ax1, x, m, strokewidth = 1, strokecolor = :black, normalization = :pdf),
ecdfplot!(ax2, x)

using Distributions

σ = √S²
nd = Normal(μ, σ)
ndpdf(x) = pdf(nd, x)
ndcdf(x) = cdf(nd, x)

xrange = Mmin:0.1:Mmax
# lines!(ax1, xrange, ndpdf.(xrange), color = :red)
lines!(ax2, xrange, ndcdf.(xrange), color = :red)

# --- Lab. 2

function get_sample_mean(x)
    return sum(x) / length(x)
end

function get_sample_variance(x)
    μ = get_sample_mean(x)
    n = length(x)
    return 1/(n-1) * sum((x -> (x - μ)^2).(x))
end

function get_mx_confidence(x, γ)
    sm = get_sample_mean(x)
    sv = get_sample_variance(x)
    n  = length(x)
    td = TDist(n - 1)
    t = quantile(td, (1 + γ) / 2)
    c = t * √sv / √n 
    return (sm - c, sm + c)
end

function get_dx_confidence(x, γ)
    sv = get_sample_variance(x)
    n  = length(x)
    hd = Chisq(n - 1)
    tlo = quantile(hd, (1 + γ) / 2)
    thi = quantile(hd, (1 - γ) / 2)
    c = n * sv
    return (c / tlo, c / thi)
end

f = Figure()
a1 = Axis(f[1, 1:2])
a21 = Axis(f[2, 1])
a22 = Axis(f[2, 2])

n = length(x)
xrange = 2:n
mx = get_sample_mean(x)
dx = get_sample_variance(x)

γ = 0.9

ms   = fill(mx, n-1)
mxcs = [get_mx_confidence(x[1:i], γ) for i in 2:n] # end
mxs  = [get_sample_mean(x[1:i]) for i in 2:n] # end
mxls = [mxc[1] for mxc in mxcs] # end
mxhs = [mxc[2] for mxc in mxcs] # end

ds   = fill(dx, n-1)
dxcs = [get_dx_confidence(x[1:i], γ) for i in 2:n] # end
dxs  = [get_sample_variance(x[1:i]) for i in 2:n] # end
dxls = [dxc[1] for dxc in dxcs] # end
dxhs = [dxc[2] for dxc in dxcs] # end

lines!(a1, xrange, ms, color = :red, linewidth = 3)
lines!(a1, xrange, mxs, linewidth = 3, linestyle = :dashdot)
lines!(a1, xrange, mxls, linewidth = 3, linestyle = :dash)
lines!(a1, xrange, mxhs, linewidth = 3, linestyle = :dot)

lines!(a21, xrange, ds, color = :red, linewidth = 3)
lines!(a21, xrange, dxs, linewidth = 3, linestyle = :dashdot)
lines!(a21, xrange, dxls, linewidth = 3, linestyle = :dash)
lines!(a21, xrange, dxhs, linewidth = 3, linestyle = :dot)

lines!(a22, xrange[4:end], ds[4:end], color = :red, linewidth = 3)
lines!(a22, xrange[4:end], dxs[4:end], linewidth = 3, linestyle = :dashdot)
lines!(a22, xrange[4:end], dxls[4:end], linewidth = 3, linestyle = :dash)
lines!(a22, xrange[4:end], dxhs[4:end], linewidth = 3, linestyle = :dot)

0

# end # module main
