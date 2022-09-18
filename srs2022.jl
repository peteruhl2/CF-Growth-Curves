### plot for srs 2022

using Plots

r = 16.
Ks = 0.25
n = 3

e = 1

x = collect(range(0,e,length=100))
m2 = r*x./(Ks .+ x)
m3 = r*x.^n./(Ks^n .+ x.^n)

plot(x->r, label = "Model 1", lw = 2)
plot!(x, m2, label = "Model 2", lw = 2)
plot!(x, m3, label = "Model 3", lw = 2)
plot!(xlims = (0, e), ylims = (0,25))
plot!(xlabel = "Nutrient - z")
plot!(ylabel = "Growth Rate - f(z)")
plot!(legendfont=font(14))
# plot!(xlabelfont = font(14))