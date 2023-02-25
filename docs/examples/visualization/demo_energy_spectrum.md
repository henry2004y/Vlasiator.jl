# ---
# title: Energy spectrum at a fixed point
# id: demo_energy_spectrum
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot energy spectrum at a fixed point
# ---

This demo shows how to plot energy spectrum of oscillations at a fixed point.
```julia
using FFTW, JLD2, CurveFit, PyPlot
using Vlasiator: RE

function main()
   file = "satellites_uniform_sampled.jld2"
   data = JLD2.load(file)

   nSatellite = length(data["t"])
   nI, nJ = size(data["rho"])[2:3]

   t = data["t"]
   # Select spatial point
   i, j = 5, 5
   var = data["rho"][:,i,j]

   dt = t[2] - t[1] # uniform sample interval [s]
   Fs = 1 / dt      # sample frequency, [Hz]
   Fn = Fs / 2      # Nyquist frequency, [Hz]

   ## Frequency calculation
   nPoints = length(var)
   nFFT = nPoints
   df = Fs / nFFT
   freq_fullrange = -Fn:df:Fn

   freq = freq_fullrange[(nPoints ÷ 2 + 1):end-1]

   var_freq = fft(var)
   var_power = abs.(fftshift(var_freq))[(nPoints ÷ 2 + 1):end]

   # k is the exponential coefficient
   a, k = @views power_fit(freq[10:end], var_power[10:end])

   figure(figsize=(6,8))
   loglog(freq, var_power)
   axis("scaled")

   min_power, max_power = extrema(@view var_power[10:end])
   xlim(freq[8], Fs)
   ylim(min_power * 0.75, max_power * 2.0)

   xlabel("Frequency [Hz]"; fontsize=14)
   ylabel("Power Density "; fontsize=14)
   title(string(round.(data["locations"][i,j]./RE, digits=1))*"Re"; fontsize=14)

   loglog(freq[10:end], a.*freq[10:end].^k, label="k = $(round(k, digits=1))")

   legend()
end

main()
```