using LaTeXStrings
using LinearAlgebra: ×, norm
using Statistics: mean

# Define units, LaTeX markup names, and LaTeX markup units for intrinsic values
const units_predefined = Dict(
   :rhom => ("kg/m3", L"$\rho_m$",  L"$\mathrm{kg}\,\mathrm{m}^{-3}$"),
   :rhoq => ("C/m3", L"$\rho_q$", L"$\mathrm{C}\,\mathrm{m}^{-3}$"),
   :rho => ("1/m3", L"$n_\mathrm{p}$", L"$\mathrm{m}^{-3}$"),
   :rhobackstream => ("1/m3",  L"$n_\mathrm{p,st}$", L"$\mathrm{m}^{-3}$"),
   :rhononbackstream => ("1/m3", L"$n_\mathrm{p,th}$", L"$\mathrm{m}^{-3}$"),
   :rho_v => ("1/m2s", L"$\Gamma_\mathrm{p}$", L"$\mathrm{m}^{-2}$s"),
   :rhovbackstream => ("1/m2s", L"$\Gamma_\mathrm{p,st}$", L"$\mathrm{m}^{-2}$s"),
   :rhovnonbackstream => ("1/m2s", L"$\Gamma_\mathrm{p,th}$", L"$\mathrm{m}^{-2}$s"),
   :v => ("m/s", L"$V$", L"$\mathrm{m}\,\mathrm{s}^{-1}$"),
   :vbackstream => ("m/s", L"$V_\mathrm{p,st}$", L"$\mathrm{m}\,\mathrm{s}^{-1}$"),
   :vNonbackstream => ("m/s", L"$V_\mathrm{p,th}$", L"$\mathrm{m}\,\mathrm{s}^{-1}$"),
   :b => ("T", L"$B$", "T"),
   :b_vol => ("T",  L"$B_\mathrm{vol}$", "T"),
   :background_b => ("T",  L"$B_\mathrm{bg}$", "T"),
   :perturbed_b => ("T", L"$B_\mathrm{pert}$", "T"),
   :bgb => ("T", L"$B_\mathrm{bg}$", "T"),
   :perb => ("T", L"B_\mathrm{pert}$", "T"),
   :perb_vol => ("T", L"B_\mathrm{vol,pert}$", "T"),
   :e => ("V/m", L"$E$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :e_vol => ("V/m",  L"$E_\mathrm{vol}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :exhall_000_100 => ("V/m", L"$E_\mathrm{Hall,000,100}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :exhall_001_101 => ("V/m", L"$E_\mathrm{Hall,001,101}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :exhall_010_110 => ("V/m", L"$E_\mathrm{Hall,010,110}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :exhall_011_111 => ("V/m", L"$E_\mathrm{Hall,011,111}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :eyhall_000_010 => ("V/m", L"$E_\mathrm{Hall,000,010}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :eyhall_001_011 => ("V/m", L"$E_\mathrm{Hall,001,011}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :eyhall_100_110 => ("V/m", L"$E_\mathrm{Hall,100,110}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :eyhall_101_111 => ("V/m", L"$E_\mathrm{Hall,101,111}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :ezhall_000_001 => ("V/m", L"$E_\mathrm{Hall,000,001}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :ezhall_010_011 => ("V/m", L"$E_\mathrm{Hall,010,011}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :ezhall_100_101 => ("V/m", L"$E_\mathrm{Hall,100,101}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :ezhall_110_111 => ("V/m", L"$E_\mathrm{Hall,110,111}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
   :pressure => ("Pa",  L"$P$", "Pa"),
   :pressure_dt2 => ("Pa", L"$P_{\mathrm{d}t/2}}$", "Pa"),
   :pressure_r => ("Pa",  L"$P_r$", "Pa"),
   :pressure_v => ("Pa", L"$P_v$", "Pa"),
   :ptensordiagonal => ("Pa", L"$\mathcal{P}_\mathrm{diag}$", "Pa"),
   :ptensoroffdiagonal => ("Pa", L"$\mathcal{P}_\mathrm{off-diag}$", "Pa"),
   :ptensorbackstreamdiagonal => ("Pa", L"$\mathcal{P}_\mathrm{st,diag}$", "Pa"),
   :ptensorbackstreamoffdiagonal => ("Pa",  L"$\mathcal{P}_\mathrm{st,off-diag}$", "Pa"),
   :ptensornonbackstreamdiagonal => ("Pa", L"$\mathcal{P}_\mathrm{th,diag}$", "Pa"),
   :ptensornonbackstreamoffdiagonal => ("Pa", L"$\mathcal{P}_\mathrm{th,off-diag}$", "Pa"),
   :t => ("K", L"$T$", "K"),
   :tpar => ("K", L"$T$", "K"),
   :tperp => ("K", L"$T$", "K"),
   :max_v_dt => ("s", L"$\Delta t_{\mathrm{max},v}$", "s"),
   :max_r_dt => ("s", L"$\Delta t_{\mathrm{max},r}$", "s"),
   :max_fields_dt => ("s", L"$\Delta t_\mathrm{max,FS}$", "s"),
   :minvalue => ("s3/m6", L"$f_\mathrm{Min}$", L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$"),
   :effectivesparsitythreshold => ("s3/m6", L"$f_\mathrm{Min}$", L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$"),
   :rho_loss_adjust => ("1/m3", L"$\Delta_\mathrm{loss} n_\mathrm{p}$", L"$\mathrm{m}^{-3}$"),
   :energydensity => ("eV/cm3", L"$\rho_{\mathrm{energy}}$", L"$\mathrm{eV}\,\mathrm{cm}^{-3}$"),
   :precipitationdiffflux => ("1/(cm2 sr s eV)", L"$'Delta F_\mathrm{precipitation}$", L"$\mathrm{cm}^{-2} \,\mathrm{sr}^{-1}\,\mathrm{s}^{-1}\,\mathrm{eV}^{-1}$"),
)

# Define derived parameters
const variables_predefined = Dict(
   :Bmag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = isempty(ids) ?
         readvariable(meta, meta.variable[rho_]) :
         readvariable(meta, meta.variable[rho_], ids)
      if hasvariable(meta, "vg_b_vol")
         Bmag = isempty(ids) ?
            sqrt.(sum(readvariable(meta, "vg_b_vol").^2, dims=1)) :
            sqrt.(sum(readvariable(meta, "vg_b_vol", ids).^2, dims=1))
      else
         @assert isempty(ids) "Do not support reading selected cells from FSGrid!"
         Bmag = sqrt.(sum(readvariable(meta, "fg_b").^2, dims=1))
      end
      @inbounds for i = eachindex(ρ) # sparsity/inner boundary
         Bmag[i] == 0.0 && (Bmag[i] = NaN)
      end
      Bmag
   end,
   :Emag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = isempty(ids) ?
         readvariable(meta, meta.variable[rho_]) :
         readvariable(meta, meta.variable[rho_], ids)
      if hasvariable(meta, "vg_e_vol")
         Emag = isempty(ids) ?
            sqrt.(sum(readvariable(meta, "vg_e_vol").^2, dims=1)) :
            sqrt.(sum(readvariable(meta, "vg_e_vol", ids).^2, dims=1))
      else
         @assert isempty(ids) "Do not support reading selected cells from FSGrid!"
         Emag = sqrt.(sum(readvariable(meta, "fg_e").^2, dims=1))
      end
      @inbounds for i = eachindex(ρ) # sparsity/inner boundary
         Emag[i] == 0.0 && (Emag[i] = NaN)
      end
      Emag
   end,
   :Vmag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = isempty(ids) ?
         readvariable(meta, meta.variable[rho_]) :
         readvariable(meta, meta.variable[rho_], ids)
      Vmag = isempty(ids) ?
         vec(sqrt.(sum(readvariable(meta, "proton/vg_v").^2, dims=1))) :
         vec(sqrt.(sum(readvariable(meta, "proton/vg_v", ids).^2, dims=1)))
      @inbounds for i = eachindex(ρ) # sparsity/inner boundary
         Vmag[i] == 0.0 && (Vmag[i] = NaN)
      end
      Vmag
   end,
   :Rhom => function (meta, ids=UInt64[])
      if hasvariable(meta, "vg_rhom")
         ρm = isempty(ids) ?
            readvariable(meta, "vg_rhom") :
            readvariable(meta, "vg_rhom", ids)
      elseif hasvariable(meta, "proton/vg_rho")
         ρm = isempty(ids) ?
            readvariable(meta, "proton/vg_rho") .* mᵢ :
            readvariable(meta, "proton/vg_rho", ids) .* mᵢ
      end
      ρm
   end,
   :P => function (meta, ids=UInt64[]) # scalar pressure
      if hasvariable(meta, "vg_pressure")
         P = isempty(ids) ?
            readvariable(meta, "vg_pressure") :
            readvariable(meta, "vg_pressure", ids)
      else
         Pdiag = isempty(ids) ?
            readvariable(meta, "proton/vg_ptensor_diagonal") :
            readvariable(meta, "proton/vg_ptensor_diagonal", ids)
         P = vec(mean(Pdiag, dims=1))
      end
      P
   end,
   :VS => function (meta, ids=UInt64[]) # sound speed
      P = readvariable(meta, "P", ids)
      ρm = readvariable(meta, "Rhom", ids)
      @inbounds for i = eachindex(ρm) # sparsity/inner boundary
         ρm[i] == 0.0 && (ρm[i] = NaN)
      end
      vs = @. √( (P*5.0f0/3.0f0) / ρm )
   end,
   :VA => function (meta, ids=UInt64[]) # Alfvén speed
      ρm = readvariable(meta, "Rhom", ids)
      @inbounds for i = eachindex(ρm) # sparsity/inner boundary
         ρm[i] == 0.0 && (ρm[i] = NaN)
      end
      Bmag = readvariable(meta, "Bmag", ids)
      VA = @. $vec(Bmag) / √(ρm*μ₀)
   end,
   :MA => function (meta, ids=UInt64[]) # Alfvén Mach number
      V = readvariable(meta, "Vmag", ids)
      @inbounds for i = eachindex(V) # sparsity/inner boundary
         V[i] == 0.0 && (V[i] = NaN)
      end
      VA = readvariable(meta, "VA", ids)
      V ./ VA
   end,
   :MS => function (meta, ids=UInt64[]) # Sonic Mach number
      V = readvariable(meta, "Vmag", ids)
      @inbounds for i = eachindex(V) # sparsity/inner boundary
         V[i] == 0.0 && (V[i] = NaN)
      end
      VS = readvariable(meta, "VS", ids)
      V ./ VS
   end,
   :Vpar => function (meta, ids=UInt64[]) # velocity ∥ B
      V = isempty(ids) ?
         readvariable(meta, "proton/vg_v") :
         readvariable(meta, "proton/vg_v", ids)
      if isempty(ids)
         b = readvariable(meta, "vg_b_vol") ./ readvariable(meta, "Bmag")
      else
         b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)
      end
      [V[:,i] ⋅ b[:,i] for i in 1:size(V,2)]
   end,
   :Vperp => function (meta, ids=UInt64[]) # velocity ⟂ B
      V = isempty(ids) ?
         readvariable(meta, "proton/vg_v") :
         readvariable(meta, "proton/vg_v", ids)
      if isempty(ids)
         b = readvariable(meta, "vg_b_vol") ./ readvariable(meta, "Bmag")
      else
         b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)
      end
      Vperp = zeros(eltype(V), size(V, 2))
      # Avoid sqrt of negative values, but does not guarantee orthogonality.
      @inbounds for i in eachindex(Vperp)
         Vpar = V[:,i] ⋅ b[:,i] .* b[:,i]
         Vperp[i] = norm(V[:,i] - Vpar)
      end
      Vperp
   end,
   :T => function (meta, ids=UInt64[]) # scalar temperature
      P = readvariable(meta, "P", ids)
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      @inbounds for i = eachindex(n) # sparsity/inner boundary
         n[i] == 0.0 && (n[i] = NaN)
      end
      T = @. P / (n*kB)
   end,
   :Ppar => function (meta, ids=UInt64[]) # P component ∥ B
      P = readvariable(meta, "Protated", ids)
      @views P[3,3,:]
   end,
   :Pperp => function (meta, ids=UInt64[]) # P component ⟂ B
      P = readvariable(meta, "Protated", ids)
      Pperp = [0.5f0(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]
   end,
   :Tpar => function (meta, ids=UInt64[]) # T component ∥ B
      P = readvariable(meta, "Protated", ids)
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      @inbounds for i = eachindex(n) # sparsity/inner boundary
         n[i] == 0.0 && (n[i] = NaN)
      end
      @. P[3,3,:] / (n*kB)
   end,
   :Tperp => function (meta, ids=UInt64[]) # scalar T component ⟂ B
      P = readvariable(meta, "Protated", ids)
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      @inbounds for i = eachindex(n) # sparsity/inner boundary
         n[i] == 0.0 && (n[i] = NaN)
      end
      Pperp = [0.5(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]
      @. Pperp / (n*kB)
   end,
   :Tanisotropy => function (meta, ids=UInt64[]) # T⟂ / T∥
      Tperp = readvariable(meta, "Tperp", ids)
      Tpar = readvariable(meta, "Tpar", ids)
      @. Tperp / Tpar
   end,
   :J => function (meta, ids=UInt64[])
      @assert isempty(ids) "Do not support current calculation for selected cells!"
      B = readvariable(meta, "vg_b_vol")
      B = reshape(B, 3, meta.ncells...)
      J = curl(meta.dcoord, B) ./ μ₀
      J = reshape(J, 3, :) # To be consistent with shape assumptions
   end,
   :Protated => function (meta, ids=UInt64[])
      # Rotate the pressure tensor to align the 3rd direction with B
      B = isempty(ids) ?
         readvariable(meta, "vg_b_vol") :
         readvariable(meta, "vg_b_vol", ids)
      Pdiag = isempty(ids) ?
         readvariable(meta, "proton/vg_ptensor_diagonal") :
         readvariable(meta, "proton/vg_ptensor_diagonal", ids)
      Podiag = isempty(ids) ?
         readvariable(meta, "proton/vg_ptensor_offdiagonal") :
         readvariable(meta, "proton/vg_ptensor_offdiagonal", ids)
      P = zeros(Float32, 3, 3, size(Pdiag, 2))
      @inbounds for i = 1:size(P, 3)
         P[1,1,i] = Pdiag[1,i]
         P[2,2,i] = Pdiag[2,i]
         P[3,3,i] = Pdiag[3,i]
         P[1,2,i] = P[2,1,i] = Podiag[1,i]
         P[2,3,i] = P[3,2,i] = Podiag[2,i]
         P[3,1,i] = P[1,3,i] = Podiag[3,i]
         @views rotateWithB!(P[:,:,i], B[:,i])
      end
      P
   end,
   :Panisotropy => function (meta, ids=UInt64[]) # P⟂ / P∥
      PR = readvariable(meta, "Protated", ids)
      @. 0.5f0*(PR[1,1,:] + PR[2,2,:]) / PR[3,3,:]
   end,
   :Pdynamic => function (meta, ids=UInt64[])
      V = readvariable(meta, "Vmag", ids)
      ρm = readvariable(meta, "Rhom", ids)
      @. ρm * V * V
   end,
   :Pb => function (meta, ids=UInt64[])
      B² = isempty(ids) ?
         vec(sum(readvariable(meta, "vg_b_vol").^2, dims=1)) :
         vec(sum(readvariable(meta, "vg_b_vol", ids).^2, dims=1))
      @inbounds for i = eachindex(B²) # sparsity/inner boundary
         B²[i] == 0.0 && (B²[i] = NaN)
      end
      mu2inv = 0.5/μ₀
      @. B² * mu2inv
   end,
   :Poynting => function (meta, ids=UInt64[])
      if hasvariable(meta, "vg_b_vol") && hasvariable(meta, "vg_e_vol")
         E = isempty(ids) ?
            readvariable(meta, "vg_e_vol") :
            readvariable(meta, "vg_e_vol", ids)
         B = isempty(ids) ?
            readvariable(meta, "vg_b_vol") :
            readvariable(meta, "vg_b_vol", ids)
      elseif hasvariable(meta, "fg_e") && hasvariable(meta, "fg_b")
         @assert isempty(ids) "Do not support selecting cells in FSGrid!"
         E = readvariable(meta, "fg_e")
         B = readvariable(meta, "fg_b")
      else
         E = readvariable(meta, "E")
         B = readvariable(meta, "B")
      end
      F = similar(E)
      Rpost = CartesianIndices(size(E)[2:end])
      @inbounds for i in Rpost
         F[:,i] = E[:,i] × B[:,i] ./ μ₀
      end
      F
   end,
   :Agyrotropy => function (meta, ids=UInt64[])
      # non-gyrotropy measure Q [Swisdak 2016]
      # Original derivation for electrons. Here we do protons first.
      if isempty(ids)
         Pdiag = readvariable(meta, "proton/vg_ptensor_diagonal")
         Podiag = readvariable(meta, "proton/vg_ptensor_offdiagonal")
         B = readvariable(meta, "vg_b_vol")
      else
         Pdiag = readvariable(meta, "proton/vg_ptensor_diagonal", ids)
         Podiag = readvariable(meta, "proton/vg_ptensor_offdiagonal", ids)
         B = readvariable(meta, "vg_b_vol", ids)
      end
      Pxx = selectdim(Pdiag, 1, 1)
      Pyy = selectdim(Pdiag, 1, 2)
      Pzz = selectdim(Pdiag, 1, 3)
      Pxy = selectdim(Podiag, 1, 1)
      Pyz = selectdim(Podiag, 1, 2) # Warning: the order may be wrong!
      Pxz = selectdim(Podiag, 1 ,3) # Warning: the order may be wrong!

      Bx = selectdim(B, 1, 1)
      By = selectdim(B, 1, 2)
      Bz = selectdim(B, 1, 3)

      I₁ = @. Pxx + Pyy + Pzz
      I₂ = @. Pxx*Pyy + Pxx*Pzz + Pyy*Pzz - Pxy*Pxy - Pyz*Pyz - Pxz*Pxz

      Ppar = @. (Bx*Bx*Pxx + By*By*Pyy + Bz*Bz*Pzz +
	      2*(Bx*By*Pxy + Bx*Bz*Pxz + By*Bz*Pyz))/B²
      Qsqr = @. √(1 - 4I₂/((I₁ - Ppar)*(I₁ + 3Ppar)))
   end,
   :Beta => function (meta, ids=UInt64[])
      P = readvariable(meta, "P", ids)
      B² = isempty(ids) ?
         vec(sum(readvariable(meta, "vg_b_vol").^2, dims=1)) :
         vec(sum(readvariable(meta, "vg_b_vol", ids).^2, dims=1))
      @inbounds for i = eachindex(B²) # sparsity/inner boundary
         B²[i] == 0.0 && (B²[i] = NaN)
      end
      @. 2 * μ₀ * P / B²
   end,
   :IonInertial => function (meta, ids=UInt64[])
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      Z = 1
      ωi = @. √(n/(mᵢ*μ₀)) * Z * qᵢ
      di = @. c / ωi
   end,
   :Larmor => function (meta, ids=UInt64[])
      if isempty(ids)
         Vperp = readvariable(meta, "Vperp")
         B = readvariable(meta, "Bmag")
      else
         Vperp = readvariable(meta, "Vperp", ids)
         B = readvariable(meta, "Bmag", ids)
      end
      rg = @. mᵢ * Vperp / (qᵢ * B)
   end,
   :Gyroperiod => function (meta, ids=UInt64[])
      B = readvariable(meta, "Bmag", ids)
      T = @. 2π * mᵢ / (qᵢ * B)
   end,
   :Gyrofrequency => function (meta, ids=UInt64[]) # [1/s]
      B = isempty(ids) ? readvariable(meta, "Bmag") : readvariable(meta, "Bmag", ids)
      f = @. qᵢ * B / (mᵢ * 2π)
   end,
   :Plasmaperiod => function (meta, ids=UInt64[])
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      T = @. 2π / (qᵢ * √(n  / (mᵢ * ϵ₀)))
   end,
   :Omegap => function (meta, ids=UInt64[]) # plasma frequency, [1/s]
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids)
      ωₚ = @. qᵢ * √(n  / (mᵢ * ϵ₀)) / 2π
   end,
)