# Predefined physical constants

const qₑ = -1.60217662e-19  # electron charge, [C]
const mₑ = 9.10938356e-31   # electron mass, [kg]
const qᵢ = 1.60217662e-19   # proton mass, [C]
const mᵢ = 1.673557546e-27  # proton mass, [kg]
const c  = 299792458.       # speed of light, [m/s]
const μ₀ = 4π*1e-7          # Vacuum permeability, [H/m]
const ϵ₀ = 1/(c^2*μ₀)       # Vacuum permittivity, [F/m]
const kB = 1.38064852e-23   # Boltzmann constant, [m²kg/(s²K)]
const RE = 6.371e6          # Earth radius, [m]

# Define units, LaTeX markup names, and LaTeX markup units for intrinsic values
const units_predefined = Dict(
   :Rhom => ("kg/m3", L"$\rho_m$",  L"$\mathrm{kg}\,\mathrm{m}^{-3}$"),
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
   :max_v_dt => ("s", L"$\Delta t_{\mathrm{max},v}$", "s"),
   :max_r_dt => ("s", L"$\Delta t_{\mathrm{max},r}$", "s"),
   :max_fields_dt => ("s", L"$\Delta t_\mathrm{max,FS}$", "s"),
   :minvalue => ("s3/m6", L"$f_\mathrm{Min}$", L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$"),
   :effectivesparsitythreshold => ("s3/m6", L"$f_\mathrm{Min}$", L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$"),
   :rho_loss_adjust => ("1/m3", L"$\Delta_\mathrm{loss} n_\mathrm{p}$", L"$\mathrm{m}^{-3}$"),
   :energydensity => ("eV/cm3", L"$\rho_{\mathrm{energy}}$", L"$\mathrm{eV}\,\mathrm{cm}^{-3}$"),
   :precipitationdiffflux => ("1/(cm2 sr s eV)", L"$\Delta F_\mathrm{precipitation}$", L"$\mathrm{cm}^{-2} \,\mathrm{sr}^{-1}\,\mathrm{s}^{-1}\,\mathrm{eV}^{-1}$"),
   :T => ("K", L"$T$", "K"),
   :Tpar => ("K", L"$T$", "K"),
   :Tperp => ("K", L"$T$", "K"),
   :Panisotropy => ("", L"$P_\perp / P_\parallel$", ""),
   :Tanisotropy => ("", L"$T_\perp / T_\parallel$", ""),
   :VS => ("m/s", L"$V_S$", L"$\mathrm{m}\,\mathrm{s}^{-1}$"),
   :VA => ("m/s", L"$V_A$", L"$\mathrm{m}\,\mathrm{s}^{-1}$"),
   :MS => ("", L"$M_S$", ""),
   :MA => ("", L"$M_A$", ""),
   :Ppar => ("Pa", L"$P_\parallel$", "Pa"),
   :Pperp => ("Pa", L"$P_\perp$", "Pa"),
   :Beta => ("",L"$\beta$", ""),
   :BetaStar => ("",L"$\beta^\ast$", ""),
   :Gyroperiod => ("s", L"$T_{gyro}$", "s"),
   :PlasmaPeriod => ("s", L"$T_{plasma}$", "s"),
   :Gyrofrequency => ("/s", L"$\omega_{g}$", L"s^{-1}"),
   :Omegap => ("/s", L"$\omega_{p}$", L"s^{-1}"),
)

# Define derived parameters
const variables_predefined = Dict(
   :Bmag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = readvariable(meta, meta.variable[rho_], ids)
      B = readvariable(meta, "vg_b_vol", ids)
      Bmag = sum(x -> x*x, B, dims=1) |> vec
      @. Bmag = √(Bmag)
      _fillinnerBC!(Bmag, ρ)
      Bmag::Vector{Float32}
   end,
   :Emag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = readvariable(meta, meta.variable[rho_], ids)
      E = readvariable(meta, "vg_e_vol", ids)
      Emag = sum(x -> x*x, E, dims=1) |> vec
      @. Emag = √(Emag)
      _fillinnerBC!(Emag, ρ)
      Emag::Vector{Float32}
   end,
   :Vmag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = readvariable(meta, meta.variable[rho_], ids)
      V = readvariable(meta, "proton/vg_v", ids)
      Vmag = sum(x -> x*x, V, dims=1) |> vec
      @. Vmag = √(Vmag)
      _fillinnerBC!(Vmag, ρ)
      Vmag::Vector{Float32}
   end,
   :Rhom => function (meta, ids=UInt64[])
      if hasvariable(meta, "vg_rhom")
         ρm = readvariable(meta, "vg_rhom", ids)
      elseif hasvariable(meta, "proton/vg_rho")
         ρm = readvariable(meta, "proton/vg_rho", ids) .* mᵢ
      end
      ρm::Union{Vector{Float32}, Vector{Float64}}
   end,
   :P => function (meta, ids=UInt64[]) # scalar pressure
      if hasvariable(meta, "vg_pressure")
         P = readvariable(meta, "vg_pressure", ids)
      else
         if hasvariable(meta, "proton/vg_ptensor_diagonal")
            Pdiag_str = "proton/vg_ptensor_diagonal"
         else
            Pdiag_str = "PTensorDiagonal"
         end
         Pdiag = readvariable(meta, Pdiag_str, ids)
         P = vec(mean(Pdiag, dims=1))
      end
      P::Vector{Float32}
   end,
   :VS => function (meta, ids=UInt64[]) # sound speed
      P = readvariable(meta, "P", ids)
      ρm = readvariable(meta, "Rhom", ids)
      _fillinnerBC!(ρm, ρm)
      VS = @. √( (P*5.0f0/3.0f0) / ρm )
      VS::Union{Vector{Float32}, Vector{Float64}}
   end,
   :VA => function (meta, ids=UInt64[]) # Alfvén speed
      ρm = readvariable(meta, "Rhom", ids)
      _fillinnerBC!(ρm, ρm)
      Bmag = readvariable(meta, "Bmag", ids)
      VA = @. Bmag / √(ρm*μ₀)
      VA::Union{Vector{Float32}, Vector{Float64}}
   end,
   :MA => function (meta, ids=UInt64[]) # Alfvén Mach number
      V = readvariable(meta, "Vmag", ids)
      VA = readvariable(meta, "VA", ids)
      MA = V ./ VA
      MA::Union{Vector{Float32}, Vector{Float64}}
   end,
   :MS => function (meta, ids=UInt64[]) # Sonic Mach number
      V = readvariable(meta, "Vmag", ids)
      VS = readvariable(meta, "VS", ids)
      MS = V ./ VS
      MS::Union{Vector{Float32}, Vector{Float64}}
   end,
   :Vpar => function (meta, ids=UInt64[]) # velocity ∥ B
      V = readvariable(meta, "proton/vg_v", ids)
      b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)'
      Vpar = @views [V[:,i] ⋅ b[:,i] for i in axes(V,2)]::Vector{Float32}
   end,
   :Vperp => function (meta, ids=UInt64[]) # velocity ⟂ B
      V = readvariable(meta, "proton/vg_v", ids)
      b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)'
      Vperp = Vector{Float32}(undef, size(V, 2))

      function _computeVperp!()
         # Avoid sqrt of negative values, but does not guarantee orthogonality.
         @inbounds @views for i in eachindex(Vperp)
            vb = V[:,i] ⋅ b[:,i]
            Vpar = SVector{3,Float32}(vb*b[1,i], vb*b[2,i], vb*b[3,i])
            Vperp[i] = norm(V[:,i] - Vpar)
         end
      end

      _computeVperp!()

      Vperp
   end,
   :T => function (meta, ids=UInt64[]) # scalar temperature
      P = readvariable(meta, "P", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      T = @. P / (n*kB)
      T::Vector{Float64}
   end,
   :Vth => function (meta, ids=UInt64[]) # thermal velocity
      T = readvariable(meta, "T", ids)
      Vth = @. √(3 * kB * T / mᵢ) # assume proton
      Vth::Vector{Float64}
   end,
   :Ppar => function (meta, ids=UInt64[]) # P component ∥ B
      P = readvariable(meta, "Protated", ids)
      Ppar = P[3,3,:]::Vector{Float32}
   end,
   :Pperp => function (meta, ids=UInt64[]) # P component ⟂ B
      P = readvariable(meta, "Protated", ids)
      Pperp = [0.5f0(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]::Vector{Float32}
   end,
   :Tpar => function (meta, ids=UInt64[]) # T component ∥ B
      P = readvariable(meta, "Protated", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      Tpar = @. P[3,3,:] / (n*kB)
      Tpar::Vector{Float64}
   end,
   :Tperp => function (meta, ids=UInt64[]) # scalar T component ⟂ B
      P = readvariable(meta, "Protated", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      Pperp = @inbounds [0.5(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]
      Tperp = @. Pperp / (n*kB)
      Tperp::Vector{Float64}
   end,
   :Tanisotropy => function (meta, ids=UInt64[]) # T⟂ / T∥
      Tperp = readvariable(meta, "Tperp", ids)
      Tpar = readvariable(meta, "Tpar", ids)
      Taniso = @. Tperp / Tpar
      Taniso::Vector{Float64}
   end,
   :J => function (meta, ids=UInt64[])
      @assert isempty(ids) "Do not support current calculation for selected cells!"
      B = readvariable(meta, "vg_b_vol")
      B = reshape(B, 3, meta.ncells...)
      J = curl(B, meta.dcoord) ./ μ₀
      J = reshape(J, 3, :)::Array{Float64, 2}
   end,
   :Protated => function (meta, ids=UInt64[])
      # Rotate the pressure tensor to align the 3rd direction with B
      B = readvariable(meta, "vg_b_vol", ids)
      Pdiag = readvariable(meta, "proton/vg_ptensor_diagonal", ids)
      Podiag = readvariable(meta, "proton/vg_ptensor_offdiagonal", ids)
      P = Array{Float32}(undef, 3, 3, size(Pdiag, 2))

      function rotate_tensor!()
         @inbounds for i in axes(P, 3)
            P[1,1,i] = Pdiag[1,i]
            P[2,2,i] = Pdiag[2,i]
            P[3,3,i] = Pdiag[3,i]
            P[1,2,i] = P[2,1,i] = Podiag[3,i]
            P[2,3,i] = P[3,2,i] = Podiag[1,i]
            P[3,1,i] = P[1,3,i] = Podiag[2,i]
            P[:,:,i] = @views rotateTensorToVectorZ(P[:,:,i], B[:,i])
         end
      end

      rotate_tensor!()
      P
   end,
   :Panisotropy => function (meta, ids=UInt64[]) # P⟂ / P∥
      PR = readvariable(meta, "Protated", ids)
      Paniso = @. 0.5f0*(PR[1,1,:] + PR[2,2,:]) / PR[3,3,:]
      Paniso::Vector{Float32}
   end,
   :Pram => function (meta, ids=UInt64[])
      V = readvariable(meta, "Vmag", ids)
      ρm = readvariable(meta, "Rhom", ids)
      Pram = @. ρm * V * V
      Pram::Union{Vector{Float32}, Vector{Float64}}
   end,
   :Pb => function (meta, ids=UInt64[])
      μ2⁻¹ = 0.5/μ₀
      B = readvariable(meta, "vg_b_vol", ids)
      Pb = sum(x -> x*x*μ2⁻¹, B, dims=1) |> vec
      _fillinnerBC!(Pb, Pb)
      Pb::Vector{Float64}
   end,
   :Epar => function (meta, ids=UInt64[])
      E = readvariable(meta, "vg_e_vol", ids)
      b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)'

      @views [E[:,i] ⋅ b[:,i] for i in axes(E,2)]::Vector{Float32}
   end,
   :Eperp => function (meta, ids=UInt64[])
      E = readvariable(meta, "vg_e_vol", ids)
      b = readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)'

      Eperp = Vector{Float32}(undef, size(E, 2))

      function _computeEperp!()
         # Avoid sqrt of negative values, but does not guarantee orthogonality.
         @inbounds @views for i in eachindex(Eperp)
            eb = E[:,i] ⋅ b[:,i]
            Epar = SVector{3,Float32}(eb*b[1,i], eb*b[2,i], eb*b[3,i])
            Eperp[i] = norm(E[:,i] - Epar)
         end
      end

      _computeEperp!()
      Eperp
   end,
   :Poynting => function (meta, ids=UInt64[]) #WARN: real Poynting flux needs perturbed vars
      if hasvariable(meta, "vg_b_vol") && hasvariable(meta, "vg_e_vol")
         E = readvariable(meta, "vg_e_vol", ids)
         B = (meta, "vg_b_vol", ids)
      elseif hasvariable(meta, "fg_e") && hasvariable(meta, "fg_b")
         @assert isempty(ids) "Do not support selecting cells in FSGrid!"
         E = readvariable(meta, "fg_e")
         B = readvariable(meta, "fg_b")
      else
         E = readvariable(meta, "E")
         B = readvariable(meta, "B")
      end
      F = similar(E)

      function computecross!()
         Rpost = CartesianIndices(size(E)[2:end])
         @inbounds @views for i in Rpost
            F[:,i] = E[:,i] × B[:,i] ./ μ₀
         end
      end

      computecross!()

      F
   end,
   :Agyrotropy => function (meta, ids=UInt64[])
      # non-gyrotropy measure Q [Swisdak 2016]
      # Original version for electrons. Here we do protons first.

      Pdiag = readvariable(meta, "proton/vg_ptensor_diagonal", ids)
      Podiag = readvariable(meta, "proton/vg_ptensor_offdiagonal", ids)
      B = readvariable(meta, "vg_b_vol", ids)

      Pxx = selectdim(Pdiag, 1, 1)
      Pyy = selectdim(Pdiag, 1, 2)
      Pzz = selectdim(Pdiag, 1, 3)
      Pxy = selectdim(Podiag, 1, 3)
      Pyz = selectdim(Podiag, 1, 1)
      Pxz = selectdim(Podiag, 1 ,2)

      Bx = selectdim(B, 1, 1)
      By = selectdim(B, 1, 2)
      Bz = selectdim(B, 1, 3)

      Qsqr = Vector{Float32}(undef, size(Pdiag,2))

      @inbounds for ic in eachindex(Qsqr)
         I₁ = Pxx[ic] + Pyy[ic] + Pzz[ic]
         I₂ = Pxx[ic]*Pyy[ic] + Pxx[ic]*Pzz[ic] + Pyy[ic]*Pzz[ic] -
            Pxy[ic]*Pxy[ic] - Pyz[ic]*Pyz[ic] - Pxz[ic]*Pxz[ic]

         Ppar = (Bx[ic]*Bx[ic]*Pxx[ic] + By[ic]*By[ic]*Pyy[ic] + Bz[ic]*Bz[ic]*Pzz[ic] +
	         2*(Bx[ic]*By[ic]*Pxy[ic] + Bx[ic]*Bz[ic]*Pxz[ic] + By[ic]*Bz[ic]*Pyz[ic])) /
            (Bx[ic]^2 + By[ic]^2 + Bz[ic]^2)

         Qsqr[ic] = √(1 - 4I₂/((I₁ - Ppar)*(I₁ + 3Ppar)))
      end
      Qsqr
   end,
   :Beta => function (meta, ids=UInt64[])
      P = readvariable(meta, "P", ids)
      B = readvariable(meta, "vg_b_vol", ids)

      Beta = Vector{Float64}(undef, size(P))
      @inbounds for i in eachindex(Beta)
         Bmag = B[1,i]^2 + B[2,i]^2 + B[3,i]^2
         Beta[i] = (Bmag != 0) ? 2 * μ₀ * P[i] / Bmag : NaN
      end
      Beta
   end,
   :BetaStar => function (meta, ids=UInt64[])
      P = readvariable(meta, "P", ids)
      Pram = readvariable(meta, "Pram", ids)
      B = readvariable(meta, "vg_b_vol", ids)

      βstar = Vector{Float64}(undef, size(P))
      @inbounds for i in eachindex(βstar)
         Bmag = B[1,i]^2 + B[2,i]^2 + B[3,i]^2
         βstar[i] = (Bmag != 0) ? 2 * μ₀ * (P[i] + Pram[i]) / Bmag : NaN
      end
      βstar
   end,
   :IonInertial => function (meta, ids=UInt64[])
      n = readvariable(meta, "n", ids)
      Z = 1
      fi = @. √(n/(mᵢ*μ₀)) * Z * qᵢ / 2π
      di = @. c / fi
      di::Vector{Float64}
   end,
   :Larmor => function (meta, ids=UInt64[])
      Vth = readvariable(meta, "Vth", ids)
      B = readvariable(meta, "Bmag", ids)
      rg = @. mᵢ * Vth / (qᵢ * B)
      rg::Vector{Float64}
   end,
   :Gyroperiod => function (meta, ids=UInt64[])
      B = readvariable(meta, "Bmag", ids)
      T = @. 2π * mᵢ / (qᵢ * B)
      T::Vector{Float64}
   end,
   :Gyrofrequency => function (meta, ids=UInt64[]) # [1/s]
      B = readvariable(meta, "Bmag", ids)
      f = @. qᵢ * B / (mᵢ * 2π)
      f::Vector{Float64}
   end,
   :PlasmaPeriod => function (meta, ids=UInt64[])
      n = readvariable(meta, "n", ids)
      T = @. 2π / (qᵢ * √(n  / (mᵢ * ϵ₀)))
      T::Vector{Float64}
   end,
   :Omegap => function (meta, ids=UInt64[]) # plasma frequency, [1/s]
      n = readvariable(meta, "n", ids)
      fₚ = @. qᵢ * √(n  / (mᵢ * ϵ₀)) / 2π
      fₚ::Vector{Float64}
   end,
   :n => function (meta, ids=UInt64[])
      n = readvariable(meta, "proton/vg_rho", ids)
      n::Vector{Float32}
   end,
   :MagneticTension => function (meta, ids=UInt64[])
      B = readvariable(meta, "vg_b_vol")
      B = reshape(B, 3, meta.ncells...)

      dx2⁻¹ = @. inv(2*meta.dcoord)
      @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

      tension = similar(B)
      if ndims(meta) == 3
         # Warning: this works only for regular 3D B field
         @inbounds for k in 2:size(B,4)-1, j in 2:size(B,3)-1, i in 2:size(B,2)-1
            ∂Bx∂x = (-Bx[i-1,j  ,k  ] + Bx[i+1,j,  k  ]) * dx2⁻¹[1]
            ∂Bx∂y = (-Bx[i  ,j-1,k  ] + Bx[i,  j+1,k  ]) * dx2⁻¹[2]
            ∂Bx∂z = (-Bx[i  ,j  ,k-1] + Bx[i,  j,  k+1]) * dx2⁻¹[3]

            ∂By∂x = (-By[i-1,j  ,k  ] + By[i+1,j,  k  ]) * dx2⁻¹[1]
            ∂By∂y = (-By[i  ,j-1,k  ] + By[i,  j+1,k  ]) * dx2⁻¹[2]
            ∂By∂z = (-By[i  ,j  ,k-1] + By[i,  j,  k+1]) * dx2⁻¹[3]

            ∂Bz∂x = (-Bz[i-1,j  ,k  ] + Bz[i+1,j,  k  ]) * dx2⁻¹[1]
            ∂Bz∂y = (-Bz[i  ,j-1,k  ] + Bz[i,  j+1,k  ]) * dx2⁻¹[2]
            ∂Bz∂z = (-Bz[i  ,j  ,k-1] + Bz[i,  j,  k+1]) * dx2⁻¹[3]

            tension[1,i,j,k] = (Bx[i,j,k]*∂Bx∂x + By[i,j,k]*∂Bx∂y + Bz[i,j,k]*∂Bx∂z) / μ₀
            tension[2,i,j,k] = (Bx[i,j,k]*∂By∂x + By[i,j,k]*∂By∂y + Bz[i,j,k]*∂By∂z) / μ₀
            tension[3,i,j,k] = (Bx[i,j,k]*∂Bz∂x + By[i,j,k]*∂Bz∂y + Bz[i,j,k]*∂Bz∂z) / μ₀
         end
      elseif ndims(meta) == 2
         if meta.ncells[3] == 1
            @inbounds for j in 2:size(B,3)-1, i in 2:size(B,2)-1
               ∂Bx∂x = (-Bx[i-1,j  ,1] + Bx[i+1,j,  1]) * dx2⁻¹[1]
               ∂Bx∂y = (-Bx[i  ,j-1,1] + Bx[i,  j+1,1]) * dx2⁻¹[2]

               ∂By∂x = (-By[i-1,j  ,1] + By[i+1,j,  1]) * dx2⁻¹[1]
               ∂By∂y = (-By[i  ,j-1,1] + By[i,  j+1,1]) * dx2⁻¹[2]

               ∂Bz∂x = (-Bz[i-1,j  ,1] + Bz[i+1,j,  1]) * dx2⁻¹[1]
               ∂Bz∂y = (-Bz[i  ,j-1,1] + Bz[i,  j+1,1]) * dx2⁻¹[2]

               tension[1,i,j,1] = (Bx[i,j,1]*∂Bx∂x + By[i,j,1]*∂Bx∂y) / μ₀
               tension[2,i,j,1] = (Bx[i,j,1]*∂By∂x + By[i,j,1]*∂By∂y) / μ₀
               tension[3,i,j,1] = (Bx[i,j,1]*∂Bz∂x + By[i,j,1]*∂Bz∂y) / μ₀
            end
         elseif meta.ncells[2] == 1
            @inbounds for k in 2:size(B,4)-1, i in 2:size(B,2)-1
               ∂Bx∂x = (-Bx[i-1,1,k  ] + Bx[i+1,1,  k  ]) * dx2⁻¹[1]
               ∂Bx∂z = (-Bx[i  ,1,k-1] + Bx[i,  1,  k+1]) * dx2⁻¹[3]

               ∂By∂x = (-By[i-1,1,k  ] + By[i+1,1,  k  ]) * dx2⁻¹[1]
               ∂By∂z = (-By[i  ,1,k-1] + By[i,  1,  k+1]) * dx2⁻¹[3]

               ∂Bz∂x = (-Bz[i-1,1,k  ] + Bz[i+1,1,  k  ]) * dx2⁻¹[1]
               ∂Bz∂z = (-Bz[i  ,1,k-1] + Bz[i,  1,  k+1]) * dx2⁻¹[3]

               tension[1,i,1,k] = (Bx[i,1,k]*∂Bx∂x + Bz[i,1,k]*∂Bx∂z) / μ₀
               tension[2,i,1,k] = (Bx[i,1,k]*∂By∂x + Bz[i,1,k]*∂By∂z) / μ₀
               tension[3,i,1,k] = (Bx[i,1,k]*∂Bz∂x + Bz[i,1,k]*∂Bz∂z) / μ₀
            end
         end
      else
         @error "Magnetic tension requires 2D/3D B field!"
      end

      tension::Array{Float32, 4}
   end,
)

function _fillinnerBC!(data::AbstractArray{T}, dataRef::AbstractArray{T}) where
   T<:AbstractFloat
   @inbounds for i = eachindex(dataRef) # sparsity/inner boundary
      dataRef[i] == 0 && (data[i] = NaN)
   end
end