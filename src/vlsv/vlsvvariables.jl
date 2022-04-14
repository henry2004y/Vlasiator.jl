# Predefined quantities

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
   :Beta => ("",L"$\beta$", "")
#   :Gyroperiod => ("s", L"$T_{gyro}$", "s"),
#   :Plasmaperiod => ("s", L"$T_{plasma}$", "s"),
#   :Gyrofrequency => ("/s", L"$\omega_{g}$", L"s^{-1}"),
#   :Omegap => ("/s", L"$\omega_{p}$", L"s^{-1}"),
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
            sqrt.(sum(x -> x*x, readvariable(meta, "vg_b_vol"), dims=1)) :
            sqrt.(sum(x -> x*x, readvariable(meta, "vg_b_vol", ids), dims=1))
      else
         @assert isempty(ids) "Do not support reading selected cells from FSGrid!"
         Bmag = sqrt.(sum(x -> x*x, readvariable(meta, "fg_b"), dims=1))
      end
      _fillinnerBC!(Bmag, ρ)
      Bmag
   end,
   :Emag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = isempty(ids) ?
         readvariable(meta, meta.variable[rho_]) :
         readvariable(meta, meta.variable[rho_], ids)
      if hasvariable(meta, "vg_e_vol")
         Emag = isempty(ids) ?
            sqrt.(sum(x -> x*x, readvariable(meta, "vg_e_vol"), dims=1)) :
            sqrt.(sum(x -> x*x, readvariable(meta, "vg_e_vol", ids), dims=1))
      else
         @assert isempty(ids) "Do not support reading selected cells from FSGrid!"
         Emag = sqrt.(sum(x -> x*x, readvariable(meta, "fg_e"), dims=1))
      end
      _fillinnerBC!(Emag, ρ)
      Emag
   end,
   :Vmag => function (meta, ids=UInt64[])
      rho_ = findfirst(endswith("rho"), meta.variable)
      ρ = isempty(ids) ?
         readvariable(meta, meta.variable[rho_]) :
         readvariable(meta, meta.variable[rho_], ids)
      Vmag = isempty(ids) ?
         vec(sqrt.(sum(x -> x*x, readvariable(meta, "proton/vg_v"), dims=1))) :
         vec(sqrt.(sum(x -> x*x, readvariable(meta, "proton/vg_v", ids), dims=1)))
      _fillinnerBC!(Vmag, ρ)
      Vmag
   end,
   :Rhom => function (meta, ids=UInt64[])
      if hasvariable(meta, "vg_rhom")
         ρm = isempty(ids) ?
            readvariable(meta, "vg_rhom") :
            readvariable(meta, "vg_rhom", ids) |> vec
      elseif hasvariable(meta, "proton/vg_rho")
         ρm = isempty(ids) ?
            readvariable(meta, "proton/vg_rho") .* mᵢ :
            readvariable(meta, "proton/vg_rho", ids) .* mᵢ |> vec
      end
      ρm
   end,
   :P => function (meta, ids=UInt64[]) # scalar pressure
      if hasvariable(meta, "vg_pressure")
         P = isempty(ids) ?
            readvariable(meta, "vg_pressure") :
            readvariable(meta, "vg_pressure", ids)
      else
         if hasvariable(meta, "proton/vg_ptensor_diagonal")
            Pdiag_str = "proton/vg_ptensor_diagonal"
         else
            Pdiag_str = "PTensorDiagonal"
         end
         Pdiag = isempty(ids) ?
            readvariable(meta, Pdiag_str) :
            readvariable(meta, Pdiag_str, ids)
         P = vec(mean(Pdiag, dims=1))
      end
      P
   end,
   :VS => function (meta, ids=UInt64[]) # sound speed
      P = readvariable(meta, "P", ids)
      ρm = readvariable(meta, "Rhom", ids)
      _fillinnerBC!(ρm, ρm)
      VS = @. √( (P*5.0f0/3.0f0) / ρm )
   end,
   :VA => function (meta, ids=UInt64[]) # Alfvén speed
      ρm = readvariable(meta, :Rhom, ids)
      _fillinnerBC!(ρm, ρm)
      Bmag = readvariable(meta, "Bmag", ids) |> vec
      VA = @. Bmag / √(ρm*μ₀)
   end,
   :MA => function (meta, ids=UInt64[]) # Alfvén Mach number
      V = readvariable(meta, "Vmag", ids)
      _fillinnerBC!(V, V)
      VA = readvariable(meta, "VA", ids)
      V ./ VA
   end,
   :MS => function (meta, ids=UInt64[]) # Sonic Mach number
      V = readvariable(meta, "Vmag", ids)
      _fillinnerBC!(V, V)
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
      [V[:,i] ⋅ b[:,i] for i in axes(V,2)]
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

      function _computeVperp!()
         # Avoid sqrt of negative values, but does not guarantee orthogonality.
         @inbounds for i in eachindex(Vperp)
            Vpar = V[:,i] ⋅ b[:,i] .* b[:,i]
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
   end,
   :Vth => function (meta, ids=UInt64[]) # thermal velocity
      T = readvariable(meta, "T", ids)
      @. √(3 * kB * T / mᵢ) # assume proton
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
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      @. P[3,3,:] / (n*kB)
   end,
   :Tperp => function (meta, ids=UInt64[]) # scalar T component ⟂ B
      P = readvariable(meta, "Protated", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
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
      J = curl(B, meta.dcoord) ./ μ₀
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

      function rotate_tensor!()
         @inbounds for i = 1:size(P, 3)
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
      @. 0.5f0*(PR[1,1,:] + PR[2,2,:]) / PR[3,3,:]
   end,
   :Pram => function (meta, ids=UInt64[])
      V = readvariable(meta, "Vmag", ids)
      ρm = readvariable(meta, "Rhom", ids)
      @. ρm * V * V
   end,
   :Pb => function (meta, ids=UInt64[])
      B² = readvariable(meta, "B²", ids)
      _fillinnerBC!(B², B²)
      mu2inv = 0.5/μ₀
      @. B² * mu2inv
   end,
   :Epar => function (meta, ids=UInt64[])
      E, b =
         if isempty(ids)
            readvariable(meta, "vg_e_vol"),
            readvariable(meta, "vg_b_vol") ./ readvariable(meta, "Bmag")
         else
            readvariable(meta, "vg_e_vol", ids),
            readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)
         end
      [E[:,i] ⋅ b[:,i] for i in axes(E,2)]
   end,
   :Eperp => function (meta, ids=UInt64[])
      E, b =
         if isempty(ids)
            readvariable(meta, "vg_e_vol"),
            readvariable(meta, "vg_b_vol") ./ readvariable(meta, "Bmag")
         else
            readvariable(meta, "vg_e_vol", ids),
            readvariable(meta, "vg_b_vol", ids) ./ readvariable(meta, "Bmag", ids)
         end
      Eperp = zeros(eltype(E), size(E, 2))

      function _computeEperp!()
         # Avoid sqrt of negative values, but does not guarantee orthogonality.
         @inbounds for i in eachindex(Eperp)
            Epar = E[:,i] ⋅ b[:,i] .* b[:,i]
            Eperp[i] = norm(E[:,i] - Epar)
         end
      end

      _computeEperp!()
      Eperp
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

      function computecross!()
         Rpost = CartesianIndices(size(E)[2:end])
         @inbounds for i in Rpost
            F[:,i] = E[:,i] × B[:,i] ./ μ₀
         end
      end

      computecross!()

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
      Pxy = selectdim(Podiag, 1, 3)
      Pyz = selectdim(Podiag, 1, 1)
      Pxz = selectdim(Podiag, 1 ,2)

      Bx = selectdim(B, 1, 1)
      By = selectdim(B, 1, 2)
      Bz = selectdim(B, 1, 3)

      Qsqr = Vector{eltype(Pdiag)}(undef, size(Pdiag,2))

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
      B² = readvariable(meta, "B²", ids)
      _fillinnerBC!(B², B²)
      @. 2 * μ₀ * P / B²
   end,
   :BetaStar => function (meta, ids=UInt64[])
      P = readvariable(meta, "P", ids)
      Pram = readvariable(meta, "Pram", ids)
      B² = readvariable(meta, "B²", ids)
      _fillinnerBC!(B², B²)
      @. 2 * μ₀ * (P + Pram) / B²
   end,
   :IonInertial => function (meta, ids=UInt64[])
      n = readvariable(meta, "n", ids)
      Z = 1
      fi = @. √(n/(mᵢ*μ₀)) * Z * qᵢ / 2π
      di = @. c / fi
   end,
   :Larmor => function (meta, ids=UInt64[])
      Vth = readvariable(meta, "Vth", ids)
      B = readvariable(meta, "Bmag", ids)
      rg = @. mᵢ * Vth / (qᵢ * B)
   end,
   :Gyroperiod => function (meta, ids=UInt64[])
      B = readvariable(meta, "Bmag", ids)
      T = @. 2π * mᵢ / (qᵢ * B)
   end,
   :Gyrofrequency => function (meta, ids=UInt64[]) # [1/s]
      B = readvariable(meta, "Bmag", ids)
      f = @. qᵢ * B / (mᵢ * 2π)
   end,
   :Plasmaperiod => function (meta, ids=UInt64[])
      n = readvariable(meta, "n", ids)
      T = @. 2π / (qᵢ * √(n  / (mᵢ * ϵ₀)))
   end,
   :Omegap => function (meta, ids=UInt64[]) # plasma frequency, [1/s]
      n = readvariable(meta, "n", ids)
      fₚ = @. qᵢ * √(n  / (mᵢ * ϵ₀)) / 2π
   end,
   :n => function (meta, ids=UInt64[])
      n = isempty(ids) ?
         readvariable(meta, "proton/vg_rho") :
         readvariable(meta, "proton/vg_rho", ids) |> vec
   end,
   :B² => function (meta, ids=UInt64[])
      B² = isempty(ids) ?
         vec(sum(x -> x*x, readvariable(meta, "vg_b_vol"), dims=1)) :
         vec(sum(x -> x*x, readvariable(meta, "vg_b_vol", ids), dims=1))
   end,
   :MagneticTension => function (meta, ids=UInt64[])
      B = readvariable(meta, "vg_b_vol")
      B = reshape(B, 3, meta.ncells...)
   
      invdx = @. inv(2*meta.dcoord)
      @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

      tension = zeros(eltype(B), size(B))
      if ndims(meta) == 3
         # Warning: this works only for regular 3D B field
         @inbounds for k in 2:size(B,4)-1, j in 2:size(B,3)-1, i in 2:size(B,2)-1
            ∂Bx∂x = (-Bx[i-1,j  ,k  ] + Bx[i+1,j,  k  ]) * invdx[1]
            ∂Bx∂y = (-Bx[i  ,j-1,k  ] + Bx[i,  j+1,k  ]) * invdx[2]
            ∂Bx∂z = (-Bx[i  ,j  ,k-1] + Bx[i,  j,  k+1]) * invdx[3]

            ∂By∂x = (-By[i-1,j  ,k  ] + By[i+1,j,  k  ]) * invdx[1]
            ∂By∂y = (-By[i  ,j-1,k  ] + By[i,  j+1,k  ]) * invdx[2]
            ∂By∂z = (-By[i  ,j  ,k-1] + By[i,  j,  k+1]) * invdx[3]

            ∂Bz∂x = (-Bz[i-1,j  ,k  ] + Bz[i+1,j,  k  ]) * invdx[1]
            ∂Bz∂y = (-Bz[i  ,j-1,k  ] + Bz[i,  j+1,k  ]) * invdx[2]
            ∂Bz∂z = (-Bz[i  ,j  ,k-1] + Bz[i,  j,  k+1]) * invdx[3]

            tension[1,i,j,k] = (Bx[i,j,k]*∂Bx∂x + By[i,j,k]*∂Bx∂y + Bz[i,j,k]*∂Bx∂z) / μ₀
            tension[2,i,j,k] = (Bx[i,j,k]*∂By∂x + By[i,j,k]*∂By∂y + Bz[i,j,k]*∂By∂z) / μ₀
            tension[3,i,j,k] = (Bx[i,j,k]*∂Bz∂x + By[i,j,k]*∂Bz∂y + Bz[i,j,k]*∂Bz∂z) / μ₀
         end
      elseif ndims(meta) == 2
         if meta.ncells[3] == 1
            @inbounds for j in 2:size(B,3)-1, i in 2:size(B,2)-1
               ∂Bx∂x = (-Bx[i-1,j  ,1] + Bx[i+1,j,  1]) * invdx[1]
               ∂Bx∂y = (-Bx[i  ,j-1,1] + Bx[i,  j+1,1]) * invdx[2]
   
               ∂By∂x = (-By[i-1,j  ,1] + By[i+1,j,  1]) * invdx[1]
               ∂By∂y = (-By[i  ,j-1,1] + By[i,  j+1,1]) * invdx[2]
   
               ∂Bz∂x = (-Bz[i-1,j  ,1] + Bz[i+1,j,  1]) * invdx[1]
               ∂Bz∂y = (-Bz[i  ,j-1,1] + Bz[i,  j+1,1]) * invdx[2]
   
               tension[1,i,j,1] = (Bx[i,j,1]*∂Bx∂x + By[i,j,1]*∂Bx∂y) / μ₀
               tension[2,i,j,1] = (Bx[i,j,1]*∂By∂x + By[i,j,1]*∂By∂y) / μ₀
               tension[3,i,j,1] = (Bx[i,j,1]*∂Bz∂x + By[i,j,1]*∂Bz∂y) / μ₀
            end
         elseif meta.ncells[2] == 1
            @inbounds for k in 2:size(B,4)-1, i in 2:size(B,2)-1
               ∂Bx∂x = (-Bx[i-1,1,k  ] + Bx[i+1,1,  k  ]) * invdx[1]
               ∂Bx∂z = (-Bx[i  ,1,k-1] + Bx[i,  1,  k+1]) * invdx[3]
   
               ∂By∂x = (-By[i-1,1,k  ] + By[i+1,1,  k  ]) * invdx[1]
               ∂By∂z = (-By[i  ,1,k-1] + By[i,  1,  k+1]) * invdx[3]
   
               ∂Bz∂x = (-Bz[i-1,1,k  ] + Bz[i+1,1,  k  ]) * invdx[1]
               ∂Bz∂z = (-Bz[i  ,1,k-1] + Bz[i,  1,  k+1]) * invdx[3]
   
               tension[1,i,1,k] = (Bx[i,1,k]*∂Bx∂x + Bz[i,1,k]*∂Bx∂z) / μ₀
               tension[2,i,1,k] = (Bx[i,1,k]*∂By∂x + Bz[i,1,k]*∂By∂z) / μ₀
               tension[3,i,1,k] = (Bx[i,1,k]*∂Bz∂x + Bz[i,1,k]*∂Bz∂z) / μ₀
            end
         end
      else
         @error "Magnetic tension requires 2D/3D B field!"
      end

      tension
   end,
)



function _fillinnerBC!(data, dataRef)
   @inbounds for i = eachindex(dataRef) # sparsity/inner boundary
      dataRef[i] == 0.0 && (data[i] = NaN)
   end
end