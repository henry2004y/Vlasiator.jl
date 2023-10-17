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
   :e_hall => ("V/m", L"$E_\mathrm{Hall}$", L"$\mathrm{V}\,\mathrm{m}^{-1}$"),
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
   :Gyrofrequency => ("rad/s", L"$\omega_{g}$", L"\mathrm{rad}/\mathrm{s}"),
   :Omegap => ("rad/s", L"$\omega_{p}$", L"\mathrm{rad}/\mathrm{s}"),
)

# Define derived parameters
const variables_predefined = Dict(
   :Bmag => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      ρ = _readrho(meta, ids)
      B = readvariable(meta, "vg_b_vol", ids)::Array{Float32,2}
      Bmag = sum(x -> x*x, B, dims=1) |> vec
      @. Bmag = √(Bmag)
      _fillinnerBC!(Bmag, ρ)
      Bmag::Vector{Float32}
   end,
   :Emag => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      ρ = _readrho(meta, ids)
      E = readvariable(meta, "vg_e_vol", ids)::Array{Float32,2}
      Emag = sum(x -> x*x, E, dims=1) |> vec
      @. Emag = √(Emag)
      _fillinnerBC!(Emag, ρ)
      Emag::Vector{Float32}
   end,
   :Vmag => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      ρ = _readrho(meta, ids)
      V = readvariable(meta, "proton/vg_v", ids)::Array{Float32,2}
      Vmag = sum(x -> x*x, V, dims=1) |> vec
      @. Vmag = √(Vmag)
      _fillinnerBC!(Vmag, ρ)
      Vmag::Vector{Float32}
   end,
   :Bhat => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      ρ = _readrho(meta, ids)
      Bhat = readvariable(meta, "vg_b_vol", ids)::Array{Float32,2}
      Bmag = sum(x -> x*x, Bhat, dims=1) |> vec
      @. Bmag = √(Bmag)
      Bhat ./= Bmag'
      _fillinnerBC!(Bhat, ρ)
      Bhat::Matrix{Float32}
   end,
   :Rhom => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      if hasvariable(meta, "vg_rhom")
         ρm = readvariable(meta, "vg_rhom", ids)
      elseif hasvariable(meta, "proton/vg_rho")
         ρm = readvariable(meta, "proton/vg_rho", ids) .* mᵢ
      end
      ρm::Union{Vector{Float32}, Vector{Float64}}
   end,
   :P => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # scalar pressure
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
   :VS => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # sound speed
      P = readvariable(meta, "P", ids)::Vector{Float32}
      ρm = readvariable(meta, "Rhom", ids)::Vector{<:Real}
      _fillinnerBC!(ρm, ρm)
      VS = @. √( (P*5/3) / ρm )
      VS::Union{Vector{Float32}, Vector{Float64}}
   end,
   :VA => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # Alfvén speed
      ρm = readvariable(meta, "Rhom", ids)::Vector{<:Real}
      _fillinnerBC!(ρm, ρm)
      Bmag = readvariable(meta, "Bmag", ids)
      VA = @. Bmag / √(ρm*μ₀)
      VA::Vector{Float64}
   end,
   :MA => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # Alfvén Mach number
      V = readvariable(meta, "Vmag", ids)::Vector{Float32}
      VA = readvariable(meta, "VA", ids)::Vector{Float64}
      MA = V ./ VA
      MA::Vector{Float64}
   end,
   :MS => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # Sonic Mach number
      V = readvariable(meta, "Vmag", ids)::Vector{Float32}
      VS = readvariable(meta, "VS", ids)
      MS = V ./ VS
      MS::Union{Vector{Float32}, Vector{Float64}}
   end,
   :Vpar => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # velocity ∥ B
      V = readvariable(meta, "proton/vg_v", ids)::Array{Float32,2}
      b = readvariable(meta, "Bhat", ids)::Array{Float32,2}
      Vpar = @views [V[:,i] ⋅ b[:,i] for i in axes(V,2)]::Vector{Float32}
   end,
   :Vperp => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # velocity ⟂ B
      V = readvariable(meta, "proton/vg_v", ids)::Array{Float32,2}
      b = readvariable(meta, "Bhat", ids)::Array{Float32,2}
      Vperp = Vector{Float32}(undef, size(V, 2))

      # Avoid sqrt of negative values, but does not guarantee orthogonality.
      @inbounds @views for i in eachindex(Vperp)
         vb = V[:,i] ⋅ b[:,i]
         Vpar = SVector{3,Float32}(vb*b[1,i], vb*b[2,i], vb*b[3,i])
         Vperp[i] = norm(V[:,i] - Vpar)
      end

      Vperp
   end,
   :T => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # scalar temperature
      P = readvariable(meta, "P", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      T = @. P / (n*Float32(kB))
      T::Vector{Float32}
   end,
   :Vth => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # thermal velocity
      T = readvariable(meta, "T", ids)
      Vth = @. √(3 * Float32(kB) * T / Float32(mᵢ)) # assume proton
      Vth::Vector{Float32}
   end,
   :Ppar => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # P component ∥ B
      P = readvariable(meta, "Protated", ids)
      Ppar = P[3,3,:]::Vector{Float32}
   end,
   :Pperp => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # P component ⟂ B
      P = readvariable(meta, "Protated", ids)
      Pperp = [0.5f0(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]::Vector{Float32}
   end,
   :Tpar => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # T component ∥ B
      P = readvariable(meta, "Protated", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      Tpar = @. P[3,3,:] / (n*Float32(kB))
      Tpar::Vector{Float32}
   end,
   :Tperp => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # scalar T component ⟂ B
      P = readvariable(meta, "Protated", ids)
      n = readvariable(meta, "n", ids)
      _fillinnerBC!(n, n)
      Pperp = @inbounds [0.5f0(P[1,1,i] + P[2,2,i]) for i in 1:size(P,3)]
      Tperp = @. Pperp / (n*Float32(kB))
      Tperp::Vector{Float32}
   end,
   :Tanisotropy => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # T⟂ / T∥
      Tperp = readvariable(meta, "Tperp", ids)
      Tpar = readvariable(meta, "Tpar", ids)
      Taniso = @. Tperp / Tpar
      Taniso::Vector{Float32}
   end,
   :J => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      @assert isempty(ids) "Do not support current calculation for selected cells!"
      B = fillmesh(meta, ["vg_b_vol"]; maxamronly=true)[1][1][1]
      J = curl(B, meta.dcoord) ./ Float32(μ₀)
      J::Array{Float32, 4}
   end,
   :Jpar => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # velocity ∥ B
      B = fillmesh(meta, ["vg_b_vol"]; maxamronly=true)[1][1][1]
      J = curl(B, meta.dcoord)
      μ₀32 = Float32(μ₀)
      Jpar = @views [J[:,i,j,k] ⋅ B[:,i,j,k] / (μ₀32 * norm(B[:,i,j,k]))
         for i in axes(J,2), j in axes(J,3), k in axes(J,4)]::Array{Float32,3}
   end,
   :Jperp => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # velocity ⟂ B
      B = fillmesh(meta, ["vg_b_vol"]; maxamronly=true)[1][1][1]
      J = curl(B, meta.dcoord)

      μ₀32 = Float32(μ₀)
      Jperp = Array{Float32, 3}(undef, size(J)[2:end])
      # Avoid sqrt of negative values, but does not guarantee orthogonality.
      @inbounds @views for k in axes(J,4), j in axes(J,3), i in axes(J,2)
         b = SVector{3,Float32}(B[:,i,j,k] ./ norm(B[:,i,j,k]))
         jb = J[:,i,j,k] ⋅ b
         Jpar = SVector{3,Float32}(jb*b[1], jb*b[2], jb*b[3])
         Jperp[i,j,k] = norm(J[:,i,j,k] - Jpar) / μ₀32
      end

      Jperp::Array{Float32,3}
   end,
   :Protated => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
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
   :Panisotropy => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # P⟂ / P∥
      PR = readvariable(meta, "Protated", ids)
      Paniso = @. 0.5f0*(PR[1,1,:] + PR[2,2,:]) / PR[3,3,:]
      Paniso::Vector{Float32}
   end,
   :Pram => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      V = readvariable(meta, "Vmag", ids)
      ρm = readvariable(meta, "Rhom", ids)
      Pram = @. ρm * V * V
      Pram::Union{Vector{Float32}, Vector{Float64}}
   end,
   :Pb => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      μ2⁻¹ = 0.5/μ₀
      B = readvariable(meta, "vg_b_vol", ids)
      Pb = sum(x -> x*x*μ2⁻¹, B, dims=1) |> vec
      _fillinnerBC!(Pb, Pb)
      Pb::Vector{Float64}
   end,
   :Epar => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      E = readvariable(meta, "vg_e_vol", ids)::Array{Float32,2}
      b = readvariable(meta, "Bhat", ids)::Array{Float32,2}

      @views [E[:,i] ⋅ b[:,i] for i in axes(E,2)]::Vector{Float32}
   end,
   :Eperp => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      E = readvariable(meta, "vg_e_vol", ids)::Array{Float32,2}
      b = readvariable(meta, "Bhat", ids)::Array{Float32,2}

      Eperp = Vector{Float32}(undef, size(E, 2))

      # Avoid sqrt of negative values, but does not guarantee orthogonality.
      @inbounds @views for i in eachindex(Eperp)
         eb = E[:,i] ⋅ b[:,i]
         Epar = SVector{3,Float32}(eb*b[1,i], eb*b[2,i], eb*b[3,i])
         Eperp[i] = norm(E[:,i] - Epar)
      end

      Eperp
   end,
   :Ehallx => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      Ex1 = readvariable(meta, "fg_e_hall_0")::Array{Float32,3}
      Ex2 = readvariable(meta, "fg_e_hall_5")::Array{Float32,3}
      Ex3 = readvariable(meta, "fg_e_hall_8")::Array{Float32,3}
      Ex4 = readvariable(meta, "fg_e_hall_11")::Array{Float32,3}

      @. 0.25 * (Ex1 + Ex2 + Ex3 + Ex4)
   end,
   :Ehally => function (meta::MetaVLSV, ids::Vector{Int}=Int[])      
      Ey1 = readvariable(meta, "fg_e_hall_1")::Array{Float32,3}
      Ey2 = readvariable(meta, "fg_e_hall_3")::Array{Float32,3}
      Ey3 = readvariable(meta, "fg_e_hall_9")::Array{Float32,3}
      Ey4 = readvariable(meta, "fg_e_hall_10")::Array{Float32,3}

      @. 0.25 * (Ey1 + Ey2 + Ey3 + Ey4)
   end,
   :Ehallz => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      Ez1 = readvariable(meta, "fg_e_hall_2")::Array{Float32,3}
      Ez2 = readvariable(meta, "fg_e_hall_4")::Array{Float32,3}      
      Ez3 = readvariable(meta, "fg_e_hall_6")::Array{Float32,3}
      Ez4 = readvariable(meta, "fg_e_hall_7")::Array{Float32,3}

      @. 0.25 * (Ez1 + Ez2 + Ez3 + Ez4)
   end,
   :Poynting => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) #WARN: real Poynting flux needs perturbed vars
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
   :Agyrotropy => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
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
   :Beta => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      P = readvariable(meta, "P", ids)
      B = readvariable(meta, "vg_b_vol", ids)

      μ₀32 = Float32(μ₀)
      Beta = Vector{Float32}(undef, size(P))
      @inbounds for i in eachindex(Beta)
         Bmag = B[1,i]^2 + B[2,i]^2 + B[3,i]^2
         Beta[i] = (Bmag != 0) ? 2 * μ₀32 * P[i] / Bmag : NaN32
      end
      Beta::Vector{Float32}
   end,
   :BetaStar => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      P = readvariable(meta, "P", ids)
      Pram = readvariable(meta, "Pram", ids)
      B = readvariable(meta, "vg_b_vol", ids)

      μ₀32 = Float32(μ₀)
      βstar = Vector{Float32}(undef, size(P))
      @inbounds for i in eachindex(βstar)
         Bmag = B[1,i]^2 + B[2,i]^2 + B[3,i]^2
         βstar[i] = (Bmag != 0) ? 2 * μ₀32 * (P[i] + Pram[i]) / Bmag : NaN32
      end
      βstar::Vector{Float32}
   end,
   :IonInertial => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      n = readvariable(meta, "n", ids)
      Z = 1
      ωᵢ = @. √(n) * Z * Float32(qᵢ / √(mᵢ*μ₀))
      di = @. Float32(c) / ωᵢ
      di::Vector{Float32}
   end,
   :Larmor => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      Vth = readvariable(meta, "Vth", ids)
      B = readvariable(meta, "Bmag", ids)
      rg = @. Float32(mᵢ) * Vth / (Float32(qᵢ) * B)
      rg::Vector{Float32}
   end,
   :Gyroperiod => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      B = readvariable(meta, "Bmag", ids)
      T = @. Float32(2π * mᵢ / qᵢ) / B
      T::Vector{Float32}
   end,
   :Gyrofrequency => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # [1/s]
      B = readvariable(meta, "Bmag", ids)
      f = @. Float32(qᵢ  / (mᵢ * 2π)) * B
      f::Vector{Float32}
   end,
   :PlasmaPeriod => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      n = readvariable(meta, "n", ids)
      T = @. Float32(2π / qᵢ * √(mᵢ * ϵ₀)) / √n
      T::Vector{Float32}
   end,
   :Omegap => function (meta::MetaVLSV, ids::Vector{Int}=Int[]) # plasma frequency, [1/s]
      n = readvariable(meta, "n", ids)
      ωₚ = @. Float32(qᵢ / √(mᵢ * ϵ₀)) * √n
      ωₚ::Vector{Float32}
   end,
   :n => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
      n = readvariable(meta, "proton/vg_rho", ids)
      n::Vector{Float32}
   end,
   :MagneticTension => function (meta::MetaVLSV, ids::Vector{Int}=Int[])
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

function _readrho(meta::MetaVLSV, ids::Vector{Int}=Int[])
   ρ_ = findfirst(endswith("rho"), meta.variable)
   ρ = readvariable(meta, meta.variable[ρ_], ids)::Vector{Float32}
end

function _fillinnerBC!(data::AbstractArray{T}, dataRef::AbstractArray{T}) where
   T<:AbstractFloat
   @inbounds for i in eachindex(dataRef) # sparsity/inner boundary
      dataRef[i] == 0 && (data[i] = NaN)
   end
end

function _fillinnerBC!(data::AbstractMatrix{T}, dataRef::AbstractVector{T}) where
   T<:AbstractFloat
   @inbounds for i in eachindex(dataRef) # sparsity/inner boundary
      dataRef[i] == 0 && (data[:,i] .= NaN)
   end
end