using LaTeXStrings
using LinearAlgebra: ×, norm

# Define units for intrinsic values
const units_predefined = Dict(
   "rhom" => "kg/m3",
   "rhoq" => "C/m3",
   "rho" => "1/m3",
   "rhobackstream" => "1/m3",
   "rhononbackstream" => "1/m3",
   "rho_v" => "1/m2s",
   "rhovbackstream" => "1/m2s",
   "rhovnonbackstream" => "1/m2s",
   "v" => "m/s",
   "vbackstream" => "m/s",
   "nNonbackstream" => "m/s",
   "b" => "T",
   "b_vol" => "T",
   "background_b" => "T",
   "perturbed_b" => "T",
   "bgb" => "T",
   "perb" => "T",
   "perb_vol" => "T",
   "e" => "V/m",
   "e_vol" => "V/m",
   "exhall_000_100" => "V/m",
   "exhall_001_101" => "V/m",
   "exhall_010_110" => "V/m",
   "exhall_011_111" => "V/m",
   "eyhall_000_010" => "V/m",
   "eyhall_001_011" => "V/m",
   "eyhall_100_110" => "V/m",
   "eyhall_101_111" => "V/m",
   "ezhall_000_001" => "V/m",
   "ezhall_010_011" => "V/m",
   "ezhall_100_101" => "V/m",
   "ezhall_110_111" => "V/m",
   "pressure" => "Pa",
   "pressure_dt2" => "Pa",
   "pressure_r" => "Pa",
   "pressure_v" => "Pa",
   "ptensordiagonal" => "Pa",
   "ptensoroffdiagonal" => "Pa",
   "ptensorbackstreamdiagonal" => "Pa",
   "ptensorbackstreamoffdiagonal" => "Pa",
   "ptensornonbackstreamdiagonal" => "Pa",
   "ptensornonbackstreamoffdiagonal" => "Pa",
   "max_v_dt" => "s",
   "max_r_dt" => "s",
   "max_fields_dt" => "s",
   "minvalue" => "s3/m6",
   "effectivesparsitythreshold" => "s3/m6",
   "rho_loss_adjust" => "1/m3",
   "energydensity" => "eV/cm3",
   "precipitationdiffflux" => "1/(cm2 sr s eV)"
)

# Define LaTeX markup names for intrinsic values
const latex_predefined = Dict(
   "rhom" => L"$\rho_m$",
   "rhoq" => L"$\rho_q$",
   "rho" => L"$n_\mathrm{p}$",
   "rhobackstream" => L"$n_\mathrm{p,st}$",
   "rhononbackstream" => L"$n_\mathrm{p,th}$",
   "rho_v" => L"$\Gamma_\mathrm{p}$",
   "rhovbackstream" => L"$\Gamma_\mathrm{p,st}$",
   "rhovnonbackstream" => L"$\Gamma_\mathrm{p,th}$",
   "v" => L"$V$",
   "vbackstream" => L"$V_\mathrm{p,st}$",
   "vnonbackstream" => L"$V_\mathrm{p,th}$",
   "b" => L"$B$",
   "b_vol" => L"$B_\mathrm{vol}$",
   "background_b" => L"$B_\mathrm{bg}$",
   "perturbed_b" => L"$B_\mathrm{pert}$",
   "bgb" => L"$B_\mathrm{bg}$",
   "perb" => L"B_\mathrm{pert}$",
   "perb_vol" => L"B_\mathrm{vol,pert}$",
   "e" => L"$E$",
   "e_vol" => L"$E_\mathrm{vol}$",
   "exhall_000_100" => L"$E_\mathrm{Hall,000,100}$",
   "exhall_001_101" => L"$E_\mathrm{Hall,001,101}$",
   "exhall_010_110" => L"$E_\mathrm{Hall,010,110}$",
   "exhall_011_111" => L"$E_\mathrm{Hall,011,111}$",
   "eyhall_000_010" => L"$E_\mathrm{Hall,000,010}$",
   "eyhall_001_011" => L"$E_\mathrm{Hall,001,011}$",
   "eyhall_100_110" => L"$E_\mathrm{Hall,100,110}$",
   "eyhall_101_111" => L"$E_\mathrm{Hall,101,111}$",
   "ezhall_000_001" => L"$E_\mathrm{Hall,000,001}$",
   "ezhall_010_011" => L"$E_\mathrm{Hall,010,011}$",
   "ezhall_100_101" => L"$E_\mathrm{Hall,100,101}$",
   "ezhall_110_111" => L"$E_\mathrm{Hall,110,111}$",
   "pressure" => L"$P$",
   "pressure_dt2" => L"$P_{\mathrm{d}t/2}}$",
   "pressure_r" => L"$P_r$",
   "pressure_v" => L"$P_v$",
   "ptensordiagonal" => L"$\mathcal{P}_\mathrm{diag}$",
   "ptensoroffdiagonal" => L"$\mathcal{P}_\mathrm{off-diag}$",
   "ptensorbackstreamdiagonal" => L"$\mathcal{P}_\mathrm{st,diag}$",
   "ptensorbackstreamoffdiagonal" => L"$\mathcal{P}_\mathrm{st,off-diag}$",
   "ptensornonbackstreamdiagonal" => L"$\mathcal{P}_\mathrm{th,diag}$",
   "ptensornonbackstreamoffdiagonal" => L"$\mathcal{P}_\mathrm{th,off-diag}$",
   "max_v_dt" => L"$\Delta t_{\mathrm{max},v}$",
   "max_r_dt" => L"$\Delta t_{\mathrm{max},r}$",
   "max_fields_dt" => L"$\Delta t_\mathrm{max,FS}$",
   "minvalue" => L"$f_\mathrm{Min}$",
   "effectivesparsitythreshold" => L"$f_\mathrm{Min}$",
   "rho_loss_adjust" => L"$\Delta_\mathrm{loss} n_\mathrm{p}$",
)


# Define LaTeX markup units for intrinsic values
const latexunits_predefined = Dict(
   "rhom" => L"$\mathrm{kg}\,\mathrm{m}^{-3}$",
   "rhoq" => L"$\mathrm{C}\,\mathrm{m}^{-3}$",
   "rho" => L"$\mathrm{m}^{-3}$",
   "rhobackstream" => L"$\mathrm{m}^{-3}$",
   "rhononbackstream" => L"$\mathrm{m}^{-3}$",
   "rho_v" => L"$\mathrm{m}^{-2}$s",
   "rhovbackstream" => L"$\mathrm{m}^{-2}$s",
   "rhovnonbackstream" => L"$\mathrm{m}^{-2}$s",
   "v" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "vbackstream" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "vnonbackstream" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "b" => L"T",
   "b_vol" => L"T",
   "background_b" => L"T",
   "perturbed_b" => L"T",
   "bgb" => L"T",
   "perb" => L"T",
   "perb_vol" => L"T",
   "e" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "e_vol" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_000_100" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_001_101" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_010_110" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_011_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_000_010" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_001_011" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_100_110" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_101_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_000_001" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_010_011" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_100_101" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_110_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "pressure" => L"Pa",
   "pressure_dt2" => L"Pa",
   "pressure_r" => L"Pa",
   "pressure_v" => L"Pa",
   "ptensordiagonal" => L"Pa",
   "ptensoroffdiagonal" => L"Pa",
   "ptensorbackstreamdiagonal" => L"Pa",
   "ptensorbackstreamoffdiagonal" => L"Pa",
   "ptensornonbackstreamdiagonal" => L"Pa",
   "ptensornonbackstreamoffdiagonal" => L"Pa",
   "max_v_dt" => L"s",
   "max_r_dt" => L"s",
   "max_fields_dt" => L"s",
   "minvalue" => L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$",
   "effectivesparsitythreshold" => L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$",
   "rho_loss_adjust" => L"$\mathrm{m}^{-3}$",
   "energydensity" => L"$\mathrm{eV}\,\mathrm{cm}^{-3}$",
   "precipitationdiffflux" => L"$\mathrm{cm}^{-2} \,\mathrm{sr}^{-1}\,\mathrm{s}^{-1}\,\mathrm{eV}^{-1}$"
)

# Define derived parameters
const variables_predefined = Dict(
   "bmag" => meta -> dropdims(sqrt.(sum(read_variable(meta, "vg_b_vol").^2, dims=1)), dims=1),
   "emag" => meta -> dropdims(sqrt.(sum(read_variable(meta, "vg_e_vol").^2, dims=1)), dims=1),
   "vmag" => meta -> dropdims(sqrt.(sum(read_variable(meta, "proton/vg_v").^2, dims=1)), dims=1),
   "vs" => function (meta) # sound speed
      Pdiag = read_variable(meta, "proton/vg_ptensor_diagonal")
      P = dropdims(sum(Pdiag, dims=1) ./ 3, dims=1)
      ρm = read_variable(meta, "proton/vg_rho") .* mᵢ
      # Handling sparse storage of the Vlasov solver
      for i = 1:length(ρm)
         ρm[i] == 0.0 && (ρm[i] = Inf)
      end
      vs = @. √( (P*5.0/3.0) / ρm )
   end,
   "va" => function (meta) # Alfvén speed
      ρm = read_variable(meta, "proton/vg_rho") .* mᵢ
      # Handling sparse storage of the Vlasov solver
      for i = 1:length(ρm)
         ρm[i] == 0.0 && (ρm[i] = Inf)
      end
      Bmag = get_variable_derived(meta, "bmag")
      VA = @. Bmag / √(ρm*μ₀)
   end,
   "MA" => function (meta) # Alfven Mach number
      V = read_variable(meta, "vmag")
      VA = get_variable_derived(meta, "va")
      VA ./ V 
   end,
   "upar" => function (meta) # Parallel velocity to the B field
      v = read_variable(meta, "proton/vg_v")
      B = read_variable(meta, "vg_b_vol")
      BmagInv = inv.(get_variable_derived(meta, "bmag"))
      [v[:,i] ⋅ (B[:,i] .* BmagInv[i]) for i in 1:size(v,2)]
   end,
   "uperp" => function (meta) # Perpendicular velocity to the B field
      v = read_variable(meta, "proton/vg_v")
      B = read_variable(meta, "vg_b_vol")
      BmagInv = inv.(get_variable_derived(meta, "bmag"))
      upar = [v[:,i] ⋅ (B[:,i] .* BmagInv[i]) for i in 1:size(v,2)]
      vmag2 = dropdims(sum(v.^2, dims=1), dims=1)
      uperp = @. √(vmag2 - upar^2) # This may be errorneous due to Float32!
   end,
   "PRotated" => function (meta)
      # Rotate the pressure tensor to align the z-component with the B field
      B = read_variable(meta, "vg_b_vol")
      Pdiag = read_variable(meta, "proton/vg_ptensor_diagonal")
      Podiag = read_variable(meta, "proton/vg_ptensor_offdiagonal")
      P = zeros(Float32, 3, 3, size(Pdiag, 2))
      @inbounds for i = 1:size(P, 2)
         P[1,1,i] = Pdiag[1,i]
         P[2,2,i] = Pdiag[2,i]
         P[3,3,i] = Pdiag[3,i]
         P[1,2,i] = P[2,1,i] = Podiag[1,i]
         P[2,3,i] = P[3,2,i] = Podiag[2,i]
         P[3,1,i] = P[1,3,i] = Podiag[3,i]
         rotateTensorToVectorZ!(P[:,:,i], B[:,i])
      end
      P
   end,
   "Anisotropy" => function (meta) # perpendicular / parallel component ratio
      P_rotated = get_variable_derived(meta, "PRotated")
      @. 0.5*(P_rotated[1,1,:] + P_rotated[2,2,:]) / P_rotated[3,3,:]
   end,
   "Pdynamic" => function (meta)
      vmag = get_variable_derived(meta, "vmag")
      ρm = read_variable(meta, "proton/vg_rho") .* mᵢ
      rhom.*Vmag.*Vmag
   end,
   "Poynting" => function (meta)
      if has_variable(meta.footer, "vg_b_vol")
         E = read_variable(meta, "vg_e_vol")
         B = read_variable(meta, "vg_b_vol")
      elseif has_variable(meta.footer, "B_vol")
         E = read_variable(meta, "E_vol")
         B = read_variable(meta, "B_vol")
      else
         E = read_variable(meta, "E")
         B = read_variable(meta, "B")
      end
      F = similar(E)
      @inbounds for i = 1:size(F,2)
         F[:,i] = E[:,i] × B[:,i] ./ μ₀
      end
   end,
   "aGyrotropy" => function (meta)
      # non-gyrotropy measure Q [Swisdak 2016]
      I₁ = @. Pxxe + Pyye + Pzze
      I₂ = @. Pxxe*Pyye + Pxxe*Pzze + Pyye*Pzze - Pxye*Pxye - Pyze*Pyze - Pxze*Pxze

      Ppar = @. (Bx*Bx*Pxxe + By*By*Pyye + Bz*Bz*Pzze +
	      2*(Bx*By*Pxye + Bx*Bz*Pxze + By*Bz*Pyze))/B²
      Qsqr = @. √(1 - 4I₂/((I₁ - Ppar)*(I₁ + 3Ppar)))
   end,
   "beta" => function (meta)
      Pressure = read_variable(meta, "pressure")
      Magneticfield = P = read_variable(meta, "B")   
      2.0 * μ₀ * Pressure / sum(Magneticfield^2)
   end,
   "ion_inertial" => function (meta)

   end,
   "larmor" => function (meta)

   end,
   "gyroperiod" => function (meta)

   end,
   "plasmaperiod" => function (meta)

   end,
)