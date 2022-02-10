# Physical parameter monitoring.

"""
    check_plasma_characteristics(n, v, T, B)

Display characteristic plasma parameters given density `n`, bulk velocity `v`, temperature
`T`, and magnetic field strength `B` in SI units.
"""
function check_plasma_characteristics(n, v, T, B)

   γ = 5/3    # adiabatic index

   p = n*kB*T # [Pa]

   # Characteristic parameters
   di  = √(mᵢ*ϵ₀/n)*c/qᵢ           # ion inertial length, [m]
   ωci = qᵢ*B/mᵢ                   # [/s]
   vA  = B / √(μ₀ * n * mᵢ)        # Alfven speed, [m/s]
   vS  = √(γ * p / (n * mᵢ))       # sonic speed, [m/s]
   vT  = √(kB * T / mᵢ)            # thermal speed, [m/s]
   rᵢ  = vT / ωci

   println("--------------------------------------------------")
   println("* Characteristic plasma properties")
   println("Density             : ", rpad(round(n/1e6; digits=2), 9), "amu/cc")
   println("Velocity            : ", rpad(round(v/1e3; digits=2), 9), "km/s")
   println("Pressure            : ", rpad(round(p*1e9; digits=3), 9), "nPa")
   println("Temperature         : ", rpad(round(T; digits=3), 9), "K")
   println("Magnetic field      : ", rpad(round(B*1e9; digits=2), 9), "nT")
   println("Ion inertial length : ", rpad(round(di/1e3; digits=2),9), "km")
   println("Ion gyroradius      : ", rpad(round(rᵢ/1e3; digits=2),9), "km")
   println("Gyrofrequency       : ", rpad(round(ωci; digits=2),   9), "Hz")
   println("Alfvén speed        : ", rpad(round(vA/1e3; digits=2),9), "km/s")
   println("Sonic speed         : ", rpad(round(vS/1e3; digits=2),9), "km/s")
   println("--------------------------------------------------")

   return true # for unit testing
end