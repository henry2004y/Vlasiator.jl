# 2D magnetic flux function.

"""
    compute_flux_function(b, Δ, nG=2)

Calculate the 2D magnetic flux function ψ from the magnetic field `b` and discrete steps
`Δ`. `nG` is the number of ghost cells along each dimension in the vector field.
ψ is defined as
``\\psi = \\int B_x dz = - \\int B_z dx`` from ``\\mathbf{B} = \\hat{y}\\times\\nabla\\psi``
in the X-Z plane and Y is the out-of-plane direction. This is strictly true if the guide
field By is zero and B is divergence-free. However, numerically there will be errors.
The current implementation calculates ψ by integrating along -z boundary first, and then
going along z.
Reference: https://doi.org/10.1063/1.3657424
"""
function compute_flux_function(b::AbstractArray{T,N}, Δ, nG=2) where {T,N}
   dx, dz = Δ
   flux = zeros(T, size(b,2), size(b,3))

   # Find ψ along the -z boundary cells
   @inbounds for i in axes(b,2)[1+nG:end-nG]
      bz = b[3,i,1+nG]
      flux[i,1+nG] = flux[i-1,1+nG] - bz * dx
   end

   # For each row, integrate in z-direction
   @inbounds for j in axes(b,3)[2+nG:end-nG], i in axes(b,2)[1+nG:end-nG]
      bx = b[1,i,j]
      flux[i,j] = flux[i,j-1] + bx * dz
   end

   flux
end

"""
    find_reconnection_points(ψ, retol=1e-3) -> indices_x, indices_o

Find X-point and O-point indices in 2D magnetic field topology from flux function `ψ`.
`retol` determines the ratio w.r.t. |∇ψ|² to accept a gradient as 0. The current
implementation does not work for the 2 layers near the boundary.
"""
function find_reconnection_points(ψ, retol=1e-3)
   ∂ψ = gradient(ψ)
   ∂²ψ = gradient(view(∂ψ,1,:,:))
   # fill boundary layers
   ∂²ψ[:,2,:] = ∂²ψ[:,1,:]
   ∂²ψ[:,end-1,:] = ∂²ψ[:,end,:]
   ∂²ψ[:,:,2] = ∂²ψ[:,:,1]
   ∂²ψ[:,:,end-1] = ∂²ψ[:,:,end]
   # 2nd derivatives
   ∂²ψ∂x² = ∂²ψ[1,:,:]
   ∂²ψ∂x∂y = ∂²ψ[2,:,:] # == ∂²ψ∂y∂x = ∂²ψ[1,:,:]

   ∂²ψ = gradient(view(∂ψ,2,:,:))
   # fill boundary layers
   ∂²ψ[:,2,:] = ∂²ψ[:,1,:]
   ∂²ψ[:,end-1,:] = ∂²ψ[:,end,:]
   ∂²ψ[:,:,2] = ∂²ψ[:,:,1]
   ∂²ψ[:,:,end-1] = ∂²ψ[:,:,end]

   ∂²ψ∂y² = ∂²ψ[2,:,:]

   indices_x = Matrix{Int64}(undef, 2, 0)
   indices_o = Matrix{Int64}(undef, 2, 0)

   ∂ψmag² = [∂ψ[1,i,j]^2 + ∂ψ[2,i,j]^2 for j in axes(∂ψ,3), i in axes(∂ψ,2)]
   ∂ψmean = mean(∂ψmag²)

   for j in axes(∂ψ,3)[2:end-1], i in axes(∂ψ,2)[2:end-1]
      ∂ψ[1,i,j]^2 + ∂ψ[2,i,j]^2 > retol*∂ψmean && continue
      # Hessian matrix det(H) < 0 => saddle point (x-point)
      if ∂²ψ∂x²[i,j] * ∂²ψ∂y²[i,j] - ∂²ψ∂x∂y[i,j]^2 < 0
         indices_x = hcat(indices_x, [i,j])
      else # o-point
         indices_o = hcat(indices_o, [i,j])
      end
   end

   return indices_x, indices_o
end