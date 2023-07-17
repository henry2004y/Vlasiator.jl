# 2D magnetic flux function.

"""
    compute_flux_function(b::AbstractArray{T,N}, Δ::Vector{T}, nG::Int=2) where {T,N}

Calculate the 2D magnetic flux function ψ from the magnetic field `b` and discrete steps
`Δ`. `nG` is the number of ghost cells along each dimension in the vector field.
ψ is defined as
``\\psi = \\int B_x dz = - \\int B_z dx`` from ``\\mathbf{B} = \\hat{y}\\times\\nabla\\psi``
in the X-Z plane and Y is the out-of-plane direction. This is strictly true if B is
divergence-free and the guide field By is constant. However, numerically there will be
errors. The current implementation calculates ψ by integrating along -z boundary first,
and then going along z.
Reference: [Flux function](https://doi.org/10.1063/1.3657424)
"""
function compute_flux_function(b::AbstractArray{T,N}, Δ::Vector{T}, nG::Int=2) where {T,N}
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
    find_reconnection_points(ψ::Array{T,2}; retol::Float64=1e-4,
       method::Int=1) -> indices_x, indices_o

Find X-point and O-point indices in 2D magnetic field topology from flux function `ψ`.
The current implementation does not work for the 3 layers near the boundary.

# Keywords
- `retol=1e-1`: determines the relative tolerance of the ratio w.r.t. |∇ψ|² to accept a gradient as 0.
- `method=1`: method 1 compute the cell-centered 1st and 2nd order derivatives and check the Hessian matrix; method 2 check the flux function at each point against its 8 neighbors, which is more deterministic.
"""
function find_reconnection_points(ψ::Array{T,2}; retol::Float64=1e-4, method::Int=1) where T
   indices_x = Matrix{Int}(undef, 2, 0)
   indices_o = Matrix{Int}(undef, 2, 0)

   if method == 1
      dx = ones(eltype(ψ), 2)
      ∂ψ = gradient(ψ, dx)
      ∂²ψ = gradient(view(∂ψ,1,:,:), dx)
      # fill boundary layers
      ∂²ψ[:,2,:] = ∂²ψ[:,1,:]
      ∂²ψ[:,end-1,:] = ∂²ψ[:,end,:]
      ∂²ψ[:,:,2] = ∂²ψ[:,:,1]
      ∂²ψ[:,:,end-1] = ∂²ψ[:,:,end]
      # 2nd derivatives
      ∂²ψ∂x² = ∂²ψ[1,:,:]
      ∂²ψ∂x∂y = ∂²ψ[2,:,:] # == ∂²ψ∂y∂x = ∂²ψ[1,:,:]

      ∂²ψ = gradient(view(∂ψ,2,:,:), dx)
      # fill boundary layers
      ∂²ψ[:,2,:] = ∂²ψ[:,1,:]
      ∂²ψ[:,end-1,:] = ∂²ψ[:,end,:]
      ∂²ψ[:,:,2] = ∂²ψ[:,:,1]
      ∂²ψ[:,:,end-1] = ∂²ψ[:,:,end]

      ∂²ψ∂y² = ∂²ψ[2,:,:]

      ∂ψmag² = [∂ψ[1,i,j]^2 + ∂ψ[2,i,j]^2 for j in axes(∂ψ,3), i in axes(∂ψ,2)]
      ∂ψmean = mean(∂ψmag²)

      for j in axes(∂ψ,3)[4:end-3], i in axes(∂ψ,2)[4:end-3]
         ∂ψ[1,i,j]^2 + ∂ψ[2,i,j]^2 > retol*∂ψmean && continue
         # Hessian matrix det(H) < 0 => saddle point (X-point)
         if ∂²ψ∂x²[i,j] * ∂²ψ∂y²[i,j] - ∂²ψ∂x∂y[i,j]^2 < 0
            indices_x = hcat(indices_x, [i,j])
         else # O-point
            indices_o = hcat(indices_o, [i,j])
         end
      end
   elseif method == 2
      for j in axes(ψ,2)[4:end-3], i in axes(ψ,1)[4:end-3]
         if issaddle(@views ψ[i-1:i+1,j-1:j+1])
            indices_x = hcat(indices_x, [i,j])
         elseif isextrema(@views ψ[i-1:i+1,j-1:j+1])
            indices_o = hcat(indices_o, [i,j])
         end
      end
   end

   return indices_x, indices_o
end

"Check if the center point in a 3x3 matrix `ψ` is a saddle point."
function issaddle(ψ)
   nmax, nmin = 0, 0
   minmax1 = extrema(@views ψ[:,2])
   minmax2 = extrema(@views ψ[2,:])
   minmax3 = extrema([ψ[1,1],ψ[2,2],ψ[3,3]])
   minmax4 = extrema([ψ[1,3],ψ[2,2],ψ[3,1]])
   if ψ[2,2] == minmax1[1]
      nmin += 1
   elseif ψ[2,2] == minmax1[2]
      nmax += 1
   end
   if ψ[2,2] == minmax2[1]
      nmin += 1
   elseif ψ[2,2] == minmax2[2]
      nmax += 1
   end
   if ψ[2,2] == minmax3[1]
      nmin += 1
   elseif ψ[2,2] == minmax3[2]
      nmax += 1
   end
   if ψ[2,2] == minmax4[1]
      nmin += 1
   elseif ψ[2,2] == minmax4[2]
      nmax += 1
   end
   if nmin > 0 && nmax > 0
      return true
   else
      return false
   end
end

"Check if the center point in a 3x3 matrix `ψ` is an extrema point."
function isextrema(ψ)
   if ψ[2,2] == maximum(ψ) || ψ[2,2] == minimum(ψ)
      return true
   else
      return false
   end
end
