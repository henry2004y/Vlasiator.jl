"""
    fg_grad(dx::AbstractVector, A::AbstractArray)

Calculate 2nd order cell-centered ∇A where `A` is a 4D array of size (nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function fg_grad(dx::AbstractVector, A::AbstractArray{T,N}) where {T,N}
   @assert N == 3 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   B = zeros(T, 3, size(A)[1], size(A)[2], size(A)[3])
   invdx = inv.(dx*2)
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂A∂x = (-A[i-1,j,1] + A[i+1,j,1]) * invdx[1]
            ∂A∂y = (-A[i,j-1,1] + A[i,j+1,1]) * invdx[2]
            ∂A∂z = 0

            Bx[i,j,1] = ∂A∂x
            By[i,j,1] = ∂A∂y
            Bz[i,j,1] = ∂A∂z
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂A∂x = (-A[i-1,1,k] + A[i+1,1,k]) * invdx[1]
            ∂A∂y = 0
            ∂A∂z = (-A[i,1,k-1] + A[i,1,k+1]) * invdx[3]

            Bx[i,1,k] = ∂A∂x
            By[i,1,k] = ∂A∂y
            Bz[i,1,k] = ∂A∂z
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,3)-1, j in 2:size(A,2)-1, i in 2:size(A,1)-1
         ∂A∂x = (-A[i-1,j,k] + A[i+1,j,k]) * invdx[1]
         ∂A∂y = (-A[i,j-1,k] + A[i,j+1,k]) * invdx[2]
         ∂A∂z = (-A[i,j,k-1] + A[i,j,k+1]) * invdx[3]

         Bx[i,j,k] = ∂A∂x
         By[i,j,k] = ∂A∂y
         Bz[i,j,k] = ∂A∂z
      end
   end
   B
end

"""
    fg_grad_component(dx::AbstractVector, Ain::AbstractArray, c::Integer)

Calculate 2nd order cell-centered ∇A where `A` is a 4D array of size (nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension, operating on component c of A.
"""
function fg_grad_component(dx::AbstractVector, Ain::AbstractArray{T,N}, c::Integer) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views A = Ain[c,:,:,:]
   #return fg_grad(dx, A) #?
   B = zeros(T, 3, size(A)[1], size(A)[2], size(A)[3])
   invdx = inv.(dx*2)
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]
   
   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂A∂x = (-A[i-1,j,1] + A[i+1,j,1]) * invdx[1]
            ∂A∂y = (-A[i,j-1,1] + A[i,j+1,1]) * invdx[2]
            ∂A∂z = 0

            Bx[i,j,1] = ∂A∂x
            By[i,j,1] = ∂A∂y
            Bz[i,j,1] = ∂A∂z
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂A∂x = (-A[i-1,1,k] + A[i+1,1,k]) * invdx[1]
            ∂A∂y = 0
            ∂A∂z = (-A[i,1,k-1] + A[i,1,k+1]) * invdx[3]

            Bx[i,1,k] = ∂A∂x
            By[i,1,k] = ∂A∂y
            Bz[i,1,k] = ∂A∂z
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,3)-1, j in 2:size(A,2)-1, i in 2:size(A,1)-1
         ∂A∂x = (-A[i-1,j,k] + A[i+1,j,k]) * invdx[1]
         ∂A∂y = (-A[i,j-1,k] + A[i,j+1,k]) * invdx[2]
         ∂A∂z = (-A[i,j,k-1] + A[i,j,k+1]) * invdx[3]

         Bx[i,j,k] = ∂A∂x
         By[i,j,k] = ∂A∂y
         Bz[i,j,k] = ∂A∂z
      end
   end
   B
end

"""
    fg_curl(dx::AbstractVector, A::AbstractArray)

Calculate 2nd order cell-centered ∇×A where `A` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function fg_curl(dx::AbstractVector, A::AbstractArray{T,N}) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]
   
   B = zeros(T, size(A))
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

   invdx = inv.(dx*2)

   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,j,1] + Ay[i+1,j,1]) * invdx[1]
            ∂Az∂x = (-Az[i-1,j,1] + Az[i+1,j,1]) * invdx[1]
            ∂Ax∂y = (-Ax[i,j-1,1] + Ax[i,j+1,1]) * invdx[2]
            ∂Az∂y = (-Az[i,j-1,1] + Az[i,j+1,1]) * invdx[2]
            ∂Ax∂z = 0
            ∂Ay∂z = 0

            Bx[i,j,1] = ∂Az∂y - ∂Ay∂z
            By[i,j,1] = ∂Ax∂z - ∂Az∂x
            Bz[i,j,1] = ∂Ay∂x - ∂Ax∂y
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,1,k] + Ay[i+1,1,k]) * invdx[1]
            ∂Az∂x = (-Az[i-1,1,k] + Az[i+1,1,k]) * invdx[1]
            ∂Ax∂y = 0
            ∂Az∂y = 0
            ∂Ax∂z = (-Ax[i,1,k-1] + Ax[i,1,k+1]) * invdx[3]
            ∂Ay∂z = (-Ay[i,1,k-1] + Ay[i,1,k+1]) * invdx[3]

            Bx[i,1,k] = ∂Az∂y - ∂Ay∂z
            By[i,1,k] = ∂Ax∂z - ∂Az∂x
            Bz[i,1,k] = ∂Ay∂x - ∂Ax∂y
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
         ∂Ay∂x = (-Ay[i-1,j,k] + Ay[i+1,j,k]) * invdx[1]
         ∂Az∂x = (-Az[i-1,j,k] + Az[i+1,j,k]) * invdx[1]
         ∂Ax∂y = (-Ax[i,j-1,k] + Ax[i,j+1,k]) * invdx[2]
         ∂Az∂y = (-Az[i,j-1,k] + Az[i,j+1,k]) * invdx[2]
         ∂Ax∂z = (-Ax[i,j,k-1] + Ax[i,j,k+1]) * invdx[3]
         ∂Ay∂z = (-Ay[i,j,k-1] + Ay[i,j,k+1]) * invdx[3]

         Bx[i,j,k] = ∂Az∂y - ∂Ay∂z
         By[i,j,k] = ∂Ax∂z - ∂Az∂x
         Bz[i,j,k] = ∂Ay∂x - ∂Ax∂y
      end
   end
   B
end


"""
    fg_div(dx::AbstractVector, A::AbstractArray)

Calculate 2nd order cell-centered ∇.A where `A` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function fg_div(dx::AbstractVector, A::AbstractArray{T,N}) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]
   
   B = zeros(T, size(A)[2:end])

   invdx = inv.(dx*2)

   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,j,1] + Ay[i+1,j,1]) * invdx[1]
            ∂Az∂x = (-Az[i-1,j,1] + Az[i+1,j,1]) * invdx[1]
            ∂Ax∂y = (-Ax[i,j-1,1] + Ax[i,j+1,1]) * invdx[2]
            ∂Az∂y = (-Az[i,j-1,1] + Az[i,j+1,1]) * invdx[2]
            ∂Ax∂z = 0
            ∂Ay∂z = 0

            Bx[i,j,1] = ∂Az∂y - ∂Ay∂z
            By[i,j,1] = ∂Ax∂z - ∂Az∂x
            Bz[i,j,1] = ∂Ay∂x - ∂Ax∂y
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,1,k] + Ay[i+1,1,k]) * invdx[1]
            ∂Az∂x = (-Az[i-1,1,k] + Az[i+1,1,k]) * invdx[1]
            ∂Ax∂y = 0
            ∂Az∂y = 0
            ∂Ax∂z = (-Ax[i,1,k-1] + Ax[i,1,k+1]) * invdx[3]
            ∂Ay∂z = (-Ay[i,1,k-1] + Ay[i,1,k+1]) * invdx[3]

            Bx[i,1,k] = ∂Az∂y - ∂Ay∂z
            By[i,1,k] = ∂Ax∂z - ∂Az∂x
            Bz[i,1,k] = ∂Ay∂x - ∂Ax∂y
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
         ∂Ax∂x = (-Ax[i-1,j,k] + Ax[i+1,j,k]) * invdx[1]
         ∂Ay∂y = (-Ay[i,j-1,k] + Ay[i,j+1,k]) * invdx[2]
         ∂Az∂z = (-Az[i,j,k-1] + Az[i,j,k+1]) * invdx[3]

         B[i,j,k] = ∂Ax∂x + ∂Ay∂y + ∂Az∂z
      end
   end
   B
end

"""
    fg_matder(dx::AbstractVector, A::AbstractArray, a::AbstractArray)

Calculate the material derivative of A along a, using 2nd order cell-centered differences.
`A` is a 4D array of size (3, nx, ny, nz),
`a` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""

function fg_matder(dx::AbstractVector, A::AbstractArray{T,N}, a::AbstractArray{T,N}) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]
   @views ax, ay, az = a[1,:,:,:], a[2,:,:,:], a[3,:,:,:]
   
   B = zeros(T, size(A))
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

   invdx = inv.(dx*2)

   if any(==(1), size(A)) # 2D
      println("Feel free to implement 2D")
      return missing
   else # 3D
      @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
         ∂Ax∂x = (-Ax[i-1,j,k] + Ax[i+1,j,k]) * invdx[1]
         ∂Ay∂x = (-Ay[i-1,j,k] + Ay[i+1,j,k]) * invdx[1]
         ∂Az∂x = (-Az[i-1,j,k] + Az[i+1,j,k]) * invdx[1]
         ∂Ax∂y = (-Ax[i,j-1,k] + Ax[i,j+1,k]) * invdx[2]
         ∂Ay∂y = (-Ay[i,j-1,k] + Ay[i,j+1,k]) * invdx[2]
         ∂Az∂y = (-Az[i,j-1,k] + Az[i,j+1,k]) * invdx[2]
         ∂Ax∂z = (-Ax[i,j,k-1] + Ax[i,j,k+1]) * invdx[3]
         ∂Ay∂z = (-Ay[i,j,k-1] + Ay[i,j,k+1]) * invdx[3]
         ∂Az∂z = (-Az[i,j,k-1] + Az[i,j,k+1]) * invdx[3]

         Bx[i,j,k] = ax[i,j,k]*∂Ax∂x + ay[i,j,k]*∂Ax∂y + az[i,j,k]*∂Ax∂z
         By[i,j,k] = ax[i,j,k]*∂Ay∂x + ay[i,j,k]*∂Ay∂y + az[i,j,k]*∂Ay∂z
         Bz[i,j,k] = ax[i,j,k]*∂Az∂x + ay[i,j,k]*∂Az∂y + az[i,j,k]*∂Az∂z
      end
   end
   B
end

function fg_curvature(dx::AbstractVector, A::AbstractArray{T,N}) where {T,N}
   An = fg_normalize(A)
   return fg_matder(dx, An, An)
end

function fg_normalize(A::AbstractArray{T,N}) where {T,N}
   norms = sqrt.(sum(A.^2, dims=1))
   An = A./norms
   return An
end

"""
   fg_kappaC: gradient scale of curvature of A along said curvature
"""
function fg_kappac(dx::AbstractVector, A::AbstractArray{T,N}) where {T,N}
   An = fg_normalize(A)
   Ad = fg_matder(dx, An, An)
   dropdims(sum(fg_normalize(Ad).*Ad,dims=1), dims=1)
end

"""
   fg_kappa: gradient scale of scalar A along curvature B
"""
function fg_kappa(dx::AbstractVector, A::AbstractArray{T,N}, B::AbstractArray{T,M}) where {T,N,M}
   Alng = fg_grad(dx, log.(A))
   dropdims(sum(fg_normalize(B).*Alng, dims=1), dims=1)
end
