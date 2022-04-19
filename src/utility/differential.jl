"""
    curl(A::AbstractArray{T,N}, dx=ones(T,3)) where {T,N}

Calculate 2nd order cell-centered ∇×A where `A` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function curl(A::AbstractArray{T,N}, dx=ones(T,3)) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]

   B = zeros(T, size(A))
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

   invdx = @. inv(2*dx)

   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,j,1] + Ay[i+1,j,1]) * invdx[1]
            ∂Az∂x = (-Az[i-1,j,1] + Az[i+1,j,1]) * invdx[1]
            ∂Ax∂y = (-Ax[i-1,j,1] + Ax[i+1,j,1]) * invdx[2]
            ∂Az∂y = (-Az[i-1,j,1] + Az[i+1,j,1]) * invdx[2]
            ∂Ax∂z = (-Ax[i-1,j,1] + Ax[i+1,j,1]) * invdx[3]
            ∂Ay∂z = (-Ay[i-1,j,1] + Ay[i+1,j,1]) * invdx[3]

            Bx[i,j,1] = ∂Az∂y - ∂Ay∂z
            By[i,j,1] = ∂Ax∂z - ∂Az∂x
            Bz[i,j,1] = ∂Ay∂x - ∂Ax∂y
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (-Ay[i-1,1,k] + Ay[i+1,1,k]) * invdx[1]
            ∂Az∂x = (-Az[i-1,1,k] + Az[i+1,1,k]) * invdx[1]
            ∂Ax∂y = (-Ax[i-1,1,k] + Ax[i+1,1,k]) * invdx[2]
            ∂Az∂y = (-Az[i-1,1,k] + Az[i+1,1,k]) * invdx[2]
            ∂Ax∂z = (-Ax[i-1,1,k] + Ax[i+1,1,k]) * invdx[3]
            ∂Ay∂z = (-Ay[i-1,1,k] + Ay[i+1,1,k]) * invdx[3]

            Bx[i,1,k] = ∂Az∂y - ∂Ay∂z
            By[i,1,k] = ∂Ax∂z - ∂Az∂x
            Bz[i,1,k] = ∂Ay∂x - ∂Ax∂y
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
         ∂Ay∂x = (-Ay[i-1,j,  k  ] + Ay[i+1,j  ,k  ]) * invdx[1]
         ∂Az∂x = (-Az[i-1,j,  k  ] + Az[i+1,j  ,k  ]) * invdx[1]
         ∂Ax∂y = (-Ax[i  ,j-1,k  ] + Ax[i  ,j+1,k  ]) * invdx[2]
         ∂Az∂y = (-Az[i  ,j-1,k  ] + Az[i  ,j+1,k  ]) * invdx[2]
         ∂Ax∂z = (-Ax[i  ,j,  k-1] + Ax[i  ,j  ,k+1]) * invdx[3]
         ∂Ay∂z = (-Ay[i  ,j,  k-1] + Ay[i  ,j  ,k+1]) * invdx[3]

         Bx[i,j,k] = ∂Az∂y - ∂Ay∂z
         By[i,j,k] = ∂Ax∂z - ∂Az∂x
         Bz[i,j,k] = ∂Ay∂x - ∂Ax∂y
      end
   end
   B
end

"""
    gradient(A::AbstractArray{T,N}, dx=ones(T, 3)) where {T,N}

Calculate 2nd order cell-centered ∇A where `A` is a scalar array and `dx` is a vector of
grid intervals in each dimension.
!!! warning
    The current implementation has issues at the boundary if gradient is taken multiple
    times.
"""
function gradient(A::AbstractArray{T,N}, dx=ones(T, N)) where {T,N}
   @assert N < 4 "$N dimension array A detected!"
   @assert N == length(dx) "The array A shall have the same dimension as dx!"
   @assert all(!=(1), size(A)) "No singular dimension is allowed!"

   if N == 3 # 3D
      B = zeros(T, 3, size(A)[1], size(A)[2], size(A)[3])
      invdx = @. inv(2*dx)
      @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

      @inbounds for k in 2:size(A,3)-1, j in 2:size(A,2)-1, i in 2:size(A,1)-1
         ∂A∂x = (-A[i-1,j  ,k  ] + A[i+1,j,  k  ]) * invdx[1]
         ∂A∂y = (-A[i  ,j-1,k  ] + A[i,  j+1,k  ]) * invdx[2]
         ∂A∂z = (-A[i  ,j  ,k-1] + A[i,  j,  k+1]) * invdx[3]

         Bx[i,j,k] = ∂A∂x
         By[i,j,k] = ∂A∂y
         Bz[i,j,k] = ∂A∂z
      end
   elseif N == 2 # 2D
      B = zeros(T, 2, size(A)[1], size(A)[2])
      invdx = @. inv(2*dx)
      @views Bx, By = B[1,:,:], B[2,:,:]

      @inbounds for j in 2:size(A,2)-1, i in 2:size(A,1)-1
         ∂A∂x = (-A[i-1,j  ] + A[i+1,j  ]) * invdx[1]
         ∂A∂y = (-A[i  ,j-1] + A[i  ,j+1]) * invdx[2]

         Bx[i,j] = ∂A∂x
         By[i,j] = ∂A∂y
      end
   else # 1D
      B = zeros(T, size(A)[1])
      invdx = @. inv(2*dx)
      @views Bx = B[:]

      @inbounds for i in 2:size(A,1)-1
         ∂A∂x = (-A[i-1] + A[i+1]) * invdx[1]

         Bx[i] = ∂A∂x
      end     
   end
   B
end

"""
    divergence(A::AbstractArray{T,N}, dx=ones(T, 3)) where {T,N}

Calculate 2nd order cell-centered ∇⋅A where `A` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function divergence(A::AbstractArray{T,N}, dx=ones(T, 3)) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"
   @assert all(!=(1), size(A)) "Input vector must be from 3D data!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]

   B = zeros(T, size(A)[2:end])

   invdx = @. inv(2*dx)

   @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
      ∂Ax∂x = (-Ax[i-1,j  ,k  ] + Ax[i+1,j  ,k  ]) * invdx[1]
      ∂Ay∂y = (-Ay[i  ,j-1,k  ] + Ay[i  ,j+1,k  ]) * invdx[2]
      ∂Az∂z = (-Az[i  ,j  ,k-1] + Az[i  ,j  ,k+1]) * invdx[3]

      B[i,j,k] = ∂Ax∂x + ∂Ay∂y + ∂Az∂z
   end
   B
end