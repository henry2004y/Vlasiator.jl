"""
    curl(dx::AbstractVector, A::AbstractArray)

Calculate 2nd order cell-centered ∇×A where `A` is a 4D array of size (3, nx, ny, nz) and
`dx` is a vector of grid intervals in each dimension.
"""
function curl(dx, A::AbstractArray{T,N}) where {T,N}
   @assert N == 4 && length(dx) == 3 "Input vector shall be indexed in 3D!"

   @views Ax, Ay, Az = A[1,:,:,:], A[2,:,:,:], A[3,:,:,:]
   
   B = zeros(T, size(A))
   @views Bx, By, Bz = B[1,:,:,:], B[2,:,:,:], B[3,:,:,:]

   invdx = inv.(dx)

   if any(==(1), size(A)) # 2D
      if size(A,4) == 1
         @inbounds for j in 2:size(A,3)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (Ay[i-1,j,1] - 2Ay[i,j,1] + Ay[i+1,j,1]) * invdx[1]
            ∂Az∂x = (Az[i-1,j,1] - 2Az[i,j,1] + Az[i+1,j,1]) * invdx[1]
            ∂Ax∂y = (Ax[i-1,j,1] - 2Ax[i,j,1] + Ax[i+1,j,1]) * invdx[2]
            ∂Az∂y = (Az[i-1,j,1] - 2Az[i,j,1] + Az[i+1,j,1]) * invdx[2]
            ∂Ax∂z = (Ax[i-1,j,1] - 2Ax[i,j,1] + Ax[i+1,j,1]) * invdx[3]
            ∂Ay∂z = (Ay[i-1,j,1] - 2Ay[i,j,1] + Ay[i+1,j,1]) * invdx[3]

            Bx[i,j,1] = ∂Az∂y - ∂Ay∂z
            By[i,j,1] = ∂Ax∂z - ∂Az∂x
            Bz[i,j,1] = ∂Ay∂x - ∂Ax∂y
         end
      elseif size(A,3) == 1
         @inbounds for k in 2:size(A,4)-1, i in 2:size(A,2)-1
            ∂Ay∂x = (Ay[i-1,1,k] - 2Ay[i,1,k] + Ay[i+1,1,k]) * invdx[1]
            ∂Az∂x = (Az[i-1,1,k] - 2Az[i,1,k] + Az[i+1,1,k]) * invdx[1]
            ∂Ax∂y = (Ax[i-1,1,k] - 2Ax[i,1,k] + Ax[i+1,1,k]) * invdx[2]
            ∂Az∂y = (Az[i-1,1,k] - 2Az[i,1,k] + Az[i+1,1,k]) * invdx[2]
            ∂Ax∂z = (Ax[i-1,1,k] - 2Ax[i,1,k] + Ax[i+1,1,k]) * invdx[3]
            ∂Ay∂z = (Ay[i-1,1,k] - 2Ay[i,1,k] + Ay[i+1,1,k]) * invdx[3]

            Bx[i,1,k] = ∂Az∂y - ∂Ay∂z
            By[i,1,k] = ∂Ax∂z - ∂Az∂x
            Bz[i,1,k] = ∂Ay∂x - ∂Ax∂y
         end
      end
   else # 3D
      @inbounds for k in 2:size(A,4)-1, j in 2:size(A,3)-1, i in 2:size(A,2)-1
         ∂Ay∂x = (Ay[i-1,j,k] - 2Ay[i,j,k] + Ay[i+1,j,k]) * invdx[1]
         ∂Az∂x = (Az[i-1,j,k] - 2Az[i,j,k] + Az[i+1,j,k]) * invdx[1]
         ∂Ax∂y = (Ax[i-1,j,k] - 2Ax[i,j,k] + Ax[i+1,j,k]) * invdx[2]
         ∂Az∂y = (Az[i-1,j,k] - 2Az[i,j,k] + Az[i+1,j,k]) * invdx[2]
         ∂Ax∂z = (Ax[i-1,j,k] - 2Ax[i,j,k] + Ax[i+1,j,k]) * invdx[3]
         ∂Ay∂z = (Ay[i-1,j,k] - 2Ay[i,j,k] + Ay[i+1,j,k]) * invdx[3]

         Bx[i,j,k] = ∂Az∂y - ∂Ay∂z
         By[i,j,k] = ∂Ax∂z - ∂Az∂x
         Bz[i,j,k] = ∂Ay∂x - ∂Ax∂y
      end
   end
   B
end