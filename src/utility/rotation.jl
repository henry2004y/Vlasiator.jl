# Utility functions for common algebraic operations.

"""
    rotateTensorToVectorZ(tensor, vector)

Rotate `tensor` with a rotation matrix that aligns the 3rd direction with `vector`, which is
equivalent to change the basis from (i,j,k) to (i′,j′,k′) where k′ ∥ vector.
Reference: https://math.stackexchange.com/questions/2303869/tensor-rotation
"""
function rotateTensorToVectorZ(tensor::AbstractMatrix{T}, v::AbstractVector{T}) where T
   k = SVector{3, T}(0.0, 0.0, 1.0)
   axis = v × k::SVector{3, T}
   if axis[1] == axis[2] == 0
      return tensor
   else
      normalize!(axis)
      angle = acos(v ⋅ k / hypot(v[1], v[2], v[3]))
      R = getRotationMatrix(axis, angle)
      return R * tensor * R'
   end
end

"""
    getRotationMatrix(axis, angle) --> SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in
radians.
Reference: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
"""
function getRotationMatrix(v::AbstractVector, θ)
   sinθ, cosθ = sincos(eltype(v)(θ))
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end

