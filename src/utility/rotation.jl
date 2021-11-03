# Utility functions for common algebraic operations.

"""
    rotateTensorToVectorZ(tensor, vector)

Rotates `tensor` with a rotation matrix that aligns the 3rd direction with `vector`.
"""
function rotateTensorToVectorZ(tensor::AbstractMatrix{T}, v::AbstractVector{T}) where T
   unitz = SVector{3, T}(0.0, 0.0, 1.0)
   vz = v × unitz::SVector{3, T}
   if vz[1] == vz[2] == 0
      return tensor
   else
      vz ./= hypot(vz[1], vz[2], vz[3])
      angle = acos(v ⋅ unitz / hypot(v[1], v[2], v[3]))
      R = getRotationMatrix(vz, angle)
      return R * tensor * R'
   end
end

"""
    getRotationMatrix(vector, angle)

Creates a rotation matrix that rotates around a unit `vector` by an `angle` in radians.
References: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
"""
function getRotationMatrix(v::AbstractVector, θ)
   sinθ, cosθ = sincos(eltype(v)(θ))
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end
