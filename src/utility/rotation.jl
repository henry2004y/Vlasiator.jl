# Utility functions for common algebraic operations.

import LinearAlgebra: ⋅, ×
using StaticArrays

"""
    rotateTensorToVector(tensor, vector)

Rotates `tensor` with a rotation matrix that aligns z-axis with `vector`.
"""
function rotateTensorToVectorZ!(T, v)
   unitz = @SVector [0.0, 0.0, 1.0]
   vz = v × unitz
   vz ./= hypot(vz...)
   angle = acos(v ⋅ unitz / hypot(v...))
   R = rotation_matrix(vz, angle)
   # Rotate Tensor
   R * T * R'
end

"""
    get_rotation_matrix(vector, angle)

Creates a rotation matrix that rotates around a unit `vector` by an `angle` in
radians.
References: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
"""
function get_rotation_matrix(v, θ)
   cosθ, sinθ = cos(θ), sin(θ)
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end

"""
    getRotationB(B) -> SMatrix

Obtain a rotation matrix with each column being a unit vector which is parallel (`v3`) and
perpendicular (`v1,v2`) to the magnetic field `B`. The two perpendicular directions are
chosen based on the reference vector of z-axis in the Cartesian coordinates.
"""
function getRotationB(B)
   # reference vector
   v0 = @SVector [0.0, 0.0, 1.0]
   b = hypot(B...)
   if B[3] / b ≈ 1.0 # B aligned with reference vector
      R = @SMatrix [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
   else
      # vector along B
      v3 = @SVector [B[1]/b, B[2]/b, B[3]/b]
      v1 = v0 × v3
      v2 = v3 × v1

      R = @SMatrix [v1[1] v2[1] v3[1]; v1[2] v2[2] v3[2]; v1[3] v2[3] v3[3]]
   end
   R
end

"""
    rotateWithB(T, B) -> Matrix

Rotate the tensor `T` with the 3rd direction aligned with `B`.
See also: [`rotateWithB!`](@ref)
"""
function rotateWithB(T, B)
   R = getRotationB(B)
   R * T * R'
end

"""
    rotateWithB!(T, B)

Rotate the tensor `T` with the 3rd direction aligned with `B`.
See also: [`rotateWithB`](@ref)
"""
function rotateWithB!(T, B)
   R = getRotationB(B)
   T[:] = R * T * R'
end