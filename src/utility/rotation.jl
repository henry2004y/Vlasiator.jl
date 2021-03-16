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
   T_rotated = R * T * R'
end

""" 
    rotation matrix(vector, angle)

Creates a rotation matrix that rotates around a unit `vector` by an `angle` in
radians.
Reference: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
"""
function rotation_matrix(v, θ)
   cosθ, sinθ = cos(θ), sin(θ)
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end