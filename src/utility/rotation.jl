# Utility functions for common algebraic operations.

"""
    rotateTensorToVectorZ(tensor, vector)

Rotates `tensor` with a rotation matrix that aligns z-axis with `vector`.
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

"""
    getRotationB(B) -> SMatrix

Obtain a rotation matrix with each column being a unit vector which is parallel (`v3`) and
perpendicular (`v1,v2`) to the magnetic field `B`. The two perpendicular directions are
chosen based on the reference vector of z-axis in the Cartesian coordinates.
"""
function getRotationB(B::AbstractVector{T}) where T
   b = hypot(B[1], B[2], B[3])
   R =
      if B[3] / b ≈ -1.0 # B aligned with reference vector
         SMatrix{3,3,T}([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0])
      elseif B[3] / b ≈ 1.0 # B aligned with reference vector
         SMatrix{3,3,T}([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0])
      else
         @warning "bug: not working correctly for general cases!"
         v0 = SVector{3,T}(0.0, 0.0, 1.0) # reference vector
         v3 = SVector{3,T}(B[1]/b, B[2]/b, B[3]/b) # unit vector along B
         v1 = v0 × v3::SVector{3, T}
         v2 = v3 × v1::SVector{3, T}

         @SMatrix [v1[1] v2[1] v3[1]; v1[2] v2[2] v3[2]; v1[3] v2[3] v3[3]]
      end
   R
end

"""
    rotateWithB(T, B) -> Matrix

Rotate the tensor `T` with the 3rd direction aligned with `B`.
See also: [`rotateWithB!`](@ref)
"""
function rotateWithB(T::AbstractMatrix, B::AbstractVector)
   R = getRotationB(B)
   R * T * R'
end

"""
    rotateWithB!(T, B)

Rotate the tensor `T` with the 3rd direction aligned with `B`.
See also: [`rotateWithB`](@ref)
"""
function rotateWithB!(T::AbstractMatrix, B::AbstractVector)
   R = getRotationB(B)
   T[:] = R * T * R'
end