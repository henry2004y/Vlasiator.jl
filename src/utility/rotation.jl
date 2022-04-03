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
    getRotationMatrix(axis::AbstractVector, angle) --> SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in
radians.
Reference: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

# Example

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = v ./ norm(v) 
angle = -74 / 180 * π
R = getRotationMatrix(v̂, angle)
```
"""
function getRotationMatrix(v::AbstractVector{<:AbstractFloat}, θ)
   sinθ, cosθ = sincos(eltype(v)(θ))
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end

"""
    getRotationMatrix(e1::Matrix, e2::Matrix) --> SMatrix{3,3}

Obtain the rotation matrix from orthgonal base vectors `e1` to `e2`, such that a vector
``\\mathbf{u}_1`` in `e1` can be expressed as ``\\mathbf{u}_1 = M\\cdot \\mathbf{u}_2``,
where ``M`` is the rotation matrix and ``\\mathbf{u}_2`` is the same vector in `e2`.

# Example

```julia
e1 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
e2 = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
R = getRotationMatrix(e1, e2)
```
"""
function getRotationMatrix(e1::Matrix, e2::Matrix)
   @views begin
      r11 = e1[:,1] ⋅ e2[:,1]
      r12 = e1[:,1] ⋅ e2[:,2]
      r13 = e1[:,1] ⋅ e2[:,3]
      r21 = e1[:,2] ⋅ e2[:,1]
      r22 = e1[:,2] ⋅ e2[:,2]
      r23 = e1[:,2] ⋅ e2[:,3]
      r31 = e1[:,3] ⋅ e2[:,1]
      r32 = e1[:,3] ⋅ e2[:,2]
      r33 = e1[:,3] ⋅ e2[:,3]
   end
   R = @SMatrix [r11 r12 r13; r21 r22 r23; r31 r32 r33]
end
