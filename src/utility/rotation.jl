# Utility functions for common algebraic operations.

"""
    rotateTensorToVectorZ(tensor::AbstractMatrix, vector::AbstractVector) -> SMatrix{3,3}

Rotate `tensor` with a rotation matrix that aligns the 3rd direction with `vector`, which is
equivalent to change the basis from (i,j,k) to (i′,j′,k′) where k′ ∥ vector.
Reference: [Tensor rotation](https://math.stackexchange.com/questions/2303869/tensor-rotation)
"""
function rotateTensorToVectorZ(tensor::AbstractMatrix{T}, v::AbstractVector{T}) where T
   k = SVector{3, T}(0.0, 0.0, 1.0)
   axis = v × k::SVector{3, T}
   if axis[1] == axis[2] == 0
      return tensor
   else
      angle = acos(v ⋅ k / √(v[1]^2 + v[2]^2 + v[3]^2))
      R = AngleAxis(angle, axis...)
      return R * tensor * R'
   end
end

"""
    getRotationMatrix(e1::AbtractMatrix) -> AngleAxis

Obtain the rotation matrix from the Cartesian coordinate to orthgonal base vectors `vbase`,
where each column represents a base vector.

# Example

```julia
e1 = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
R = getRotationMatrix(e1)
```
"""
function getRotationMatrix(vbase::AbstractMatrix{T}) where T
   k = SVector{3, T}(0, 0, 1)
   vref = SVector{3, T}(vbase[1,3], vbase[2,3], vbase[3,3])

   if vref == k
      vref = SVector{3, T}(vbase[1,2], vbase[2,2], vbase[3,2])
      k = SVector{3, T}(0, 1, 0)
   end
   if vref == k
      axis = SVector{3, T}(1, 0, 0)
      angle = 0.0
   else
      axis = k × vref
      angle = acos(k ⋅ vref)
   end

   R = AngleAxis(angle, axis...)
end