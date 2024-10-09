###########################################################################
## Misc utils

"""
    matrix_tensor_product(A,B)

TBW
"""
matrix_tensor_product(A,B) = reshape(reshape(A, :, size(A)[end]) * reshape(B, size(A)[end], :), size(A)[1:end-1]..., size(B)[2:end]...)

