###########################################################################
## Misc utils

"""
    matrix_tensor_product(A,B)

TBW
"""
matrix_tensor_product(A,B) = reshape(A * reshape(B, size(A,2), :), size(A,1), size(B)[2:end]...)

