module LyapTest

#=
* Symmetry is not exploited when solving diagonal entries
* Perhaps it would be better with separate methods for the real and
  complex versions of the block solve algorithms
* What convention for signs and tranposes?
* Should error checking be done? In that case where?
=#


using LinearAlgebra
using StaticArrays

"""
    `b, d, nblocks = _schurstructure(R::AbstractMatrix)`

Return the block strucutre of a quasi-traingular Schur matrix `R`.

`d` contains the block sizees of each diagonal block (`1` or `2`)

`b` contains the indices of the blocks

`nblocks` is the number of blocks

"""
function _schurstructure(R::AbstractMatrix, ul=Val(:U)::Union{Val{:U}, Val{:L}})
    n = size(R,1)

    d = Vector{Int}(undef, n) # block sizes
    b = Vector{UnitRange{Int64}}(undef, n) # block indices

    j = 1 # column if ul=:u, row if ul=:l
    k = 0 # block number
    while j <= n
        k += 1
        if j == n
            d[k] = 1
        else
            if ul == Val(:U)
                d[k] = iszero(R[j+1, j]) ? 1 : 2
            else
                d[k] = iszero(R[j, j+1]) ? 1 : 2
            end
        end
        b[k] = j:j+d[k]-1
        j += d[k]
    end
    resize!(d, k)
    resize!(b, k)
    return d, b, k
end

# FIXME: better handling of uniform scaling?!
issquare(A::Number) = true
issquare(A::AbstractMatrix) = size(A,1) == size(A,2)
function _check_lyap_inputs(A, Q)
    if !issquare(A); error("The A matrix must be square"); end
    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end
    if size(Q, 1) != size(A, 1); error("The A and Q matrices must have the same dimensions"); end
end

function _check_sylv_inputs(A, B, C)
    if !issquare(A); error("The A matrix must be square"); end
    if !issquare(B); error("The B matrix must be square"); end
    if size(C, 1) != size(A, 1); error("The A and C matrices have inconsistent dimensions"); end
    if size(C, 2) != size(B, 2); error("The B and C matrices have inconsistent dimensions"); end
end

# Should preferably be fixed in LinearAlgebra
LinearAlgebra.schur(A::AbstractMatrix{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))


sylvcsoltype(A, B, C) = Base.promote_op(sylvc, eltype(A), eltype(B), eltype(C))
sylvdsoltype(A, B, C) = Base.promote_op(sylvd, eltype(A), eltype(B), eltype(C))

"""
    _sylvc!(A, B, C)

Solve the continuous-time sylvester equation

`AX + XB = C`

for small matrices (1x1, 1x2, 2x1, 2x2)
"""
function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvc!(C, A, B, C, Val{2}(), Val{2}())
    elseif M == 2 && N == 1
        _sylvc!(C, A, B, C, Val{2}(), Val{1}())
    elseif M == 1 && N == 2
        _sylvc!(C, A, B, C, Val{1}(), Val{2}())
    elseif M == 1 && N == 1
        _sylvc!(C, A, B, C, Val{1}(), Val{1}())
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvc!(X::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    A′ = SMatrix{M,M}(A)
    B′ = SMatrix{N,N}(B)
    Cv = SVector{M*N}(C[:]) # vectorization of C

    Xv = lu(kron(SMatrix{N,N}(I), A′) + kron(transpose(B′), SMatrix{M,M}(I))) \ Cv # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X) (with A = I or B = I)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    X .= reshape(Xv, M, N)
end


"""
    _sylvd!(X, A, B, C)

Solve the discrete-time sylvester equation

`AXB - X = C`

for small matrices (1x1, 1x2, 2x1, 2x2).
"""
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvd!(C, A, B, C, Val{2}(), Val{2}())
    elseif M == 2 && N == 1
        _sylvd!(C, A, B, C, Val{2}(), Val{1}())
    elseif M == 1 && N == 2
        _sylvd!(C, A, B, C, Val{1}(), Val{2}())
    elseif M == 1 && N == 1
        _sylvd!(C, A, B, C, Val{1}(), Val{1}())
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvd!(X::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    A′ = SMatrix{M,M}(A)
    B′ = SMatrix{N,N}(B)
    Cv = SVector{M*N}(C[:]) # vectorization of C

    Xv = (kron(transpose(B′), A′) - I) \ Cv # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    X .= reshape(Xv, M, N)
end

"""
    sylvc(A, B, C)

Solve the continuous-time Sylvester equation

`AX + XB = C`

A solution `X` exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λ + μ = 0.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvc(A, B, C)
    _check_sylv_inputs(A, B, C)

    At2, UA = schur(A')
    B2, UB = schur(B)

    C2 = UA'*C*UB # C2 should have the right type

    Y = _sylvc_schur!(Matrix(At2'), B2, C2, Val(:sylv))

    X = UA*Y*UB'
end
@inline function sylvc(a::Number, b::Number, c::Number)
    x = c / (a + b)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    x
end

"""
    sylvd(A, B, C)

Solve the discrete-time Sylvester equation

`AXB - X = C`

A solution `X` exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λμ = 1.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvd(A, B, C)
    _check_sylv_inputs(A, B, C)

    At2, UA = schur(A')
    B2, UB = schur(B)

    C2 = UA'*C*UB

    Y = _sylvd_schur!(Matrix(At2'), B2, C2, Val(:sylv))

    X = UA*Y*UB'
end
@inline function sylvd(a::Number, b::Number, c::Number)
    x = c / (a * b - 1)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    x
end

"""
    lyapc(A, Q)

Computes the solution `X` of the continuous-time Lyapunov equation

`AX + XA' + Q = 0`

A solution exists unless `A` has an eigenvalue λ = ±1 or an eigenvalue pair λ₁λ₂ = 1.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function lyapc(A, Q)

     _check_lyap_inputs(A, Q)

    At2, U = schur(A')

    Q2 = U'*Q*U # Could use in-place, Small savings could be possible here, see [1]

    Y = _sylvc_schur!(Matrix(At2'), At2, lmul!(-1, Q2), Val(:lyap))

    X = U*Y*U' # Small savings could be possible here, see [1]
end
"""
    `X = lyapd(A, Q)`

Find the solution `X` to the discrete-time Lyapunov equation

`AXA' - X + Q = 0`

A solution `X` unless `A` has an eigenvalue λ = ±1 or an eigenvalue pair λ₁λ₂ = 1 .


[1] **Barraud, A.** (1977) "A numerical algorithm to solve A'XA - X = Q"
    IEEE Transactions on Automatic Control

[2] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.

"""
function lyapd(A, Q)

    _check_lyap_inputs(A, Q)

    At2, U = schur(A')

    #tmp = similar(Q, sylvdsoltype(A,A,Q))
    #Q2 = U'*mul!(tmp, Q, U) # Some savings could be possible here, see [1]
    Q2 = U'*Q*U # Some savings could be possible here, see [1]

    Y = _sylvd_schur!(Matrix(At2'), At2, lmul!(-1, Q2), Val(:lyap))

    X = U*Y*U' # mul!(Y.data, U, mul!(tmp, Y, U')) # Y*U*U' # Some savings could be possible here, see [1]
end


"""
Solve the continuous-time Sylvester equation

`AX + XB = C`

where
`A` is assumed to have lower Schur form (quasi-triangular, 1x1 & 2x2 blocks on the diagonal)
`B` is assumed to have upper Schur form

If `alg == Val(:lyap)` then `C` should be Hermitian

The input matrix `C` is overwritten with the solution `X`,
so this matrix should have the same type as the solution, typically `typeof(C[1] / (A[1] + B[1]))`

See also `sylvc`
"""
function sylvc_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    _check_sylv_inputs(A, B, C)
    T = sylvcsoltype(A, B, C)
    _sylvc_schur!(convert(Matrix, A), convert(Matrix, B), convert(Matrix{T}, C), Val(:sylv))
end
function lyapc_schur!(A::AbstractMatrix, Q::AbstractMatrix)
    _check_lyap_inputs(A, Q)
    T = sylvcsoltype(A, A, Q)
    _sylvc_schur!(convert(Matrix, A), convert(Matrix, A'), lmul!(-1, Matrix{T}(Q)), Val(:lyap))
end
function _sylvc_schur!(A::Matrix, B::Matrix, C::Matrix, alg::Union{Val{:sylv},Val{:lyap}},
    schurtype::Union{Val{:real},Val{:complex}} = isreal(A) || isreal(B) ? Val(:real) : Val(:complex)) where {T <: Number}
    # The user should preferably use sylvc_schur! and lyapc_schur!
    # I.e., this method does not check whether C is hermitian
    # The matrix C is successively replaced with the solution X
    # if alg == Val(:lyap), only the lower triangle of C is computed
    # after which an Hermitian view is applied

    # get block indices and nbr of blocks
    if schurtype == Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:L)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:U))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    for j=1:nblocksb
        i0 = (alg == Val(:lyap) ? j : 1)
        for i=i0:nblocksa
            if schurtype == Val(:complex)
                if i > 1; C[i,j] -= sum(A[i, k] * C[k, j] for k=1:i-1); end
                if j > 1; C[i,j] -= sum(C[i, k] * B[k, j] for k=1:j-1); end

                C[i,j] = sylvc(A[i, i], B[j, j], C[i, j]) # C[i,j] now contains  solution Y[i,j]

                if alg == Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end
            else
                Cij = view(C, ba[i], bb[j])

                if i > 1; @views mul!(Cij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], -1, 1); end
                if j > 1; @views mul!(Cij, C[ba[i], 1:bb[j-1][end]], B[1:bb[j-1][end], bb[j]], -1, 1); end

                @views _sylvc!(A[ba[i], ba[i]], B[bb[j], bb[j]], Cij) # Cij now contains the solution Yij

                if alg == Val(:lyap) && i > j
                    #copyto_noalias!(view(C, ba[j], bb[i]), Cij')
                    view(C, ba[j], bb[i]) .= Cij'
                end
            end
        end
    end
    return C
end

"""
Solve the discrete-time Sylvester equation

`AXB - X = C`

where `A` is assumed to have lower Schur form (quasi-triangular, 1x1 & 2x2 blocks on the diagonal)
`B` is assumed to have upper Schur form

If `alg == Val(:lyap)` then `C` should be Hermitian

The input matrix `C` is overwritten with the solution `X`,
so this matrix should have the same type as the solution, typically `typeof(C[1] / (A[1] * B[1] - 1)`

See also `sylvd`
"""
function sylvd_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    _check_sylv_inputs(A, B, C)
    T = sylvdsoltype(A, B, C)
    _sylvd_schur!(convert(Matrix, A), convert(Matrix, B), convert(Matrix{T}, C), Val(:sylv))
end
function lyapd_schur!(A::AbstractMatrix, Q::AbstractMatrix)
    _check_lyap_inputs(A, Q)
    T = sylvdsoltype(A, A, Q)
    _sylvd_schur!(convert(Matrix, A), convert(Matrix, A'), lmul!(-1, Matrix{T}(Q)), Val(:lyap))
end
function _sylvd_schur!(A::Matrix, B::Matrix, C::Matrix, alg::Union{Val{:sylv},Val{:lyap}},
    schurtype::Union{Val{:real},Val{:complex}} = isreal(A) || isreal(B) ? Val(:real) : Val(:complex)) where {T <: Number}

    G = zeros(eltype(C), size(A,1), size(B, 1)) # Keep track of A*X for improved performance

    # get block dimensions, block indices, nbr of blocks
    if schurtype == Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:L)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:U))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    for j=1:nblocksb
        i0 = (alg == Val(:lyap) ? j : 1)
        for i=i0:nblocksa
            if schurtype == Val(:complex)
                # Compute Gij up to the contribution from Aii*Yij which is added at the end of each iteration
                if i > 1; G[i,j] += sum(A[i,k] * C[k,j] for k=1:i-1); end

                C[i,j] -= sum(G[i,k] * B[k,j] for k=1:j)

                C[i,j] = sylvd(A[i,i], B[j,j], C[i,j]) # C[i,j] now contains  solution Y[i,j]

                if alg == Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end

                G[i,j] += A[i, i] * C[i, j]
            else
                Aii = A[ba[i], ba[i]]
                Bjj = B[bb[j], bb[j]]

                Cij = view(C, ba[i], bb[j])
                Gij = view(G, ba[i], bb[j])

                if i > 1
                    @views mul!(Gij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], 1, 1)
                end

                @views mul!(Cij, G[ba[i], 1:bb[j][end]], B[1:bb[j][end], bb[j]], -1, 1)

                _sylvd!(Aii, Bjj, Cij) # Cij now contains the solution Yij

                if alg == Val(:lyap) && i > j
                    #copyto_noalias!(view(C, ba[j], bb[i]), Cij')
                    view(C, ba[j], bb[i]) .= Cij'
                end

                mul!(Gij, Aii, Cij, 1, 1)
            end
        end
    end
    return C
end




function copyto_noalias!(dest::AbstractArray{T1,N}, src::AbstractArray{T2,N}) where {T1,T2,N}
     for I in eachindex(IndexStyle(src,dest), src)
         @inbounds dest[I] = src[I]
     end
     dest
end


end
