sys = ss([-1 0; 0 -2], [1; 0], [1 0], 0)

P = gram(sys, :c)
Q = gram(sys, :o)

S = cholesky(Hermitian(P), Val(true), check=false).U
R = cholesky(Hermitian(Q), Val(true), check=false).U
U,Σ,V = svd(S*R') # Σ are the Hankel singular values

tol = 1e-8
k = findlast(>(tol), Σ)

Σ1 = diagm(sqrt.(Σ[1:k]))


T = Σ1 \ (V[:, 1:k]' * R)
Tp = (S'*U[:, 1:k]) / Σ1

Pz = T*P*Tp
Qz = T*Q*Tp
if norm(Pz-Qz) > sqrt(eps())
    @warn("balreal: Result may be inaccurate")
    println("Controllability gramian before transform")
    display(P)
    println("Controllability gramian after transform")
    display(Pz)
    println("Observability gramian before transform")
    display(Q)
    println("Observability gramian after transform")
    display(Qz)
    println("Singular values of PQ")
    display(Σ)
end

sysr = ss(T*sys.A*Tp, T*sys.B, sys.C*Tp, sys.D, sys.Ts)
