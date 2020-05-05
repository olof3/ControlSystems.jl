@testset "test_simplification" begin
## SMINREAL ##
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])
@test sminreal(G) == G
@test sminreal(G[1, 1]) == ss(0)
@test sminreal(G[1, 2]) == ss([-1 -2.5; 4 0], [1; 0], [3 0], [0])
@test sminreal(G[2, 1]) == ss([-5], [2], [-2], [1])
@test sminreal(G[2, 2]) == ss([-6], [2], [1], [0])


## MINREAL ##

# From transfer functions

s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sysmin = minreal(sys)

@test size(sysmin.A,1) == 3 # Test that the reduction of sys worked

@test hinfnorm(sys - sysmin)[1] < 1e-15 # And that the answer is correct

@test_broken balreal(sys-sysmin)

@test all(sigma(sys-sysmin, [0.0, 1.0, 2.0])[1] .< 1e-15)  # Previously crashed because of zero dimensions in tzero

t = 0:0.1:10
y1,x1 = step(sys,t)[[1,3]]
y2,x2 = step(sysmin,t)[[1,3]]
@test sum(abs2,y1.-y2) < 1e-6 # Test that the output from the two systems are the same



# minreal of statespace realizations

@test minreal(ss([-1 0; 0 -2], [1; 0], [1 0], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss([-1 0; 0 -2], [1; 0], [1 1], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss([-1 0; 0 -2], [1; 1], [1 0], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss([-2 0; 0 -1], [0; 1], [0 1], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss([-2 0; 0 -1], [0; 1], [1 1], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss([-2 0; 0 -1], [1; 1], [0 1], 0)) == ss(-1, 1, 1, 0)

@test minreal(ss(diagm([-1, -2, -3, -4]), [1; zeros(3)], [1 zeros(1,3)], 0)) == ss(-1, 1, 1, 0)
@test minreal(ss(diagm([-1, -2, -3, -4]), [zeros(3); 1], [zeros(1,3) 1], 0)) == ss(-4, 1, 1, 0)


Ω = [1, 2, 3]
sys0 = ss(-1, 1, 1, 0)

@test freqresp(minreal(0.2*sys0 + 0.3*sys0 + 0.5*sys0), Ω) ≈ freqresp(sys0, Ω)
@test freqresp(minreal(sum([0.1*sys0 for _=1:10])), Ω) ≈ freqresp(sys0, Ω)



Ω = [1, 2, 3]

sys0 = ss([-1 1; -2 0], [0; 1], [1 0], 0)
#sys0 = ssrand(2,2,2,stable=true, proper=true)

eigvals(gram(sys0, :o))
eigvals(gram(sys0, :c))


freqresp(minreal(sys0/3 + sys0/3 + sys0/3), Ω) - freqresp(sys0, Ω)


@test freqresp(minreal(sum([0.2*sys0 for _=1:5])), Ω) ≈ freqresp(sys0, Ω)


freqresp(minreal(sum([0.2*sys0 for _=1:5])), Ω) - freqresp(sys0, Ω)

C = cholesky(Hermitian(P), Val(true), check=false)

Juno.@enter minreal(ss(diagm([-1, -2, -3, -4]), [zeros(3); 1], [zeros(1,3) 1], 0))

# Discrete time
@test minreal(ss([0.5 0; 0 0.2], [1; 0], [1 0], 0, 1)) == ss(0.5, 1, 1, 0, 1)
@test minreal(ss([0.5 0; 0 0.2], [1; 0], [1 1], 0, 1)) == ss(0.5, 1, 1, 0, 1)
@test minreal(ss([0.5 0; 0 0.2], [1; 1], [1 0], 0, 1)) == ss(0.5, 1, 1, 0, 1)



end
