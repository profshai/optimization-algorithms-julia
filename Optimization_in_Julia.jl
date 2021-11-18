import Pkg
Pkg.add("Ipopt")
Pkg.add("JuMP")

using JuMP
using Ipopt

model =  Model(Ipopt.Optimizer)

gamma1, gamma2, gamma3 =   2.0, 10.0, 20.0
w1, w2, w3 = 1.0,2.0,3.0
p1, p2, p3 = 1.0,2.0,6.0
I = 10.0


c1_start = 10.0/3
c2_start  = (I-p1*c1_start )/ p2
c3_start  = (I-p1*c1_start )/ p3
@variable(model, c1, start = c1_start)
@variable(model, c2, start = c2_start)
@variable(model, c3, start = c3_start)

@NLobjective(model, Max, w1*c1^(1-gamma1)/(1-gamma1) + w2*c2^(1-gamma2)/(1-gamma2) + w3*c3^(1-gamma3)/(1-gamma3))
@constraint(model, p1*c1 + p2*c2 + p3*c3== I)
@constraint(model, c1 >=0)
@constraint(model, c2 >=0)
@constraint(model, c3 >=0)

optimize!(model)

# Optimal consumption values:
println("c1 = ", value(c1), " c2 = ", value(c2), " c3 = ", value(c3))

function bisection(f,a,b,tol)
    sa = sign(f(a))
    sb = sign(f(b))
    if sa == sb
        error("Interval is not valid")
    end
    while abs(b-a) > tol
        c = (a+b)/2.0
        sc = sign(f(c))
        if sc == 0.0
            a,b = c,c
        elseif sa == sc
            a,b = c,b
        else
            a,b = a,c
        end
    end
    return (a+b)/2.0
end

gamma = [2.0; 10.0; 20.0]
w = [1.0; 2.0; 3.0]
p = [1.0;2.0;6.0]
I = 10.0

struct consumer
    gamma
    w
    p
    I
end
con = consumer(gamma, w, p, I)
function nc(e::consumer)
    return size(e.gamma,1)
end

function choice(e::consumer, lambda, i)
    d = lambda *e.p[i]
    return (e.w[i]/d)^(1.0/e.gamma[i])
end

function choices(e::consumer, lambda)
    n = nc(e)
    c = Vector{Float64}(undef, n)
    for i=1:n
        c[i] = choice(e,lambda,i)
    end
    return c
end

function foc(e::consumer, lambda)
    s = 0.0
    n = nc(e)
    for i=1:n
        s = s + choice(e,lambda,i)*e.p[i]
    end
    return e.I - s
end

function optimal_bisection(e,a,b)
    tol = 1e-16
    function goal(lambda)
        return foc(e, lambda)
    end
    lambda_star = bisection(goal, a,b,tol)
    return choices(e, lambda_star)
end
a = 0.000000001
b =1.0
c_values = optimal_bisection(con,a,b)
display(c_values)


