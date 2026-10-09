import Pkg; Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("LaTeXStrings")

using JuMP, Ipopt
 
function solve_mayer_free(; N::Int=500, eps_::Float64=1e-6,
                            warm_start::Union{Nothing,NamedTuple}=nothing)

    T = π
    hstep = T / N
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    # variables 
    @variable(model, x1[0:N])
    @variable(model, x2[0:N])
    @variable(model, x3[0:N])
    @variable(model, x4[0:N])
    @variable(model, 0 <= u[0:N] <= 1)
    
    @variable(model, eps_ <= x0 <= 1)
    @variable(model, -1.0 <= y0 <= 1.0)

    
    # avoid the solver to approach 0
    for k in 0:N
        set_lower_bound(x3[k], eps_)
        set_lower_bound(x1[k], eps_)
        set_upper_bound(x1[k], 1 - eps_)
    end

    if warm_start === nothing
        set_start_value(x0, 0.0) 
        set_start_value(y0, 0.0)
        for k in 0:N
            set_start_value(x1[k], 0.5 + 0 * k / N) 
            set_start_value(x2[k], 0.0)
            set_start_value(x3[k], 0.0)
            set_start_value(x4[k], 0.0)
            set_start_value(u[k], 1)
        end
    else
        set_start_value(x0, warm_start.x0)
        set_start_value(y0, warm_start.y0)
        for k in 0:N
            set_start_value(x1[k], warm_start.x1[k+1])
            set_start_value(x2[k], warm_start.x2[k+1])
            set_start_value(x3[k], warm_start.x3[k+1])
            set_start_value(x4[k], warm_start.x4[k+1])
            set_start_value(u[k],  warm_start.u[k+1])
        end
    end

    
    # initial-terminal constraints 
    @constraint(model, x1[0] == x0)
    @constraint(model, x2[0] == y0)
    fix(x3[0], 0.0; force = true)
    fix(x4[0], 0.0; force = true)

    @constraint(model, x1[N] == 1 - x0)
    @constraint(model, x2[N] == -y0)

    # dynamics def  
    f1(x1,x2,x3,x4,u) = x2
    f2(x1,x2,x3,x4,u) = -x1 + u
    f3(x1,x2,x3,x4,u) = x1*(u - 1.0) + 0.5
    f4(x1,x2,x3,x4,u) = 0.5/(x1^2) + 0.5/(x1 - 1.0)^2

    # discretization control system 
    for k in 0:N-1
        @constraint(model, x1[k+1] == x1[k] + hstep/2 * (
            f1(x1[k],x2[k],x3[k],x4[k],u[k]) + f1(x1[k+1],x2[k+1],x3[k+1],x4[k+1],u[k+1])))
        @constraint(model, x2[k+1] == x2[k] + hstep/2 * (
            f2(x1[k],x2[k],x3[k],x4[k],u[k]) + f2(x1[k+1],x2[k+1],x3[k+1],x4[k+1],u[k+1])))
        @constraint(model, x3[k+1] == x3[k] + hstep/2 * (
            f3(x1[k],x2[k],x3[k],x4[k],u[k]) + f3(x1[k+1],x2[k+1],x3[k+1],x4[k+1],u[k+1])))
        @constraint(model, x4[k+1] == x4[k] + hstep/2 * (
            f4(x1[k],x2[k],x3[k],x4[k],u[k]) + f4(x1[k+1],x2[k+1],x3[k+1],x4[k+1],u[k+1])))
    end

    # objective function 
    @objective(model, Min, (x3[N]+x2[0]) * x4[N])

    # direct optimization via ipopt
    optimize!(model)

    println("Status: ", termination_status(model))
    println("Objective = ", objective_value(model))
    println("x0 = h(0)  = ", value(x0))
    println("y0 = h'(0) = ", value(y0))
    println("Number of iterations",MOI.get(model,MOI.BarrierIterations()))

    return model, value.(x1), value.(x2), value.(x3), value.(x4), value.(u), value(x0), value(y0)
end

# execution
model, x1, x2, x3, x4, u, x0, y0 = solve_mayer_free(N=500)
