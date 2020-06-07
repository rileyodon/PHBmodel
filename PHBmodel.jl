using Plots
using DifferentialEquations

plotlyjs()

mutable struct stream
    """State properties
    """
    q
    Ss
    Snh
    Xb
    Xphb
end

function kinetics(s)
    """Returns the rate constants for the growth (r1) and accumulation (r2) reaction, given the state.
    """
    
    μₘₐₓ_g = 4.56;  # d-1	Maximum biomass growth rate
    μₘₐₓ_p = 3.6;  # d-1	Maximum PHA production rate
    KSₛ_g = 2.98;  # kg/m3	half sat coeff for growth
    KSₛ_p = 9.3;  # kg/m3	half sat coeff for PHA production
    KSₙₕ_g = 0.021;  # gN/m3	Ammonia sat coeff for growth
    KIₙₕ_p = 0.021;  # gN/m3	Ammonia inhib coeff for PHA uptake
    KSᵢ = 2.68; # kg/m3 substate inhibition coef for microbe activity
    PHAₘₐₓ = 0.8;  # max fraction of PHA content

    V, Ss, Snh, Xb, Xphb = s; # extract variables from state
    
    Xₚₕₐ = Xphb / (Xb + Xphb);

    r1 = μₘₐₓ_g * Xb * Ss / (KSₛ_g + Ss) * Snh / (KSₙₕ_g + Snh); # [kg/m3.d] growth
    r2 = μₘₐₓ_p * Xb * Ss / (KSₛ_p + Ss) * KIₙₕ_p / (KIₙₕ_p + Snh) * PHAₘₐₓ / (PHAₘₐₓ + Xₚₕₐ) * KSᵢ / (KSᵢ + Ss); # [kg/m3.d] accumulation

    r1, r2
end


function pha_reactions!(ds, s, p, t)
    V, Ss, Snh, Xb, Xphb = s
    IN, OUT = p
    
    r1, r2 = kinetics(s)

    # mass stoichiometry
    v = [
    -2 -0.71740413 -0.150442478 1 0 0.986430678; # growth
    -1.666666667 -0.103359173 0 0 1 0.397932817; # accumulation
    ];

    ds[1] = dV = IN.q - OUT.q # reactor volume
    ds[2] = dSs = IN.q / V * (IN.Ss - Ss) + v[1, 1] * r1 + v[2, 1] * r2 # Substate
    ds[3] = dSnh = IN.q / V * (IN.Snh - Snh) + v[1, 3] * r1 + v[2, 3] * r2 # Ammonia
    ds[4] = dXb = IN.q / V * (IN.Xb - Xb) + v[1, 4] * r1 + v[2, 4] * r2 # bacteria
    ds[5] = dXphb = IN.q / V * (IN.Xphb - Xphb) + v[1, 5] * r1 + v[2, 5] * r2 # PHB

    ds

end


function batch_feed_policy(t, F_Ss, Q)
    """
    t, time [days]
    F_Ss, Fresh feed glucose concentration [kg/m3]
    Q, fresh feed rate [L/d ]
    """
    Ss = F_Ss 
    Snh = 0  
    if t < 1
        # first day conditions for growth
        Snh = 10.0 # Ammonia for first day
        Ss = 140.00  # g/L in both feeds.
    end
    stream(Q, Ss, Snh, 0, 0)
end

function no_flow(t)
    """Outlet condition of zero flow
    """
    stream(0, 0, 0, 0, 0)
end


function fedBatchConfig(S0, tspan, F_Ss, Q = 2.0)
    inlet = t->batch_feed_policy(t, F_Ss, Q)
    outlet = t->no_flow(t)
    p = (inlet, outlet)
    prob = ODEProblem(batchSystem!, S0, tspan, p) 
    sol = solve(prob, Rodas4())
end

function batchConfig(S0, tspan)
    inlet = t->no_flow(t)
    outlet = t->no_flow(t)
    p = (inlet, outlet)
    prob = ODEProblem(batchSystem!, S0, tspan, p) 
    sol = solve(prob, Rodas4())
end

function batchSystem!(ds, S, P, t)
    p = (P[1](t), P[2](t))
    ds = pha_reactions!(ds, S, p, t)
end

function chemostat_reactor(ds, s, p, t)
    
    ds = pha_reactions!(ds, s, p, t)
    
    in, out = p
    
    out = stream(out.q, s[2], s[3], s[4], s[5])

    ds, out
end

function chemostat_system!(DS, S, P, t)
    stream_states = 5
    IN, OUT = P
    num_reactors = length(IN)
    s = S[1:stream_states]
    ds = DS[1:stream_states]
    p = (IN[1], stream(IN[1].q, 0, 0, 0, 0))
    ds, out = chemostat_reactor(ds, s, p, t)
    DS[1:stream_states] = ds   

    for i = 2:num_reactors
        in = mixer(IN[i], out)
        s = S[stream_states * (i - 1) + 1:stream_states * i]
        ds = DS[stream_states * (i - 1) + 1:stream_states * i]
        p = (in, stream(in.q, 0, 0, 0, 0))
        ds, out = chemostat_reactor(ds, s, p, t)
        DS[stream_states * (i - 1) + 1:stream_states * i] = ds
    end
    DS
end

function mixer(m1, m2)
    q = m1.q + m2.q

    Ss = (m1.q * m1.Ss + m2.q * m2.Ss) / q
    Snh = (m1.q * m1.Snh + m2.q * m2.Snh) / q
    Xb = (m1.q * m1.Xb + m2.q * m2.Xb) / q
    Xphb = (m1.q * m1.Xphb + m2.q * m2.Xphb) / q

    stream(q, Ss, Snh, Xb, Xphb)
end

function gen_Qfresh_Array(QF, R, n)
    Rs = [R^(i - 2) for i = 2:n]
    base_R = sum(Rs)
    Qfreshs = QF .* Rs ./ base_R
end

function process_flow_config_n(V, tspan, S0, n, QF = 1.0, R = 1.0)

    Qfresh = QF / (n - 1)
    IN = [stream(2.0, 140.00, 10.0, 0.0, 0.0)]

    Qfreshs = gen_Qfresh_Array(QF, R, n)

    for Qfresh in Qfreshs
        push!(IN, stream(Qfresh, 140.0, 0.0, 0.0, 0.0))
    end
    OUT = stream(Qfresh, 0, 0, 0, 0)
    P = (IN, OUT)
    prob = ODEProblem(chemostat_system!, S0, tspan, P)
    sol = solve(prob, Rodas4())
    # gui()
    # plot(sol)
    return sol
end

function make_S0(V, s)
    base_S0 = [V, 1e-6, 1e-6, 1e-6, 1e-6];
    S0 = [];
    for i = 1:s
        append!(S0, base_S0);
    end
    S0;
end

function initial_HRT_plot(V_span, S0s, labels, step = 0.1)
    rs = length(S0s)
    Cs = Array{Float64}(undef, 0, rs);

    Vs = [V for V = V_span[1]:step:V_span[2]];

    QF = 0.75; # baseline fresh feed rate
    R = 1.0; # baseline fresh feed distribution ratio

    for V in Vs # calculate for each system volume
        tspan = (0.0, 100.0);
        C = Array{Float64}(undef, 1, rs); # array to store steady state PHB conc in

        for (i, S0) in enumerate(S0s) # calculate for each number of accumulation reactors
            num_r = length(S0) / 5.0
            Vr = V / num_r; # accumulation reactor volume
            for k = 6:5:length(S0)
                S0[k] = Vr # specify accumulation reactor volums
            end
            sol = process_flow_config_n(V, tspan, S0, num_r, QF, R);
            
            # extract steady state PHB concentration              
            Xphb = sol.u[end][end];
            C[i] = Xphb;

        end

        Cs = [Cs; C];
    end

    sysVol = Vs .+ S0s[1][1];
    sysInletQ = 2 + QF;
    HRT = sysVol / sysInletQ;
    yield = (Cs .* sysInletQ) ./ (sysInletQ * 140.00)

    plot(HRT, yield, xlabel = "Hydraulic Retention Time (d)", ylabel = "Yield (g PHB / g Glucose)", label = labels, lw = 3)
    # savefig("chemostat_trains.svg")

    # add fed batch to plot
    tspan = (0.0, 8.0)
    S0 = [0.5, 140.0, 10.0, 5.0, 0.0]
    F_Ss = 140.0;
    fedbatch_sol = fedBatchConfig(S0, tspan, F_Ss)
    HRT = (fedbatch_sol[1,:] .- 0.5) ./ 2;
    yield = (fedbatch_sol[5,:] .* fedbatch_sol[1,:]) ./ (S0[1] * S0[2] .+ F_Ss .* fedbatch_sol[1,:])
    plot!(HRT, yield, label = "FB(Qf=2)", lw = 3, linestyle = :dot)

    # basic batch
    tspan = (0.0, 8.0)
    S0 = [10, 140.0, 10.0, 5.0, 0.0]
    batch_sol = batchConfig(S0, tspan)
    yield = (batch_sol[5,:] .* batch_sol[1,:]) ./ (S0[2] * S0[1])

    plot!(batch_sol.t, yield, label = "B", lw = 3, linestyle = :dot)
    # savefig("Reactor Configuration Comparison.svg")

    # concentration profiles for batch reactors
    gui()
    plot(batch_sol, title = "Batch Reactor")
    gui()
    plot(fedbatch_sol, title = "Fed Batch")


end



# Chemostat initial conditions
V0 = 2;
S01 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0];
S02 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0];
S03 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0];
S04 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0];
S05 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0];

labels = [ "C(n=1, Qf=0.75, R=1)" "C(n=2, Qf=0.75, R=1)" "C(n=3, Qf=0.75, R=1)" "C(n=4, Qf=0.75, R=1)" "C(n=5, Qf=0.75, R=1)" ];


S0s = [S01, S02, S03, S04, S05]

######################
# initial assessment #
######################
initial_HRT_plot([0.1, 20.00], S0s, labels)


################
# optimisation #
################

# fed batch conc comparison:
tspan = (0.0, 11.0)
S0 = [0.5, 140.0, 10.0, 5.0, 0.0]

function FB_glucose_feed_comp()
    # gui()
    # plot()
    for F_Ss in [10.0, 20.0, 40.0, 70.0, 100.0, 140.0]
        # gui()
        gui()
        fedbatch_sol = fedBatchConfig(S0, tspan, F_Ss);
        HRT = (fedbatch_sol[1,:] .- 0.5) ./ 2;
        yield = (fedbatch_sol[5,:] .* fedbatch_sol[1,:]) ./ (S0[1] * S0[2] .+ F_Ss .* (fedbatch_sol[1,:] .- S0[1]));
        plot!(HRT, fedbatch_sol[5,:], label = "Feed Glucose $F_Ss g/L", lw = 3, palette = :Spectral_8)
    end
    plot!(xlabel = "Hydraulic Retention Time (d)", ylabel = "PHB Concentration (g / L)")
end
gui()
plot()
FB_glucose_feed_comp()


function FB_feed_rate()
    # gui()
    # plot()
    F_Ss = 140.0
    for Q in [0.0, 0.05, 0.10, 0.15, 0.25,  0.5, 1.0, 2.0]
        
        # gui()
        gui()
        fedbatch_sol = fedBatchConfig(S0, tspan, F_Ss, Q);
        HRT = (fedbatch_sol[1,:] .- 0.5) ./ Q;
        yield = (fedbatch_sol[5,:] .* fedbatch_sol[1,:]) ./ (S0[1] * S0[2] .+ F_Ss .* (fedbatch_sol[1,:] .- S0[1]));
        plot!(HRT, yield, label = "FB(Qf=$Q L/d)", lw = 3, palette = :Spectral_8)
    end
    plot!(xlabel = "Hydraulic Retention Time (d)", ylabel = "Yield (g PHB / g Glucose)")
    # savefig("FB_feed_rate.svg")
end
gui()
plot()
FB_feed_rate()


function C_fresh_feed_rate()
    # gui()
    # plot()
    step = 0.1

    V0 = 2.0;
    S0 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0]
    

    V_span = [0.1, 40.00];

    Vs = [V for V = V_span[1]:step:V_span[2]];

    # QFs = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    QFs = [0.85, 1.0, 1.15, 1.25, 1.35, 1.5]

    rs = length(QFs)
    Cs = Array{Float64}(undef, 0, rs);
    yields = Array{Float64}(undef, 0, rs);
    HRTs = Array{Float64}(undef, 0, rs);

    for V in Vs
        tspan = (0.0, 100.0);
        C = Array{Float64}(undef, 1, rs);
        yield = Array{Float64}(undef, 1, rs);
        HRT = Array{Float64}(undef, 1, rs);

        for (i, QF) in enumerate(QFs)
            # S0 = make_S0(V, length(S0s[i]));
            # println(S0)
            num_r = length(S0) / 5.0
            Vr = V / num_r;
            for k = 6:5:length(S0)
                S0[k] = Vr
            end
            sol = process_flow_config_n(V, tspan, S0, num_r, QF);
            
            Xphb = sol.u[end][end];
            C[i] = Xphb;
            yield[i] = Xphb * (2.0 + QF) / ((2.0 + QF) * 140.0)
            sysVol = V + S0[1];
            HRT[i] = sysVol / (2.0 + QF);
            # S0s[i] = sol.u[end]

        end

        Cs = [Cs; C];
        yields = [yields; yield]
        HRTs = [HRTs; HRT]
    end

    labels = Array{String}(undef, 1, rs);
    for (i, QF) in enumerate(QFs)
        labels[i] = "C(n=3, Qf=$QF, R=1)"
    end

    plot(HRTs, yields, xlabel = "Hydraulic Retention Time (d)", ylabel = "Yield (g PHB / g Glucose)", xlims = (0, 7.5), label = labels, lw = 3, palette = :Spectral_6)
    # savefig("chemostat_feed_rate.svg")
end
gui()
plot()
C_fresh_feed_rate()


function C_fresh_feed_distribution()
    # gui()
    # plot()
    step = 0.1

    V0 = 2.0;
    S0 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0]
    

    V_span = [0.1, 40.00];

    Vs = [V for V = V_span[1]:step:V_span[2]];

    Rs = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

    rs = length(Rs)
    Cs = Array{Float64}(undef, 0, rs);
    yields = Array{Float64}(undef, 0, rs);
    HRTs = Array{Float64}(undef, 0, rs);

    QF = 1.25;

    for V in Vs
        tspan = (0.0, 100.0);
        C = Array{Float64}(undef, 1, rs);
        yield = Array{Float64}(undef, 1, rs);
        HRT = Array{Float64}(undef, 1, rs);

        for (i, R) in enumerate(Rs)
            # S0 = make_S0(V, length(S0s[i]));
            # println(S0)
            num_r = length(S0) / 5.0
            Vr = V / num_r;
            for k = 6:5:length(S0)
                S0[k] = Vr
            end
            sol = process_flow_config_n(V, tspan, S0, num_r, QF, R);
                        
            Xphb = sol.u[end][end];
            C[i] = Xphb;
            yield[i] = Xphb * (2.0 + QF) / ((2.0 + QF) * 140.0)
            sysVol = V + S0[1];
            HRT[i] = sysVol / (2.0 + QF);
            # S0s[i] = sol.u[end]

        end

        Cs = [Cs; C];
        yields = [yields; yield]
        HRTs = [HRTs; HRT]
    end

    labels = Array{String}(undef, 1, rs);
    for (i, R) in enumerate(Rs)
        Qfreshs = gen_Qfresh_Array(QF, R, length(S0) / 5)
        Q1 = Qfreshs[1]

        labels[i] = "C(n=3, Qf=$QF, R=$R)"
    end

    plot(HRTs, yields, xlabel = "Hydraulic Retention Time (d)", ylabel = "Yield (g PHB / g Glucose)", xlims = (0, 7.5), label = labels, lw = 3, palette = :Spectral_6)
    # savefig("chemostat_feed_ratio.svg")
end
gui()
plot()
C_fresh_feed_distribution()

###
# Best against original #
###

function chemstat_result_com(QF = 1.25, R = 1.0)

    step = 0.1

    V0 = 2.0;
    S0 = [V0, 140.0, 10.0, 100.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0, V0, 40.0, 0.0, 10.0, 0.0]

    V_span = [0.1, 40.00];

    Vs = [V for V = V_span[1]:step:V_span[2]];

    Rs = [R]

    rs = length(Rs)
    Cs = Array{Float64}(undef, 0, rs);
    yields = Array{Float64}(undef, 0, rs);
    HRTs = Array{Float64}(undef, 0, rs);


    for V in Vs
        tspan = (0.0, 100.0);
        C = Array{Float64}(undef, 1, rs);
        yield = Array{Float64}(undef, 1, rs);
        HRT = Array{Float64}(undef, 1, rs);

        for (i, R) in enumerate(Rs)

            num_r = length(S0) / 5.0
            Vr = V / num_r;
            for k = 6:5:length(S0)
                S0[k] = Vr
            end
            sol = process_flow_config_n(V, tspan, S0, num_r, QF, R);
            
            
            Xphb = sol.u[end][end];
            C[i] = Xphb;
            yield[i] = Xphb * (2.0 + QF) / ((2.0 + QF) * 140.0)
            sysVol = V + S0[1];
            HRT[i] = sysVol / (2.0 + QF);
            # S0s[i] = sol.u[end]

        end

        Cs = [Cs; C];
        yields = [yields; yield]
        HRTs = [HRTs; HRT]
    end

    labels = Array{String}(undef, 1, rs);
    for (i, R) in enumerate(Rs)
        Qfreshs = gen_Qfresh_Array(QF, R, length(S0) / 5)
        Q1 = Qfreshs[1]

        labels[i] = "C(n=3, Qf=$QF, R=$R)"
    end

    return HRTs, yields, labels
end

gui()
F_Ss = 140.0
Q = 2.0
fedbatch_sol = fedBatchConfig(S0, tspan, F_Ss, Q);
HRT = (fedbatch_sol[1,:] .- 0.5) ./ Q;
yield = (fedbatch_sol[5,:] .* fedbatch_sol[1,:]) ./ (S0[1] * S0[2] .+ F_Ss .* (fedbatch_sol[1,:] .- S0[1]));
plot(HRT, yield, label = "FB(Qf=$Q)", lw = 3, linestyle = :dot, linecolor = RGB(50 / 255, 136 / 255, 189 / 255))

Q = 0.25
fedbatch_sol = fedBatchConfig(S0, tspan, F_Ss, Q);
HRT = (fedbatch_sol[1,:] .- 0.5) ./ Q;
yield = (fedbatch_sol[5,:] .* fedbatch_sol[1,:]) ./ (S0[1] * S0[2] .+ F_Ss .* (fedbatch_sol[1,:] .- S0[1]));
plot!(HRT, yield, label = "FB(Qf=$Q)", lw = 3, linecolor = RGB(50 / 255, 136 / 255, 189 / 255))

HRTs, yields, labels = chemstat_result_com(0.75, 1.0)
plot!(HRTs, yields, xlims = (0, 7.5), label = labels, lw = 3, linestyle = :dot, linecolor = RGB(213 / 255, 62 / 255, 79 / 255))

HRTs, yields, labels = chemstat_result_com(1.25, 0.8)
plot!(HRTs, yields, xlims = (0, 7.5), label = labels, lw = 3, linecolor = RGB(213 / 255, 62 / 255, 79 / 255))


plot!(xlabel = "Hydraulic Retention Time (d)", ylabel = "Yield (g PHB / g Glucose)")


# savefig("optimised_improvement.svg")