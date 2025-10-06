
# ------------------- READ ME---------------------
# This is the Julia code for the Dynamic Coupling Rotor (DCR) model of muscle contraction, as described in the manuscript.
# The model is replete with discontinuities, non-linearities and stiff equations, making it a challenging numerical problem.
# There are several parameters that can be adjusted to help the solver along, but these should not be changed if the solver is working fine.
# The key parameters to adjust are: 

# 1) switch - this changes how the rotors' force is calculated to make them less backdrivable. This *can not* be turned on when looking for oscillatory motion, only unidirectional.
#
# 2) events - this turns on the event framework that detects when a rotor contacts a tooth and enforces the coupling constraint. 
#             This will make the solver slower, but is effective.
# 3) solver tolerances and dtmax - Making these smaller will increase accuracy, but also typically increase solve time.


# Indications the solver is struggling are;
#      - solve times that seem to keep increasing (check the progress bar in the bottom left of the REPL)
#      - warnings about the maximum number of iterations being reached (convergence failures)
#      - results that show a system that synchronises, then quickly desynchronises and remains disordered (synchronisation should be stable).

# Typically parameter choices with low alpha values (viscous resistance) and high gamma values (isotonic force) are more challenging for the solver.

# ------------------- END READ ME---------------------

# These packages allow the use of the required solvers - the SUNDIALS suite of solver packages is particularly useful for solving non-linear and discontinuous equations.
# Makie is a better visualisation package than Plots, the base julia data visualisation package.


using DifferentialEquations
using GLMakie
using Sundials
using ProgressLogging
using FFTW
GLMakie.activate!();

# alpha is viscous resistance, beta is auxotonix boundary constraint, gamma is isotonic force.

alpha = 2;
beta = 0.00;
gamma = 2;

# If you're struggling to get the solver to work and you're only interested in unidirectional motion, try changing switch from 0 to 1. 
# It changes how the rotors' force is calculated to make them less backdrivable, which creates fewer runaway events. This *can not* be 
# turned on when looking for oscillatory motion. Switch can be increased to up to 10 or so for very low alpha, high gamma, effectively 
# preventing the backbone moving in reverse.

switch = 0;

# Events may make the system solve, but it will take much much longer. Leave this as false if you can. 

events = false;

# StS is the same as `c' in the paper. This must be between 0 and 1 and represents the range of oscillator contact. Probably don't make it 1, as you'll end up dividing by zero.

StS = 0.14;

# Rotor and tooth period is set as required. The preset is physiologically determined from 38nm and 43nm

rot_per = 8.6*StS
tooth_per = 7.8*StS

# Number of rotors/motors and the time (time is dimensionless and measured as a function of the natural period of an oscillator)

rot_num = 100;
time = 200;



# Heavyside function allowing detatchment and determining if a motor *could* be contacted

function H(theta, StS)

    if mod(theta, 2*pi) > pi - asin(StS) && mod(theta, 2*pi) < pi + asin(StS)

        heavy = 1

    else

        heavy = 0

    end

    return heavy

end



# f2 is the non-linear function to be solved. The IDA DAE will try to minimise the residuals, kept in the out array

function f2(out, du, u, p, t)

    alpha = p[1]
    beta = p[2]
    StS = p[4]
    rot_per = p[5]
    tooth_per = p[6]
    rot_num = Int(p[7])
    gamma = p[8]

    cont = zeros(rot_num)

    tooth_num = Int(ceil((rot_num*rot_per)/tooth_per));
    tooth_orig_x = zeros(tooth_num);
    tooth_loc_x = zeros(tooth_num);
    tooth_loc_ind = zeros(tooth_num);
    @inbounds for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end


    rot_forces = zeros(rot_num)

    @inbounds for i = 1:rot_num
        if H(u[i+1], StS) == 1  # Check H first to avoid unnecessary loops
            @inbounds for j = 1:tooth_num
                if abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1]))) < 0.02
                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]              
                    cont[i] = 1
                    break  # Exit once contact found
                end
            end
        end
    end


    @inbounds for i = 1:rot_num
        if cont[i] ==0
            out[i+1] = 2*pi - du[i+1]
        end
    end






    # This is the normal version of rot forces.
  
    @inbounds for i = 1:1:rot_num
        
        if switch == 0
            atten =  (du[i+1]/(2*pi))*(0.5*tanh(-40*du[i+1])+0.5);
        else
            atten = -(du[i+1]/(2*pi))*(0.5*tanh(-40*du[i+1]-2)+0.5)*(switch-1);
        end

        rot_forces[i] = ((1 - du[i+1]/(2*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS);

    end

    # This version of rot forces is for when attenuating the maximum motor force to better mimic electric motors.

    # for i = 1:1:rot_num

    #     velad = du[i+1]-1.5*pi;
        
    #     atten =  (velad/(0.5*pi))*(0.5*tanh(-80*velad)+0.5);

    #     rot_forces[i] = ((1 - velad/(0.5*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS);

    # end


      
       spring = beta


    out[1] = sum(rot_forces)  - alpha*du[1] - spring*(u[1])- gamma*(0.5*tanh(2*(u[1]-1.3))+0.5)


    

    



end

#f2_smooth is a modified version of f2 that applies a damping effect to rotors as they approach teeth, rather than enforcing an immediate velocity constraint. This can help with numerical stability and convergence in some scenarios.

function f2_smooth(out, du, u, p, t)
    alpha = p[1]
    beta = p[2] 
    StS = p[4]
    rot_per = p[5]
    tooth_per = p[6]
    rot_num = Int(p[7])
    gamma = p[8]

    tooth_num = Int(ceil((rot_num*rot_per)/tooth_per));
    tooth_orig_x = zeros(tooth_num);
    tooth_loc_x = zeros(tooth_num);
    tooth_loc_ind = zeros(tooth_num);
    @inbounds for j = 1:tooth_num
        tooth_orig_x[j] = j*tooth_per - tooth_per;
        tooth_loc_x[j] = tooth_orig_x[j] + u[1];
        tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num);
    end

    
    rot_forces = zeros(rot_num)
    
    @inbounds for i = 1:rot_num
        # Calculate proximity to nearest tooth
        min_distance = Inf
        in_contact = false
        
        if H(u[i+1], StS) == 1
            @inbounds for j = 1:tooth_num
                distance = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1])))
                min_distance = min(min_distance, distance)
                
                if distance < 0.02
                    in_contact = true
                    # Rotor is in contact - enforce coupling constraint
                    out[i+1] = du[1]/abs(cos(u[i+1])) - du[i+1]
                    break
                end
            end
        end
        
        if !in_contact
            # Apply proximity-based damping for approaching rotors
            target_velocity = 2*pi
            if H(u[i+1], StS) == 1 && min_distance < 0.05  # Buffer zone
                damping_factor = min_distance / 0.1  # 0 at contact, 1 at buffer edge
                target_velocity = 2*pi * (0.3 + 0.7 * damping_factor)  # Slow to 30% of normal
            end
            
            out[i+1] = target_velocity - du[i+1]
        end
        
        # Calculate forces
        atten = (du[i+1]/(2*pi))*(0.5*tanh(-40*du[i+1])+0.5)
        rot_forces[i] = ((1 - du[i+1]/(2*pi))+atten)*abs(cos(u[i+1]))*H(u[i+1], StS)
    end
    
    # Backbone equation
    spring = beta
    out[1] = sum(rot_forces) - alpha*du[1] - spring*(u[1]) - gamma*(0.5*tanh(2*(u[1]-1.3))+0.5)
end

function condition1_safe(out, u, t, integrator)
    try
        rot_per = integrator.p[5]
        tooth_per = integrator.p[6]
        rot_num = Int(integrator.p[7])
        StS = integrator.p[4]
       
        tooth_num = Int(ceil((rot_num*rot_per)/tooth_per))
        tooth_orig_x = zeros(tooth_num)
        tooth_loc_x = zeros(tooth_num)
        tooth_loc_ind = zeros(tooth_num)
        
        @inbounds for j = 1:tooth_num
            tooth_orig_x[j] = j*tooth_per - tooth_per
            tooth_loc_x[j] = tooth_orig_x[j] + u[1]
            tooth_loc_ind[j] = mod(tooth_loc_x[j], tooth_per*tooth_num)
        end

        @inbounds for i = 1:rot_num
            out[i] = Inf
            
            if H(u[i+1], StS) == 1  && integrator.du[i+1] > 0.9*2*pi
                @inbounds for j = 1:tooth_num
                    distance = abs(tooth_loc_ind[j] - ((i-1)*rot_per + sin(-u[i+1])))
                    out[i] = min(out[i], distance - 0.023)
                end
            else
                out[i] = 1.0
            end
            
            # Safety check for problematic values
            if !isfinite(out[i]) || isnan(out[i])
                println("WARNING: Non-finite value in condition at t=$(integrator.t), i=$i")
                out[i] = 1.0
            end
        end
        
    catch e
        println("ERROR in condition1 at t=$(integrator.t): $e")
        fill!(out, 1.0)  # Safe fallback
    end
end


function affect_safe!(integrator, idx)

    StS = integrator.p[4]
    rot_num = Int(integrator.p[7])
    spring = integrator.p[2]
    alpha = integrator.p[1]
    gamma = integrator.p[8]  

    try
        
        @inbounds for i = 1:1:rot_num
            if idx == i
                old_du = integrator.du[i+1]
                cos_val = cos(integrator.u[i+1])
                
                # Safety check for division by zero
                if abs(cos_val) < 1e-12
                    println("WARNING: Near-zero cos value at t=$(integrator.t), rotor $i")
                    cos_val = sign(cos_val) * 1e-12
                end
                
                new_du = 0.5*(abs(old_du - (integrator.du[1]/abs(cos_val))))
                
                
                integrator.du[i+1] = new_du
                println("Event: t=$(integrator.t), rotor $i, du: $old_du -> $new_du")
            end
        end

  
    catch e
        println("ERROR in affect! at t=$(integrator.t): $e")
    end
end



#VCC is the event framework. It exists to detect when a rotor contacts a tooth and enforce the coupling constraint. This can often cause more trouble than it's worth. Try 
# running f2_smooth without the callback if you're having trouble getting the solver to work.

#cb = VectorContinuousCallback(condition1, affect!, rot_num)#, terminate_integrator = false);
cb = VectorContinuousCallback(condition1_safe, affect_safe!, rot_num,
                             rootfind=true);

Force_Vel_Data = zeros(20, 3)


c = 1;
  
# This conditional allows you to start the system pre organised, if preffered  

if rot_per > tooth_per
    phase_diff =  mod(((tooth_per - rot_per)/tooth_per)*2*pi, 2*pi);
elseif tooth_per > rot_per
    phase_diff =   mod(((tooth_per -  rot_per)/tooth_per)*2*pi, 2*pi);
else
    phase_diff = 0;
end

# These are the rotor initial conditions the options for u0 are zerod, pre organised and random

u0 = zeros(rot_num+1);
du0 = zeros(rot_num+1);

#phase_diff = -0.2*2*pi

for i = 1:rot_num
    du0[i+1] = 2*pi;

    #u0[i+1] = 0;
    #u0[i+1] = mod(i*phase_diff, 2*pi);
    u0[i+1] = rand()*2*pi
end




function debug_callback(integrator)
    if integrator.t > 0.1  # Every 0.1 time units, print status
        println("t = $(integrator.t), max |u| = $(maximum(abs.(integrator.u)))")
        println("max |du| = $(maximum(abs.(integrator.du)))")
    end
    return false  # Don't terminate
end





# The numerical framework

tspan = (0.0, time)
debug_cb = DiscreteCallback((u,t,integrator) -> mod(t, 0.1) < integrator.dt, 
                           debug_callback, save_positions=(true,true))

differential_vars = zeros(rot_num+1)

@inbounds for i = 1:rot_num+1
    differential_vars[i] = true;
end

println("Starting integration with ", rot_num, " rotors")
println("Time span: ", tspan)

# P is the parameter tuple
p = (alpha, beta, switch, StS, rot_per, tooth_per, rot_num, gamma);


    

prob_smooth = DAEProblem(f2_smooth, du0, u0, tspan, p, differential_vars = differential_vars)
prob = DAEProblem(f2, du0, u0, tspan, p, differential_vars = differential_vars)
combined_cb = CallbackSet(cb, debug_cb)

if smooth == true
    prob = prob_smooth
end

if events == false
    sol = solve(prob, IDA(linear_solver  = :LapackDense), maxiters = 10^7,  dtmax = 1e-4, reltol = 1e-7, abstol = 1e-7, progress=true, progress_steps=1)
else
    sol = solve(prob, IDA(linear_solver  = :LapackDense), maxiters = 10^7, callback=combined_cb,  dtmax = 1e-4, reltol = 1e-7, abstol = 1e-7, progress=true, progress_steps=1)
end



# The rest of this is just various ways to manipulate the data and show some of the results.


tvals = (time/10000:time/10000:time);
uvals = sol.(tvals);
uvals = hcat(uvals...);


#displacements = [u[:] for u in sol.u];
#displacements = hcat(displacements...);

vel = zeros(length(tvals), length(uvals[:, 1]));
bound = zeros(length(tvals), length(uvals[1, :]));

@inbounds for j = 1:1:length(uvals[:, 1])
    @inbounds for i = 2:1:length(tvals)
        vel[i, j] = (uvals[j, i] - uvals[j, i-1])/(tvals[i] - tvals[i-1]);
        if vel[i, j] <  6.2
            bound[i, j] = 1;
        end
    end
end 
tot_bound = zeros(length(tvals), 1);
tot_mot_frac = zeros(length(tvals), 1);
@inbounds for i = 1:1:length(tvals)

    tot_bound[i] = sum(bound[i, 2:end]);
    tot_mot_frac[i] = tot_bound[i]/(rot_num);

end
av_vel = zeros(length(tvals))
av_mot_frac = zeros(length(tvals))
window = 100;
for i = 1:1:length(tvals)-window
    av_vel[i] = sum(vel[i:i+window, 1])/window;
    av_mot_frac[i] = sum(tot_mot_frac[i:i+window])/window;
end





println("Motor fraction = ", sum(tot_mot_frac[8001:10000])/2000);
println("Average velocity = ", sum(vel[8001:10000])/2000); 



phase_diff = zeros(rot_num-1, 10000);
for j = 1:10000
    for i = 1:1:rot_num-1
        phase_diff[i, j] = mod(uvals[i+2, j] - uvals[i+1, j], 2*pi)
    end
end





av_phase_diff = sum(phase_diff[8000:end])/length(phase_diff[8000:end])
println("Average phase difference = ", av_phase_diff)




#phase_error = abs.(phase_diff .- phase_diff_cal)
phase_error = abs.(phase_diff .- av_phase_diff)

G_phase_error = zeros(length(phase_error[1, :]), 1)

for i = 1:1:length(phase_error[1, :])
    G_phase_error[i] = sum((phase_error[:, i]))/(rot_num-1)
end




G_phase_error_avg = sum(G_phase_error[8001:10000])/2000
display(G_phase_error_avg)

error = G_phase_error_avg
println("Error for gamma = ", gamma, " is ", G_phase_error_avg)


fig1 = Figure()
ax1 = Axis(fig1[1, 1], xlabel = "Time", ylabel = "Displacement", title = "Backbone Displacement")
ax2 = Axis(fig1[1, 2], xlabel = "Time", ylabel = "Motor Fraction", title = "Motor Fraction Over Time")
ax3 = Axis(fig1[2, 1], xlabel = "Time", ylabel = "Velocity", title = "Backbone Velocity")
ax4 = Axis(fig1[2, 2], xlabel = "Displacement", ylabel = "Phase Error", title = "Phase Error vs Displacement")


lines!(ax1, tvals[2:2:end], uvals[1, 2:2:end]) 

# lines!(ax1, tvals[2:2:end],  cos.(uvals[2, 2:2:end]))
lines!(ax2, tvals[10:10:end], tot_mot_frac[10:10:end, 1])
lines!(ax2, tvals[10:10:end], av_mot_frac[10:10:end, 1])
lines!(ax3, tvals[1:1:end], vel[1:1:end, 1])
lines!(ax4, uvals[1, 1:1:end], G_phase_error[1:1:end])
# #lines!(ax4, 1:60, error, color = :red)


 display(fig1)



