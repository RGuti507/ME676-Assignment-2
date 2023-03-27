#Course: Robotics
#Student: Ricardo Gutierrez

#Instructions:
#First: Run startup.jl
#Second: Run v1.jl

using GLM
using RigidBodyDynamics
using LinearAlgebra
using Plots
using PyPlot



function traj(t::Float64)
    # compute the desired joint angle at time t
    Vmax=(qd-q0)/8;
    T1=2;   
    T2=8;
    Tmax=10;
    y1= (t)->  q0+0.5*Vmax*t^2/T1
    q1=y1(T1);    
    y2=(t)-> q1+Vmax*(t-T1)
    q2=y2(T2);
    y3=(t)-> q2+Vmax*(t-T2)*Tmax/(Tmax-T2)-0.5*Vmax*(t^2-T2^2)/(Tmax-T2)
    q3=y3(Tmax)
    if t <=T1
        return q0+(0.5*Vmax*t^2)/T1, Vmax*t/T1, Vmax/T1
    elseif t > T1 && t <= T2
        return q1+Vmax*(t-T1), Vmax, zeros(9,1)
    elseif t > T2 && t <= Tmax
        return q2+Vmax*(t-T2)*Tmax/(Tmax-T2)-0.5*Vmax*(t^2-T2^2)/(Tmax-T2), Vmax*t/(T2-Tmax)-(Tmax*Vmax)/(T2-Tmax), Vmax/(T2-Tmax)
    elseif t > Tmax
        return q3, zeros(9,1), zeros(9,1)
    end
end

function control!(τ, t, state)
    # Set-point Regulation using PD
    qdes, qdes_dot, qdes_ddot = traj(t)
    τ .= -100*(velocity(state)-qdes_dot) - 500*(configuration(state) - qdes)
    act_sat = 50; # Actuator limits
    τ .= map( x -> x > act_sat ? act_sat : x,τ)
    τ .= map( x -> x < -act_sat ? -act_sat : x,τ)
    return τ
end

function control_CTC!(τ2, t, state)
    # Compute a value for τ
    qdes, qdes_dot, qdes_ddot = traj(t)
    aq=similar(velocity(state))    
    aq .=-100*(velocity(state)-qdes_dot) -500*(configuration(state) - qdes)
    τ2 .=inverse_dynamics(state,aq)
    # Saturate
    act_sat = 50; # Actuator limits
    τ2 .= map( x -> x > act_sat ? act_sat : x,τ2)
    τ2 .= map( x -> x < -act_sat ? -act_sat : x,τ2)
    return τ2
end

# Load mechanism info
delete!(vis)
#display_urdf("panda.urdf",vis)
urdfPath = "panda.urdf"
mvis, mechanism = display_urdf(urdfPath,vis)
mechanism = parse_urdf(Float64, urdfPath)
state = MechanismState(mechanism)

global q0= [0.01;-0.5;-0.0;-2.0;-0.3;1.5;-0.7;0.1;0.1];
global qd = [0.0;0.0;0.0;0.0;0.0;pi;0.01;0.01;0.01];
set_configuration!(state,q0)
zero_velocity!(state)
set_configuration!(mvis, configuration(state))



#-------------------------------------------------------------------------------------------------------
problem2 = ODEProblem(Dynamics(mechanism,control_CTC!), state, (0.00, 11.00));
solution2 = transpose(solve(problem2, Vern7()))*(180/3.14159);
fig1, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(8, 50))
plt.show()
time2=solve(problem2, Vern7()).t
ax1[1].plot(time2,solution2[:,1:end])
ax1[1].set_xlabel("Time-sec")
ax1[1].set_ylabel("Angle-Degrees")
ax1[1].text(0,150,"Time response using CTC Controller")
qend2=solution2[end,1:9];
ess2=norm(qend2*pi/180-qd);
ax1[1].text(8,100,"ess= "*string(round(ess2,digits=10)))
print("The norm of the Steady-State Error with CTC controller is: $ess2")

#-------------------------------------------------------------------------------------------------------
problem = ODEProblem(Dynamics(mechanism,control!), state, (0.00, 11.00));
solution = transpose(solve(problem, Vern7()))*(180/3.14159);
time1=solve(problem, Vern7()).t
ax1[2].plot(time1,solution[:,1:end])
ax1[2].set_ylabel("Angle-Degrees")
ax1[2].text(0,150,"Time response using PD Controller")
qend=solution[end,1:9];
ess=norm(qend*pi/180-qd);
ax1[2].text(8,100,"ess= "*string(round(ess,digits=5)))
print("The norm of the Steady-State Error with PD controller is: $ess %f")

