% We solve the problem of stability of a time invariant linear (discrete)
% system using Lyapunov theory. The problem is formulated as a LMI system
% A'*P + P*A < 0 , P > I (continuous time)
% or
% A'*P*A - P < 0 , P > I (discrete time)
% We minimize the trace of P.

function penm = test_Lyapunov()

A=[0 1 0;0 0 1;-1 -2 -3];
penm = lyapu(A);
