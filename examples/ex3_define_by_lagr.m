function [penm] = ex3_define_by_lagr()
% define Example 3 as is necessary to handle by PenLab
% same as ex3_define() but this time using Hessian of the Lagrangian
% call:
%   penm=ex3_define_by_lagr();    % to define structure 
%   prob=penlab(penm);    % to convert the structure and initialize the problem
%   prob.opts....=...;    % to change option settings if desired
%   prob.solve();         % to start the solver
%   prob.x                % to retrieve the final point (solution)

  penm = [];

  % [optional] set problem name/comment (for log files)
  penm.probname = 'examples/ex3';
  penm.comment = 'Source: user external definition of functions, using Hessian of the Lagrangian';

  penm.Nx = 2;
  penm.lbx = [-0.5; -Inf];
  penm.ubx = [0.5; 1];

  penm.NgNLN = 2;
  %penm.NgLIN = 0;          % optional, 0 by default
  penm.lbg = [0; 0];
  %penm.ubg = [Inf, Inf];   % optional, +Inf for upper bounds by default

  % possible even as anonymous?? Try it
  penm.objfun = @ex3_objfun;
  penm.confun = @ex3_confun;
  penm.objgrad = @ex3_objgrad;
  penm.congrad = @ex3_congrad;
  %penm.objhess = @ex3_objhess;
  %penm.conhess = @ex3_conhess;
  penm.lagrhess = @ex3_lagrhess;

  % [optional] set starting point
  %penm.x = zeros(2,1);
  penm.x = [-2; 1];        % suggested starting point
  
  end

function [f,userdata] = ex3_objfun(x,Y,userdata)
% objective function values for Example 3

  f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
  
end

function [g,userdata] = ex3_confun(x,Y,userdata)
% function values for Example 3

  g = zeros(2,1);
  g(1) = x(1)^2 + x(2);
  g(2) = x(1)   + x(2)^2;
  
end

function [df, userdata]=ex3_objgrad(x,Y,userdata)
% Gradients for Example 3, note that they are stored in columns!

  df = [-400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1)); 200*(x(2)-x(1)^2)];
  
end

function [dg, userdata]=ex3_congrad(x,Y,userdata)
% Gradients for Example 3, note that they are stored in columns!

  dg = [ 2*x(1), 1; 1, 2*x(2)];
  
end

function [ddL, userdata] = ex3_lagrhess(x,Y,v,userdata)
% Hessian of the lagrangian, Example 3
%   H=nabla^2 f + sum v_i nabla^2 g_i
% v should be of size NgNLN, thus in this example NgNLN=2

  hess_f = [1200*x(1)^2-400*x(2)+2, -400*x(1); -400*x(1), 200];
  hess_g1 = [2, 0; 0, 0];
  hess_g2 = [0, 0; 0, 2];
  ddL = hess_f + v(1)*hess_g1 + v(2)*hess_g2;
  
end



