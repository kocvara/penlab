function [penm] = ex1_define()
% define Example 1 as is necessary to handle by PenLab
% call:
%   penm=ex1_define();    % to define structure
%   prob=penlab(penm);    % to convert the structure and initialize the problem
%   prob.opts....=...;    % to change option settings if desired
%   prob.solve();         % to start the solver
%   prob.x                % to retrieve the final point (solution)

penm = [];
% [optional] set problem name/comment (for log files)
penm.probname = 'examples/ex1';
penm.comment = 'Source: user external definition of functions';

penm.Nx = 3;
penm.lbx = zeros(3,1);   % allow 0 as well

penm.NgNLN = 1;
penm.NgLIN = 1;
penm.lbg = [4, 24];   % would column vector work?
penm.ubg = [Inf, 24];

%
penm.objfun = @ex1_objfun;
penm.confun = @ex1_confun;
penm.objgrad = @ex1_objgrad;
penm.congrad = @ex1_congrad;
penm.objhess = @ex1_objhess;
penm.conhess = @ex1_conhess;
%penm.lagrhess = @ex1_lagrhess;

end

function [f,userdata] = ex1_objfun(x,Y,userdata)
% objective function values for Example 1

f = x(1)^2 + 4*x(2)^2 - x(3)^2 + x(1)*x(2) - 2*x(1)*x(3);

end

function [g,userdata] = ex1_confun(x,Y,userdata)
% function values for Example 1

g = zeros(2,1);
g(1) = x(1)^2 + x(2)^2 + x(3)^2;
g(2) = 2*x(1) + 6*x(2) + 4*x(3);

end

function [df, userdata]=ex1_objgrad(x,Y,userdata)
% Gradients for Example 1, note that they are stored in columns!

df = [ 2.*x(1) + x(2) - 2*x(3); 8.*x(2) + x(1); -2.*x(3) - 2.*x(1)];

end

function [dg, userdata]=ex1_congrad(x,Y,userdata)
% Gradients for Example 1, note that they are stored in columns!

dg = [ [ 2*x(1); 2*x(2); 2*x(3)] , [2; 6; 4]];

end

function [ddf, userdata] = ex1_objhess(x,Y,userdata)
% Hessians of the objective function, Example 1

ddf = [2, 1, -2; 1, 8, 0; -2, 0, -2];

end

function [ddgk, userdata] = ex1_conhess(x,Y,k,userdata)
% Hessians of the constraints, Example 1

ddgk = 2.*eye(3,3);

end






