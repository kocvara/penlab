function [penm] = ncm1y_define()
% define penm structure for Nearest Correlation Matrix Example 2
% without the constrained on the condition number
% (section 7.2 of the original Pennon manual exactly as it is computed
% there, i.e., just one triangle in the objective)
% using matrix variables Y:
%
%    min_Y sum_{i<=j} (y_ij-h_ij)^2
%    s.t.  Y_ii=1
%          Y >= 0 (pos. semidefinite)
%
%


  penm = [];

  penm.probname= 'NCM Example 2 using matrix variables';
  penm.comment = 'matrix 6x6 as a matrix variable';

  % matrix H
  H = [ 1.00, -0.44, -0.20,  0.81, -0.46, -0.05;
       -0.44,  1.00,  0.87, -0.38,  0.81, -0.58;
       -0.20,  0.87,  1.00, -0.17,  0.65, -0.56;
        0.81, -0.38, -0.17,  1.00, -0.37, -0.15;
       -0.46,  0.81,  0.65, -0.37,  1.00,  0.08;
       -0.05, -0.58, -0.56, -0.15,  0.08,  1.00];

  % keep the whole structure
  penm.userdata=H;

  % no 'normal' variables
  penm.Nx=0;
  % one matrix variable of the given 3-diagonal structure (same as H)
  penm.NY=1;
  penm.Y{1}=H;

  %penm.Yinit{1}=zeros(6,6);
  penm.Yinit{1}= eye(6);

  % box constraints on matrix variables
  penm.lbY = [0];
  penm.ubY = [Inf];
  %penm.ubY = [10];

  % treat the box constraint by a log-barrier
  penm.lbYbar = [1];
  %penm.ubYbar = [1];

  penm.NgLIN=6;
  penm.lbg = ones(6,1);
  penm.ubg = ones(6,1);

  penm.objfun = @ncm2y_objfun;
  penm.objgrad = @ncm2y_objgrad;
  penm.objhess = @ncm2y_objhess;

  penm.confun = @ncm2y_confun;
  penm.congrad = @ncm2y_congrad;
  %penm.conhess = @sdp_conhess;  not needed because all linear
  
end

function [f,userdata] = ncm2y_objfun(x,Y,userdata)
  % f = sum_ij  (y_ij - h_ij)^2
  % matrix H is stored in userdata

  YH = tril(Y{1}-userdata);
  f = YH(:)'*YH(:);

end

function [df, userdata]=ncm2y_objgrad(x,Y,userdata)
% return gradient of the objective w.r.t. all variables (even matrix)
% return NYnnz x 1

  % obj = sum(x_ij - h_ij)^2
  % d/dx_ij obj = 2*(x_ij-h_ij)
  % and we go only over lower triangle --> get indices which belong only to it
  idx = find(tril(ones(size(Y{1}))));
  YH = full(Y{1}-userdata);
  df = 2*YH(idx);
  
end

function [ddf, userdata] = ncm2y_objhess(x,Y,userdata)
% Hessians of the objective function and constraints

  [n m] = size(Y{1});
  dim = n*(n+1)/2;
  ddf = 2*speye(dim,dim);

end

function [g,userdata] = ncm2y_confun(x,Y,userdata)
% only one constraints: tr(Y)=1, return value of the body of the cunstraint
% thus tr(Y)
% vector Ng x 1

  g=diag(Y{1});

end

function [dg, userdata]=ncm2y_congrad(x,Y,userdata)
% returns all constraints at once
% vector Ng x Nx --> 1xNYnnz TODO no ... transpose! (Nx+NYnnz) x Ng
% this needs to be improved to be automatic from size of Y
% Y is organized as lower triangle column oriented

  %dg=[1;0;0;0;0;0; 1;0;0;0;0; 1;0;0;0; 1;0;0; 1;0; 1];

  % only the elements belonging to the diagonal of Y have derivative
  % find their indices

  % this is not the best way but for now ... :-(
  [n m] = size(Y{1});
  dim = n*(n+1)/2;
  M = tril(ones(n)+eye(n));  % my variables are nonzero and these of diag are =2
  idx_allvar = find(M);
  idx_diagvar = find(M(idx_allvar)==2);

  dg = sparse(idx_diagvar, [1:n], ones(n,1), dim,n);
  
end

