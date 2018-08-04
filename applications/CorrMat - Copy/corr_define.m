function [penm] = corr_define(H)
% CORR_DEFINE defines PenLab structure for the
% nearest correlation matrix problem.
% Both input arguments are optional to redefine the default values.
%
% Call:
% >> penm = corr_define;
% >> problem = penlab(penm);
% >> problem.solve();
% >> eig(problem.Y{1}*problem.x)

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 19 August 2016

  penm = [];

  penm.probname = 'NCM, Example 7.1 from PENLAB paper, no constraints';
  penm.comment  = 'a scalar and a full matrix as variables';

  % matrix H
  if (nargin<1)
    H = [1 -0.44 -.2 .81 -.46 -0.05;
      -.44 1 .87 -.38 .81 -.58;
      -.2 .87 1 -.17 .65 -.56;
      .81 -.38 -.17 1 -.37 -.15;
      -.46 .81 .65 -.37 1 .08;
      -.05 -.58 -.56 -.15 .08 1];
  end

  n = size(H,1);

  % keep the whole structure
  penm.userdata.H = H;
  [dum,idiag] = svec2(H); clear dum
  penm.userdata.idiag = idiag;

  % one matrix variable
  penm.NY = 1;
  penm.Y{1} = ones(n,n); % to define sparsity structure of Y

  % box constraints on matrix variables
  penm.lbY = [0];
  
  % box constraints on matrix elements
  mlb = -Inf(n,n); mub = Inf(n,n); 
  % diag elements equal to one
  for i=1:n; mlb(i,i) = 1; mub(i,i) = 1; end;
  
  penm.lbYx=cell(1,1); 
  penm.ubYx=cell(1,1); 
  penm.lbYx{1} = mlb;
  penm.ubYx{1} = mub;  

  penm.objfun = @corr_objfun;
  penm.objgrad = @corr_objgrad;
  penm.objhess = @corr_objhess;

  % starting point
  penm.Yinit{1} = eye(n,n);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [f,userdata] = corr_objfun(x,Y,userdata)

  YH = packmat(Y{1}-userdata.H);
  f = YH(:)'*YH(:);
  end

  
 function [df, userdata] = corr_objgrad(x,Y,userdata)

  YH=svec2(Y{1}-userdata.H);
  
  df(1) = sum(2*svec2(Y{1}).*YH);
  df(1:length(YH)) = 2.*YH;
  
  df = df';
 end

  
function [ddf, userdata] = corr_objhess(x,Y,userdata)
  
  YH = packmat(Y{1}-userdata.H);
  yy = packmat(Y{1});
  n = length(yy);
  ddf = zeros(n,n);
  
  for i= 1:n
      ddf(i,i) = 2;
  end
end

