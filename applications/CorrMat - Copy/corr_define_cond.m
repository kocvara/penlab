function [penm] = corr_define_cond(kappa,H)
% CORR_DEFINE defines PenLab structure for Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.
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

  penm.probname = 'NCM, Example 7.1 from PENLAB paper';
  penm.comment  = 'a scalar and a full matrix as variables';

  % matrix H
  if (nargin<2)
    H = [1 -0.44 -.2 .81 -.46 -0.05;
      -.44 1 .87 -.38 .81 -.58;
      -.2 .87 1 -.17 .65 -.56;
      .81 -.38 -.17 1 -.37 -.15;
      -.46 .81 .65 -.37 1 .08;
      -.05 -.58 -.56 -.15 .08 1];
  end

  if (nargin<1)
    kappa = 10;
  end

  n = size(H,1);

  % keep the whole structure
  penm.userdata.H = H;
  [dum,idiag] = svec2(H); clear dum
  penm.userdata.idiag = idiag;

  % one 'normal' variable
  penm.Nx = 1;
  % one matrix variable
  penm.NY = 1;
  penm.Y{1} = ones(n,n); % to define sparsity structure of Y

  % box constraints on matrix variables
  penm.lbY = [1];
  penm.ubY = [kappa];
  
  % nonlinear constraints - diagonal 
  penm.NgNLN = n;
  penm.lbg = [ones(n,1)] ;
  penm.ubg = [ones(n,1)] ;

  penm.objfun = @corr_objfun;
  penm.objgrad = @corr_objgrad;
  penm.objhess = @corr_objhess;

  penm.confun = @corr_confun;
  penm.congrad = @corr_congrad;
  penm.conhess = @corr_conhess; 
  
  % starting point
  penm.Yinit{1} = eye(n,n);
  penm.xinit = 1;
  
  end
  
  function [g,userdata] = corr_confun(x,Y,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

  n=length(Y{1});
  for i=1:n
      g(i,1)=x(1)*Y{1}(i,i);
  end
  end
  
  function [dg, userdata] = corr_congrad(x,Y,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.
% It returns all gradients of (standard) constraints at once
% as a rectangular matrix (Nx+NYnnz) x Ng.

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

  n = length(Y{1});
  nn = n*(n+1)/2;
  idiag = userdata.idiag;
  
  dg = sparse((nn+1),n); % #variables x #constraints
  dg(1,:) = [diag(Y{1})];
  
  for i=1:n
      dg(idiag(i)+1,i) = x(1);
  end
  end
  
  
  function [ddgk, userdata] = corr_conhess(x,Y,k,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

  n = length(Y{1});
  nn = n*(n+1)/2;
  idiag = userdata.idiag;
  ddgk = zeros((nn+1),(nn+1));
  ddgk(1,1) = 1;
  ddgk(idiag(k)+1,idiag(k)+1) = 1;
  end
  
  function [f,userdata] = corr_objfun(x,Y,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

  YH = packmat(x(1).*Y{1}-userdata.H);
  f = YH(:)'*YH(:);
  end

  
 function [df, userdata] = corr_objgrad(x,Y,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.
% return gradient of the objective w.r.t. all variables

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

  YH=svec2(x(1).*Y{1}-userdata.H);
  
  df(1) = sum(2*svec2(Y{1}).*YH);
  df(2:length(YH)+1) = 2*x(1).*YH;
  
  df = df';
 end

  
function [ddf, userdata] = corr_objhess(x,Y,userdata)
% Example 7.1 from the PENLAB paper,
% nearest correlation matrix with the constrained condition number.
% Hessians of the objective function.

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013
  
  YH = packmat(x(1).*Y{1}-userdata.H);
  yy = packmat(Y{1});
  n = length(yy);
  ddf = zeros(n+1,n+1);
  
  ddf(1,1) = 2*sum(yy.^2);
  ddf(1,2:n+1) = 2.*(x(1).*yy+YH);
  ddf(2:n+1,1) = 2.*(x(1).*yy'+YH');
  for i= 1:n
      ddf(i+1,i+1) = 2*x(1)^2;
  end
end

