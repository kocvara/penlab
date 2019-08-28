function [penm] = ex72_define()
% Example 7.2, problem (17) from PENNON user's guide
%

  penm = [];

  penm.probname= 'Example 7.2 from PENNON User''s Guide';
  penm.comment = 'a scalar and a 6x6 full matrix as variables';

  % matrix H
  H = [1 -0.44 -.2 .81 -.46 -0.05;
    -.44 1 .87 -.38 .81 -.58;
    -.2 .87 1 -.17 .65 -.56;
    .81 -.38 -.17 1 -.37 -.15;
    -.46 .81 .65 -.37 1 .08;
    -.05 -.58 -.56 -.15 .08 1];

  % keep the whole structure
  penm.userdata.H = H;
  [dum,idiag] = packmat(H); clear dum
  penm.userdata.idiag = idiag;

  % one 'normal' variable
  penm.Nx=1;
  % one matrix variable
  penm.NY=1;
  penm.Y{1}=H; %to define sparsity structure of Y

  % box constraints on matrix variables
  penm.lbY = [1];
  penm.ubY = [10];
  
  % nonlinear constraints
  penm.NgNLN=6;
  penm.lbg = [1;1;1;1;1;1];
  penm.ubg = [1;1;1;1;1;1];

  penm.objfun = @ex72_objfun;
  penm.objgrad = @ex72_objgrad;
  penm.objhess = @ex72_objhess;

  penm.confun = @ex72_confun;
  penm.congrad = @ex72_congrad;
  penm.conhess = @ex72_conhess; 
  
  penm.Yinit{1}=eye(6);
  penm.xinit = 1;

end

function [f,userdata] = ex72_objfun(x,Y,userdata)
% Example 7.2, problem (17) from PENNON user's guide
  % matrix H is stored in userdata

  YH = svec2(x(1).*Y{1}-userdata.H);
  f = YH(:)'*YH(:);
  
end

function [df, userdata]=ex72_objgrad(x,Y,userdata)
% Example 7.2, problem (17) from PENNON user's guide
% return gradient of the objective w.r.t. all variables

  % 
  YH=svec2(x(1).*Y{1}-userdata.H);
  
  df(1) = sum(2*svec2(Y{1}).*YH);
  df(2:length(YH)+1) = 2*x(1).*YH;
  
  df = df';
  
end

function [ddf, userdata] = ex72_objhess(x,Y,userdata)
% Example 7.2, problem (17) from PENNON user's guide
% Hessians of the objective function and constraints
  
  YH=packmat(x(1).*Y{1}-userdata.H);
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

function [g,userdata] = ex72_confun(x,Y,userdata)
% Example 7.2, problem (17) from PENNON user's guide

  % 
  for i=1:length(Y{1})
      g(i,1)=x(1)*Y{1}(i,i);
  end

end

function [dg, userdata]=ex72_congrad(x,Y,userdata)
% Example 7.2, problem (17) from PENNON user's guide
% returns all constraints at once
% rectangule matrix (Nx+NYnnz) x Ng

  %
  n =length(Y{1});
  nn = n*(n+1)/2;
  idiag = userdata.idiag;
  
  dg = sparse((nn+1),n); % #variables x #constraints
  dg(1,:) = diag(Y{1});
  
  for i=1:n
      dg(idiag(i)+1,i) = x(1);
  end
  
end

function [ddgk, userdata]=ex72_conhess(x,Y,k,userdata)
% Example 7.2, problem (17) from PENNON user's guide

  %
  n =length(Y{1});
  nn = n*(n+1)/2;
  idiag = userdata.idiag;
  
  ddgk = zeros((nn+1),(nn+1));
  ddgk(1,1) = 1;
  ddgk(idiag(k)+1,idiag(k)+1) = 1;
  
end