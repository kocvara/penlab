function [penm] = VTSdual_define(par)
% define penm structure for the dual VTS problem with lower and upper bounds

nx = par.nx; ny = par.ny;
A = par.A; IA = par.IA; nelem=par.nelem; nnod=par.nnod;
nloads = par.nloads;
ff = par.ff{1};
ny2 = ny/2; nelem2 = nelem/2; nnod2 = 2*nx*ny2; nnodi = 2*nx;
n = nnod ; m = nelem;
volume = nelem/3;
lbound = 1e-4;
ubound = 1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = 1;
for ie=1:nelem
    len = IA(ie);
    %clear ss;
    ss = sparse(A(lb:lb+len-1),A(lb+len:lb+2*len-1),A(lb+2*len:lb+3*len-1),nnod,nnod);
    hel = triu(ss) + tril(ss)' - diag(diag(ss));
    ss = hel +hel' - diag(diag(hel)) ;
    AA{1,ie} = ss;
    lb = lb + 3*len;
end

sdpdata.nelem = m;
sdpdata.nnod = n;
sdpdata.rhs = ff;
sdpdata.lbound = lbound;
sdpdata.ubound = ubound;
sdpdata.A = AA;
sdpdata.V = volume;
clear A;

penm = [];

penm.userdata=sdpdata;

penm.Nx=n+1+m+m;
penm.lbx=[-Inf.*ones(n+1,1);zeros(m+m,1)];

penm.NgNLN = m;
penm.ubg = zeros(m,1);

penm.objfun = @VTSd_lu_objfun;
penm.objgrad = @VTSd_lu_objgrad;
penm.objhess = @VTSd_lu_objhess;

penm.confun = @VTSd_lu_confun;
penm.congrad = @VTSd_lu_congrad;
penm.conhess = @VTSd_lu_conhess;

end

%%%%% functions
function [g,userdata] = VTSd_lu_confun(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  lbound = userdata.lbound;
  ubound = userdata.ubound;
  A = userdata.A;
  
  u = x(1:nnod);
  alpha = x(nnod+1);
  rho_l = x(nnod+2:nnod+2+nelem);
  rho_u = x(nnod+2+nelem:end);
  
  for i = 1:nelem
      g(i,1) = 0.5.*u'*A{1,i}*u - alpha + rho_l(i) - rho_u(i);
  end
end

function [dg, userdata]=VTSd_lu_congrad(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  lbound = userdata.lbound;
  ubound = userdata.ubound;
  A = userdata.A;
  
  u = x(1:nnod);
  alpha = x(nnod+1);
  rho_l = x(nnod+2:nnod+2+nelem);
  rho_u = x(nnod+2+nelem:end);
  
  for i = 1:nelem
      ei = sparse(nelem,1); ei(i) = 1;
      dg(:,i) = [A{1,i}*u; -1; ei; -ei];
  end
end

function [ddf, userdata]=VTSd_lu_conhess(x,Y,k,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  A = userdata.A;
  
  ddf = sparse(nelem+nelem+1+nnod,nelem+nelem+1+nnod);
  ddf(1:nnod,1:nnod) = A{1,k};
end

function [f,userdata] = VTSd_lu_objfun(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  lbound = userdata.lbound;
  ubound = userdata.ubound;
  V = userdata.V;
  
  alpha = x(nnod+1);
  
  f = V*alpha - rhs'*x(1:nnod) - lbound*sum(x(nnod+2:nnod+2+nelem)) + ubound*sum(x(nnod+2+nelem:end));
end

function [df, userdata]=VTSd_lu_objgrad(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  lbound = userdata.lbound;
  ubound = userdata.ubound;
  V = userdata.V;
  
  alpha = x(nnod+1);
  
  df = [-rhs; V; -lbound.*ones(nelem,1); ubound.*ones(nelem,1)];
end

function [ddf, userdata] = VTSd_lu_objhess(x,Y,userdata)
  ddf = [];
end