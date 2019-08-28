function [dg, userdata]=VTSd_congrad(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  ubound = userdata.ubound;
  A = userdata.A;
  
  u = x(1:nnod);
  alpha = x(nnod+1);
  rho = x(nnod+2:end);
  
  for i = 1:nelem
      ei = sparse(nelem,1); ei(i) = 1;
      dg(:,i) = [A{1,i}*u; -1; -ei];
  end
