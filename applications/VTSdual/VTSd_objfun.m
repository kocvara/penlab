function [f,userdata] = VTSd_objfun(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  ubound = userdata.ubound;
  V = userdata.V;
  
  alpha = x(nnod+1);
  
  f = V*alpha - rhs'*x(1:nnod) + ubound*sum(x(nnod+2:end));

