function [g,userdata] = VTSd_confun(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  ubound = userdata.ubound;
  A = userdata.A;
  
  u = x(1:nnod);
  alpha = x(nnod+1);
  rho = x(nnod+2:end);
  
  for i = 1:nelem
      g(i,1) = 0.5.*u'*A{1,i}*u - alpha - rho(i);
  end



