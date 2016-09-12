function [Akx, userdata] = bmi_mconfun(x,Y,k,userdata)
% evaluate A_k(x) based on sdpdata, k denotes a block number
% sdpdata is a structure as obtained from readsdpa.m
% note that we aim for A_k(x)<=0 thus A_k(x) is 'reversed' than usual:
%   A_k(x) = F_0 - sum x_i*F_i

  Akx=[];
  if (k<=0 || k>userdata.Na)
    return;
  end

  Akx=userdata.A{k,1};
  for i=userdata.Adep{k}
    % won't work if F{} is [] ... <== dims must match
    Akx = Akx - x(i).*userdata.A{k,i+1};
    % or - ?
  end
  
  indK=userdata.indK{k};
  for i=1:length(x)-1
      li = indK(i,1);
      if li>0
          for jj=1:li
              j = indK(1,jj+1);
              Akx = Akx - x(i)*x(j).*userdata.K{k,i,j};
          end
      end
  end
              
              
      

