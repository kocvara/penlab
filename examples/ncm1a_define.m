function [penm] = ncm1a_define()
% define penm structure for Nearest Correlation Matrix Example 1
% using matrix inequalities A(x):
%    min_x sum (x_ij-h_ij)^2
%    s.t.  tr(X)=1
%          X >= 0 (pos. semidefinite)
%
  penm = [];

  penm.probname= 'NCM Example 1 using normal variables & matrix constraints';
  penm.comment = '3 diagonal matrix 3x3 as a vector of variables';

  % matrix H
  H = [2.2, -1.1, 0; -1.1, 1.9, -1.15; 0, -1.15, 2.1]./6;
  %H = [2.2, -1.1, 0; -1.1, 1.9, -1.1; 0, -1.1, 2.1]./6;
  %H = [2, -1, 0; -1, 2, -1; 0, -1, 2]./6;
  %H = [2.2, -1.1, 0; -1.1, 1.9, -1.1; 0, -1.1, 2.1];

  % matrix X will be expressed as:
  %   X = [ x_1, x_2,   0;
  %         x_2, x_3, x_4;
  %           0, x_4, x_5]  
  %
  %   XX = [ x(1), x(2), 0; x(2), x(3), x(4); 0, x(4), x(5)]  
  % by using matrices A{i} which just point where x(i) will be projected
  %   X = A{1}*x(1) + A{2}*x(2) + ... + A{5}*x(5)
  A = cell(5,1);
  A{1} = [1,0,0;
          0,0,0;
          0,0,0];
  A{2} = [0,1,0;
          1,0,0;
          0,0,0];
  A{3} = [0,0,0;
          0,1,0;
          0,0,0];
  A{4} = [0,0,0;
          0,0,1;
          0,1,0];
  A{5} = [0,0,0;
          0,0,0;
          0,0,1];

  userdata.H=H;
  userdata.A=A;

  % keep the whole structure
  penm.userdata=userdata;

  penm.Nx = 5;

  % starting point
  penm.xinit=[1;0;1;0;1];

  penm.NgLIN = 1;
  penm.lbg = [1];
  penm.ubg = [1];

  penm.NALIN=1;
  % and it should be positive semidefinite
  penm.lbA=[0];

  penm.objfun = @ncm1a_objfun;
  penm.objgrad = @ncm1a_objgrad;
  penm.objhess = @ncm1a_objhess;

  penm.confun = @(x,Y,userdata) deal(x(1)+x(3)+x(5), userdata);
  penm.congrad = @(x,Y,userdata) deal([1;0;1;0;1], userdata);

  penm.mconfun = @ncm1a_mconfun;
  penm.mcongrad = @(x,Y,k,i,userdata) deal(sparse(userdata.A{i}), userdata);
  
end

function [f,userdata] = ncm1a_objfun(x,Y,userdata)
  % f = sum_ij  (x_ij - h_ij)^2
  % matrix H is stored in userdata

  % I could build Y and do this, but let's write it directly
  %YH = Y-userdata;
  %f = YH(:)'*YH(:);

  H=userdata.H;
  f=(x(1)-H(1,1))^2 + 2*(x(2)-H(2,1))^2 + (x(3)-H(2,2))^2 + 2*(x(4)-H(3,2))^2 + (x(5)-H(3,3))^2;

end
  
function [df, userdata]=ncm1a_objgrad(x,Y,userdata)
% return gradient of the objective w.r.t. all variables (even matrix)
% return Nx x 1

  H=userdata.H;
  df=[ 2*(x(1)-H(1,1)); 
       4*(x(2)-H(2,1)); 
       2*(x(3)-H(2,2)); 
       4*(x(4)-H(3,2));
       2*(x(5)-H(3,3)) ];
   
end

function [ddf, userdata] = ncm1a_objhess(x,Y,userdata)
% Hessians of the objective function and constraints

  ddf = diag([2,4,2,4,2]);  
  
end
  
function [Akx, userdata] = ncm1a_mconfun(x,Y,k,userdata)
% There is only one matrix variable A(x) = sum x(i)*userdata.A{i}

  A=userdata.A;
  Akx=x(1).*A{1} + x(2).*A{2} + x(3).*A{3} + x(4).*A{4} + x(5).*A{5};

end


