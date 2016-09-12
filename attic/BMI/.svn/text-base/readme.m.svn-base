% *Bilinear matrix inequalities*
% 
% We want to solve an optimization problem with constraints in the form of
% bilinear matrix inequalities. The problem is store in a file
% |bmidata.mat| (obviously, the name of the structure is optional). The
% structure of the datafile is explained in Appendix BMI of PENlab Manual. 
% For instance, you can create a sample structure by calling
%
bmi_ex3
%
% or by loading the data from a file, e.g.,
%
load bmidata
%
% Once the structure is in the memory, all the user has to
% do to solve the problem is to call the following sequence of commands:
%
penm=bmi_define(bmidata);
prob=penlab(penm);
prob.solve();