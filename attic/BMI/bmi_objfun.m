function [f,userdata] = bmi_objfun(x,Y,userdata)

  f = userdata.c'*x;

