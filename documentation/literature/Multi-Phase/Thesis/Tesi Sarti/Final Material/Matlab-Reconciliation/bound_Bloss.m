function [c,ceq]=bound_Bloss(x)
c=[];
rad_s=x(12)*2*pi/60;
W=rad_s*(x(1)+0.33);
ceq=W-x(2)*(x(10)-x(11))+x(8);