function y = forwardmodel(b1,b2,t_sampled,dSeq,x0)
%FORWARDMODEL simulate impulsive system with given parameters
if nargin<5
    x0=[0;0];
end
A=[-b1 0; 1 -b2];
x = xSol(A, dSeq,t_sampled,x0);
y=x(2,:);
end