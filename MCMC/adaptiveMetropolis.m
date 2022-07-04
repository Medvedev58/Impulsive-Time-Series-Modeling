function [thetas, sl, naccept, yEst] = adaptiveMetropolis(i0, thetas, thetaUL, thetaLL, t_sampled, y, mnoisestd, sl, sigma, usex0)
%ADAPTIVEMETROPOLIS Adaptive metropolis chain

if nargin<10
    usex0=true;
end
nMCMC = size(thetas,2);
ntheta = size(thetas,1);

if usex0
    n_imps = (ntheta-4)/2;
else
    n_imps = (ntheta-2)/2;
end
theta = thetas(:,1)';
L = sl(:,1);

xcov = eye(ntheta)*0.1;
xmean = theta;
Sdpar = 0.02;
epar = Sdpar/100;

naccept = 0;
yEst = [];

for i=2:nMCMC
    if i > i0
        [xcov,xmean]=covRec(thetas(:,i-1)',xcov,xmean,i-1);
        sigma = Sdpar*xcov + Sdpar*epar*eye(numel(xmean));
    end
    
    step = mvnrnd(zeros(ntheta,1), sigma);
    theta_new = theta+step;
    
    if all(theta_new<thetaUL) && all(theta_new>thetaLL) && theta_new(1)<theta_new(2)
        b1_next = theta_new(1);
        b2_next = theta_new(2);
        d_next = theta_new(3:3+n_imps-1);
        t_next = theta_new(3+n_imps:3+2*n_imps-1);
        if usex0
            x0_next = theta_new(end-1:end)';
        else
            x0_next = [0;0];
        end
        
        [t_sorted,perm]=sort(t_next);
        d_sorted = d_next(perm);
        dSeq_next=[t_sorted;d_sorted];
        ytemp=forwardmodel(b1_next,b2_next,t_sampled,dSeq_next,x0_next);
        L_new = -1/(2*mnoisestd^2)*sum((y'-ytemp).^2) - numel(ytemp)/2*log(2*pi) - numel(ytemp)/2*log(mnoisestd^2);
        ys=rand();
        
        a = min(1, exp(L_new - L));
        if ys<a
            theta = theta_new;
            L = L_new;
            naccept = naccept+1;
            yEst = [yEst ytemp'];
        end
    end
    
    sl(:,i) = L;
    thetas(:,i) = theta;
end