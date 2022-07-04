%B1B2ESTIMATIONDEMO demonstrates least squares based and MCMC estimation 
%   of b1 and b2 parameters on synthetic data


%% Data generation parameters
n_datasets=30; % number of synthetic datasets

n_impulses=3; % number of impulses in generated data
min_impulse=0.4; % smallest impulse weight

timescale = 5;
min_time = 1.5; % distance between impulses uniformly distributed between
% min_time and timescale

min_b1=0.4; % lowest b1 value. b1 uniformly distributed between min_b1 and
% min_b1+1
min_bdiff = 0.3; % lowest difference between b1 and b2. b2 uniformly
% distributed between b1+min_bdiff and b1+min_bdiff+1
impulsescale = 4; % Impulse amplitudes uniformly distributed between
% min_impulse and impulsescale
endtime = 5; % fixed time from last impulse to end of horizon
delta_t_sampled=0.5; % sampling time

mnoisestd=0.004; % noise std

%% Estimation parameters

minimplimit=0.05; % minimum limit for counting impulse

plotandpause = true; % pause and generate plots after each dataset

nMCMC = 3e4; % MCMC chain length
i0 = 2;
iBurnin = 1e4; % MCMC burnin length

%%

dist_MCMC_gammaP=zeros(n_datasets,nMCMC-iBurnin); % distance between MCMC sample and LS gammaP estimate
b2s=zeros(n_datasets,nMCMC-iBurnin); % LS estimated gamma_P b2 value
b1s=zeros(n_datasets,nMCMC-iBurnin); % LS estimated gamma_P b1 value
trueb=zeros(n_datasets,2); % true b1 and b2 values

R=cell(1,n_datasets); % Gelman-Rubin statistic for each dataset

t_coarse = zeros(1,n_datasets); % computation time coarse grid
t_fine = zeros(1,n_datasets); % computation time fine grid
t_MCMC = zeros(1,n_datasets); % computation time MCMC

for k=1:n_datasets
    disp(['Run ' num2str(k)])
    
    %% Generate sampled noisy data
    rng((k-1)*5)
    b1= rand(1)+min_b1;
    b2= b1+min_bdiff+rand(1);
    trueb(k,:)=[b1 b2];
    A=[-b1 0; 1 -b2];
    dSeq = [cumsum([1 ones(1,n_impulses-1)*(timescale-min_time)].*rand(1,n_impulses)+min_time);...
        (impulsescale-min_impulse)*rand(1,n_impulses)+min_impulse];
    t_end=dSeq(1,end)+endtime;
    
    t=0:0.02:t_end;
    x = xSol(A, dSeq,t,[0;0]);

    t_sampled=0:delta_t_sampled:t(end);
    x_sampled=interp1(t,x(2,:),t_sampled);
    x_noisy = x_sampled+mnoisestd*randn(size(x_sampled));

    y=x_noisy;
    
    %% Setup optimization parameters
    tk = t_sampled(2:end-1);
    m = length(y);
    n = length(tk);
    options = optimset('lsqlin');
    options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);
    y=y';
    
    b1_min=b1*0.5;
    b1_max=(b1+b2)/2;
    b_delta_coarse=0.02;
    b_delta_fine=0.002;
    b2_min=b1_max+b_delta_coarse;
    b2_max=b2*1.5;
    b1_range_coarse=b1_min:b_delta_coarse:b1_max;
    b2_range_coarse=b2_min:b_delta_coarse:b2_max;
    b2_range_fine=b2_min:b_delta_fine:b2_max;
    num_b1_fine=numel(-3*b_delta_coarse:b_delta_fine:3*b_delta_coarse);
    
    ressum_coarse=zeros(numel(b1_range_coarse),numel(b2_range_coarse));
    beta_coarse=zeros(numel(b1_range_coarse),numel(b2_range_coarse),n+2);
    ressum_fine=zeros(num_b1_fine,numel(b2_range_fine));
    beta_fine=zeros(num_b1_fine,numel(b2_range_fine),n+2);
    
    %% Gridding for LS estimation
    tic
    for k1=1:numel(b1_range_coarse)
        b1_=b1_range_coarse(k1);
        parfor k2=1:numel(b2_range_coarse)
            b2_=b2_range_coarse(k2);
            pars=[b1_ b2_];
            [beta0,ressum_coarse(k1,k2)] = constrImpLS(pars,y,t_sampled,options);
            beta_coarse(k1,k2,:)=beta0;
        end
    end
    t_coarse(k)= toc;
    
    %% Coarse gamma_P estimation
    maximpulses = numel(t_sampled)*0.5;
    sumimpulses = sum(beta_coarse(:,:,2:end)>minimplimit*mean(beta_coarse(:,:,2:end),3),3);
    
    ressumRed = ressum_coarse;
    ressumRed(sumimpulses>maximpulses)=NaN;
    
    ressumgrad = -ressumRed(2:end-1,2:end-1)./((ressumRed(3:end,2:end-1)-ressumRed(1:end-2,2:end-1))/(2*b_delta_coarse));
    
    [m1val,k1range]=min(ressumgrad,[],1);
    k1range=k1range+1;
    b1est=b1_range_coarse(k1range);
    %b1estadj=b1est+m1val;
    
    %% fine gamma_P estimation    
    tic
    b1_range_fine=zeros(num_b1_fine,numel(b2_range_fine));
    for k2=1:numel(b2_range_fine)
        b2_=b2_range_fine(k2);
        [~,coarseIx]=min(abs(b2_range_coarse(2:end-1)-b2_));
        b1_range_fine(:,k2) = b1est(coarseIx)-3*b_delta_coarse:b_delta_fine:b1est(coarseIx)+3*b_delta_coarse;
        parfor k1=1:num_b1_fine
            b1_=b1_range_fine(k1,k2);
            pars=[b1_ b2_];
            [beta0,ressum_fine(k1,k2)] = constrImpLS(pars,y,t_sampled,options);
            beta_fine(k1,k2,:)=beta0;
        end
    end
    t_fine(k)= toc;
    
    sumimpulses = sum(beta_fine(:,:,2:end)>minimplimit*mean(beta_fine(:,:,2:end),3),3);
    
    ressumRed = ressum_fine;
    ressumRed(sumimpulses>maximpulses)=NaN;
    
    ressumgrad = -ressumRed(2:end-1,:)./((ressumRed(3:end,:)-ressumRed(1:end-2,:))/(2*b_delta_fine));
    
    [m1val,k1range]=min(ressumgrad,[],1);
    k1range=k1range+1;
    idx=sub2ind(size(b1_range_fine),k1range,1:size(b2_range_fine,2));
    b1est=b1_range_fine(idx);
    b1estadj=b1est+m1val;
    
    %% MCMC estimation
    
    disp('MCMC')
    n_imps=n_impulses+2; % additional impulses for robustness
    sigma=diag([0.05 0.05 0.15*ones(1,n_imps) 0.5*ones(1,n_imps)].^2);
    
    b1_=rand(1,1)*(b1_max-b1_min)+b1_min;
    b2_=rand(1,1)*(b2_max-b2_min)+b2_min;
    d_=(impulsescale-min_impulse)*rand(1,n_imps)+min_impulse;
    t_=linspace(0,t_end,n_imps+2);
    t_=t_(2:end-1);
    
    dSeq=[t_;d_];
    y0=forwardmodel(b1_,b2_,t_sampled,dSeq);
    L0 = -1/(2*mnoisestd^2)*sum((y'-y0).^2) - numel(y0)/2*log(2*pi) - numel(y0)/2*log(mnoisestd^2);
    sl_ = zeros(1,nMCMC);
    sl_(1) = L0;
    
    thetas_ = zeros(2+n_imps*2,nMCMC);
    thetas_(:,1) = [b1_ b2_ d_ t_]';
    
    thetaUL = inf*ones(1,size(thetas_,1));
    thetaUL(end-n_imps+1:end) = t_end; % max time
    thetaUL(3:3+n_imps) = impulsescale*2; % max impulse
    thetaUL(1:2) =[b1_max b2_max];
    thetaLL = zeros*ones(1,size(thetas_,1));
    thetaLL(1:2) =[b1_min b2_min];
    
    npar=4; % number of parallel chains
    thetas =cell(1,npar);
    sl = cell(1,npar);
    naccept = cell(1,npar);
    yEst = cell(1,npar);
    tic
    parfor n=1:npar
        rng(n*121+k-18)
        [thetas{n}, sl{n}, naccept{n}, yEst{n}] = adaptiveMetropolis(i0, thetas_, thetaUL, thetaLL, t_sampled, y, mnoisestd, sl_, sigma, false);
    end
    t_MCMC(k)= toc;
    
    
    % Impulse post-processing
    parfor n=1:npar
        thetasred{n}=thetas{n}(:,iBurnin+1:end);
        [thetasred{n},redtype{n}]=mergethetas(thetasred{n},t_sampled,0.1,0.1,false);
    end
    
    % Gelman-Rubin statistic
    X_all = cell2mat(thetasred);
    X = reshape(X_all,[size(X_all,1),size(X_all,2)/numel(thetas),numel(thetas)]);
    
    % chain mean
    mx = squeeze(mean(X,2));
    L = size(X,2);
    % grand mean
    mmx = mean(mx,2);
    J = size(X,3);
    % between chain variance
    B = L/(J-1) * sum( (mx - mmx).^2,2);
    
    % within chain variance
    s = squeeze(var(X,[],2));
    W = mean(s,2);
    
    % GR-statistic
    R{k} = (W*(L-1)/L + B/L) ./ W;
    
    disp('GR-statistic')
    disp(num2str(R{k}))
    
    for k1=1:size(X,2)
        [minsquare,minix]=min((X(1,k1,1)-b1estadj).^2 + (X(2,k1,1)-b2_range_fine).^2);
        if minix==1 % outside line
            dist_MCMC_gammaP(k,k1)=NaN;
            b2s(k,k1)=NaN;
            b1s(k,k1)=NaN;
        else
            dist_MCMC_gammaP(k,k1) = sign(X(1,k1,1)-b1estadj(minix))* sqrt(minsquare);
            b2s(k,k1) = b2_range_fine(minix);
            b1s(k,k1) = b1estadj(minix);
        end
    end
    if plotandpause

        figure(1)
        plot(t,x(2,:))
        hold on
        plot(t_sampled,x_noisy,'.')
        hold off
        xlabel('Time')
        ylabel('y')

        for n=1:npar
            figure(1+n)
            histogram2(thetas{n}(1,iBurnin+1:end),thetas{n}(2,iBurnin+1:end),'FaceColor','b')
            zl=zlim;
            hold on
            surf([b1estadj;b1estadj],[b2_range_fine;b2_range_fine],...
                [zeros(size(b1estadj));zl(2)/2*ones(size(b1estadj))],'FaceColor','r')
            alpha 0.5
            zl(1)=zl(1)-2;
            plot3([b1 b1],[b2 b2],zl,'r','LineWidth',2)
            hold off
            xlabel('b_1')
            ylabel('b_2')
        end
        [hc,b1edges,b2edges] = histcounts2(X_all(1,:),X_all(2,:));
        figure(npar+1)
        imagesc(b2edges, b1edges, hc)
        colormap(flipud(gray))
        h = colorbar;
        ylabel(h, 'Count')
        set(gca,'YDir','normal')
        hold on
        plot(b2_range_fine,b1estadj,'b','LineWidth',1)
        plot(b2,b1,'r*','MarkerSize',12)
        hold off
        xlabel('b_2')
        ylabel('b_1')

        pause
    end
end
%% Quantiles for estimated posterior
a1q=zeros(2,n_datasets);
a2q=zeros(2,n_datasets);
d25=zeros(1,n_datasets);
d75=zeros(1,n_datasets);
d5_95=zeros(2,n_datasets);

for k=1:n_datasets
    if nanmean(R{k})<1.2
        a1q(:,k)=quantile(b1s(k,:),[0.25 0.75]);
        a2q(:,k)=quantile(b2s(k,:),[0.25 0.75]);
        d25(k)=quantile(dist_MCMC_gammaP(k,:),0.25);
        d75(k)=quantile(dist_MCMC_gammaP(k,:),0.75);
        d5_95(:,k)=quantile(dist_MCMC_gammaP(k,:),[0.05 0.95]);
    else
        a1q(:,k)=NaN;
        a2q(:,k)=NaN;
        d25(k)=NaN;
        d75(k)=NaN;
        d5_95(:,k)=NaN;
    end
end
