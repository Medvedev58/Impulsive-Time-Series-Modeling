%B1B2INPUTESTIMATIONDEMO demonstrates estimation of both parameters and
%   input impulses from synthetic data

%
% b1=1.1777 b2=1.7152 dSeq = [2.3243    7.2044; 3.9014    2.0324]
%% Data
y = [0.0015 0.0047 0.0015 0.0022 -0.0068 0.5278 0.9905 0.8557 0.5973 0.3880...
    0.2338 0.1367 0.0756 0.0464 0.0342 0.4013 0.5218 0.4222 0.2786 0.1744...
    0.1113 0.0670 0.0296 0.0249 0.0149]; % synthetic data
t_delta=0.5;
t_sampled = 0:t_delta:(numel(y)-1)*t_delta;

%% Parameters
options = optimset('lsqlin');
options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);

% parameter ranges
b1_min=0.7;
b1_max=1.4;

b2_min=1.41;
b2_max=2.4;

% parameter gridding
b1_delta=0.01;
b2_delta=0.01;
b1_range=b1_min:b1_delta:b1_max;
b2_range=b2_min:b2_delta:b2_max;

minimplimit=0.05; % minimum limit for counting impulse

dmin = 0.02; % limit for removing impulses

%% Setup optimization
tk = t_sampled(2:end-1);
m = numel(y);
n = numel(tk);
ressum=zeros(numel(b1_range),numel(b2_range));
betaVals=zeros(numel(b1_range),numel(b2_range),n+2);
    
%% Solve LS over grid of b1 and b2
for k1=1:numel(b1_range)
    b1_=b1_range(k1);
    for k2=1:numel(b2_range)
        b2_=b2_range(k2);
        pars=[b1_ b2_];
        [beta0,ressum(k1,k2)] = constrImpLS(pars,y,t_sampled,options);
        betaVals(k1,k2,:)=beta0;
    end
end

%% Gamma_P estimation
maximpulses = numel(t_sampled)*0.5;
sumimpulses = sum(betaVals(:,:,2:end)>minimplimit*mean(betaVals(:,:,2:end),3),3);

ressumRed = ressum;
ressumRed(sumimpulses>maximpulses)=NaN;

ressumgrad = -ressumRed(2:end-1,:)./((ressumRed(3:end,:)-ressumRed(1:end-2,:))/(2*b1_delta));

[m1val,k1range]=min(ressumgrad,[],1);
k1range=k1range+1;

b1est=b1_range(k1range);
b1estadj=b1est+m1val;

figure(4)
plot(b2_range,b1estadj)
xlabel('b2')
ylabel('b1')

%% Sum impulses along gamma_P

figure(1)
subplot(2,1,1)
plot(t_sampled,y,'b.')

ressumB=zeros(1,numel(b2_range));
dSeqB=cell(1,numel(b2_range));
impt=ones(10,numel(b2_range))*-1;
for k2=1:numel(b2_range)
    b2_=b2_range(k2);
    b1_=b1estadj(k2);
    pars=[b1_ b2_];
    [beta0,ressum0] = constrImpLS(pars,y,t_sampled,options);
    %% Remove small impulses
    useIx = beta0(2:end)>dmin*max(beta0(2:end));
    tk = t_sampled(1:end-1);
    tk=tk(useIx);
    tOpt=tk(1:end);
    [betaSparse,ressumB(k2)]=constrImpLS(pars,y,t_sampled,options,tOpt);
    dSeq = [[0 tOpt];betaSparse(2:end)'];
    dSeqB{k2} = mergeImpulses(dSeq,b1_,b2_,t_sampled);
    
    impt(1:numel(dSeqB{k2}(1,:)),k2)=dSeqB{k2}(1,:);
    
    % Plot for some example values
    if b2_==1.5
        if size(dSeqB{k2},2)==1
            dSeqB{k2}=[dSeqB{k2} [0;0]];
        end
        ASim = [-b1_ 0; 1 -b2_];
        t=t_sampled(1):0.1:t_sampled(end);
        xSim = xSol(ASim, dSeqB{k2}(:,2:end),t,[betaSparse(2);betaSparse(1)]);
        subplot(2,1,1)
        hold on
        plot(t,xSim(2,:),'r')
        hold off
        subplot(2,1,2)
        stem(dSeqB{k2}(1,:),dSeqB{k2}(2,:),'r')
        xlabel('Time')
        ylabel('Impulse weight')
    end
    if b2_==1.81
        if size(dSeqB{k2},2)==1
            dSeqB{k2}=[dSeqB{k2} [0;0]];
        end
        ASim = [-b1_ 0; 1 -b2_];
        t=t_sampled(1):0.1:t_sampled(end);
        xSim = xSol(ASim, dSeqB{k2}(:,2:end),t,[betaSparse(2);betaSparse(1)]);
        subplot(2,1,1)
        hold on
        plot(t,xSim(2,:),'g')
        hold off
        subplot(2,1,2)
        hold on
        stem(dSeqB{k2}(1,:),dSeqB{k2}(2,:),'g')
        hold off
    end
    if b2_==2.1
        if size(dSeqB{k2},2)==1
            dSeqB{k2}=[dSeqB{k2} [0;0]];
        end
        ASim = [-b1_ 0; 1 -b2_];
        t=t_sampled(1):0.1:t_sampled(end);
        xSim = xSol(ASim, dSeqB{k2},t,[0;betaSparse(1)]);
        subplot(2,1,1)
        hold on
        plot(t,xSim(2,:),'m')
        xl=xlim();
        hold off
        subplot(2,1,2)
        hold on
        stem(dSeqB{k2}(1,:),dSeqB{k2}(2,:),'m')
        hold off
        xlim(xl);
    end
end
impt(impt==-1)=NaN;
