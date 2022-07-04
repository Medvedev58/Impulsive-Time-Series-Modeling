function [thetas,redtype] = mergethetas(thetas,t_sampled,dlim,trelax,usex0)
%MERGETHETAS merge or remove impulses based on small separation, low
%   weights or being outside the simulation horizon
nMCMC = size(thetas,2);
ntheta = size(thetas,1);

if usex0
    n_imps = (ntheta-4)/2;
else
    n_imps = (ntheta-2)/2;
end

d = thetas(3:3+n_imps-1,:);
t = thetas(3+n_imps:3+2*n_imps-1,:);

[t_,perm]=sort(t);

d_ = zeros(size(d));
redtype=zeros(1,4);

for k=1:nMCMC
    d_(:,k) = d(perm(:,k),k);
    for k1=1:n_imps
        if k1<n_imps
            condT = t_sampled>t_(k1,k);
            condT2 = t_sampled>t_(k1+1,k);
            condT2rel = t_sampled>t_(k1+1,k)-trelax;
            if (any(condT2) && find(condT,1)==find(condT2,1))...
                    ||((any(condT) && any(condT2rel)) && (find(condT,1)==find(condT2rel,1)))
                b1=thetas(1,k);
                b2=thetas(2,k);
                delta=t_(k1+1,k)-t_(k1,k);
                tauk = 1/(b2-b1)* log(d_(k1+1,k)*(exp(b2*delta)-exp(b1*delta))/(d_(k1,k) + exp(b1*delta)*d_(k1+1,k)) + 1) + t_(k1,k);
                gk = exp(-b1*(tauk-t_(k1,k)))*(exp(b1*delta)*d_(k1+1,k)+d_(k1,k));
                d_(:,k)=[d_(1:k1-1,k); gk; d_(k1+2:end,k); NaN];
                t_(:,k)=[t_(1:k1-1,k); tauk;  t_(k1+2:end,k); NaN];
                redtype(1)=redtype(1)+1;
                if k1+1<n_imps
                    kf1 = find(condT,1);
                    condT3 = t_sampled>t_(k1+1,k);
                    condT3rel = t_sampled>t_(k1+1,k)-trelax;
                    if any(condT3) && (find(condT3,1)==kf1 || find(condT3rel,1)==kf1)
                        delta=t_(k1+1,k)-tauk;
                        tauk = 1/(b2-b1)* log(d_(k1+1,k)*(exp(b2*delta)-exp(b1*delta))/(d_(k1,k) + exp(b1*delta)*d_(k1+1,k)) + 1) + t_(k1,k);
                        gk = exp(-b1*(tauk-t_(k1,k)))*(exp(b1*delta)*d_(k1+1,k)+d_(k1,k));
                        d_(:,k)=[d_(1:k1-1,k); gk; d_(k1+2:end,k); NaN];
                        t_(:,k)=[t_(1:k1-1,k); tauk;  t_(k1+2:end,k); NaN];
                        redtype(4)=redtype(4)+1;
                    end
                end
            end
        end
        if d_(k1,k)<dlim*mean(d(:,k))
            d_(:,k)=[d_(1:k1-1,k); NaN; d_(k1+1:end,k)];
            t_(:,k)=[t_(1:k1-1,k); NaN; t_(k1+1:end,k)];
            redtype(2)=redtype(2)+1;
        end
        if t_(k1,k)+trelax>t_sampled(end)
            d_(:,k)=[d_(1:k1-1,k); NaN; d_(k1+1:end,k)];
            t_(:,k)=[t_(1:k1-1,k); NaN; t_(k1+1:end,k)];
            redtype(3)=redtype(3)+1;
        end
    end
    d_(:,k)=[d_(~isnan(d_(:,k)),k);d_(isnan(d_(:,k)),k)];
    t_(:,k)=[t_(~isnan(t_(:,k)),k);t_(isnan(t_(:,k)),k)];
end

thetas = [thetas(1:2,:); d_; t_; thetas(3+2*n_imps:end,:)];