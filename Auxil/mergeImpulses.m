function dSeq_m = mergeImpulses(dSeq,b1,b2,t_sampled)
% MERGEIMPULSES iterates through impulse vector dSeq and replaces
% consecutive impulses with a single impulse in between
dSeq_m=[];
startedPair=false;
if size(dSeq,2)==1
    dSeq_m=dSeq;
    return;
end
for k=2:size(dSeq,2)
    if find(t_sampled==dSeq(1,k))-find(t_sampled==dSeq(1,k-1))==1&& ~startedPair
        delta=dSeq(1,k)-dSeq(1,k-1);
        tauk = 1/(b2-b1)* log(dSeq(2,k)*(exp(b2*delta)-exp(b1*delta))...
            /(dSeq(2,k-1) + exp(b1*delta)*dSeq(2,k)) + 1) + dSeq(1,k-1);
        gk = exp(-b1*(tauk-dSeq(1,k-1)))*(exp(b1*delta)*dSeq(2,k)+dSeq(2,k-1));
        dSeq_m=[dSeq_m [tauk;gk]];
        startedPair=true;
    elseif startedPair
        startedPair=false;
    else
        dSeq_m=[dSeq_m dSeq(:,k-1)];
    end
end
if ~startedPair
    dSeq_m=[dSeq_m dSeq(:,k)];
end