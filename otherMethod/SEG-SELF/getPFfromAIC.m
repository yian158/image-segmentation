function [PF] = getPFfromAIC(stats,k)

PF = zeros(length(stats),k);

for i=1:length(stats),
    AIC = stats(i).AIC;
    
    if length(AIC) < k
        AIC(length(AIC):k) = min(AIC);
    end
    
    %AIC = nCompl*log(1-TotalPerf)+2*NUMEllipses;
    for j=1:k,
        PF(i,j) = 1-exp((AIC(j)-2*j) / stats(i).NumCom);
    end
end

