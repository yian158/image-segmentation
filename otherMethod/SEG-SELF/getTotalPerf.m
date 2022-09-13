%Gets the total performance of EL set and the AIC and BIC measures

function [TotalPerf,AIC,BIC] = getTotalPerf(EL,Crit)

if Crit == 1,
    perf = [EL.tomh_enwsh];
    [min_perf, pos] = min(perf);
elseif Crit == 2,
    perf = [EL.outPixels];
    [min_perf, pos] = max(perf);
end
TotalPerf = sum([EL.InArea]) / area;

