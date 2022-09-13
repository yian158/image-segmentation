function [ GT ] = correctGT( GT)

% vals = GT(:);
% id1 = max(vals);
% Len = zeros(1,id1);
% for i=1:length(vals),
%     x = vals(i);
%     if x > 0
%         Len(x) = Len(x)+1;
%     end
% end
% S = sum(Len);

% for i=1:length(Len),
%     if Len(i) < 0.05*S,
%         GT(GT == i) = 0;
%     end
% end
% GT(GT ~= 0) = 1;

GT(GT >= 1) = 1;
GT(GT < 1) = 0;


end

