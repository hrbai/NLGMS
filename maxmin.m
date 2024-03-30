function [nmi_max,acc_max] = maxmin(X,class_num,Y,repeat_times)
    nmis = zeros(repeat_times,1);
    accs = zeros(repeat_times,1);
    for i=1:repeat_times
      label = litekmeans(X,class_num,'Replicates',1);
      % NMI and acc
      nmis(i) = nmi(Y,label);
      bestLabel = bestMap(Y,label);
      accs(i) = sum(bestLabel==Y)/length(Y);
    end
    % max 
    nmi_max = max(nmis);
    acc_max = max(accs);
end