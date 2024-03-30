function [nmi_mean,nmi_std,acc_mean,acc_std] = repeatkmeans(X,class_num,Y,repeat_times)
    nmis = zeros(repeat_times,1);
    accs = zeros(repeat_times,1);
    for i=1:repeat_times
      label = litekmeans(X,class_num,'Replicates',1);
      % NMI and acc
      nmis(i) = nmi(Y,label);
      bestLabel = bestMap(Y,label);
      accs(i) = sum(bestLabel==Y)/length(Y);
    end
    % mean and std
    nmi_mean = mean(nmis);
    nmi_std = std(nmis);
    acc_mean = mean(accs);
    acc_std = std(accs);
end