function [nmi_mean,nmi_std,acc_mean,acc_std] = write_baseline_result(X,class_num,Y,fpath)
    fs = fopen(fpath, 'a+');
    fprintf(fs, ['\r\n***** tuning parameters for ',fpath,'*****\r\n']);
    % cluster
    [nmi_mean,nmi_std,acc_mean,acc_std] = repeatkmeans(X,class_num,Y,20);
    fprintf(fs,['Clustering using all the ',num2str(size(X,2)),' features. Clustering MIhat: ',num2str(nmi_mean),...
        ', std: ',num2str(nmi_std),'\r\n']);
    fprintf(fs,['Selected feature num: ',num2str(size(X,2)),', Clustering ACC: ',num2str(acc_mean),...
        ', std: ',num2str(acc_std),'\r\n']);
   
end