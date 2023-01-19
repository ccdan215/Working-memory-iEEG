% Subtract trial-wise mean (per condition).
%
% E. L. Johnson, PhD
% Copyright (c) 2017
% UC Berkeley
% eljohnson@berkeley.edu

function [cond1, cond2, cond3] = subtracterp(data, ncond)

cfg = [];

cfga = [];
cfga.avgoverrpt = 'yes';

if ncond == 1
    cfg.trials = find(data.trialinfo(:,2)==1);
    
    cond1 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond1);
    for r = 1:length(cond1.trial)
        cond1.trial{1,r} = cond1.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
    
elseif ncond == 2
    cfg.trials = find(data.trialinfo(:,2)==1 & data.trialinfo(:,3)==1);
    
    cond1 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond1);    
    for r = 1:length(cond1.trial)
        cond1.trial{1,r} = cond1.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
    
    cfg.trials = find(data.trialinfo(:,2)==1 & data.trialinfo(:,3)==2);
    
    cond2 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond2);    
    for r = 1:length(cond2.trial)
        cond2.trial{1,r} = cond2.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
    
    
elseif ncond == 3
    cfg.trials = find(data.trialinfo(:,2)==1 & data.trialinfo(:,4)==1);
    
    cond1 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond1);    
    for r = 1:length(cond1.trial)
        cond1.trial{1,r} = cond1.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
    
    cfg.trials = find(data.trialinfo(:,2)==1 & data.trialinfo(:,4)==2);
    
    cond2 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond2);    
    for r = 1:length(cond2.trial)
        cond2.trial{1,r} = cond2.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
    
    cfg.trials = find(data.trialinfo(:,2)==1 & data.trialinfo(:,5)==2);
    
    cond3 = ft_selectdata(cfg, data);
    cond_mean = ft_selectdata(cfga, cond3);    
    for r = 1:length(cond3.trial)
        cond3.trial{1,r} = cond3.trial{1,r}-cond_mean.trial{1,1};
    end
    clear cond_mean
end

clearvars -except cond*

end
