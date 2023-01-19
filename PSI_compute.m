% Compute phase slope index (psi), null distribution of psi, and psi-zscore. 
% Requires fieldtrip toolbox.
% subfunctions: subtracterp, from the function shared by Johnson EL, et al.,2018.
% Key:
% subid = 's1' (required)
% elecs = {'Hipp1','Hipp2','Amy1','Amy2'} (required)
% electrodes from the hippocampus and the amygdala within the same hemisphere

function [] = PSI_compute(subid, elecs)
% PSI_compute('s1', {'Hipp1','Hipp2','Amy1','Amy2'});

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%

dataDir = [path '/data_indiv/'];
saveDir = [path '/psi/'];
mkdir(saveDir);

% load data
load([dataDir subid '/data_derived']);

% initialize fieldtrip
ft_defaults;

% raw data were first converted into the format used in EEGlab toolbox, 
% and then it was converted to the format used in filedtrip toolbox
data=eeglab2fieldtrip(EEG,'preprocessing');

% select elecs
cfgc = [];
cfgc.channel = elecs;
data_select = ft_selectdata(cfgc, data);

% compute encoding PSI
cfgt=[];
cfgt.latency = [1 3];
data_select = subtracterp(data_select, 1); % subtract ERP
data_enc = ft_selectdata(cfgt, data_select);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 1:1:40; % 1-40 Hz
cfg.tapsmofrq = 1; % 1-Hz half-band
cfg.pad = 4; % pad to 4 s
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
data_enc_cs = ft_freqanalysis(cfg, data_enc); % spectral decomposition

cfgp = [];
cfgp.method = 'psi';
cfgp.bandwidth = 5;
psi_enc = ft_connectivityanalysis(cfgp, data_enc_cs);

% compute maintenance PSI
cfgt=[];
cfgt.latency = [3 6];
data_maint = ft_selectdata(cfgt, data_select);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 1:1:40; % 1-40 Hz
cfg.tapsmofrq = 1; % 1-Hz half-band
cfg.pad = 6; % pad to 6 s
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
data_maint_cs = ft_freqanalysis(cfg, data_maint); % spectral decomposition

cfgp = [];
cfgp.method = 'psi';
cfgp.bandwidth = 5;
psi_maint = ft_connectivityanalysis(cfgp, data_maint_cs);

% save by subject
save([saveDir subid '_psi'], 'psi_enc', 'psi_maint');

%%% Compute null distribution of PSI using statistical bootstrapping.
% By shuffling the order of trials at each channel

% define number of permutation
num = 200;
EEG_null = EEG;
for n = 1:num
    disp(n)
    % shuffle trials at each chan
    for c=1:size(EEG.data,1)
        EEG_null.data(c,:,:)=EEG.data(c,:,randperm(size(EEG.data,3)));
    end
    
    % initialize fieldtrip
    ft_defaults;
    data_null=eeglab2fieldtrip(EEG_null,'preprocessing');

    % select elecs
    cfgc = [];
    cfgc.channel = elecs;
    data_select = ft_selectdata(cfgc, data);

    % compute encoding PSI
    cfgt=[];
    cfgt.latency = [1 3];
    data_select = subtracterp(data_select, 1); % subtract ERP
    data_enc = ft_selectdata(cfgt, data_select);

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foi = 1:1:40; % 1-40 Hz
    cfg.tapsmofrq = 1; % 1-Hz half-band
    cfg.pad = 4; % pad to 4 s
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    data_enc_cs = ft_freqanalysis(cfg, data_enc); % spectral decomposition

    cfgp = [];
    cfgp.method = 'psi';
    cfgp.bandwidth = 5;
    psi_enc_null{n} = ft_connectivityanalysis(cfgp, data_enc_cs);

    % compute maintenance PSI
    cfgt=[];
    cfgt.latency = [3 6];
    data_maint = ft_selectdata(cfgt, data_select);

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foi = 1:1:40; % 1-40 Hz
    cfg.tapsmofrq = 1; % 1-Hz half-band
    cfg.pad = 6; % pad to 6 s
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    data_maint_cs = ft_freqanalysis(cfg, data_maint); % spectral decomposition

    cfgp = [];
    cfgp.method = 'psi';
    cfgp.bandwidth = 5;
    psi_maint_null{n} = ft_connectivityanalysis(cfgp, data_maint_cs);
    clear data_*    
end

% save by subject
save([saveDir subid '_psi_null'], 'psi_enc_null', 'psi_maint_null');

%%% get PSI-zscore based on all subjects
% get raw PSI cycle through subjects
npp = 14;
psi_enc_subjects = [];
psi_maint_subjects = [];
for pp = 1:npp
    subid = strcat('s',num2str(pp),'_psi.mat');
    load([dataDir subid]);    
    psi_enc_subjects(:,:,:,pp)=psi_enc.psispctrm; 
    psi_maint_subjects(:,:,:,pp)=psi_maint.psispctrm;
end

% reshape raw psi by channel-pairs * frequency * subjects
psi_enc_subjects = reshape(psi_enc_subjects,[size(psi_enc_subjects,1)*size(psi_enc_subjects,2), size(psi_enc_subjects,3), size(psi_enc_subjects,4)]);
psi_maint_subjects = reshape(psi_maint_subjects,[size(psi_maint_subjects,1)*size(psi_maint_subjects,2), size(psi_maint_subjects,3), size(psi_maint_subjects,4)]);

% get null PSI cycle through subjects
psi_enc_subjects_null = [];
psi_maint_subjects_null = [];
for pp = 1:npp
    subid = strcat('s',num2str(pp),'_psi_null.mat');
    load([dataDir subid]);
    for n = 1:num            
        psi_enc_subjects_null(:,:,:,n,pp)=psi_enc_null{n}.psispctrm; 
        psi_maint_subjects_null(:,:,:,n,pp)=psi_maint_null{n}.psispctrm;
    end
end

% reshape null psi by channel-pairs * frequency * permutations * subjects
psi_enc_subjects_null = reshape(psi_enc_subjects_null,[size(psi_enc_subjects_null,1)*size(psi_enc_subjects_null,2), size(psi_enc_subjects_null,3), size(psi_enc_subjects_null,4), size(psi_enc_subjects_null,5)]);
psi_maint_subjects_null = reshape(psi_maint_subjects_null,[size(psi_maint_subjects_null,1)*size(psi_maint_subjects_null,2), size(psi_maint_subjects_null,3), size(psi_maint_subjects_null,4), size(psi_maint_subjects_null,5)]);

% get PSI-zscore based on all subjects
psi_enc_subjects_avg = mean(mean(psi_enc_subjects,3));
psi_enc_subjects_nullavg = mean(mean(psi_enc_subjects_null,4));
psi_enc_subjects_zscore = (psi_enc_subjects_avg-mean(psi_enc_subjects_nullavg,3))./std(psi_enc_subjects_nullavg,1,3);

psi_maint_subjects_avg = mean(mean(psi_maint_subjects,3));
psi_maint_subjects_nullavg = mean(mean(psi_maint_subjects_null,4));
psi_maint_subjects_zscore = (psi_maint_subjects_avg-mean(psi_maint_subjects_nullavg,3))./std(psi_maint_subjects_nullavg,1,3);

% save zscore
save([saveDir 'PSI_zscore.mat'], 'psi_enc_subjects_zscore','psi_maint_subjects_zscore');

end


