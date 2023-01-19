%% Compute encoding-encoding dissimilarity (EED) matrix based on the power.
% Key:
% subid = 's1' (required)

function [] = EED_compute(subid)
% EED_compute('s1');

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%

dataDir = [path '/power/'];
saveDir = [path '/EED/'];
mkdir(saveDir);

% load data
load([dataDir subid '/power']);

% prepare power matrix
zpower=zpower.powspctrm;
num_trial=size(zpower,1);

% pick trial-pairs from the number of trials
s=nchoosek([1:num_trial],2);
ss=randperm(length(s),length(s));

% encoding power
T=1001:3200;
zpower_enc=zpower(:,:,:,T);

% 200ms time window with step of 10ms
wt=size(T,2);
wn(1,:)=[1:10:wt-199];
wn(2,:)=wn(1,:)+199;

% 1-40hz
fre=1:100;
ind=[1:40];
EED=[];
for n=1:100 % we use 100 trial-pairs
    disp(n)
    for t1=1:length(wn)
        for t2=1:length(wn)
            a=mean(zpower_enc(s(ss(n),1),:,ind,wn(1,t1):wn(2,t1)),4);
            a=a(:);
            b=mean(zpower_enc(s(ss(n),2),:,ind,wn(1,t2):wn(2,t2)),4);
            b=b(:);
            r=corr(a,b,'type','spearman');            
            z=0.5*log((1+r)/(1-r));            
            EED(t1,t2,n)=1-z;% dissimilarity
            clear a b r z
        end
    end
end

% get averaged EED matrix across trial-pairs
EED_avg=mean(EED,3);

% save by subject
save([saveDir subid '_EED'], 'EED_avg');

end

%% Plot EED matrix across subjects in hippocampus or amygdala.

function [] = EED_plot()

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%
dataDir = [path '/EED/'];
saveDir = [path '/EED/'];

% define participants
npp = 14;

% load EED matrix cycle through participants
EED_subjects = [];
for pp = 1:npp
    subid = strcat('s',num2str(pp),'_EED.mat');
    load([dataDir subid]);
    EED_subjects(:,:,pp)=EED_avg;    
end

% save EED_subjects in hippocampus or amygdala
save([saveDir 'XX_EED_subjects.mat'], 'EED_subjects');

% average EED matrix across subjects
EED_avg = mean(EED_subjects,3);

% plot EED matrix
wn(1,:)=[1:10:2200-199];
wn(2,:)=wn(1,:)+199;

T=1:2200;
t=T(wn(1,:));

figure
colormap('jet')
imagesc(t,t,EED_avg);
xlabel('Encoding Time (ms)');
ylabel('Encoding Time (ms)');
xlim([t(1) t(end)]);
ylim([t(1) t(end)]);
set(gca,'YDir','normal');
caxis([0.8,1]); 
colorbar('location','EastOutside');
set(gca,'fontsize',12,'fontweight','bold');
titlename=['XX-EED']; % hippocampus or amygdala
title(titlename,'fontsize',12,'fontweight','bold');
name=[titlename,'.fig'];
savefig(name)

end

%% Compare EED matrix between the hippocampus and the amygdala

% subfunction: permutest,written by Edden M. Gerber, lab of Leon Y. Deouell, 2014
% Key: input EED matrix from the hippocampus as well as the amygdala

function [] = EED_compare(EED_hippocampus, EED_amygdala)

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%
dataDir = [path '/EED/'];
saveDir = [path '/EED/'];

% load EED matrix
hipp = load([dataDir 'hippocampus_EED_subjects.mat']);
amy = load([dataDir 'amygdala_EED_subjects.mat']);
EED_hipp = hipp.EED_subjects;
EED_amy = amy.EED_subjects;

% compare EED matrix using cluster-based permutation test
[clusters, p_values, t_sums, ~ ] = permutest( EED_hipp, EED_amy, true, ...
    0.05, 1000, true );

% extract significant clusters
sig=zeros(200,200);
for i=1:length(p_values)
    if p_values(i)<0.05
        tmp=[];
        tmp=clusters{i};
        for k=1:size(tmp,1)
            a=tmp(k,:);
            sig(a(1),a(2))=1;
        end
    end
end

% define number of subjects
npp = 14;

% extract significance at each subject
EED_hippsig = [];
EED_amysig = [];
for i = 1:npp
    % hipp
    tmp = [];
    tmp = EED_hipp(:,:,i);
    tmp = tmp.*sig;
    EED_hippsig(i,1) = sum(tmp(:))/ sum(sum(tmp~=0));
    
    % amy
    tmp = [];
    tmp = EED_amy(:,:,i);
    tmp = tmp.*sig;
    EED_amysig(i,1) = sum(tmp(:))/ sum(sum(tmp~=0));   
end

% compare EED_sig between the hippocampus and the amygdala using paired t-test
[h,p,ci,stats] = ttest(EED_amysig,EED_hippsig);

% save statistics
save([saveDir 'EED_comparision.mat'], 'h','p','ci','stats');

end
