%% Compute encoding-maintenance similarity (EMS) matrix based on the power.
% Key: subid = 's1' (required)

function [] = EMS_compute(subid)
% EMS_compute('s1');

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%

dataDir = [path '/power/'];
saveDir = [path '/EMS/'];
mkdir(saveDir);

% load data
load([dataDir subid '/power']);

% prepare power matrix
zpower=zpower.powspctrm;
num_trial=size(zpower,1);

% power
T=1001:8000;
zpower=zpower(:,:,:,T);

% 200ms time window with step of 10ms
wt=size(T,2);
wn(1,:)=[1:10:wt-199];
wn(2,:)=wn(1,:)+199;

% 1-40hz
fre=1:100;
ind=find(fre>=1 & fre<=40);
RSA=[];
for n=1:num_trial
    disp(n)
    for t1=1:length(wn)
        for t2=1:length(wn)
            a=mean(zpower(n,:,ind,wn(1,t1):wn(2,t1)),4);
            a=a(:);
            b=mean(zpower(n,:,ind,wn(1,t2):wn(2,t2)),4);
            b=b(:);
            r=corr(a,b,'type','spearman');
            z=0.5*log((1+r)/(1-r));
            RSA(t1,t2,n)=z;
            clear a b r z
        end
    end
end

% get EMS matrix across trials
inde=find(T<=3000);% encoding
indm=find(T>3000 & T<=6000); % maintenance
EMS_avg=mean(RSA(inde,indm,:),3);

% save by subject
save([saveDir subid '_EMS'], 'EMS_avg');

end

%% Plot EMS matrix across subjects in hippocampus or amygdala.

function [] = EMS_plot()

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%
dataDir = [path '/EMS/'];
saveDir = [path '/EMS/'];

% define participants
npp = 14;

% load EMS matrix cycle through participants
EMS_subjects = [];
for pp = 1:npp
    subid = strcat('s',num2str(pp),'_EMS.mat');
    load([dataDir subid]);
    EMS_subjects(:,:,pp)=EMS_avg;    
end

% save EMS_subjects in hippocampus or amygdala
save([saveDir 'XX_EMS_subjects.mat'], 'EMS_subjects');

% average EMS matrix across subjects
EMS_avg = mean(EMS_subjects,3);

% plot EMS matrix
wn(1,:)=[1:10:7000-199];
wn(2,:)=wn(1,:)+199;

T=1:7000;
t=T(wn(1,:));

inde=find(t<=2000);
indm=find(t>2000 & t<=5000);

figure
colormap('jet')
imagesc(t(indm),t(inde),EMS_avg);
xlabel('Maintenance Time (ms)');
ylabel('Encoding Time (ms)');
xlim([2000 5000]);
ylim([0 2000]);
set(gca,'YDir','normal');
caxis([0,0.35]); 
colorbar('location','EastOutside');
set(gca,'fontsize',12,'fontweight','bold');
titlename=['XX-EMS']; % hippocampus or amygdala
title(titlename,'fontsize',12,'fontweight','bold');
name=[titlename,'.fig'];
savefig(name)

end

%% Compare EMS matrix between the hippocampus and the amygdala

% subfunction: permutest,written by Edden M. Gerber, lab of Leon Y. Deouell, 2014
% Key: input EMS matrix from the hippocampus as well as the amygdala

function [] = EMS_compare(EMS_hippocampus, EMS_amygdala)

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%
dataDir = [path '/EMS/'];
saveDir = [path '/EMS/'];

% load EMS matrix
hipp = load([dataDir 'hippocampus_EMS_subjects.mat']);
amy = load([dataDir 'amygdala_EMS_subjects.mat']);
EMS_hipp = hipp.EMS_subjects;
EMS_amy = amy.EMS_subjects;

% compare EMS matrix using cluster-based permutation test
[clusters, p_values, t_sums, ~ ] = permutest( EMS_hipp, EMS_amy, true, ...
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
EMS_hippsig = [];
EMS_amysig = [];
for i = 1:npp
    % hipp
    tmp = [];
    tmp = EMS_hipp(:,:,i);
    tmp = tmp.*sig;
    EMS_hippsig(i,1) = sum(tmp(:))/ sum(sum(tmp~=0));
    
    % amy
    tmp = [];
    tmp = EMS_amy(:,:,i);
    tmp = tmp.*sig;
    EMS_amysig(i,1) = sum(tmp(:))/ sum(sum(tmp~=0));   
end

% compare EMS_sig between the hippocampus and the amygdala using paired t-test
[h,p,ci,stats] = ttest(EMS_amysig,EMS_hippsig);

% save statistics
save([saveDir 'EMS_comparision.mat'], 'h','p','ci','stats');

end
