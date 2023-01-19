%% Compute task-induced power. 

% subfunctions: 
% BOSC_tf_power,from the Better OSCillation detection (BOSC) library.
% zbaseline, from the function shared by Johnson EL, et al.,2018.

% Key:
% subid = 's1' (required)
% elecs_index = [index of hippocampus or amygdala in channels] (required)

function [] = power_compute(subid, elecs_index)
% power_compute('s1', [1 2 9 10])

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%

dataDir = [path '/data_indiv/'];
saveDir = [path '/power/'];
mkdir(saveDir);

% load data, format of data: channels * timepoints * trials
data=load([dataDir subid '/data_derived']);
data_ind=data(elecs_index,:,:); 
[nchan,npts,ntrials]=size(data_ind);

% compute power
% define parameters
Fre = [1:100];
Fsample=1000;
wavenum=6;

% compute power at each trial and each channel
power=zeros(ntrials,nchan,length(Fre),npts);
for i=1:nchan
    for j=1:ntrials
        tmp=[];
        tmp=data_ind(i,:,j);
        [B,T,F]=BOSC_tf_power(tmp,Fre,Fsample,wavenum);
        power(j,i,:,:)=B;
        clear B
    end
end
clear data_ind

% zscore on baseline
% -0.5~0s as baseline
T=T-1;
wm.powspctrm=power;
bt=find(T>=-0.5 & T<0);
base.powspctrm=power(:,:,:,bt);
zpower = zbaseline(wm, base);
clear power wm bt base

% save by subject
save([saveDir subid '_power'], 'zpower');

end

%% Plot task-induced power across subjects in hippocampus or amygdala. 

function [] = power_plot()

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%
dataDir = [path '/power/'];

% define participants
npp = 14;

% load power matrix cycle through participants
power_subjects = [];
for pp = 1:npp
    subid = strcat('s',num2str(pp),'_power.mat');
    load([dataDir subid]);
    % average power across trials and channels
    power_subjects(:,:,pp)=squeeze(mean(mean(zpower,2)));    
end

% average power across subjects
power_avg=mean(power_subjects,3);

% plot power
t=-1000:6999;
f=1:100;
figure
colormap(jet)
imagesc(t,f,power_avg);
xlabel('Times(ms)','fontsize',12,'fontweight','bold');    
ylabel('Frequency(Hz)','fontsize',12,'fontweight','bold');
set(gca,'YDir','normal');
set(gca,'XLim',[-500 7000]);
set(gca,'YLim',[1 100]);
line([0 0],[1 100],'LineWidth',1,'Color','w'); %stim onset;
line([2000 2000],[1 100],'LineWidth',1,'Color','w'); %stim onset;
line([5000 5000],[1 100],'LineWidth',1,'Color','w'); %stim onset;
caxis([-2,2]); 
colorbar('location','EastOutside');
set(gca,'fontsize',14,'fontweight','bold');
titlename=['XX-zpower'];% hippocampus or amygdala 
title(titlename,'fontsize',14,'fontweight','bold');
saveas(gca,titlename); 

end


