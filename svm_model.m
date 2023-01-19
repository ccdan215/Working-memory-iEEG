% Support vector machine (SVM) classifier for WM load. 
% Requires LIBSVM package proposed by Chang et al.,2011.
% Take EED features at trial level as example.
% This classification was performed for hippocampus and amygdala,separately.
% The procedure is very similar to use EMS features as well as PSI feature to predict WM load.

function [] = svm_model(EED_features)

% INPUT PATH TO DATA DIRECTORY
path = %%INSERT%%

dataDir = [path '/EED/'];
saveDir = [path '/SVM/'];
mkdir(saveDir);

% load EED matrix with WM load at trial level from any region
EED_set4=load([dataDir 'EED_set4_trials.mat'],'EED_set4');
EED_set6=load([dataDir 'EED_set6_trials.mat'],'EED_set6');
EED_set8=load([dataDir 'EED_set8_trials.mat'],'EED_set8');

EED_set4=EED_set4.EED_set4;
EED_set6=EED_set6.EED_set6;
EED_set8=EED_set8.EED_set8;

EED_set4=reshape(EED_set4,[size(EED_set4,1)*size(EED_set4,2) size(EED_set4,3)]);
EED_set4=EED_set4';
EED_set6=reshape(EED_set6,[size(EED_set6,1)*size(EED_set6,2) size(EED_set6,3)]);
EED_set6=EED_set6';
EED_set8=reshape(EED_set8,[size(EED_set8,1)*size(EED_set8,2) size(EED_set8,3)]);
EED_set8=EED_set8';

% define number of cross-validation
num = 100;

% SVM training and testing
for i=1:num
    disp(i)
    A=EED_set4;
    B=EED_set6;
    C=EED_set8;

    % split all subjects' dataset into 70% as training and 30% as testing
    [train_set4_idx, ~, test_set4_idx] = dividerand(size(A,1), 0.7, 0, 0.3);
    test_set4=A(test_set4_idx,:);
    train_set4=A(train_set4_idx,:);
    [train_set6_idx, ~, test_set6_idx] = dividerand(size(B,1), 0.7, 0, 0.3);
    test_set6=B(test_set6_idx,:);
    train_set6=B(train_set6_idx,:);
    [train_set8_idx, ~, test_set8_idx] = dividerand(size(C,1), 0.7, 0, 0.3);
    test_set8=A(test_set8_idx,:);
    train_set8=A(train_set8_idx,:);
    
    % pool training dataset as well as testing dataset at each load together
    train_set4_label=ones(size(train_set4,1),1);            
    train_set6_label=2.*ones(size(train_set6,1),1);            
    train_set8_label=3.*ones(size(train_set8,1),1);
    train=[train_set4;train_set6;train_set8];
    train_label=[train_set4_label;train_set6_label;train_set8_label];
    test=[test_set4;test_set6;test_set8];
    test_label=[1*ones(size(test_set4,1),1);2*ones(size(test_set6,1),1);3*ones(size(test_set8,1),1)];

    % normalization training and testing dataset
    [Train,PS]=mapminmax(train');
    train_zscore=Train';
    Test=mapminmax('apply',test',PS);
    test_zscore=Test';

    % principle component analysis (PCA)
    [pc,score,latent,tsquare] = pca(train_zscore,'Algorithm','svd','Centered',false);
    tmp=cumsum(latent)./sum(latent);
    indd=find(tmp>=0.99);% keep principle components that explain 99% of the data
    train_pc=train_zscore*pc(:,1:indd(1));           
    test_pc=test_zscore*pc(:,1:indd(1));

    % training and testing model
    parameter=['-s ' num2str(1) ' -t ' num2str(0) ' -c ' num2str(1) ' -g ' num2str(0.1) ' -b 1'];
    model = svmtrain(train_label,train_pc,parameter);             
    [~,acc,dec_values] = svmpredict(test_label,test_pc,model);    
    predict_accall(i)=acc(1);% get classification accuracy for each cross-validation
    clear A B C test train pc score latent tsquare tmp indd PS acc dec_values          
end

% save accuracy
save([saveDir 'EED_svm_accuracy.mat'], 'predict_accall');

end
