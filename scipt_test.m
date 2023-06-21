clear all;
close all; 
clc;
%% loading simtb sources
sp.M =4;
sp.nT = 300; 
sp.nV = 50;
%rho = 6; % sources generated using rho = 6; 
load simtb_sources
 
%% parameters
nCC = 8; % common components
tstd  = sqrt(0.3); %% noise vector generation for time sources
sstd  = sqrt(0.01); %0.0025
srcs = nCC+1; % 8 common + 1 unique source
K = srcs; %+9
K2 = srcs;
nIter = 15;
Dp = dctbases(sp.nT,sp.nT);
Dp = Dp(:,2:end);
D_tol = 1e-8;


%% data generation
iM = 1; % first subject
SM_gw = reshape(iSM(iM,1:nCC,:),nCC,sp.nV*sp.nV);
TC_gw = reshape(iTC(iM,:,1:nCC),sp.nT,nCC);
for sub=1:sp.M
    iSM_sw(sub,:,:) = reshape(iSM(sub,[1:nCC nCC+sub],:),nCC+1,sp.nV*sp.nV);
    iTC_sw(sub,:,:) = reshape(zscore(iTC(sub,:,[1:nCC-2 (nCC-1)+3*(sub-1):(nCC-1)+3*(sub-1)+2])),sp.nT,nCC+1);
    SM_gw = [SM_gw; reshape(iSM_sw(sub,nCC+1,:),1,sp.nV*sp.nV)];  % Pulling out the 9th source from sw-SMs
    TC_gw = [TC_gw  reshape(iTC_sw(sub,:,nCC+1),sp.nT,1)];
end


rng('default')
rng(5,'twister') 
for sub=1:sp.M
    TC_sw{sub} = reshape(iTC_sw(sub,:,:),sp.nT,srcs);
    SM_sw{sub} = reshape(iSM_sw(sub,:,:),srcs,sp.nV*sp.nV);
    Y{sub} = (TC_sw{sub}+tstd(1)*randn(sp.nT,nCC+1))*(SM_sw{sub}+sstd(1)*randn(nCC+1,sp.nV*sp.nV));
    Y{sub} = Y{sub}-repmat(mean(Y{sub}),size(Y{sub},1),1);
end


%% ACSD
tic
for i =1:sp.M
    [Zt{1}(:,:,i),Zs{1}(:,:,i),E{1}(:,i),C{1}(:,i)]= my_ACSD(Y{i},Dp(:,1:K),12,nIter,TC_sw{i},SM_sw{i}); % and 2 %=
end
toc


%% ssBSS
params1.K = K;
params1.P = K; %8
params1.lam1 = 6; %6
params1.zeta1 = 30;
params1.Kp = 150;
params1.nIter = nIter;
params1.alpha = 10^-8; %1e-9
for i=1:sp.M
    [Zt{2}(:,:,i),Zs{2}(:,:,i),~,~,E{2}(:,i),C{2}(:,i)]=ssBSS_pre(Y{i},Dp,params1,TC_sw{i},SM_sw{i}); %_mod
end
Dq = [Zt{2}(:,:,1) Zt{2}(:,:,2) Zt{2}(:,:,3) Zt{2}(:,:,4)];
Xq = [Zs{2}(:,:,1); Zs{2}(:,:,2); Zs{2}(:,:,3); Zs{2}(:,:,4)];


%% swbDL
tic
for i =1:sp.M
    [Zt{3}(:,:,i),Zs{3}(:,:,i),E{3}(:,i),A1,B1,C{3}(:,i)]= swbDL(Y{i},Dq,Xq,nIter,K,8,TC_sw{i},SM_sw{i}); %Kc+Ks
end
toc


%% swsDL
tic
for sub =1:sp.M
    [Zt{4}(:,:,sub),Zs{4}(:,:,sub),E{4}(:,sub),A2,B2,C{4}(:,i)]= swsDL(Y{sub},Dq,Xq,nIter,K,24,16,TC_sw{sub},SM_sw{sub}); %Kc+Ks 12/4 32/64
end
toc


%% SW correlation coomputation
clear SMcorr_sw22; clear TCcorr_sw22
nA = 4;
for a =1:nA
    for sub = 1:sp.M
        [sD{a,sub},sX{a,sub},ind]=sort_TSandSM_spatial(TC_sw{sub},SM_sw{sub},Zt{a}(:,:,sub),Zs{a}(:,:,sub),srcs);
        for i =1:srcs
            TCcorr_sw(1,i,sub,a) =abs(corr(TC_sw{sub}(:,i),Zt{a}(:,ind(i),sub)));
            SMcorr_sw(1,i,sub,a) =abs(corr(SM_sw{sub}(i,:)',Zs{a}(ind(i),:,sub)'));
        end
    end
end
reshape(sum(mean(TCcorr_sw,3)),1,nA)
reshape(sum(mean(SMcorr_sw,3)),1,nA)



%% Plots
sDD{1} = TC_gw; sXX{1} = SM_gw;
tmptmpTC = TCcorr_sw; tmptmpSM = SMcorr_sw;
TCcorr_sw(:,:,1:4,1)=0;
SMcorr_sw(:,:,1:4,1)=0;
TCcorr_sw(:,:,1:4,2:nA+1)= tmptmpTC;
SMcorr_sw(:,:,1:4,2:nA+1)= tmptmpSM;
TCcorr_sw2(:,:,1:4,1)= TCcorr_sw(:,:,1:4,1);  TCcorr_sw2(:,10:12,1:4,1)=0;
SMcorr_sw2(:,:,1:4,1)= SMcorr_sw(:,:,1:4,1);  SMcorr_sw2(:,10:12,1:4,1)=0;

for i =2:nA+1
    sDD{i}= sD{i-1,1};  sDD{i}(:,10)= sD{i-1,2}(:,9); sDD{i}(:,11)= sD{i-1,3}(:,9); sDD{i}(:,12)= sD{i-1,4}(:,9);
    sXX{i}= sX{i-1,1};  sXX{i}(10,:)= sX{i-1,2}(9,:); sXX{i}(11,:)= sX{i-1,3}(9,:); sXX{i}(12,:)= sX{i-1,4}(9,:);
    TCcorr_sw2(:,1:9,1,i)= TCcorr_sw(:,:,1,i);
    TCcorr_sw2(:,10:12,1,i)= TCcorr_sw(:,9,2:4,i);
    SMcorr_sw2(:,1:9,1,i)= SMcorr_sw(:,:,1,i);
    SMcorr_sw2(:,10:12,1,i)= SMcorr_sw(:,9,2:4,i);
end

TC_corr= reshape(TCcorr_sw2(1,:,1,:),srcs+3,nA+1)';
SM_corr= reshape(SMcorr_sw2(1,:,1,:),srcs+3,nA+1)';

toprint = ['TC correlations are: [', repmat('%2.4f, ', 1, numel(sum(TC_corr',1))-1), '%2.4f]\n'];
fprintf(toprint, sum(TC_corr',1))
toprint = ['SM correlations are: [', repmat('%2.4f, ', 1, numel(sum(SM_corr',1))-1), '%2.4f]\n'];
fprintf(toprint, sum(SM_corr',1))

f = figure; 
f.Position = [100 100 1250 850]; 
my_subplots_horz(nA+1,srcs+3,sp.nV,sp.nV,TC_corr,SM_corr,sDD,sXX);

