% Matlab script accompanying the paper 
% "Optimal design of experiments by combining coarse and fine measurements"
% arXiv:1702.06001
%
% Alpha Lee, Michael Brenner and Lucy Colwell 
% 24/05/2017
%
% This script produces Figure 1 of the paper. 

clear all 

% use the glmnet package for ridge regression (http://web.stanford.edu/~hastie/glmnet_matlab/) 
addpath(strcat(pwd,'/glmnet_matlab'))
options.alpha =0 ; %ridge rather than lasso 
options.lambda_min = 0.000001; 

%Import data
fp_avalon = importdata('datMat_avalon_solubility'); %avalon fingerprint 
fp_maccs = importdata('datMat_MACCS_solubility'); %MACCS fingerprint 
fp_morgan = importdata('datMat_Morgan3_solubility_r6_1024bit'); %Morgan 6 fingerprint 

%the solubility data taken from ESOL (J. Delaney, J. Chem. Inf. Comput. Sci., 2004, 44 (3), pp 1000?1005)
sol_dat = importdata('sol_data'); 

%Combine the fingerprints
fp = [fp_avalon fp_morgan fp_maccs]; 

%Create the coarse dataset 
soluable_fp = fp(find(sol_dat>-2),:); 
soluable_data = sol_dat(find(sol_dat>-2));  

insoluable_fp = fp(find(sol_dat<-4),:); 
insoluable_data = sol_dat(find(sol_dat<-4)); 

%z-score and remove variables that are 0 for more than 99% of the sample
[soluable_fp, mu_s, var_s]  = zscore(soluable_fp); 

ind_s = find(mu_s < 0.01);
soluable_fp(:,ind_s) = [];
mu_s(ind_s) = []; 
var_s(ind_s) = []; 

C_sol = soluable_fp'*soluable_fp/size(soluable_fp,1); 
[v1,d1] = eig(C_sol);


[insoluable_fp, mu_is, var_is]  = zscore(insoluable_fp); 
ind_is = find(mu_is < 0.01); 

insoluable_fp(:,ind_is) = [];
mu_is(ind_is) = []; 
var_is(ind_is) = []; 

C_insol = insoluable_fp'*insoluable_fp/size(insoluable_fp,1); 
[v2,d2] = eig(C_insol);

% Remove zero modes as the correlation matrix is rank deficient
v1 = v1(:,(end-length(soluable_data)-1):end);  
v2 = v2(:,(end-length(insoluable_data)-1):end);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pred = []; 
solval = []; 
        for ii =1:10 
            
            %randomly partition data into training and verification sets 
            rind = randperm(size(fp,1)); 
            nsam = round(size(fp,1)/10); 
            fpver= fp(rind(1:nsam),:); %verification set
            fptrain= fp(rind((nsam+1):end),:); %training set 

            %Compute the independent variables for the training set 
            fp1 = fptrain; 
            fp1(:,ind_s) = [] ; 

            fp2 = fptrain; 
            fp2(:,ind_is) = [] ; 


            coeff = [fptrain (fp1*v1).^2 (fp2*v2).^2];
            
            %Compute the ridge parameter using 10-fold cross validation on training set
            fit = cvglmnet(coeff , sol_dat(rind((nsam+1):end)),[],options,'mae',10) ;

            %Compute out-of-sample prediction 
            fp1 = fpver; 
            fp1(:,ind_s) = [] ; 
            
            fp2 = fpver; 
            fp2(:,ind_is) = [] ; 

            coeff = [fpver (fp1*v1).^2 (fp2*v2).^2];

            pred = [pred ; cvglmnetPredict(fit, coeff)]; 
            solval = [solval ; sol_dat(rind(1:nsam))]; 

        end 


figure(1) 
plot(solval,pred,'o')
xlabel('log(S / mol L^{-1})') 
ylabel('predicted log(S / mol L^{-1})')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1B 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of quantitative measurements that we include 
N = 100:50:size(fp,1); 

parfor kk = 1:length(N)   
    
for mm = 1:10 
      
    %randomly shuffle the samples 
    rr = randperm(size(fp,1)); 
    fp_cv = fp(rr,:); 
    sol_dat_cv = sol_dat(rr);
    
    %take only N(kk) quantitative measurements   
    fp_N = fp_cv(1:N(kk),:); 
    sol_dat_N = sol_dat_cv(1:N(kk)); 
    pred = []; 
    solval=[];
    
        for ii =1:10 %this loop estimates the out-of-sample error 
            
            %randomly partition data into training and verification sets 
            rind = randperm(size(fp_N,1)); 
            nsam = round(N(kk)/10); 
            fpver= fp_N(rind(1:nsam),:); %verification set
            fptrain= fp_N(rind((nsam+1):end),:); %training set 

            %Compute the independent variables for the training set 
            fp1 = fptrain; 
            fp1(:,ind_s) = [] ; 

            fp2 = fptrain; 
            fp2(:,ind_is) = [] ; 


            coeff = [fptrain (fp1*v1).^2 (fp2*v2).^2];
            
            %Compute the ridge parameter using 10-fold cross validation on training set
            fit = cvglmnet(coeff , sol_dat_N(rind((nsam+1):end)),[],options,'mae',10) ;

            %Compute out-of-sample prediction 
            fp1 = fpver; 
            fp1(:,ind_s) = [] ; 
            
            fp2 = fpver; 
            fp2(:,ind_is) = [] ; 

            coeff = [fpver (fp1*v1).^2 (fp2*v2).^2];

            pred = [pred ; cvglmnetPredict(fit, coeff)]; 
            solval = [solval ; sol_dat_N(rind(1:nsam))]; 

        end 
        
        %compute mean absolute error 
        acc(kk,mm) = mean(abs(pred-solval)); 
end
end 

%compute the standard error of the mean 
errors = std(acc,[],2)/sqrt(size(acc,1)); 

figure(2)
errorbar(N,mean(acc,2),errors) 
xlabel('Number of quantitative measurements') 
ylabel('Mean absolute error') 
