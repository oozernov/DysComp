%% 
%   GANONG DATA -- LOGISTIC REG FIT 
% 
%   - General analyses of the Ganong effects in four participant groups 
%   - Finding the optimal logistic fit parameters (slope and inflection)
%   - Figure plots for logistic fit + raw prop. resp data
%   
%   author: Sung-Joo Lim (sungjoo@bu.edu)
%
%% set up

clear all
close all

addpath('~/Dropbox/matlab_code_reposit/')

% raw data -- single-trial data cleaned out (after excluding subjects who did not
% complete the testing)
Ganong_masterdataT = readtable('~/Dropbox/Doc/Projects_at_BU/other/Ganong-Ola/Ola_data/masterdata_MATLAB_FINAL.csv',...
    'Delimiter',',','ReadVariableNames',1);


Ganong_masterdata = table2array(Ganong_masterdataT);

allSubjects = unique(Ganong_masterdataT.PartID);
SubjectGroupList = zeros(length(allSubjects),3);


% 1 = g/k, 2 = d/t
isPhonemeG = unique(Ganong_masterdataT.word_category_gk); 

% voiced word-nw (1), voiced nw-word (2), nonword pair (3)
wordPairs = unique(Ganong_masterdataT.wordpair_type); 
% GK: 1 = gift, 2 = giss, 3 = gith
% DT: 1 = dash, 2 = dask, 3 = dath


% 1-7 continnum steps
continuum_steps = unique(Ganong_masterdataT.n_step); 



%% get average prop. voice resp for each subject 

% keep track of prob. resp. voiced
meanPropData = zeros(length(allSubjects),3+length(continuum_steps),...
    length(wordPairs), length(isPhonemeG));

% keep the count of no resp trials
nanCounts = zeros(length(allSubjects),3+length(continuum_steps),...
    length(wordPairs), length(isPhonemeG)); 


% get prop. resp of G on each continuum step
for ss_ind = 1:length(allSubjects)
    
    ss = allSubjects(ss_ind);
    
    ind = (Ganong_masterdata(:,1) == ss);
    adultGroup =  Ganong_masterdata(ind == 1,3); % 1 = adult, 0 = child
    dysGroup =  Ganong_masterdata(ind == 1,4); % 1 = dys, 0 = typ
    
    SubjectGroupList(ss_ind,:) = [ss, adultGroup(1), dysGroup(1)];
    
    
    for isG = 1:length(isPhonemeG) % is GK or DT condition
        
        for wordpair_ind = 1:length(wordPairs) 
            
            for step = 1:length(continuum_steps)
                
                index = (Ganong_masterdata(:,1) == ss) .* (Ganong_masterdata(:,5) == isG) .* ...
                    (Ganong_masterdata(:,7) == wordpair_ind) .* (Ganong_masterdata(:,6) == step);
                
                respVoiced = nanmean(Ganong_masterdata(index == 1,8)); % skip nans -- no resp.
                
                
                % update the arrays
                meanPropData(ss_ind, 1:3, wordpair_ind, isG) = [ss, adultGroup(1), dysGroup(1)]; % condition
                meanPropData(ss_ind, step+3, wordpair_ind, isG) = respVoiced; % 
                
                nanCounts(ss_ind, 1:3, wordpair_ind, isG) = [ss, adultGroup(1), dysGroup(1)];
                nanCounts(ss_ind, step+3, wordpair_ind, isG) = sum(isnan(Ganong_masterdata(index == 1,8)));
                
                clear respVoiced;
            end
        end
    end
    
end    
    


%% logistic fitting 
%
% find the best logistic fitting parameters (given the bounds) using
% non-linear least square fitting optimization


% logstic parameters slope and inflection
params = [0 0];
lb = [0 1]; % lower bounds
ub = [inf 7]; % upper bounds

x = 1:7;



% logistic fit optimization setting
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'Display','off');


% keep the logistic fitting parameters
logistic_fitted_params = zeros(length(allSubjects), 2+3, 3, 2);

% save the model fitted data 
model_fitted_vals = zeros(length(allSubjects),length(continuum_steps),...
    length(wordPairs), length(isPhonemeG));


% loop through each subject, each condition to find the logistic fit.
for ss_ind = 1:length(allSubjects)
    
    ss = allSubjects(ss_ind);
    
    ind = (Ganong_masterdata(:,1) == ss);
    adultGroup =  Ganong_masterdata(ind == 1,3); % 1 = adult, 0 = child
    dysGroup =  Ganong_masterdata(ind == 1,4); % 1 = dys, 0 = typ
    
    for isG = 1:length(isPhonemeG)
        
        for wordpair_ind = 1:length(wordPairs)
            

            % prop. resp voiced data is flipped in order to find the
            % positive slope of the logistic function fit
            prop_voiced_data = flip(meanPropData(ss_ind, 4:10, wordpair_ind, isG));
            
            % non-linear least square fitting
            [coeff resnorm] = lsqcurvefit(@logistic,params,x,prop_voiced_data,lb,ub,options);
            
            % save the fitted parameters in the matrix
            logistic_fitted_params(ss_ind,1:3,wordpair_ind, isG) = [ss, adultGroup(1), dysGroup(1)];
            logistic_fitted_params(ss_ind,4:5,wordpair_ind, isG) = coeff;
            
            % model fitted output using the parameters
            model_fitted_vals(ss_ind,:,wordpair_ind, isG) = logistic(coeff, x);
            
        end
    end
end



%% fitted data save to R

logistic_params_GK = [logistic_fitted_params(:,:,1,1), ones(length(logistic_fitted_params(:,:,1,1)),1);...
    logistic_fitted_params(:,:,2,1),2*ones(length(logistic_fitted_params(:,:,1,1)),1);...
    logistic_fitted_params(:,:,3,1),3*ones(length(logistic_fitted_params(:,:,1,1)),1)];

logistic_params_DT = [logistic_fitted_params(:,:,1,2), ones(length(logistic_fitted_params(:,:,1,1)),1);...
    logistic_fitted_params(:,:,2,2), 2*ones(length(logistic_fitted_params(:,:,1,1)),1);...
    logistic_fitted_params(:,:,3,2), 3*ones(length(logistic_fitted_params(:,:,1,1)),1)];


logistic_params_GK = [logistic_params_GK, ones(length(logistic_params_GK),1)];
logistic_params_DT = [logistic_params_DT, zeros(length(logistic_params_GK),1)];


variable_names = {'subj_id','isAdult','isDys','slope','inflection','wordPairType','isGK'};
logistic_params_GKTable = array2table(logistic_params_GK,'VariableNames',variable_names);
logistic_params_DTTable = array2table(logistic_params_DT,'VariableNames',variable_names);

% combined dataset
variable_names = {'subj_id','isAdult','isDys','slope','inflection','wordPairType','isGK'};
logistic_params_GKTable = array2table(logistic_params_GK,'VariableNames',variable_names);
logistic_params_DTTable = array2table(logistic_params_DT,'VariableNames',variable_names);
logistic_params_AllDataTable = array2table([logistic_params_GK; logistic_params_DT],...
    'VariableNames',variable_names);


% SAVE DATA to csv -- this is the same dataset that I provided.
% writetable(logistic_params_GKTable, '~/Dropbox/Doc/Projects_at_BU/other/Ganong-Ola/Ola_data/MATLAB-logistic/GK_matlab_logistic_fit_wBounds_flip_FINAL.csv',...
%     'Delimiter',',');
% writetable(logistic_params_DTTable, '~/Dropbox/Doc/Projects_at_BU/other/Ganong-Ola/Ola_data/MATLAB-logistic/DT_matlab_logistic_fit_wBounds_flip_FINAL.csv',...
%     'Delimiter',',');
% writetable(logistic_params_AllDataTable, '~/Dropbox/Doc/Projects_at_BU/other/Ganong-Ola/Ola_data/MATLAB-logistic/allDATA_matlab_logistic_fit_wBounds_flip_FINAL.csv',...
%     'Delimiter',',');



%% Figures -- logistic curves for each subject / each condition


% SubjectGroupList

isG = 1; 


for isAdult = 0:1
    for isDys = 0:1
        
        index = (SubjectGroupList(:,2) == isAdult) .* (SubjectGroupList(:,3) == isDys);
        
        subject_list = SubjectGroupList(index == 1,1);
        
        
        figure;
        
        for ss_ind = 1:length(subject_list)
            ss = subject_list(ss_ind);
            subj_ind = (allSubjects == ss);
            
            
            subplot(4,8,ss_ind); hold on;
            
            for wordpair_ind = 1:2 %length(wordPairs)
                
                % GK: 1 = gift, 2 = giss, 3 = gith
                % DT: 1 = dash, 2 = dask, 3 = dath

                prop_voiced_data = flip(meanPropData(subj_ind, 4:10, wordpair_ind, isG));
                
                if(wordpair_ind == 1)
                    plot(1:7, prop_voiced_data,'bo');
                    plot(1:7, model_fitted_vals(subj_ind,:, wordpair_ind, isG),'b');
                elseif(wordpair_ind == 2)
                    plot(1:7, prop_voiced_data,'ro');
                    plot(1:7, model_fitted_vals(subj_ind,:, wordpair_ind, isG),'r');
                else
                    plot(1:7, prop_voiced_data,'ko');
                    plot(1:7, model_fitted_vals(subj_ind,:, wordpair_ind, isG),'k');
                end
            end
            
            xlim([1 7]); ylim([0 1]);
            title(int2str(ss));
            if(ss_ind == 1)
                if(isG == 1)
                    legend('gift',' ', 'giss',' ', 'gith',' ');
                else
                     legend('dash',' ', 'dask',' ', 'dath',' ');
                end
                
            end
            
        end
    end
end




