% This function was used for single-trial decoding analysis
function single_trial_decoding

% listed are example files for illustration only
folder_name_hm4d={
    '20180605#3',[5,9,13,17,21,25,29],[]; %cno
    '20180912#3',[-1],[];
    };

folder_name_control={
    '2020101302', [1,5,9,13,21,25,29] [];
    '2021042802', [1,5,9,13,17,21,25,29] [];
    };

dir_names_compared={[pwd '\data\hm4d\produce_single_cell_data\'],...
    [pwd '\data\control\produce_single_cell_data\']};
folder_name_compared={folder_name_hm4d,folder_name_control};
file_prefix_compared={'PFCspike_','PFCspike_'};
dir_cate={'hm4d','control'};

bin_length_compared=[50,100,150]; % ms
pop_cell_dim_compared=[1:10];

pop_cell_dim_max=pop_cell_dim_compared(end);

decode_percent_all_compared=[];
decode_percent_overall_compared=[];
decode_percent_overall_valid_session_compared=[];
decode_percent_overall_timecurve_compared=[];
decode_percent_overall_timecurve_ts_compared=[];
decode_ssname=[];
decode_valid_ssname=[];
decode_info_bits_compared=[];
decode_pcs_compared=[];
decode_cell_firing_compared=[];
decode_stimulusID_compared=[];

result_file_sfx=[]; 

stimulusID= {'socialA';'socialB';'novel';'empty'};

stimulusIDGroup={...                 % 2-class comparison
    {'socialA';'novel'}...
    };


%% Section I: do single-trial decoding
if false % set as TRUE to re-run the decoding using the illustrative files

for bshuffleTest=[false,true]

for istimgrp=1:length(stimulusIDGroup)
stimulusID=stimulusIDGroup{istimgrp};

for idim=1:length(pop_cell_dim_compared)
pop_cell_dim=pop_cell_dim_compared(idim);

for idir=1:length(dir_names_compared)
dir_name=dir_names_compared{idir};
folder_name=folder_name_compared{idir};
file_prefix=file_prefix_compared{idir};

[folderLength,~]=size(folder_name);    
datalengthMili=3000;
rvAll_pairs=cell(folderLength,6);
rvAll_names=cell(folderLength,2);

vectorAll_during=cell(folderLength,4);
vectorAll_pre=cell(folderLength,4);

vectorAll_during_norm=cell(folderLength,4);
vectorAll_during_normss=cell(folderLength,1);

for ibin=1:length(bin_length_compared)
bin_length=bin_length_compared(ibin);

label_all_test=[];
label_all_judge=[];
ssname_all=[];
label_timecurve_test=[];
label_timecurve_judge=[];
score_timecurve_judge=[];
timestamp_timecurve=[];
pcs_all_sessions=[];
cell_firing_all_sessions=[];
stimulusID_all_sessions=[];

kk_used = [];
 for kk=1:folderLength
    
    clear raster_data raster_datab raster_labels dataVector2 databVector2 dataVector_pre2 dataVector_during2
    name_now=folder_name{kk,1};
    name_a=[file_prefix name_now(1:8) '_m' name_now(10:end) '*.mat'];  
    
    name_cell_a=dir([dir_name name_a]);
    
    if isempty(name_cell_a) 
        warning(['''' name_a '''not found!']);
        continue; 
    end
    
    load([name_cell_a(1).folder '\' name_cell_a(1).name],'raster_labels');
    
    TrialLength=length(raster_labels.stimulus_ID);
    stimulusID_now = raster_labels.stimulus_ID;
    
    if exist('stimulusIDGroup','var')
        for istimCur=1:length(stimulusID_now)
            for istimType=1:length(stimulusID)
                if ~isempty(strfind(stimulusID{istimType},stimulusID_now{istimCur}))
                    stimulusID_now{istimCur}=stimulusID{istimType};
                end
            end
        end
    end

    % only use sessions with more than 8 cells
    if length(name_cell_a)>=pop_cell_dim_max
        kk_used =[kk_used kk];
        % construct vector averaged by pre and during time
        dataVector_pre=zeros(length(name_cell_a),TrialLength);
        dataVector_during=zeros(length(name_cell_a),TrialLength);
        databVector_pre=zeros(length(name_cell_a),TrialLength);
        databVector_during=zeros(length(name_cell_a),TrialLength);
        % vector for all cells
        dataVector=zeros(length(name_cell_a),datalengthMili*2/bin_length,TrialLength);
        databVector=zeros(length(name_cell_a),datalengthMili*2/bin_length,TrialLength);

        dataVector_mean=zeros(length(name_cell_a),1);
        firing_raster=zeros(length(name_cell_a),datalengthMili*2,TrialLength);

        for mm=1:length(name_cell_a)
            clear raster_data raster_datab raster_labels raster_during_now raster_pre_now vector_length_now covVector_now covVector_now2 
            cellName_now=name_cell_a(mm).name;
            
            load([name_cell_a(mm).folder '\' cellName_now])
   
            data_count=zeros(TrialLength,datalengthMili/bin_length*2);
            datab_count=zeros(TrialLength,datalengthMili/bin_length*2);
            data_fr=zeros(TrialLength,datalengthMili/bin_length*2);
            datab_fr=zeros(TrialLength,datalengthMili/bin_length*2);
            for jj=1:TrialLength
                % use 100 ms window to count
                data_count(jj,:) = histcounts(find(raster_data(jj,:)==1),0:bin_length:datalengthMili*2) *1000/bin_length;
                datab_count(jj,:) = histcounts(find(raster_datab(jj,:)==1),0:bin_length:datalengthMili*2) *1000/bin_length;

                data_fr(jj,:)=data_count(jj,:)*1000/bin_length; % firing rate (Hz)
                datab_fr(jj,:)=datab_count(jj,:)*1000/bin_length; % firing rate (Hz)
            end
            data_pre_now = mean( data_fr(:,1:datalengthMili/bin_length),2);
            data_during_now = mean( data_fr(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2),2);
            datab_pre_now = mean( datab_fr(:,1:datalengthMili/bin_length),2);
            datab_during_now = mean( datab_fr(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2),2);
            
            dataVector_pre(mm,:)=data_pre_now;
            dataVector_during(mm,:)=data_during_now;
            databVector_pre(mm,:)=datab_pre_now;
            databVector_during(mm,:)=datab_during_now;
            
            dataVector(mm,:,:)=data_fr';
            databVector(mm,:,:)=datab_fr';
            
            dataVector_mean(mm,1)=mean(mean(data_fr));

            firing_raster(mm,:,:)=raster_data';
        end

        idx=dataVector_mean>=0.5;
        dataVector2 = dataVector(idx,:,:);
        databVector2 = databVector(idx,:,:);
        dataVector_during2 = dataVector_during(idx,:);
        dataVector_pre2 = dataVector_pre(idx,:);
        firing_raster = firing_raster(idx,:,:);
        
        if sum(dataVector_mean>=0.5)>=pop_cell_dim_max
            
            result_file_sfx=[result_file_sfx {'_high_var_cell'}];
            cat_firing=dataVector_during2;
            [~,idx_sorted]=sort(var(cat_firing,0,2),'descend');
            dataVector2 = dataVector2(idx_sorted(1:pop_cell_dim),:,:);
            databVector2 = databVector2(idx_sorted(1:pop_cell_dim),:,:);
            dataVector_during2 = dataVector_during2(idx_sorted(1:pop_cell_dim),:);
            dataVector_pre2 = dataVector_pre2(idx_sorted(1:pop_cell_dim),:);
            firing_raster = firing_raster(idx_sorted(1:pop_cell_dim),:,:);

            idx=contains(stimulusID_now,stimulusID,'IgnoreCase',true);
            stimulusID_now=stimulusID_now(idx);
            dataVector2=dataVector2(:,:,idx);
            firing_raster=firing_raster(:,:,idx);
            stimulus_prior=zeros(1,length(stimulusID));
            for idx=1:length(stimulusID)
                stimulus_prior(idx)=mean(strcmpi(stimulusID_now,stimulusID{idx}));
            end

            % limiting trial number (more than 8 trials required for each category)
            num_lim_trial=8;
            bskip=false;
            for idx=1:length(stimulusID)
                ist=find(strcmpi(stimulusID_now,stimulusID{idx}));
                if length(ist)<num_lim_trial
                    bskip=true;
                    break;
                end
            end
            if bskip
                continue;
            end
            result_file_sfx=[result_file_sfx {['_triallim' num2str(num_lim_trial)]}];

            session_test_labels=[];
            session_judge_labels=[];
            session_judge_scores=[];
            session_test_labels_timecurve=[];
            session_judge_labels_timecurve=[];
            session_judge_scores_timecurve=[];
            session_ts_timecurve=[];

            TrialLength=length(stimulusID_now);
            [ncell,nfea,nsmpl]=size(dataVector2);

            if bshuffleTest   % # shuffle labels
                stimulusID_now=stimulusID_now(randperm(length(stimulusID_now)));
                result_file_sfx=[result_file_sfx {'_shuffle'}];
            end

            % # save the original cell firing data and pcs after spatial filtering
            session_pcs=[];          % pc x feature x trial
            session_cell_firing=[]; % cell x feature x trial
            session_stimulusID=[];   % trial x 1

            % # leave one trial out
            dataSamples=dataVector2;
            idx_all_trials=1:nsmpl;
            for idx_test_trial=1:nsmpl
                idx_train_trial=setdiff(idx_all_trials,idx_test_trial);                

            % # do spatial filtering if no less than 5 cells are available         
            npcs=1;
            if length(stimulusID)==2
                datalabels=2*strcmpi(stimulusID_now,stimulusID{1})-1;
                sratio=sum(datalabels==1)/sum(datalabels==-1);
                if sratio>=1/5 && sratio<=5 ...
                        && sum(datalabels==1)>=2 && sum(datalabels==-1)>=2
                    if ncell>=5
                        pcs_all=[];
                        [mdl,pcs,~]=SpRayleigh(dataVector2(:,:,idx_train_trial),datalabels(idx_train_trial),npcs,'fisher','nozscore');
                        pcs_all(:,:,idx_train_trial)=pcs;
                        pcs_all(:,:,idx_test_trial)=TransformFea(dataVector2(:,:,idx_test_trial),mdl);

                        session_cell_firing(:,:,idx_test_trial)=dataVector2(:,:,idx_test_trial);

                        session_pcs(:,:,idx_test_trial)=pcs_all(:,:,idx_test_trial);
                        session_stimulusID{idx_test_trial}=stimulusID_now{idx_test_trial};
                        dataSamples=pcs_all;
                    end
                else
                    continue;
                end
            end
            result_file_sfx=[result_file_sfx {['_pc' num2str(npcs)]}];

            vector_all=[];
            vector_label_all=[];
            timestamp_all=[];
            vector_pre_all=[];
            vector_pre_label_all=[];
            timestamp_pre_all=[];
            vector_all_trial_index=[]; 
            vector_pre_trial_index=[]; 
            
            for ss=1:length(stimulusID)
                
                raster_index=[];
                
                for jj=1:TrialLength
                    if strcmp(stimulusID_now{1,jj},stimulusID{ss})
                        raster_index = [ raster_index jj];
                    end
                end

                vector_now = dataSamples(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2,raster_index);
                vector_pre = dataSamples(:,1:datalengthMili/bin_length,raster_index);

                if isempty(vector_now)
                    continue;
                end

                [c_a,c_b,c_c] = size(vector_now);
                
                vector_combine = reshape(vector_now,c_a,c_b*c_c);                
                vector_label = ones(c_a,c_b*c_c)*ss;
                
                vector_all=[vector_all vector_combine];
                vector_label_all=[vector_label_all vector_label(1,:)];
                timestamp_all = [timestamp_all,repmat(c_b*bin_length+[1:c_b]*bin_length,1,c_c)];

                vector_pre_all=[vector_pre_all reshape(vector_pre,c_a,numel(vector_pre)/c_a)];
                vector_pre_label_all=[vector_pre_label_all ones(1,numel(vector_pre)/c_a)*ss];
                timestamp_pre_all = [timestamp_pre_all,repmat([1:size(vector_pre,2)]*bin_length,1,c_c)];

                vector_all_trial_index=[vector_all_trial_index,reshape(ones(c_b,1)*raster_index,1,c_b*c_c)];
                vector_pre_trial_index=[vector_pre_trial_index,reshape(ones(size(vector_pre,2),1)*raster_index,1,numel(vector_pre)/c_a)];
            end
            
            % decoding
            features=vector_all';
            labels=vector_label_all';
            label_classes=unique(labels);            

            test_index=vector_all_trial_index==idx_test_trial;
            train_index=~test_index;
            
            test_features_pre=vector_pre_all(:,test_index)';
            test_labels_pre=vector_pre_label_all(test_index)';
            test_timetamp_pre=timestamp_pre_all(test_index)';
            train_features_pre=vector_pre_all(:,train_index)';
            train_labels_pre=vector_pre_label_all(train_index)';
            train_timetamp_pre=timestamp_pre_all(train_index)';

            % # use -1 ~ 0s as baseline
            result_file_sfx=[result_file_sfx {'_bl23'}];
            test_features_pre=test_features_pre(test_timetamp_pre>=2000,:);
            test_labels_pre=test_labels_pre(test_timetamp_pre>=2000,:);
            test_timetamp_pre=test_timetamp_pre(test_timetamp_pre>=2000,:);
            train_features_pre=train_features_pre(train_timetamp_pre>=2000,:);
            train_labels_pre=train_labels_pre(train_timetamp_pre>=2000,:);
            train_timetamp_pre=train_timetamp_pre(train_timetamp_pre>=2000,:);

            test_features_during=features(test_index,:);
            test_labels_during=labels(test_index,:);
            test_timestamp_during=timestamp_all(test_index)';
            train_features_during=features(train_index,:);
            train_labels_during=labels(train_index,:);
            train_timestamp_during=timestamp_all(train_index)';

            test_features=cat(1,test_features_pre,test_features_during);
            test_labels=cat(1,test_labels_pre,test_labels_during);
            test_timetamp=cat(1,test_timetamp_pre,test_timestamp_during);

            % # Decoding: Naive Bayes classifier (temporally accumulated & time-dependent modeling)
            % p(s|f(t0),f(t1),...,f(tm)) ¡Ø p(s) ¦° p(f(ti)|s)
            %                            =  p(s) ¦° p(f1(ti)|s)p(f2(ti)|s)...p(fn(ti)|s)
            %                            =  p(s) ¦° ¦° p(fi(tj)|s)
            % s--stimulus, f--feature
            % sc_accum(t0 ~ tm) = log(p(s,f(t0),...,f(tm)))
            %                   = log(p(s)) + ¦² ¦² log(p(fi(tj)|s)
            method_sfx='_nbayes_accum_timemodel';
            train_features=cat(1,train_features_pre,train_features_during);
            train_labels=cat(1,train_labels_pre,train_labels_during);
            train_timetamp=cat(1,train_timetamp_pre,train_timestamp_during);
            [test_timetamp,idx_time_sorted]=sort(test_timetamp,'ascend');
            test_features=test_features(idx_time_sorted);
            test_labels=test_labels(idx_time_sorted);
            sc_accum=[];
            for iclass=1:length(label_classes)
                m=[]; % ¦Ó x n
                v=[];
                for its=1:length(test_timetamp)
                    fea=train_features(train_labels==label_classes(iclass)&...
                        train_timetamp==test_timetamp(its),:);
                    m(its,:)=mean(fea,1);
                    v(its,:)=var(fea,0,1);
                end
                p_s=stimulus_prior(iclass);  

                scc=[];
                for its=1:size(test_timetamp,1)
                    x=test_features(its,:);
                    n=size(x,1);
                    p_fs=exp(-(x-repmat(m,n,1)).^2./repmat(v,n,1)/2)./sqrt(2*pi*repmat(v,n,1));
                    scc(its)=sum(log(p_fs(:)),'omitnan');
                end
                sc_accum(:,iclass)=log(p_s)+cumsum(scc');
            end
            [~,s_max]=max(sc_accum,[],2);
            judge_labels=label_classes(s_max);
            sc=sc_accum;
            
            session_test_labels=[session_test_labels;test_labels];
            session_judge_labels=[session_judge_labels;judge_labels];
            session_judge_scores=[session_judge_scores;sc];

            session_test_labels_timecurve=[session_test_labels_timecurve;test_labels];
            session_judge_labels_timecurve=[session_judge_labels_timecurve;judge_labels];
            session_judge_scores_timecurve=[session_judge_scores_timecurve;sc];
            session_ts_timecurve=[session_ts_timecurve,test_timetamp];

            result_file_sfx=unique([result_file_sfx {method_sfx}]);
            end

            label_all_test{end+1,1}=session_test_labels;
            label_all_judge{end+1,1}=session_judge_labels;           
            ssname_all{end+1,1}=name_now;
            label_timecurve_test{end+1,1}=session_test_labels_timecurve;
            label_timecurve_judge{end+1,1}=session_judge_labels_timecurve;
            score_timecurve_judge{end+1,1}=session_judge_scores_timecurve;
            timestamp_timecurve{end+1,1}=session_ts_timecurve;
            pcs_all_sessions{end+1,1}=session_pcs;
            cell_firing_all_sessions{end+1,1}=session_cell_firing;            
            stimulusID_all_sessions{end+1,1}=session_stimulusID; 
            
        end
    end
 end
 
valid_PFC_fileNo=length(label_all_test);
 
decode_percent_all = zeros(length(stimulusID),length(stimulusID),valid_PFC_fileNo);
decode_class_all = zeros(length(stimulusID),length(stimulusID),valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    label_test_now = label_all_test{nn,1};
    label_judge_now = label_all_judge{nn,1};
    [Ntest,Nrandom] = size(label_test_now);
    label_test_now2 = reshape(label_test_now,Ntest*Nrandom,1);
    label_judge_now2 = reshape(label_judge_now,Ntest*Nrandom,1);
    for kk=1:length(stimulusID)
         label_test_now11 = label_test_now2(label_test_now2 == kk);
          label_judge_now11 = label_judge_now2(label_test_now2 == kk);
         for jj=1:length(stimulusID)
             decode_percent_all(kk,jj,nn) = length( find( label_judge_now11 ==jj) ) /length(label_test_now11);
             decode_class_all(kk,jj,nn) = length( find( label_judge_now11 ==jj) );
         end
    end
end

validIdx=true(1,valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    if sum(sum(isnan(decode_percent_all(:,:,nn)))) ...
            || mean(decode_percent_all(:,1,nn)==1)==1
        validIdx(nn)=false;
    end
end

decode_percent_mean=mean(decode_percent_all(:,:,validIdx),3);

% overall accuracy
decode_percent_overall = zeros(1,valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    label_test_now = label_all_test{nn,1};
    label_judge_now = label_all_judge{nn,1};    
    decode_percent_overall(nn) = mean(label_test_now==label_judge_now);
end

% time curve of accuracy
decode_percent_overall_timecurve = zeros(valid_PFC_fileNo,1);
decode_percent_overall_timecurve_ts = zeros(valid_PFC_fileNo,1);
for nn=1:valid_PFC_fileNo    
    label_test_now = label_timecurve_test{nn,1};
    label_judge_now = label_timecurve_judge{nn,1};
    score_judge_now = score_timecurve_judge{nn,1};
    ts = timestamp_timecurve{nn};
    uts=unique(ts);
    for tt=1:length(uts)
        idxcal=ts==uts(tt);
        ltest=label_test_now(idxcal);
        ljudge=label_judge_now(idxcal);
        sc=score_judge_now(idxcal,:);

        idxvalid=sum(isnan(sc),2)==0;
        ltest=ltest(idxvalid);
        ljudge=ljudge(idxvalid);
        if ~isempty(ljudge)
            ucat=unique(ltest);
            acccat=zeros(1,length(ucat));
            for icat=1:length(ucat)
                truecat=ucat(icat);
                acccat(icat)=mean(ljudge(ltest==truecat)==truecat);
            end
            decode_percent_overall_timecurve(nn,tt)=mean(acccat);

        else
            decode_percent_overall_timecurve(nn,tt)=nan;
        end
    end
    decode_percent_overall_timecurve_ts(nn,1:length(uts))=uts';
end

decode_percent_all_compared{idim,idir,ibin}=decode_percent_all(:,:,validIdx);
decode_percent_overall_compared{idim,idir,ibin}=decode_percent_overall;
decode_percent_overall_valid_session_compared{idim,idir,ibin}=decode_percent_overall(validIdx);
decode_ssname{idim,idir,ibin}=ssname_all;
decode_valid_ssname{idim,idir,ibin}=ssname_all(validIdx);
decode_percent_overall_timecurve_compared{idim,idir,ibin}=decode_percent_overall_timecurve(validIdx,:);
decode_percent_overall_timecurve_ts_compared{idim,idir,ibin}=decode_percent_overall_timecurve_ts(validIdx,:);
decode_pcs_compared{idim,idir,ibin}=pcs_all_sessions(validIdx);
decode_cell_firing_compared{idim,idir,ibin}=cell_firing_all_sessions(validIdx);
decode_stimulusID_compared{idim,idir,ibin}=stimulusID_all_sessions(validIdx);

end
end
end

dslist={};
result_file_sfx=unique(result_file_sfx);
for igrp=1:length(stimulusIDGroup)
    stim_suffix=['_' strjoin(stimulusIDGroup{igrp},'_')];
    dslist{igrp}=['decode_res' stim_suffix strcat(result_file_sfx{:}) '.mat'];
end

save([pwd '\data\' dslist{istimgrp}]);
disp(['Current results saved to ' dslist{istimgrp}]);
end
end

end


%% Section II: plot using saved results

dslist_plot={...
    'decode_res_socialA_novel_bl23_high_var_cell_nbayes_accum_timemodel_pc1_triallim8.mat';...
    'decode_res_socialA_novel_bl23_high_var_cell_nbayes_accum_timemodel_pc1_shuffle_triallim8.mat';...
    };
dslist_title={'','shuffled'};

bin_compared=50;
dim_compared=7;
hfig_time_curve_all=figure;
for ids=1:length(dslist_plot)
    load([pwd '\data\' dslist_plot{ids}]);

% check sample availability
for idir=1:length(dir_cate)
    v=decode_valid_ssname(:,idir,:);
    sz=reshape(arrayfun(@(x) length(x{1}), v(:)),...
        length(pop_cell_dim_compared),length(bin_length_compared));
    disp(['Number of samples for # ' dir_cate{idir} ':']);
    disp('( Dimension x Bin length )');
    disp(sz);
    idx=find(sz==0);
    if ~isempty(idx)
        disp(['Lack of samples for ' dir_cate{idir} ...
            ' at indicies: ']);
        error(num2str(idx'));
    end
end

stimulusID_abbr=replace(stimulusID,{'socialA','socialB','novel','empty'},{'S1','S2','N','E'});

% plot time curve of decoding accuracies
dir_color={'b','r'};
ibin=find(bin_length_compared==bin_compared);
for idim=1:length(pop_cell_dim_compared)    
    for idir=1:length(dir_names_compared)       
        t=decode_percent_overall_timecurve_ts_compared{idim,idir,ibin}(1,:)/1000-3;
        v=decode_percent_overall_timecurve_compared{idim,idir,ibin};
        m=mean(v,1,'omitnan');
        s=std(v,0,1,'omitnan')/sqrt(size(v,1));
        if idim==dim_compared
            figure(hfig_time_curve_all);
            if length(dslist_plot)<=3; subplot(1,length(dslist_plot),ids);
            else; subplot(ceil(length(dslist_plot)/3),3,ids); 
            end
            htcp_all(idir)=plot(t,m,dir_color{idir});
            hold on;
            plotshaded(t,[m-s;m+s],dir_color{idir});
            title([stimulusID_abbr{1} ' vs. ' stimulusID_abbr{2} '(' dslist_title{ids} ')']);
            xlabel('time(s)');
            ylabel('accuracy');
            ylim([0,1]);
            if ids==length(dslist_plot) && idir==length(dir_names_compared)
                legend(htcp_all,dir_cate);
            end
        end
    end    
end
end
end


%% Subfunction #1: SpRayleigh
% Calculate projection vectors of spatial filtering based on Rayleigh
% quotient for multi-channel features.

function [mdl,pcs,ffeatures]=SpRayleigh(features,labels,varargin)
% [mdl,pcs,ffeatures]=SpRayleigh(features,labels,npcs,method,zopt)

npcs=size(features,1);
if length(varargin)>=1 && ~isempty(varargin{1}) ...
        && npcs>varargin{1}
    npcs=varargin{1};
end
method='fisher';
if length(varargin)>=2
    method=varargin{2};
end
zopt='zscore';
if length(varargin)>=3
    zopt=varargin{3};
end

switch zopt
    case 'zscore'
        [features,mu,sigma]=zscore(features,[],3); % nChannel x tFeature x nSample
    case 'chlzscore'
        [nc,nt,ns]=size(features);
        [features,mu,sigma]=zscore(reshape(features,nc,nt*ns),[],2); % nChannel x tFeature x nSample
        features=reshape(features,nc,nt,ns);
        mu=repmat(mu,1,nt);
        sigma=repmat(sigma,1,nt);
    case 'nozscore'
        mu=0;
        sigma=1;
    otherwise
        error('Unsupported option!');
end

switch method        
    case 'fisher' % # Fisher's rule (relevance space projection)
        S=permute(features(:,:,labels==1),[2,1,3]);% tFeature X nChannel x nSample
        N=permute(features(:,:,labels==-1),[2,1,3]);
        A=mean(S,3); % tFeature X nChannel
        B=mean(N,3); % tFeature X nChannel
        Swa=zeros(size(S,2),size(S,2));
        Swb=zeros(size(S,2),size(S,2));
        for isample=1:size(S,3)
            Swa=Swa+(S(:,:,isample)-A)'*(S(:,:,isample)-A);
        end
        Swa=Swa/size(S,3);
        for isample=1:size(N,3)
            Swb=Swb+(N(:,:,isample)-B)'*(N(:,:,isample)-B);
        end
        Swb=Swb/size(N,3);
        Sw=(Swa+Swb)/2;
        Sb=(A-B)'*(A-B);
        try
            L=chol(Sw);
        catch ME
            disp(['Sw is not symmetric positive definite, and a small perturbation was added to enable the Cholesky decompostition!'])
            p=0.01;
            L=chol((1-p)*Sw+p*eye(size(Sw,1)));
        end
        invL=inv(L);
        [V,D]=eig(invL'*Sb*invL);
        U=invL*V;
        Ua=V'*L;
        [~,I]=sort(diag(D),'descend');
        u=U(:,I);
        iu=Ua(I,:);
    otherwise
        error('Unsupported method!');
end
u=u(:,1:npcs);
iu=iu(1:npcs,:);
mdl.mu=mu;
mdl.sigma=sigma;
mdl.u=u;
mdl.iu=iu;
if nargout>1
    pcs=zeros(npcs,size(features,2),size(features,3));
    for isample=1:size(pcs,3)
        pcs(:,:,isample)=u'*features(:,:,isample);
    end
end
if nargout>2
    ffeatures=zeros(size(features));
    for isample=1:size(pcs,3)
        ffeatures(:,:,isample)=iu'*pcs(:,:,isample).*sigma+mu;
    end
end
end


%% Subfunction #2: TransformFea
% transform features using the linear model MDL

function [pcs,ffea]=TransformFea(fea,mdl)
if isempty(mdl)||~isfield(mdl,'u')||isempty(mdl.u)
    pcs=fea;
    ffea=fea;
    return; 
end

% For vector input of mu and sigma, normalization will be done only along
% the channel dimension, and tFeature x nSamples will be seen as a whole.
if size(mdl.mu,2)==1
    mdl.mu=repmat(mdl.mu,1,size(fea,2));
    mdl.sigma=repmat(mdl.sigma,1,size(fea,2));
end

mdl.sigma(mdl.sigma==0)=1;
fea=(fea-repmat(mdl.mu,1,1,size(fea,3)))./...
    repmat(mdl.sigma,1,1,size(fea,3));

npcs=size(mdl.u,2);
pcs=zeros(npcs,size(fea,2),size(fea,3));
for isample=1:size(pcs,3)
    pcs(:,:,isample)=mdl.u'*fea(:,:,isample);
end

if nargout>1
    ffea=zeros(size(fea));
    for isample=1:size(fea,3)
        ffea(:,:,isample)=mdl.iu'*pcs(:,:,isample).*mdl.sigma+mdl.mu;
    end
end
end


%% Subfunction #3: plotshaded

function varargout = plotshaded(x,y,fstr)
% x: x coordinates
% y: either just one y vector, or 2xN or 3xN matrix of y-data
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');
 
if size(y,1)>size(y,2)
    y=y';
end;
 
if size(y,1)==1 % just plot one line
    plot(x,y,fstr);
end;
 
if size(y,1)==2 %plot shaded area
    px=[x,fliplr(x)]; % make closed patch
    py=[y(1,:), fliplr(y(2,:))];
    patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
end;
 
if size(y,1)==3 % also draw mean
    px=[x,fliplr(x)];
    py=[y(1,:), fliplr(y(3,:))];
    patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
    plot(x,y(2,:),fstr);
end;
 
alpha(.2); % make patch transparent
end

