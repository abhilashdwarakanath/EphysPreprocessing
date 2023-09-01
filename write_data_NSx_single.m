function [lengthPFC,nParts]=write_data_NSx_single(params,filenames,whichFiles,directories,nevfiles)

% An alternative to the parse_data_NSx from the Rodrigo lab which can deal
% with synchronised NSPs and handle large files.
%
% Abhi. MPIBC. August 2016

%% Read the two file headers and get the timestamps

NSP1 = openNSx(['./' filenames{whichFiles(1)}], 'noread','report');

nchans = NSP1.MetaTags.ChannelCount;
sr = NSP1.MetaTags.SamplingFreq;

%% Get recording lengths and save timestamps file
pauses = 0;
if length(NSP1.MetaTags.DataPoints)==2
    lengthPFC = NSP1.MetaTags.DataPoints(2);
    nParts.n = length(NSP1.MetaTags.DataPoints);
    nParts.first = NSP1.MetaTags.DataPoints(2);
    
elseif length(NSP1.MetaTags.DataPoints)==3
    lengthPFC = NSP1.MetaTags.DataPoints(2)+NSP1.MetaTags.DataPoints(3);
    nParts.n = length(NSP1.MetaTags.DataPoints);
    nParts.first = NSP1.MetaTags.DataPoints(2);
    nParts.second = NSP1.MetaTags.DataPoints(3);
elseif length(NSP1.MetaTags.DataPoints)==1
    lengthPFC = NSP1.MetaTags.DataPoints(1);
    nParts.n = length(NSP1.MetaTags.DataPoints);
elseif length(NSP1.MetaTags.DataPoints)>3
    pauses = 1;
end

if pauses == 0
    
    lts = lengthPFC;
    
    if ~exist('NSX_TimeStamps.mat', 'file') == 2
        if sr == 1000
             TimeStamps=linspace(0,(lengthPFC-1)*1e6/sr,ceil(lengthPFC));%TimeStamps in microsec, with 0 corresponding to the first sample
            cd(directories.taskdirPFC)
            save('NSX_TimeStamps','TimeStamps','lengthPFC','lts','nchans','sr','-v7.3');
        else
             TimeStamps=linspace(0,(lengthPFC-1)*1e6/sr,lengthPFC);%TimeStamps in microsec, with 0 corresponding to the first sample
            cd(directories.PFC)
            save('NSX_TimeStamps','TimeStamps','lengthPFC','lts','nchans','sr','-v7.3');
        end
    end
    
    % Correct the timestamps in the NEV file
    cd(directories.recording)
    NEV = openNEV(['./' nevfiles{1}]);
    
    if nParts.n==1 || nParts.n==2
        save([nevfiles{1} '.mat'],'NEV')
    elseif nParts.n==3
        timeTOAdd = NSP1.MetaTags.DataPoints(nParts.n-1);
        
        % Fix pulses timings
        %     pulseTimings = double(NEV.Data.SerialDigitalIO.TimeStamp);
        %     IPI = diff(pulseTimings);
        %     [idx idx] = find(IPI<=0);
        %     if ~isempty(idx)
        %         pulsesToCorrect = idx+1:length(pulseTimings);
        %         NEV.Data.SerialDigitalIO.TimeStamp(pulsesToCorrect) = NEV.Data.SerialDigitalIO.TimeStamp(pulsesToCorrect)+timeTOAdd;
        %     end
        %
        %     % Fix comments timings
        %     commentsTimings = double(NEV.Data.Comments.TimeStamp);
        %     ICI = diff(commentsTimings);
        %     ICI = ICI./params.fs;
        %     [idx idx] = find(ICI<=5);
        %     if ~isempty(idx)
        %         commentsToCorrect = idx+1:length(commentsTimings);
        %         NEV.Data.Comments.TimeStamp(commentsToCorrect) = NEV.Data.Comments.TimeStamp(commentsToCorrect)+timeTOAdd;
        %     end
        PFCevents = NEV;
        PFCevents.timeCorr = timeTOAdd;
        save([nevfiles{1}(1:end-4) '.mat'],'PFCevents')
    end
    
    %% Read the files channel by channel and write them into the NC5 files
    
    for i = 1:NSP1.MetaTags.ChannelCount
        
        tic;
        
        cd(directories.recording)
        fprintf('Parsing Channel %d from the PFC array....\n ',i)
        datNSP1 = openNSx(['./' filenames{whichFiles(1)}],'read','report',['c:' num2str(i)]);
        
        if iscell(datNSP1.Data)==1 && length(datNSP1.Data)==2 % A condition needs to be put in here. Anton's data has 2 cells but Hayo's has 3 and Hayo's 3rd cell for the ns2 is what matters
            datNSP1 = datNSP1.Data{2};
        elseif iscell(datNSP1.Data)==1 && length(datNSP1.Data)==3
            %datNSP1 = [datNSP1.Data{2} datNSP1.Data{3}]; %% Add a condition
            %for Hayo. THIS NEEDS TO BE MANUALLY CHANGED BECAUSE HAYO's PFC
            %HAS TWO RESETS DUE TO SYNC (BECAUSE PPC IS THE MASTER!!!)
            %The parsing until now has been fine because we only used the
            %PFC for the paper, and I have a check delay function, but it
            %would be nice to fix it right here at the Parsing itself. It
            %especially matters for the ns2 file that records the eye
            %movements.
            datNSP1 = [datNSP1.Data{2} datNSP1.Data{3}];% If there was a break. This is weird.
        else
            datNSP1 = datNSP1.Data;
        end
        cd(directories.PFC)
        outfile_handle = fopen(['NSX' num2str(NSP1.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PFC NC5 file for channel %d\n',i)
        fwrite(outfile_handle,datNSP1,'int16');
        fclose(outfile_handle);
        
        clear datNSP1;
        
        el = toc;
        
        fprintf('Channel %d of %d processed. Elapsed time is %d seconds',i,params.chans,el)
        
    end
    
else
    fprintf('This recording has pauses! Checking for data consistency and parsing if consistent...\n')
    
    % Calculate delays between different segments
    
    if sr == 30000
        for i = 3:length(NSP1.MetaTags.Timestamp)
            pauseDelays(i) = NSP1.MetaTags.Timestamp(i) - (NSP1.MetaTags.Timestamp(i-1) + NSP1.MetaTags.DataPoints(i-1));
        end
    else
        for i = 3:length(NSP1.MetaTags.Timestamp)
            pauseDelays(i) = ceil(NSP1.MetaTags.Timestamp(i)/30) - (ceil(NSP1.MetaTags.Timestamp(i-1)/30) + NSP1.MetaTags.DataPoints(i-1));
        end
    end
    
    
    for i = 1:NSP1.MetaTags.ChannelCount
        cd(directories.recording)
        fprintf('Parsing Channel %d from the PFC array....\n ',i)
        datNSP1 = openNSx(['./' filenames{whichFiles(1)}],'read','report',['c:' num2str(i)]);
        
        for iPieces = 3:length(datNSP1.Data)+1
            datPiece{iPieces-1} = [datNSP1.Data{iPieces-1} int16(randn(1,pauseDelays(iPieces-1)-1))];
        end
        
        if i == 1
            tp = 0;
            for p = 2:length(datNSP1.Data)+1
                tp = tp+length(datPiece{p-1});
            end
            lengthPFC = tp;
            lts = tp;
            cd(directories.PFC)
            if ~exist('NSX_TimeStamps.mat', 'file') == 1
                TimeStamps=linspace(0,(lengthPFC-1)*1e6/params.fs,lengthPFC);%TimeStamps in microsec, with 0 corresponding to the first sample
                if sr == 1000
                    cd(directories.taskdirPFC)
                    save('NSX_TimeStamps','TimeStamps','lengthPFC','lts','nchans','sr','-v7.3');
                else
                    cd(directories.PFC)
                    save('NSX_TimeStamps','TimeStamps','lengthPFC','lts','nchans','sr','-v7.3');
                end
            end
            
        end
        
        fixedData = [];
        
        for iPieces = 3:length(datNSP1.Data)+1
            fixedData = [fixedData datPiece{iPieces-1}];
        end
        
        cd(directories.PFC)
        outfile_handle = fopen(['NSX' num2str(NSP1.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PFC NC5 file for channel %d\n',i)
        fwrite(outfile_handle,fixedData,'int16');
        fclose(outfile_handle);
        
        el = toc;
        
        fprintf('Channel %d of %d processed from a Paused Recording. Elapsed time is %d seconds',i,params.chans,el)
        nParts = 0;
        
    end
end
cd(directories.recording)