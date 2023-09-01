function [lengthPFC,lengthPPC]=write_data_NSx(params,filenames,whichFiles,directories)

% An alternative to the parse_data_NSx from the Rodrigo lab which can deal
% with synchronised NSPs and handle large files.
%
% Abhi. MPIBC. August 2016

%% Read the two file headers and get the timestamps

NSP1 = openNSx(['./' filenames{whichFiles(1)}], 'noread','report');
NSP2 = openNSx(['./' filenames{whichFiles(2)}], 'noread','report');

%% Get sync delays and align the data together

if length(NSP1.MetaTags.Timestamp) == 2 % In a clean ideal recording, there should only be two data cells for each NSP. Cell 1 is before timer concurrent reset so we don't need that. Cell 2 is with the real data. 
    
    tsPFC = NSP1.MetaTags.Timestamp(2); % Lagging NSP for us
    tsPPC = NSP2.MetaTags.Timestamp(2); % Leading NSP for us
    
    syncShift = abs(tsPFC-tsPPC);
    
    nchans = NSP1.MetaTags.ChannelCount;
    sr = NSP1.MetaTags.SamplingFreq;
    
    %% Get recording lengths and save timestamps file
    
    if tsPFC > tsPPC % check which one is lagging or leading (generalised)
        lengthPFC = NSP1.MetaTags.DataPoints(2);
        lengthPPC = NSP2.MetaTags.DataPoints(2)-syncShift; % Correct for it
    else
        lengthPPC = NSP2.MetaTags.DataPoints(2);
        lengthPFC = NSP1.MetaTags.DataPoints(2)-syncShift;
    end
    
    % Now generate the timestamps
	
    TimeStampsPFC=linspace(0,(lengthPFC-1)*1e6/sr,lengthPFC);%TimeStamps in microsec, with 0 corresponding to the first sample
    lts = lengthPFC;
    TimeStamps = TimeStampsPFC;
    cd(directories.PFC)
    save('NSX_TimeStamps','TimeStamps','lengthPFC','nchans','lts','sr','-v7.3');
    
    TimeStampsPPC=linspace(0,(lengthPPC-1)*1e6/sr,lengthPPC);%TimeStamps in microsec, with 0 corresponding to the first sample
    lts = lengthPPC;
    TimeStamps = TimeStampsPPC;
    cd(directories.PPC)
    save('NSX_TimeStamps','TimeStamps','lengthPPC','nchans','lts','sr','-v7.3');
    
    clear TimeStampsPFC; clear TimeStampsPPC; clear TimeStamps;
    
    %% Read the files channel by channel and write them into the NC5 files
    
    for i=1:128 % Replace this with the number of channels
        
        tic;
        
        cd(directories.recording)
        fprintf('Parsing Channel %d from the PFC array....\n ',i)
        datNSP1 = openNSx(['./' filenames{whichFiles(1)}],'read','report',['c:' num2str(i)]);
        fprintf('Parsing Channel %d from the PPC array....\n',i)
        datNSP2 = openNSx(['./' filenames{whichFiles(2)}],'read','report',['c:' num2str(i)]);
        
        datNSP1 = datNSP1.Data{2};
        datNSP2 = datNSP2.Data{2};
        
        remInds = 1:syncShift;
        if tsPFC > tsPPC
            datNSP2(remInds) = [];
        else
            datNSP1(remInds) = [];
        end
        
        cd(directories.PFC)
        outfile_handle = fopen(['NSX' num2str(NSP1.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PFC NC5 file for channel %d\n',i)
        fwrite(outfile_handle,datNSP1,'int16');
        fclose(outfile_handle);
        cd(directories.PPC)
        outfile_handle2 = fopen(['NSX' num2str(NSP2.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PPC NC5 file for channel %d\n',i)
        fwrite(outfile_handle2,datNSP2,'int16');
        fclose(outfile_handle2);
        
        %clear datNSP1;
        %clear datNSP2;
        
        el = toc;
        
        fprintf('Channel %d of %d of both arrays processed. Elapsed time is %d seconds\n',i,params.chans,el)
        
    end
    
elseif length(NSP1.MetaTags.Timestamp) == 3 % Dirty hack for our experiments with Hayo since Vishal waits for the timer to reset TWICE before starting the other NSP. So PFC has 2 resets, therefore the 3rd cell is with the real data
    
    tsPFC = NSP1.MetaTags.Timestamp(3);
    tsPPC = NSP2.MetaTags.Timestamp(2);
    
    syncShift = abs(tsPFC-tsPPC);
    
    nchans = NSP1.MetaTags.ChannelCount;
    sr = NSP1.MetaTags.SamplingFreq;
    
    %% Get recording lengths and save timestamps file
    
    if tsPFC > tsPPC
        lengthPFC = NSP1.MetaTags.DataPoints(3);
        lengthPPC = NSP2.MetaTags.DataPoints(2)-syncShift; % Correct the delay
    else
        lengthPPC = NSP2.MetaTags.DataPoints(2);
        lengthPFC = NSP1.MetaTags.DataPoints(3)-syncShift; % Correct the delay
    end
    
    % Now generate the timestamps
	
    TimeStampsPFC=linspace(0,(lengthPFC-1)*1e6/sr,lengthPFC);%TimeStamps in microsec, with 0 corresponding to the first sample
    lts = lengthPFC;
    TimeStamps = TimeStampsPFC;
    cd(directories.PFC)
    save('NSX_TimeStamps','TimeStamps','lengthPFC','nchans','lts','sr','-v7.3');
    
    TimeStampsPPC=linspace(0,(lengthPPC-1)*1e6/sr,lengthPPC);%TimeStamps in microsec, with 0 corresponding to the first sample
    lts = lengthPPC;
    TimeStamps = TimeStampsPPC;
    cd(directories.PPC)
    save('NSX_TimeStamps','TimeStamps','lengthPPC','nchans','lts','sr','-v7.3');
    
    clear TimeStampsPFC; clear TimeStampsPPC; clear TimeStamps;
    
    % Read the files channel by channel and write them into the NC5 files
    
    for i=1:128 % Replace this number with the number of electrodes
        
        tic;
        
        cd(directories.recording)
        fprintf('Parsing Channel %d from the PFC array....\n ',i)
        datNSP1 = openNSx(['./' filenames{whichFiles(1)}],'read','report',['c:' num2str(i)]);
        fprintf('Parsing Channel %d from the PPC array....\n',i)
        datNSP2 = openNSx(['./' filenames{whichFiles(2)}],'read','report',['c:' num2str(i)]);
        
        datNSP1 = datNSP1.Data{3};
        datNSP2 = datNSP2.Data{2};
        
        remInds = 1:syncShift;
        if tsPFC > tsPPC
            datNSP2(remInds) = [];
        else
            datNSP1(remInds) = [];
        end
        
        cd(directories.PFC)
        outfile_handle = fopen(['NSX' num2str(NSP1.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PFC NC5 file for channel %d\n',i)
        fwrite(outfile_handle,datNSP1,'int16');
        fclose(outfile_handle);
        cd(directories.PPC)
        outfile_handle2 = fopen(['NSX' num2str(NSP2.MetaTags.ChannelID(i)) '.NC5'],'w');
        fprintf('Writing the PPC NC5 file for channel %d\n',i)
        fwrite(outfile_handle2,datNSP2,'int16');
        fclose(outfile_handle2);
        
        %clear datNSP1;
        %clear datNSP2;
        
        el = toc;
        
        fprintf('Channel %d of %d of both arrays processed. Elapsed time is %d seconds\n',i,params.chans,el)
        
    end
    
elseif length(NSP1.MetaTags.Timestamp) > 3 % Fucked up recording. Pauses most likely. 
    
    sprintf('\nThis recording contains pauses! Use the appropriate parsing function\n!')
    lengthPFC = []; lengthPPC = [];
    
elseif length(NSP1.MetaTags.Timestamp) < 2 % Use with Multicentral. If this is true, then this parsing function defaults to Rodrigo's parsing function
    
    sprintf('\nThis data comes from a weirdly synchronised recording. Parsing using the old method...\n')
    [lengthPFC,sr]=parse_data_NSx(filenames{whichFiles(1)},24,directories.PFC);
    [lengthPPC,sr]=parse_data_NSx(filenames{whichFiles(2)},24,directories.PPC);
    
end

cd(directories.recording)