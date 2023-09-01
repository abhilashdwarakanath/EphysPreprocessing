function [spikes] = getSpikesGeneral(filename,params,ss,es,negorpos,filterType);

% This is based heavily on peakseek.m from MIT. peakseek is FAR faster than
% findpeaks and is much much better and more comprehensive. findpeaks.m
% should be ashamed of itself!
%
% This function also can pick up positive going spikes if at all any turn
% up, like they suprisingly did in our Utah array recordings on some days!
% It corrects for events detected which are less apart than the specified
% refractory period (refer to the global param structure). It also
% specifies a maximum spike height (I've set this quite large in our param
% structure because the "intracellular/juxtacellular" spikes are serious
% fuckin big.
%
% Potential improvements -
%
% 1. Unspecified intervals to pick up minimas, i.e. pick up a first and a
% second local minima instead of chopping the spike into two.
% 2. Fit a laplace or a gaussian or a beta distribution to the mean
% spike-shape and use the shape and scale (i.e. M and S for gaussian)
% parameters as features.
%
% Abhilash D, MPIBC May 2016

%% Read file and filter the timeseries for spikes

x = read_NC5(filename,ss,es);

fprintf('Filtering via an IIR 2nd Order Butterworth Filter ...\n')

switch filterType
    
    case{'butter','Butter','butterworth','Butterworth'}
        params.passband = [600 3000]/15000;
        [b,a] = butter(params.filterorder,params.passband);
        
        samples = filtfilt(b,a,x);
        
    case{'Tolias','tolias','Ecker','ecker','both','Both'}
        
        fprintf('Creating a minimum order FIR Bandpass filter...\n')
        
        ampPass = (1 - db2mag(-params.passripple)) / (1 + db2mag(-params.passripple));
        ampStop = db2mag(-params.stopatten);
        
        fl = [params.stopband(1) params.passband(1)]*(params.fs/2);
        fh = [params.passband(2) params.stopband(2)]*(params.fs/2);
        
        [n, fo, ao, w] = firpmord([fl fh], [0 1 0], [ampStop ampPass ampStop], params.fs);
        b = firpm(n, fo, ao, w);
        
        samples = filtfilt(b,1,x);
        
end

nChunk = 100;

for i = 1:nChunk
    r = randi([ss es-params.fs],1,1);
    noiseforms = samples(r:r+params.fs);
    thr(i) = median(abs(noiseforms-mean(noiseforms))/0.6745);
end

thr = mean(thr);

minpeakdist = (params.refractory*params.fs)/1000; % This can be changed. This is ideally the refractory period

minpeakh = params.stdmin*thr;
maxpeakh = params.stdmax*thr;

tminus = ceil((params.beforeSpike*params.fs)/1000);
tplus = ceil((params.afterSpike*params.fs)/1000);

%% Detect spikes
fprintf('Detecting spikes...\n')
if size(samples,2)==1
    samples=samples';
end

if negorpos==1 % If this argument is 1, I assume that this channel has not only POS spikes but ALSO NEG spikes
    
    % Collect positive and negative going spikes i.e. % Find all maxima/minima and ties
    poslocs=find(samples(2:end-1)>=samples(1:end-2) & samples(2:end-1)>=samples(3:end))+1;
    neglocs=find(samples(2:end-1)<=samples(1:end-2) & samples(2:end-1)<=samples(3:end))+1;
    
    % Correct for thresholds
    neglocs(samples(neglocs)>=minpeakh)=[];
    neglocs(samples(neglocs)<=maxpeakh)=[];
    
    poslocs(samples(poslocs)<=-1*minpeakh)=[];
    poslocs(samples(poslocs)>=-1*maxpeakh)=[];
    
    locs = unique(sort([neglocs poslocs]));
    
else
    
    locs=find(samples(2:end-1)<=samples(1:end-2) & samples(2:end-1)<=samples(3:end))+1;
    
    % Correct for spikes that don't really match the chosen thresholds
    locs(samples(locs)>=minpeakh)=[];
    locs(samples(locs)<=maxpeakh)=[];
    
end

% Correct and align spikes smaller than the refractory period

fprintf('Correcting and aligning spikes...\n')
if minpeakdist>1
    while 1
        
        del=diff(locs)<minpeakdist;
        
        if ~any(del), break; end % There must be a more efficient way of doing this instead of using break?
        
        pks=samples(locs);
        
        [~,mins]=min([pks(del) ; pks([false del])]);
        
        deln=find(del);
        
        deln=[deln(mins==1) deln(mins==2)+1];
        
        locs(deln)=[];
        
    end
end

%% Collect waveforms

% Maybe we should integrate getting the positive and negative spikeforms
% separately?

fprintf('Collecting waveforms...\n');

spikeForms = zeros(length(locs),tminus+tplus);

switch filterType
    
    case{'butter','Butter','butterworth','Butterworth','Tolias','tolias','Ecker','ecker'}
        
        for i = 1:length(locs)
            if locs(i) > 2*tminus && locs(i) < locs(end)-2*tminus
                spikeForms(i,:) = [samples(locs(i)-tminus+1:locs(i)+tplus)];
            end
            
        end
        
    case{'Both','both'}
        
        [b,a] = butter(params.filterorder,params.passband);
        
        samples = filtfilt(b,a,x);
        
        for i = 1:length(locs)
            if locs(i) > 2*tminus && locs(i) < locs(end)-2*tminus
                spikeForms(i,:) = [samples(locs(i)-tminus+1:locs(i)+tplus)];
            end
            
        end
end

% If there any 0 rows, i.e. any spike before tminus, remove it

spikes.waveforms = spikeForms;
spikes.times = locs;

end
