absvelocityThreshold = 0.5;
absaccelerationThreshold = 0.5;
absamplitudeThreshold = 700;
okntracePreprocessed = [];

tmp_dOKNtimeSeries = diff(filtSignal);
        % finding the second derivative
        tmp_d2OKNtimeSeries = diff(tmp_dOKNtimeSeries);
        % insert 1 NaN at the begining of 1 dev
        tmp_dOKNtimeSeries = [NaN tmp_dOKNtimeSeries];
        % insert 2 NaN at the begining of 2 dev
        tmp_d2OKNtimeSeries = [NaN NaN tmp_d2OKNtimeSeries];
        % find points where we pass the thresholds (velocity and acceleration threshold)
        tmpIndex = find(abs(filtSignal) >= absamplitudeThreshold);
        tmpIndex_d = find(abs(tmp_dOKNtimeSeries) >= absvelocityThreshold);
        tmpIndex_d = [tmpIndex_d find(isnan(tmp_dOKNtimeSeries))]; % get rid of the NaNs 
        tmpIndex_d = [tmpIndex_d find(isinf(tmp_dOKNtimeSeries))]; % get rid of the INFs
        tmpIndex_d2 = find(abs(tmp_d2OKNtimeSeries) >= absaccelerationThreshold);
        tmpIndex_d2 = [tmpIndex_d2 find(isnan(tmp_d2OKNtimeSeries))]; % get rid of the NaNs
        tmpIndex_d2 = [tmpIndex_d2 find(isinf(tmp_d2OKNtimeSeries))]; % get rid of the INFs
        % find the indices that cross all thresholds
        indicesAmpVel = intersect(tmpIndex,tmpIndex_d); %indices where amplitude and velocity are both beyond threshold
        indicesAmpAcc = intersect(tmpIndex,tmpIndex_d2); %indices where amplitude and acceleration are both beyond threshold
        indicesVelAcc = intersect(tmpIndex_d,tmpIndex_d2); %indices where acceleration and velocity are both beyond threshold
        indices = intersect(indicesAmpVel,tmpIndex_d2); % for when we want the points where amplitude, vlocity and acceleration all pass the threshold at the same time
        filtSignal(indices) = 0;
        filtSignal(indicesAmpVel) = 0;
        filtSignal(indicesAmpAcc) = 0;