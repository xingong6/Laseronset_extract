# Laseronset_extract
%% JF_loadExperiment

% you need in workspace:
% animal
% experiment
% date
% recording
% site
% load_parts
% verbose

myPaths;

%% Display progress or not
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

%% Define what to load

% Site (multiple probes) is optional
if ~exist('site', 'var') || isempty(site)
    site = [];
end

% Recording (multiple recordings on one probe) is optional
if ~exist('recording', 'var') || isempty(recording)
    recording = [];
end

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts', 'var') || isempty(load_parts)
    load_parts.cam = true;
    load_parts.imaging = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts, 'cam')
        load_parts.cam = false;
    end
    if ~isfield(load_parts, 'imaging')
        load_parts.imaging = false;
    end
    if ~isfield(load_parts, 'ephys')
        load_parts.ephys = false;
    end
end

%% Load timeline and associated inputs

[timeline_filename, timeline_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'timeline');
if ~timeline_exists
    error([animal, ' ', day, ' ', num2str(experiment), ': no timeline']);
end

if timeline_exists
    if verbose
        disp('Loading timeline...');
    end

    load(timeline_filename);

        % Get laser times
    laser_name = 'laserCmd';
    timeline_laser_idx = strcmp({Timeline.hw.inputs.name}, laser_name);
    laser_trace = Timeline.rawDAQData(:, timeline_laser_idx) ;
    laser_trace = Timeline.rawDAQData(:, timeline_laser_idx) > 0.05;  % arbitrary threshold 0.05 is too small, have noise going above it leading to misalignment, set as 0.15 to get the laser trace into a binary trace. >0.15 becomes 1, <0.15 becomes 0. Need to change this number if using lower light intensities. You need to find the right threshold to get the onset accurate even when the light is ramping. % laser_trace get a logical matrix according to whether each element are above 0.15 or not.

    % get rid of those wrong onsets later with diff_onidx, diff_offidx.

        % get the sequence/block times
        diff_onidx=[0;diff(laser_flip_onidx)];
        diff_offidx=[diff(laser_flip_offidx);0];
        laser_block_flip_start = [laser_flip_onidx(1);laser_flip_onidx(diff_onidx > 300)];
        laser_block_flip_end=[laser_flip_offidx(diff_offidx > 300);laser_flip_offidx(end)];

%     JF old lines with bugs
%     laser_flip = [find(~laser_trace(1:end-1) & ...
%         laser_trace(2:end)); find(laser_trace(1:end-1) & ...
%         ~laser_trace(2:end))];  
% get laser onset and offset timepoints % I = find(X) returns the linear indices corresponding to the nonzero entries of the array X. %
% I suspect this is wrong because in laser_flip, laser onset timepoints are
% before laser offset, which is,
% laseronset1,laseronset2,...laseronsetlast,laseroffset1,laseroffset2,...laseroffsetlast;
% and this sequence is kept in laser_flip_times, better to separate onset
% and offset.
%         laser_flip_times = Timeline.rawDAQTimestamps(laser_flip)';

     laser_flip_onidx = find(~laser_trace(1:end-1)&laser_trace(2:end));
     laser_flip_offidx= find(laser_trace(1:end-1)&~laser_trace(2:end));
     laser_flip_on_times = Timeline.rawDAQTimestamps(laser_flip_onidx);
     laser_flip_off_times = Timeline.rawDAQTimestamps(laser_flip_offidx);

% plot the laser traces.
                figure();
                clf;
                valu = 3254000;
                title('Laser');
                hold on;
                plot(Timeline.rawDAQData(1:valu, timeline_laser_idx))
                plot(-laser_trace(1:valu))
                scatter(laser_flip(find(laser_flip <= valu)), ones(size(find(laser_flip <= valu), 1), 1)-1)
                plot(medfilt1(Timeline.rawDAQData(1:valu, ...
                    timeline_laser_idx), 3))
                scatter(laser_seq_flip,ones(size(laser_seq_flip,1),1)*0.05)


    % Get camera times
    cam_name = 'pcoExposure';
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);

    cam_expose_starts = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1, timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end, timeline_cam_idx) > 2)+1);
    cam_expose_stops = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1, timeline_cam_idx) >= 2 & ...
        Timeline.rawDAQData(2:end, timeline_cam_idx) < 2)+1);

    cam_time = cam_expose_starts;
    cam_expose_times = cam_expose_stops - cam_expose_starts;

    % Get acqLive signal
    acqLive_name = 'acqLive';
    acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
    thresh = max(Timeline.rawDAQData(:, acqLive_idx)) / 2;
    acqLive_trace = Timeline.rawDAQData(:, acqLive_idx) > thresh;
    acqLive_timeline = Timeline.rawDAQTimestamps( ...
        [find(acqLive_trace, 1), find(acqLive_trace, 1, 'last') + 1]);

    % Get wheel position
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    % (this is a very strange hack to overcome a problem in the rotary
    % encoder that's known in the lab and was put on the wiki)
    wheel_position = Timeline.rawDAQData(:, rotaryEncoder_idx);
    wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;

    % Get wheel velocity by smoothing the wheel trace and taking deriv
    wheel_smooth_t = 0.05; % seconds
    wheel_smooth_samples = wheel_smooth_t / Timeline.hw.samplingInterval;
    wheel_velocity = interp1(conv(Timeline.rawDAQTimestamps, [1, 1]/2, 'valid'), ...
        diff(smooth(wheel_position, wheel_smooth_samples)), Timeline.rawDAQTimestamps)';

    % Get whether stim was flickering
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:, stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:, stimScreen_idx)) > 2;
    end

    % Get photodiode flips (compensate for screen flicker)
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    stimScreen_on = Timeline.rawDAQData(:, photodiode_idx) > 0.2;
    stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
    photodiode_trace = Timeline.rawDAQData(stimScreen_on, photodiode_idx) > photodiode_thresh;

    photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx), 3);
    photodiode_diff_thresh = range(Timeline.rawDAQData(:, photodiode_idx)) * 0.2;
    photodiode_diff_t = 20; % time (in ms) to get delayed differential
    photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
    photodiode_diff_filt = [1, zeros(1, photodiode_diff_samples), -1];
    photodiode_trace_diff = abs(conv(photodiode_trace_medfilt, photodiode_diff_filt, 'valid')) > ...
        photodiode_diff_thresh;
    photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
        photodiode_trace_diff(2:end)) + photodiode_diff_samples + 1;
    photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
    %
    %         figure();
    %         clf;
    %         valu = 10000;
    %         title('Photodiode');
    %         hold on;
    %         plot(photodiode_trace(1:valu));
    %         hold on;
    %         scatter(photodiode_flip(find(photodiode_flip <= valu)), ones(size(find(photodiode_flip <= valu), 1), 1))
    %         hold on;
    %         plot(Timeline.rawDAQData(1:valu, photodiode_idx))
    %         hold on;
    %         plot(medfilt1(Timeline.rawDAQData(1:valu, ...
    %             photodiode_idx), 3))
    %         hold on;
    %         pp = photodiode_flip(find(photodiode_flip <= valu));

    % Get flipper signal (this was added late, might not be present)
    flipper_name = 'flipper';
    flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
    flipper_thresh = 2; % TTL threshold
    flipper_trace = Timeline.rawDAQData(:, flipper_idx) > flipper_thresh;
    flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
        (flipper_trace(1:end-1) & ~flipper_trace(2:end))) + 1;
    flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';

    %     figure();
    %     title('Flipper channel');
    %     hold on;
    %     plot(flipper_trace(1:5000));
    %     hold on;
    %     scatter(flipper_flip(find(flipper_flip <= 5000)), ones(size(find(flipper_flip <= 5000), 1), 1))
    %     hold on;
    %     plot(Timeline.rawDAQData(1:5000, flipper_idx))
end

%% Load mpep protocol

[protocol_filename, protocol_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'protocol');

if protocol_exists

    if verbose
        disp('Loading mpep protocol...');
    end

    load(protocol_filename);

    % Load in hardware info
    %hwinfo_filename = AP_cortexlab_filenameJF(animal, date, experiment, 'hardware');
    %load(hwinfo_filename);

    % Stim times should just be odd (on) and even (off)
    if mod(length(photodiode_flip_times), 2) == 0
        photodiode_onsets = photodiode_flip_times(1:2:end);
        photodiode_offsets = photodiode_flip_times(2:2:end);
    else
        error('Odd number of photodiode flips');
    end

    % Get flicker/steady photodiode mode
    %photodiode_type = lower(myScreenInfo.SyncSquare.Type);

   


    % Get stim on times
    if strcmp(Protocol.xfile, 'stimOptiWaveOutput_wLaserEnable.x')
        % Xin laser protocol

        % where blocks are
        % laserFreq_idx = strcmp(Protocol.parnames, 'freq');
        % laserTotalSeqDur_idx = strcmp(Protocol.parnames, 'dur');
        % numberOfFlips = sum(ones(1,size(Protocol.pars,2)) .* (Protocol.pars(laserTotalSeqDur_idx,:) ./ 10 - 1) .*...
        %     (Protocol.pars(laserFreq_idx,:)./10)) .* Protocol.nrepeats * 2;

        % for some reason, I can't get the right number of flips from the
        % protocol paramaters alone (is this because not set # of repeats
        % in protocol, but depends on freq + duration ?). Changed to hacky
        % way of estimating - but that seems to work.


        % use the big time between blocks(conditions), 300 samples between
        % two consecutive blocks (conditions) to know which condition each
        % onset/offset (laser_time_flip elements) falls into

        % in JulieF's comments below, sequence/block means one condition

        % get the sequence/block times
        % This is from JF, but may be wrong due to the problem of laser_flip,middle element
        % of diff_start =-5410xxx due to offset following onset.
%         diff_start = [0; diff(laser_flip)];
%         laser_block_flip_start = [laser_flip(1); laser_flip(diff_start > 300)];
%         laser_block_flip_stop = [laser_flip(diff(laser_flip) > 300); laser_flip(end)];



        % assign each laser flip to one sequence
        % use the big time between blocks(conditions), 300 samples between
        % two consecutive blocks (conditions) to know which condition each
        % onset/offset (laser_time_flip elements) falls into, 
        % laser onset after the previous interval and before the next
        % interval.
        interval = [laser_block_flip_start, laser_block_flip_end];
        laser_OnBlock = arrayfun(@(x) find(laser_flip_onidx(x) >= interval(:, 1) & laser_flip_offidx(x) <= interval(:, 2)), 1:numel(laser_flip_onidx));
        % assign each sequence to sequence type
        laserParams = table('Size', [size(interval, 1), 3], 'VariableTypes', {'double', 'double', 'double'}, ...
            'VariableNames', {'Amp', 'Freq', 'Ramp'}); % XG: Add any other you might want too look at here
        
        % get the seqnums into 1-35 for each repeats and put the number
        % this is wrong due to misunderstanding of seqnums, into seqNums_real
%         for iRepeat = 1:size(Protocol.seqnums, 2)
%             seqNums_real(:, iRepeat) = Protocol.seqnums(:, iRepeat) - ((iRepeat - 1) * Protocol.npfilestimuli);
%         end
         for iRepeat = 1:Protocol.nrepeats
             for isti=1:Protocol.npfilestimuli
             seqNums_real(isti,iRepeat)=find(Protocol.seqnums(:,iRepeat)==isti+Protocol.npfilestimuli*(iRepeat-1));
             end
         end

        sequenceNames = reshape(seqNums_real, size(seqNums_real, 1)*size(seqNums_real, 2), 1);
        sequenceType_idx = arrayfun(@(x) find(sequenceNames(x) == 1:Protocol.npfilestimuli), 1:size(sequenceNames, 1));
        % get parameters for each sequence type
        % laserParams is the spreadsheet of all conditions (number equals to No.condition multiplied by No.repeats
        % (315 for this recording)) with their amp, freq, ramp,
        % dur.
        % 
        % you can add any other parameters you like, for instance,
        %         laserdcyc_idx = strcmp(Protocol.parnames, 'dcyc');
        %         laserParams.Dcyc=Protocol.pars(laserdcyc_idx, sequenceType_idx)';
        laserFreq_idx = strcmp(Protocol.parnames, 'freq');
        laserAmp_idx = strcmp(Protocol.parnames, 'amp');
        laserRamp_idx = strcmp(Protocol.parnames, 'ramp');
        laserDur_idx = strcmp(Protocol.parnames, 'dur');
        laserParams.Amp = Protocol.pars(laserAmp_idx, sequenceType_idx)';
        laserParams.Freq = Protocol.pars(laserFreq_idx, sequenceType_idx)';
        laserParams.Ramp = Protocol.pars(laserRamp_idx, sequenceType_idx)';
        laserParams.Dur = Protocol.pars(laserDur_idx, sequenceType_idx)';

        % get the amp, freq, ramp, etc for all pulses (according to the block number each laseron falls into, which is stored in laser_onBlock)
        % it is weird that we have total 21066 laser on while the repeat
        % no. is 9...
       
        laserParamsAllLaserOn = table('Size', [size(laser_OnBlock, 1), 3], 'VariableTypes', {'double', 'double', 'double'}, ...
            'VariableNames', {'Amp', 'Freq', 'Ramp'}); % XG: Add any other parameters, for instance, duration you might want too look at here
        laserParamsAllLaserOn.Amp = laserParams.Amp(laser_OnBlock)';
        laserParamsAllLaserOn.Freq = laserParams.Freq(laser_OnBlock)';
        laserParamsAllLaserOn.Ramp = laserParams.Ramp(laser_OnBlock)';
            
        % for example, amplitude of laser pulse # 2 is laserParamsAllLaserOn.Amp(2)
       
    end

end
%find

%% Load face/eyecam and processing

% Don't load if no timeline
if exist('Timeline', 'var') && load_parts.cam

    % Get cam sync from timeline
    camSync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
    camSync_thresh = max(Timeline.rawDAQData(:, camSync_idx)) / 2;
    camSync = Timeline.rawDAQData(:, camSync_idx) > camSync_thresh;
    camSync_flip = find((camSync(1:end-1) ~= camSync(2:end))) + 1;
    if length(camSync_flip) ~= 4
        error('camSync flip number ~= 4')
    end

    % EYECAM
    [eyecam_dir, eyecam_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'eyecam');

    if eyecam_exists
        if verbose
            disp('Loading eyecam...');
        end

        % Load camera processed data
        [eyecam_processed_filename, eyecam_processed_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'eyecam_processed');
        if eyecam_processed_exists
            eyecam = load(eyecam_processed_filename);
        end

        % Get camera times
        eyecam_fn = AP_cortexlab_filenameJF(animal, day, experiment, 'eyecam');
        eyecam_dir = fileparts(eyecam_fn);
        eyecam_t_savefile = [eyecam_dir, filesep, 'eyecam_t.mat'];

        if exist(eyecam_fn, 'file') && ~exist(eyecam_t_savefile, 'file')
            % Get facecam strobes
            eyeCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'eyeCameraStrobe') | ...
                strcmp({Timeline.hw.inputs.name}, 'eyeCamStrobe');
            eyeCamStrobe_thresh = max(Timeline.rawDAQData(:, eyeCamStrobe_idx)) / 5;
            eyeCamStrobe = Timeline.rawDAQData(:, eyeCamStrobe_idx) > eyeCamStrobe_thresh;
            eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end))) + 1;
            eyeCamStrobe_up_t = Timeline.rawDAQTimestamps(eyeCamStrobe_up);

            % Get sync times for cameras (or load if already done)
            [eyecam_sync_frames, n_eyecam_frames] = AP_get_cam_sync_framesJF(eyecam_fn);

            if ~isempty(eyecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~, eyecam_strobe_sync] = min(abs(camSync_flip(1)-eyeCamStrobe_up));
                eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
                eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;

                % Check that the number of frames between synchs matches
                % video and timeline
                n_eyecam_frames_syncd_movie = diff(eyecam_sync_frames) + 1;
                [~, eyecam_strobe_sync_end] = min(abs(camSync_flip(3)-eyeCamStrobe_up));
                n_eyecam_frames_syncd_timeline = eyecam_strobe_sync_end - eyecam_strobe_sync;
                if abs(n_eyecam_frames_syncd_movie-n_eyecam_frames_syncd_timeline) > 2
                    warning('Eyecam: different n frames video vs timeline');
                end

                % Get times of cam frames in timeline
                eyecam_t = nan(n_eyecam_frames, 1);
                eyecam_t(eyecam_frame_idx(eyecam_frame_idx > 0)) = eyeCamStrobe_up_t(eyecam_frame_idx > 0);

                save(eyecam_t_savefile, 'eyecam_t');
            end
        elseif exist(eyecam_fn, 'file') && exist(eyecam_t_savefile, 'file')
            load(eyecam_t_savefile);
        end

    end

    % FACECAM
    [facecam_dir, facecam_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'facecam');

    if facecam_exists
        if verbose;
            disp('Loading facecam...');
        end

        % Get camera times
        facecam_fn = AP_cortexlab_filenameJF(animal, day, experiment, 'facecam');
        facecam_dir = fileparts(facecam_fn);
        facecam_t_savefile = [facecam_dir, filesep, 'facecam_t.mat'];

        if exist(facecam_fn, 'file') && ~exist(facecam_t_savefile, 'file')
            % Get facecam strobes
            faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
            faceCamStrobe_thresh = max(Timeline.rawDAQData(:, faceCamStrobe_idx)) / 5;
            faceCamStrobe = Timeline.rawDAQData(:, faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end))) + 1;
            faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);

            % Get sync times for cameras (or load if already done)
            [facecam_sync_frames, n_facecam_frames] = AP_get_cam_sync_framesJF(facecam_fn);

            if ~isempty(facecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~, facecam_strobe_sync] = min(abs(camSync_flip(1)-faceCamStrobe_up));
                facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
                facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;

                % Check that the number of frames between syncs matches
                % video and timeline
                n_facecam_frames_syncd_movie = diff(facecam_sync_frames) + 1;
                [~, facecam_strobe_sync_end] = min(abs(camSync_flip(3)-faceCamStrobe_up));
                n_facecam_frames_syncd_timeline = facecam_strobe_sync_end - facecam_strobe_sync;
                if abs(n_facecam_frames_syncd_movie-n_facecam_frames_syncd_timeline) > 2
                    warning('Facecam: different n frames video vs timeline');
                end

                % Get times of cam frames in timeline
                facecam_t = nan(n_facecam_frames, 1);
                facecam_t(facecam_frame_idx(facecam_frame_idx > 0)) = faceCamStrobe_up_t(facecam_frame_idx > 0);

                save(facecam_t_savefile, 'facecam_t');
            end
        elseif exist(facecam_fn, 'file') && exist(facecam_t_savefile, 'file')
            load(facecam_t_savefile);
        end

        % (old/unused: etGUI and facemap)
        [facecam_processed_filename, facecam_processed_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'facecam_processed');
        if facecam_processed_exists
            facecam = load(facecam_processed_filename);
        end

        % (output from AP_mouse_movie_movement)
        [facecam_movement_filename, facecam_movement_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'facecam_movement');
        if facecam_movement_exists
            load(facecam_movement_filename);
        end

    end

end

%% Load ephys data (single long recording)

% Pick kilosort version (2 by default, 1 old if selected)
if ~exist('kilosort_version', 'var') || kilosort_version == 2
    [ephys_path, ephys_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
elseif exist('kilosort_version', 'var') && kilosort_version == 1
    [ephys_path, ephys_exists] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ks1', site, recording);
end

if ephys_exists && load_parts.ephys

    if verbose
        disp('Loading ephys...');
    end


    % Load phy sorting if it exists
    % (old = cluster_groups.csv, new = cluster_group.tsv because fuck me)
    cluster_filepattern = [ephys_path, filesep, 'cluster_group*'];
    cluster_filedir = dir(cluster_filepattern);
    if ~isempty(cluster_filedir)
        cluster_filename = [ephys_path, filesep, cluster_filedir.name];
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid, '%d%s', 'HeaderLines', 1);
        cg = cluster_groups;
        fclose(fid);
    end

    % Load sync/photodiode
    load(([ephys_path, filesep, 'sync.mat']));

    % Read header information

    header.n_channels = 385;
    ephys_sample_rate = 30000;
    header.lfp_sample_rate = 30000;
    header.filter_cutoff = 300; %Hz, defined in JF_computeLFP

    spike_times = double(readNPY([ephys_path, filesep, 'spike_times.npy'])) ./ ephys_sample_rate;
    spike_templates_0idx = readNPY([ephys_path, filesep, 'spike_templates.npy']);
    templates_whitened = readNPY([ephys_path, filesep, 'templates.npy']);
    channel_positions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
    channel_map = readNPY([ephys_path, filesep, 'channel_map.npy']);
    winv = readNPY([ephys_path, filesep, 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path, filesep, 'amplitudes.npy']);

    % Default channel map/positions are from end: make from surface
    % (hardcode this: kilosort2 drops channels)
    if isSpikeGlx
        max_depth = 2880;
        if any(max_depth-channel_positions(:, 2) < 0) %1.0
            max_depth = 3840;
            channel_positions(:, 2) = max_depth - channel_positions(:, 2);
        else
            max_depth = 2880;
            channel_positions(:, 2) = max_depth - channel_positions(:, 2);
        end

    else
        max_depth = 3840;
        channel_positions(:, 2) = max_depth - channel_positions(:, 2); % 0 = tip in NP1s, 3840 = top, reorder here
    end


    % Unwhiten templates
    templates = zeros(size(templates_whitened));
    for t = 1:size(templates_whitened, 1)
        templates(t, :, :) = squeeze(templates_whitened(t, :, :)) * winv;
    end

    % Get the waveform of all templates (channel with largest amplitude)
    [~, max_site] = max(max(abs(templates), [], 2), [], 3);
    templates_max = nan(size(templates, 1), size(templates, 2));
    for curr_template = 1:size(templates, 1)
        templates_max(curr_template, :) = ...
            templates(curr_template, :, max_site(curr_template));
    end
    waveforms = templates_max;

    % Get depth of each template
    % (get min-max range for each channel)
    template_chan_amp = squeeze(range(templates, 2));
    % (zero-out low amplitude channels)
    template_chan_amp_thresh = max(template_chan_amp, [], 2) * 0.5;
    template_chan_amp_overthresh = template_chan_amp .* (template_chan_amp >= template_chan_amp_thresh);
    % (get center-of-mass on thresholded channel amplitudes)
    template_depths = sum(template_chan_amp_overthresh.*channel_positions(:, 2)', 2) ./ sum(template_chan_amp_overthresh, 2);
    template_xdepths = channel_positions(max_site, 1);

    % Get the depth of each spike (templates are zero-indexed)
    spike_depths = template_depths(spike_templates_0idx+1);
    spike_xdepths = template_xdepths(spike_templates_0idx+1);

    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times, templates_max, ...
        sign(abs(min(templates_max, [], 2))-abs(max(templates_max, [], 2))));

    [~, waveform_trough] = min(templates_max, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(templates_max, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    % Get sync points for alignment

    % Get experiment index by finding numbered folders
    protocols_list = AP_list_experimentsJF(animal, date);
    experiment_idx = experiment == [protocols_list.experiment];
    ops.recording_software = 'SpikeGLX';
    ops.ephys_folder = [ephysAPfile, '/..'];
    [expInfo, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'expInfo', site);
    try
        [co] = mainprobe_to_timeline(ephys_path, ...
            Timeline, ops, expInfo);
        spike_times_timeline = spike_times * co(2) + co(1);
    catch
        warning('probe aligning error, using un-aligned spikes_times')
        spike_times_timeline = spike_times;
    end

    % Get "good" templates from labels
    if exist('cluster_groups', 'var') && loadClusters
        % If there's a manual classification


        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique, uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2}, {'good', 'mua', 'noise'}))
            disp('Removing manually labelled noise units...');
            good_templates_idx = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'unsorted') | strcmp(cluster_groups{2}, 'mua') | ...
                strcmp(cluster_groups{2}, 'good')));
            template_label_mua = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'mua')));
            template_label_good = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good')));
            template_labelM = good_templates_idx(ismember(good_templates_idx, template_label_mua));
            template_labelG = good_templates_idx(ismember(good_templates_idx, template_label_good));

            [good_templates, ii] = ismember(0:size(templates, 1)-1, good_templates_idx);

            [good_templatesG, ii] = ismember(0:size(templates, 1)-1, template_labelG);
            [good_templatesM, ii] = ismember(0:size(templates, 1)-1, template_labelM);


        else
            disp('Keeping manually labelled good units...');
            good_templates_idx = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good') | strcmp(cluster_groups{2}, 'mua')));
            template_label_mua = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'mua')));
            template_label_good = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good')));
            template_labelM = good_templates_idx(ismember(good_templates_idx, template_label_mua));
            template_labelG = good_templates_idx(ismember(good_templates_idx, template_label_good));

            [good_templates, ii] = ismember(0:size(templates, 1)-1, good_templates_idx);

            [good_templatesG, ii] = ismember(0:size(templates, 1)-1, template_labelG);
            [good_templatesM, ii] = ismember(0:size(templates, 1)-1, template_labelM);

            template_label = good_templatesG .* 1 + (good_templatesM) .* 2;
        end


    elseif exist('unitType', 'var')
        % If no manual but qualityMetrics are available
        if verbose;
            disp('Keeping quality metrics good units...');
        end

        % Load triage labels

        %triage_good_templates = goodUnits;

        good_templates = ...
            unitType == 1;
        good_templates_idx = find(unitType == 1) - 1;

        if exist('locationKeep', 'var')
            if verbose;
                disp('Keeping location data...');
            end
            myPaths;
            allenAt = loadStructureTreeJF([allenAtlasPath, filesep, 'allenCCF/structure_tree_safe_2017.csv']);
            probeccf = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
            load(probeccf)
            this_ccf = probe_ccf(probes(iDataset));
            theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
            theseLocationsInterest = contains(theseLocations, locationKeep);
            theseDepths = this_ccf.probe_depths(theseLocationsInterest);
            template_exists = ismember(1:max(spikeTemplates), unique(spikeTemplates));
            theseTemplates = template_depths(template_exists) >= min(theseDepths) & template_depths(template_exists) <= max(theseDepths); %correct depth units
            good_templates = goodUnits & theseTemplates';
            good_templates_idx = find(good_templates) - 1;
        end
    else
        % If no cluster groups at all, keep all
        warning([animal, ' ', day, ' - no cluster groups']);
        if verbose
            disp('No manual labeling, keeping all and re-indexing');
        end
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates, 1)-1, good_templates_idx);
    end
    % Throw out all non-good template data
    templates = templates(good_templates, :, :);
    template_depths = template_depths(good_templates);
    template_xdepths = template_xdepths(good_templates);
    waveforms = waveforms(good_templates, :);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);
    %template_label = template_label(good_templates);
    % Throw out all non-good spike data
    good_spike_idx = ismember(spike_templates_0idx, good_templates_idx);
    spike_times = spike_times(good_spike_idx);
    spike_times_full = spike_times_timeline;
    spike_templates_full = spike_templates_0idx + 1;
    spike_templates_0idx = spike_templates_0idx(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_depths = spike_depths(good_spike_idx);
    spike_xdepths = spike_xdepths(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);

    % Rename the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spike_templates_0idx)+1, 1);
    new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
    spike_templates = new_spike_idx(spike_templates_0idx+1);

end

%% Finished
if verbose
    disp('Finished loading experiment.');
end
