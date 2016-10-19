function [autocorrelation, autotime, lifetimedata, headdata, res_final] = part_crosscorrelation(inputff, tg1, tg2, channel, package_size, overlap, calc_reverse, options, calc_options)

% [autocorrelation, autotime] = part_crosscorrelation(inputff, tg1, tg2, channel, package_size, calc_reverse, overlap, options)
%
% inputff           => cell array with filenames to process
% tg1, tg2, channel => see filter_channel.m (tg1 may be a 1x2 matrix, 
%                                            tg2 has to be a 16x2 matrix,
%                                            channel has to be 16x1 matrix)
%
% unit of tg1 is expected to be seconds
% unit of tg2 is expected to be slots
%
% package_size      => size of dataset to process in one step 
%                      (min. 1 Million - max. 4/5 Million => depends on RAM)
% overlap           => size of dataset to take into each package
%                      additionally (has to be smaller than package_size!)
% calc_reverse      => if greater 0 a calculation with reversed time will
%                      be done (should smooth the resulting curve and avoid
%                      nummerical instability)
% options           => structure holding information for additional options:
%                      .statusfield => Can be a guide handle for a text
%                                      field -> A status will be printed
%                                      there
%                      .dc          => structure containing information
%                                      necessary for distributed computing
%                                      currently hold only the filed
%                                      lookupurl which must be suitable to
%                                      be used as LookupURL argument for
%                                      findResource
%
% calc_options      => options for the calculation
%                      nsub + max_t => see below
%
% NOTE: You need the cluster or more than enough memory!!!!!!
%
% Read the header data to get the overall number of photons received
% this is needed to calculate the number of packages we have to process.
% This is the overall number deviced by the package_size

% Create a temp directory we use to store our temporary files
tempdir = 'temp';
mkdir(tempdir);

% How are the 16 filters combined to form the 4 8x8 Matrix?
filter_combination = {1:8 1:8;1:8 9:16;9:16 1:8;9:16 9:16};

% Save the names of the files we can remove at the end of the run
filelist = {};

% Randompart for the filename
random_name = round(rand * 10000000);

filecount = size(inputff, 2);

% autocorrelation holds the result
% (one row per filename, calculated upwards, calculated downwards, average)
autocorrelation = cell(filecount,3);
% Auto holds the results per quadrant => with multiplefiles this has to get
% three dimensional
auto = cell(filecount,4,2);
autotime = cell(1,filecount);
lifetimedata = cell(1,filecount);
headdata = cell(1,filecount);

if exist('options', 'var')
    if isfield(options, 'statusfield')
        statusfield = options.statusfield;
    end
    if isfield(options, 'dc')
        fprintf('Distributed Computing Options found!\n');
        if isfield(options.dc, 'config')
            fprintf('Using Config to get JobManager: %s\n', options.dc.config);
            jm = findResource('scheduler','configuration', options.dc.config);
            set(jm,'configuration', options.dc.config);
        else
            fprintf('Using LookupURL to get JobManager: %s\n', options.dc.lookupurl);
            jm = findResource('scheduler', 'type', 'jobmanager', 'Name', options.dc.jmname, 'LookupURL', options.dc.lookupurl);
        end
        jobs = {};
    end
end

Nsub = 10;
max_t = 6;
if exist('calc_options', 'var')
    if isfield(calc_options, 'max_t')
        max_t = calc_options.max_t;
    end
    if isfield(calc_options, 'Nsub')
        Nsub = calc_options.Nsub;
    end
end


for f=1:filecount
    % Calculate Ncasc and Nsub
    head = ht3v2read_cl(inputff{1}, [0 0]);
%   xxx head.Sync = 1E9/head.SyncRate;
%    Timeunit=head.Sync * 1e-9;
    Timeunit=1/head.SyncRate;

    Ncasc = ceil(log2(max_t/Timeunit/Nsub));
    autotime{f} = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
    autotime{f} = autotime{f} * Timeunit;

    fprintf('Nsub: %i, Ncasc: %i\n', Nsub, Ncasc);
    
    % Make sure the result matrix is initialized - we should do it
    % here - consider an empty photon set -- this would finally crash else
    for l=1:4
        auto{f,l,1} = zeros(Ncasc * Nsub, 8, 8);
        auto{f,l,2} = zeros(Ncasc * Nsub, 8, 8);
    end
    
    ff = inputff{f};
    head = ht3v2read_cl(ff, [0 0]);
    count_packages = ceil(head.Records / package_size);

    overcount = 0;
    old_overcount = 0;

    for i=1:count_packages
        filedeps = cell(1,3);
        filedeps{1} = 'tttr2xfcs_cl.m';
        filedeps{2} = 'calc_crosscorr.m';
        name = sprintf('crosscorr_data_%i-%i-%i.mat', f, i, random_name);
        filedeps{3} = fullfile(tempdir, name);
        filelist{size(filelist,2)+1} = fullfile(tempdir, name);

        fprintf('Correlation: %03i/%03i - File %i\n', i, count_packages, f);
        if exist('statusfield', 'var')
            set(statusfield, 'String', sprintf('Correlation: %03i/%03i - File %i', i, count_packages, f));
            drawnow();
        end
        si = (i-1)*package_size - overlap;
        if si < 0
            si = 0;
        end
        ei = i*package_size;
        if ei > head.Records
            ei = head.Records;
        end
%        [head, z, tcspc, chan, markers, num, overcount] = pt3v2read(ff, [si ei-si]);
        [y, tcspc, chan, markers, num, head] = ht3v2read_cl(ff, [si ei-si+1]);
        %chan = (chan * (-1)) + 4;
        chan = chan +1;
        fprintf('Read %i photons\n', num);
        old_overcount = overcount;
        z = (overcount + y);% * 1000 / head.SyncRate;
        overcount = y(size(y,1),1) + old_overcount;

        % If we only get one mcs/tg1 timegate (this is probable) fill up the list
        % to match the size of tg2
        % 
        % use_tg1 has to be calculated, as tg1 is in seconds and z in
        % Timeunits so -- beware this is different from part_lifetime, as there
        % tg1 is expected in s (there we have the head data available, not present
        % in filter_channel):
        use_tg1 = tg1 ./ Timeunit;
        for j=size(tg1,1):size(tg2,1)
            use_tg1(j,:) = use_tg1(1,:);
        end
        
        if use_tg1(1,1) > max(z) || use_tg1(1,2) <= min(z)
            fprintf('Discarding Package -- Outside MCS Timegates\n');
            continue
        end
        
        filter_matrix = filter_channel(z, tcspc, chan, use_tg1, tg2, channel);
        
        % name contains every data, the crosscorrelation will need
        % so only the name of this file is passed to the worker and
        % everything not saved here has to be considered not available
        save(fullfile(tempdir, name), 'z', 'filter_matrix', 'Ncasc', 'Nsub', 'calc_reverse');

        clear z tcspc chan markers num filterset1 filterset2 filter_matrix si ei;
        % Do the real work, after the unneccessary data is cleared

        if exist('jm', 'var')
            fprintf('creating job %i,%i\n', f, i);
            job = createJob(jm);
            % We missuse the JobData Property to store the location for
            % the result
            set(job, 'JobData', [f,i]);
            set(job, 'FileDependencies', filedeps);
        end
        for l=1:4
            if exist('job', 'var')
                fprintf('creating task %i\n', l);
                createTask(job, @calc_crosscorr, 3, {name, filter_combination{l,1}, filter_combination{l,2}});
            else
                [t1 t2] = feval(@calc_crosscorr, fullfile(tempdir, name), filter_combination{l,1}, filter_combination{l,2});
                auto{f,l,1} = auto{f,l,1} + t1;
                auto{f,l,2} = auto{f,l,2} + t2;
            end
        end
        if exist('job', 'var')
            submit(job);
            jobs{size(jobs,2)+1} = job;
        end
    end
end

% Here we calculate the lifetime data - we can expect, that the time would
% else be spend waiting for the workers to finish processing the data
%
% Maybe put this work out into distributed computing - this has to be
% benchmarked
for f=1:filecount
    fprintf('Calculating Lifetimedata %i/%i\n', f, filecount);
    if exist('statusfield', 'var')
        set(statusfield, 'String', sprintf('Calculating Lifetimedata %i/%i', f, filecount));
        drawnow();
    end
    %[lifetimedata{f} headdata{f}] = part_lifetimedata(inputff{f}, package_size, tg1);
    % as well as the MCS
    tg.tg1 = tg1;
    tg.tg2 = tg2;
    [mcs lifetimedata{f} headdata{f} res1{f}] = part_mcs(inputff{f}, package_size, tg);
end

if exist('jm', 'var')
    fprintf('Waiting for Jobdata\n');
    if exist('statusfield', 'var')
        set(statusfield, 'String', 'Waiting for Jobdata');
        drawnow();
    end
    % Check if the jobs are ready
    finished_jobs = 0;
    while finished_jobs ~= size(jobs, 2)
        % Wait 20 seconds to give time to calculate
        finished_jobs = 0;
        for i=1:size(jobs, 2)
            job = jobs{i};
            if strcmp(get(job, 'State'), 'failed')
                % Seems to be necessary as the RZ Aachen has the brain-dead
                % config of running more workers than licences available
                jobs{i} = resubmitjob(jm, jobs{i});
                continue;
            end
            if strcmp(get(job, 'State'), 'finished')
                [pending running finished] = findTask(job);
                without_error = 1;
                for j=1:size(finished, 2)
                    if ~ isempty(get(finished(j), 'ErrorIdentifier'))
                        fprintf('Will resubmit job %i because of an Error in Task %i:\n%s\n', i, j, get(finished(j), 'ErrorMessage'));
                        without_error = 0;
                    end
                end
                if without_error == 1 || size(finished, 2) == 0
                    finished_jobs = finished_jobs + 1;
                else
                    jobs{i} = resubmitjob(jm, jobs{i});
                    continue;
                end
            end
        end
        
        fprintf('%i Jobs of %i finished\n', finished_jobs, size(jobs,2));
        if exist('statusfield', 'var')
            set(statusfield, 'String', sprintf('%i Jobs of %i finished\n', finished_jobs, size(jobs,2)));
            drawnow();
        end
        if finished_jobs ~= size(jobs,2)
           pause(20) 
        end
    end
    
    
    % Fetch the data - when all jobs have finished
    for i=1:size(jobs, 2)
        jobdata = get(jobs{i}, 'JobData');
        res = getAllOutputArguments(jobs{i});
        % Jobdata 1=>File, 2=>Package
        for j=1:size(res,1)
            auto{jobdata(1),j,1} = auto{jobdata(1),j,1} + res{j,1};
            auto{jobdata(1),j,2} = auto{jobdata(1),j,2} + res{j,2};
        end
        resb(:,1:8,1:8,i) = res{1,1};
        resb(:,1:8,9:16,i) = res{2,1};
        resb(:,9:16,1:8,i) = res{3,1};
        resb(:,9:16,9:16,i) = res{4,1};
        
        ind = [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
        
        for k = 1:16
            for l = 1:16
                resa.auto(:,k,l,i) = resb(:,ind(k),ind(l),i);
            end
        end
        
        res2{f} = resa;
        fprintf('Got result for job %i of %i\n', i, size(jobs,2));
        if exist('statusfield', 'var')
            set(statusfield, 'String', sprintf('Got result for job %i of %i', i, size(jobs,2)));
            drawnow();
        end
        destroy(jobs{i});
    end
end
clear res;
for f=1:filecount
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    res_a = res2{f};
    res.auto = res_a.auto;
    res_a = res1{f};
    res.bin = res_a.bin;
    res.tau = res_a.tau * 10e-9;
    res.tcspc = res_a.tcspc;
    res.time = res_a.time;
    res.rate = res_a.rate;
    res.tcspcdata = res_a.tcspcdata;
    res.mcs = res_a.mcs;
    res_final{f} = res;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i=1:2
        tempmatrix = zeros(size(auto{f,1,i},1),16,16);
        tempmatrix(:,1:8,1:8) = auto{f,1,i};
        tempmatrix(:,1:8,9:16) = auto{f,2,i};
        tempmatrix(:,9:16,1:8) = auto{f,3,i};
        tempmatrix(:,9:16,9:16) = auto{f,4,i};
        if i==2
            for j=1:size(tempmatrix,1)
                tempmatrix(j,:,:) = shiftdim(tempmatrix(j,:,:))';
            end
        end
        autocorrelation{f,i} = tempmatrix;
    end
    if calc_reverse > 0
        autocorrelation{f,3} = (autocorrelation{f,2} + autocorrelation{f,1}) / 2;
    else
        autocorrelation{f,3} = autocorrelation{f,1};
    end
end

% Finally remove the intermediatly created .mat files
for i=1:size(filelist,2)
    delete(filelist{i});
end

function [job] = resubmitjob(jm, old_job)
    fprintf('BEWARE -- WE RESUBMIT A JOB - THIS IS SHOULD NOT HAPPEN!\n');
    fprintf('YOU HAVE BEEN WARNED!\n');
    % Recreate the job state of the old_job in the new job
    % This is plain evil and against everything good!
    job = createJob(jm);
    set(job, 'JobData', get(old_job, 'JobData'));
    set(job, 'FileDependencies', get(old_job, 'FileDependencies'));
    old_tasks = get(old_job, 'Tasks');
    for i=1:size(old_tasks)
        old_task = old_tasks(i);
        createTask(job, get(old_task, 'Function'), get(old_task, 'NumberOfOutputArguments'), get(old_task, 'InputArguments'));
    end
    destroy(old_job);
    submit(job);
