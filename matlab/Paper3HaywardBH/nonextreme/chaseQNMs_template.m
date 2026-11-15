% QNM Computation Script - Template
% This file should be copied to each L* subdirectory as chaseQNMs.m
% Parameters will be passed via a parameter file created by Python
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)

close all
clearvars  % Clears all variables but preserves breakpoints and function permissions
format longE

% Folder of the Advanpix toolbox, used when the parameter file gives none
defaultAdvanpixPath = '/path/to/advanpix/';

% No additional paths needed - script is self-contained

% Read parameters from the parameter file created by Python
paramFile = 'matlab_params.mat';
if ~exist(paramFile, 'file')
    error('Parameter file %s not found. This script should be run by the Python orchestrator.', paramFile);
end

% Load parameters
load(paramFile, 's_value', 'L_value', 'g_value', 'resolutions', 'job_timestamp', 'advanpix_path', 'use_parallel');

% Add AdvanPix toolbox path from configuration
fprintf('Configuring AdvanPix multiprecision toolbox...\n');
if exist('advanpix_path', 'var') && ~isempty(advanpix_path)
    fprintf('  Adding AdvanPix path: %s\n', advanpix_path);
    addpath(advanpix_path);
else
    % Fallback to default path if not provided
    fprintf('  Using default AdvanPix path: %s\n', defaultAdvanpixPath);
    addpath(defaultAdvanpixPath);
end
fprintf('  AdvanPix toolbox loaded successfully\n');

% Check if parallel flag exists (for backward compatibility)
if ~exist('use_parallel', 'var')
    use_parallel = false;
end

% Display parameters
fprintf('\n========================================\n');
fprintf('QNM Computation Started\n');
fprintf('========================================\n');
fprintf('Timestamp: %s\n', job_timestamp);
fprintf('Parameters: s = %d, L = %d\n', s_value, L_value);
fprintf('Resolutions: [%s]\n', num2str(resolutions));
if use_parallel
    fprintf('Parallel Mode: Enabled\n');
else
    fprintf('Parallel Mode: Disabled\n');
end
fprintf('========================================\n\n');

% For compatibility with existing code
s_maple = s_value;
L_maple = L_value;
g_maple = g_value;  % This will be passed from Python

list = resolutions;
L = length(list);

% Job name for identification
jobName = sprintf('QNM_s%d_L%d_g%s', s_maple, L_maple, strrep(num2str(g_maple), '/', '_'));

% Timer for overall computation
overall_timer = tic;

qnms = cell(1, L);
legs = cell(1, L);

% Create a simple file-based communication for Python to read
function write_progress(message)
    % Use append mode to avoid overwriting messages
    fid = fopen('matlab_progress.txt', 'a');
    if fid ~= -1
        fprintf(fid, '%s\n', message);  % Add newline separator
        fclose(fid);
    end
end

% Write progress with timestamp for better tracking
function write_progress_with_status(resolution, status, time_taken)
    if nargin < 3
        time_taken = 0;
    end
    
    timestamp = datestr(now, 'HH:MM:SS');
    if strcmp(status, 'START')
        message = sprintf('RESOLUTION_START:%d:%s', resolution, timestamp);
        % Also write to console for immediate feedback
        fprintf('[%s] Starting eigenvalue computation for resolution %d...\n', timestamp, resolution);
    elseif strcmp(status, 'COMPLETE')
        message = sprintf('RESOLUTION_COMPLETE:%d:%.2f:%s', resolution, time_taken, timestamp);
        % Also write to console for immediate feedback
        fprintf('[%s] Completed resolution %d in %.2f seconds\n', timestamp, resolution, time_taken);
    else
        message = sprintf('RESOLUTION_UPDATE:%d:%s:%s', resolution, status, timestamp);
        fprintf('\n[%s] Resolution %d: %s\n', timestamp, resolution, status);
    end
    
    write_progress(message);
end

% Compute QNMs for each resolution with optimized approach
resolution_times = zeros(1, L);

% Clear any existing progress file
if exist('matlab_progress.txt', 'file')
    delete('matlab_progress.txt');
end

% Notify Python that we're starting MATLAB computations
write_progress('MATLAB_STARTED');

% Initialize pool_started to false. It will be set to true if the parallel pool starts successfully.
pool_started = false;

if use_parallel && L > 1
    % Parallel computation mode
    fprintf('\n--- Parallel Computation Configuration ---\n');
    fprintf('Attempting to use PARALLEL computation (%d resolutions)...\n', L);
    
    % Robust parallel pool initialization with error handling
    pool = gcp('nocreate');
    
    if isempty(pool)
        try
            % First, delete any existing pool that might be corrupted
            delete(gcp('nocreate'));
            
            % Try to start a new parallel pool with explicit configuration
            num_workers = min(L, feature('numcores'));
            fprintf('Starting parallel pool with %d workers...\n', num_workers);
            
            % Create pool with timeout and error handling
            pool = parpool('local', num_workers, 'SpmdEnabled', false, 'IdleTimeout', Inf);
            pool_started = true;
            fprintf('Parallel pool started successfully with %d workers.\n', pool.NumWorkers);
        catch ME
            fprintf('\nWARNING: Failed to start parallel pool. Error: %s\n', ME.message);
            fprintf('Detailed error information:\n');
            disp(ME);
            
            % Additional diagnostics
            fprintf('\nChecking parallel computing environment:\n');
            try
                ver_info = ver('parallel');
                if isempty(ver_info)
                    fprintf('  - Parallel Computing Toolbox not found\n');
                else
                    fprintf('  - Parallel Computing Toolbox version: %s\n', ver_info.Version);
                end
            catch
                fprintf('  - Unable to check Parallel Computing Toolbox status\n');
            end
            
            % Check for common issues
            if contains(ME.message, 'license', 'IgnoreCase', true)
                fprintf('  - This appears to be a licensing issue\n');
            elseif contains(ME.message, 'memory', 'IgnoreCase', true)
                fprintf('  - This appears to be a memory issue\n');
            end
            
            pool_started = false;
        end
    else
        pool_started = true;
        fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
    end
    
    % If parallel pool failed to start, fall back to sequential mode
    if ~pool_started
        fprintf('\nFalling back to SEQUENTIAL computation due to parallel pool failure.\n');
        write_progress('PARALLEL_FAILED_FALLBACK_SEQUENTIAL');
        use_parallel = false;
    end
end

if use_parallel && L > 1 && pool_started
    try
        % Pre-allocate cell arrays for parallel results
        parfor idx = 1:L
            res_timer = tic;
            n = list(idx);
            nstr = num2str(n);
            
            % Report start of this resolution
            fprintf('[Worker %d] Starting resolution n = %d ...\n', idx, n);
            
            mp.Digits(n);
            
            % Load matrices
            m0name = strcat('data/M0_', nstr, '.mat');
            m1name = strcat('data/M1_', nstr, '.mat');
            m2name = strcat('data/M2_', nstr, '.mat');
            
            if ~exist(m0name, 'file') || ~exist(m1name, 'file') || ~exist(m2name, 'file')
                error('Matrix files not found. Make sure Maple computation completed successfully.');
            end
            
            M0 = mp.read(m0name);
            M1 = mp.read(m1name);
            M2 = mp.read(m2name);
            
            % Compute eigenvalues
            e = polyeig(M0, mp('1i')*M1, M2);
            
            qnms{idx} = e;
            legs{idx} = strcat(' ', nstr);
            
            resolution_times(idx) = toc(res_timer);
            time_str = format_time_matlab(resolution_times(idx));
            fprintf('[Worker %d] Completed resolution n = %d. (%s)\n', idx, n, time_str);
        end
        
        % Report completion of all parallel tasks
        for idx = 1:L
            write_progress_with_status(list(idx), 'COMPLETE', resolution_times(idx));
        end
    catch ME
        fprintf('\nERROR during parallel computation: %s\n', ME.message);
        fprintf('Attempting to fall back to sequential mode...\n');
        write_progress('PARALLEL_RUNTIME_ERROR_FALLBACK');
        
        % Clean up parallel pool
        delete(gcp('nocreate'));
        
        % Force sequential mode
        use_parallel = false;
    end
end

% Check if we should use sequential mode (either by choice or fallback)
if ~exist('pool_started', 'var')
    pool_started = false;
end

if ~use_parallel || L == 1 || ~pool_started
    % Sequential computation mode
    fprintf('\n--- Sequential Computation Mode ---\n');
    if ~use_parallel
        fprintf('Using SEQUENTIAL computation (%d resolutions)...\n', L);
        fprintf('Processing resolutions one by one for detailed progress tracking\n');
    end
    
    % Sequential computation with enhanced progress reporting
    for idx = 1:L
        res_timer = tic;
        n = list(idx);
        nstr = num2str(n);
        
        % Report start of this resolution
        write_progress_with_status(n, 'START');
        
        fprintf('[%d/%d] Computing QNMs for resolution n = %d ... ', idx, L, n);
        
        mp.Digits(n);
        fprintf('Setting precision to %d digits\n', n);
        
        % Load matrices
        m0name = strcat('data/M0_', nstr, '.mat');
        m1name = strcat('data/M1_', nstr, '.mat');
        m2name = strcat('data/M2_', nstr, '.mat');
        
        if ~exist(m0name, 'file') || ~exist(m1name, 'file') || ~exist(m2name, 'file')
            error('Matrix files not found. Make sure Maple computation completed successfully.');
        end
        
        write_progress_with_status(n, 'LOADING_MATRICES');
        fprintf('  Loading matrices from data/ directory...\n');
        fprintf('    Loading M0_%s.mat...', nstr);
        M0 = mp.read(m0name);
        fprintf(' done. Size: %dx%d\n', size(M0,1), size(M0,2));
        
        fprintf('    Loading M1_%s.mat...', nstr);
        M1 = mp.read(m1name);
        fprintf(' done. Size: %dx%d\n', size(M1,1), size(M1,2));
        
        fprintf('    Loading M2_%s.mat...', nstr);
        M2 = mp.read(m2name);
        fprintf(' done. Size: %dx%d\n', size(M2,1), size(M2,2));
        
        write_progress_with_status(n, 'COMPUTING_EIGENVALUES');
        fprintf('  Computing eigenvalues using polyeig...\n');
        fprintf('    Building polynomial eigenvalue problem...\n');
        fprintf('    Matrix dimensions: %dx%d\n', size(M0,1), size(M0,2));
        fprintf('    This may take a while for large matrices...\n');
        
        % Compute eigenvalues
        eig_start = tic;
        e = polyeig(M0, mp('1i')*M1, M2);
        eig_time = toc(eig_start);
        fprintf('    Eigenvalue computation completed in %.2f seconds\n', eig_time);
        fprintf('  Found %d eigenvalues\n', length(e));
        
        qnms{idx} = e;
        legs{idx} = strcat(' ', nstr);
        
        resolution_times(idx) = toc(res_timer);
        time_str = format_time_matlab(resolution_times(idx));
        fprintf('done. (%s)\n', time_str);
        
        % Report completion with timing (keep raw seconds for Python)
        write_progress_with_status(n, 'COMPLETE', resolution_times(idx));
    end
end

% Save the computed QNM vectors
fprintf('\n--- Saving Results ---\n');
if ~exist('results', 'dir')
    fprintf('Creating results directory...\n');
    mkdir('results');
end
fprintf('Saving QNM data to results/qnms.mat...');
save('results/qnms.mat', 'qnms');
fprintf(' done\n');
fprintf('QNM data contains %d resolution sets\n', length(qnms));

% Analysis phase
analysis_timer = tic;
fprintf('\n========================================\n');
fprintf('Starting QNM Analysis\n');
fprintf('========================================\n');
write_progress('ANALYSIS_START');

% Find Common QNMs in All Three Vectors
if L ~= 3
    error('This analysis requires exactly three resolutions (L = 3).');
end

fprintf('\nSearching for common modes across all resolutions...\n');
tol = mp('1e-4');  % Tolerance for matching QNMs
fprintf('  Using tolerance: %s\n', num2str(tol));

e1 = qnms{1};
e2 = qnms{2};
e3 = qnms{3};
fprintf('  Resolution %d: %d modes\n', list(1), length(e1));
fprintf('  Resolution %d: %d modes\n', list(2), length(e2));
fprintf('  Resolution %d: %d modes\n', list(3), length(e3));

fprintf('\nComparing modes from resolution %d against others...\n', list(1));
commonQNMs = [];  % Will store modes present in all three vectors
for j = 1:length(e1)
    currentMode = e1(j);
    if (min(abs(e2 - currentMode)) < tol && min(abs(e3 - currentMode)) < tol)
        commonQNMs = [commonQNMs; currentMode];
    end
end
fprintf('  Found %d matching modes\n', length(commonQNMs));

% Remove duplicates
commonQNMs = unique(commonQNMs);

% Categorize and analyze QNMs
fprintf('\nCategorizing %d common QNMs...\n', length(commonQNMs));
write_progress(sprintf('ANALYSIS_COMMON_MODES:%d', length(commonQNMs)));
imagThreshold = 1e-16;
fprintf('  Using imaginary threshold: %e\n', imagThreshold);
purelyImaginary = [];
positiveReal = [];

for j = 1:length(commonQNMs)
    mode = commonQNMs(j);
    if abs(real(mode)) < imagThreshold
        purelyImaginary = [purelyImaginary; mode];
    elseif real(mode) > 0
        positiveReal = [positiveReal; mode];
    end
end

% Sort modes
if ~isempty(positiveReal)
    [~, sortIdx] = sort(imag(positiveReal), 'descend');
    positiveReal = positiveReal(sortIdx);
end
if ~isempty(purelyImaginary)
    [~, sortIdx] = sort(imag(purelyImaginary), 'descend');
    purelyImaginary = purelyImaginary(sortIdx);
end

uniqueModeCount = length(positiveReal) + length(purelyImaginary);

% Display results with actual mode values
fprintf('\n========================================\n');
fprintf('QNM Analysis Results\n');
fprintf('========================================\n');
fprintf('Common modes found: %d (in all resolutions)\n', length(commonQNMs));
fprintf('  - Oscillatory modes: %d\n', length(positiveReal));
fprintf('  - Overdamped modes: %d\n', length(purelyImaginary));
fprintf('========================================\n');

% Send analysis summary to Python for Telegram
analysis_msg = sprintf('ANALYSIS_SUMMARY:%d:%d:%d', length(commonQNMs), length(positiveReal), length(purelyImaginary));
write_progress(analysis_msg);

% Display actual mode values like the old version
if ~isempty(positiveReal)
    fprintf('\nOscillatory modes (positive real part):\n');
    for k = 1:length(positiveReal)
        mode_str = sprintf('%.6f%+.6f*I', real(positiveReal(k)), imag(positiveReal(k)));
        fprintf('  Mode #%-3d: %s\n', k, mode_str);
        % Send each oscillatory mode to Python
        if k <= 20  % Send first 20 oscillatory modes
            write_progress(sprintf('OSCILLATORY_MODE:%d:%s', k, mode_str));
        end
    end
end

if ~isempty(purelyImaginary)
    fprintf('\nOverdamped modes (purely imaginary):\n');
    % Show first 20 modes
    numToShow = min(20, length(purelyImaginary));
    for k = 1:numToShow
        current_imag = imag(purelyImaginary(k));
        delta_str = '      N/A';
        if k < length(purelyImaginary)
            next_imag = imag(purelyImaginary(k+1));
            delta_val = abs(current_imag - next_imag);
            delta_str = sprintf('%9.6f', delta_val);
        end
        mode_str = sprintf('%+11.6f*I', current_imag);
        fprintf('  Mode #%-3d: %s   Delta: %s\n', k, mode_str, delta_str);
        % Send overdamped modes with delta to Python
        write_progress(sprintf('OVERDAMPED_MODE:%d:%s:%s', k, mode_str, delta_str));
    end
    if length(purelyImaginary) > 20
        fprintf('  ... and %d more overdamped modes\n', length(purelyImaginary) - 20);
        write_progress(sprintf('OVERDAMPED_MORE:%d', length(purelyImaginary) - 20));
    end
end

fprintf('\nNote: %d total common QNMs found (including negative real part counterparts)\n', length(commonQNMs));

% Check for instabilities among oscillatory and overdamped modes
oscillatoryAndOverdamped = [];
for j = 1:length(commonQNMs)
    mode = commonQNMs(j);
    if abs(real(mode)) < imagThreshold || real(mode) > 0
        oscillatoryAndOverdamped = [oscillatoryAndOverdamped; mode];
    end
end

hasPosImag = any(imag(oscillatoryAndOverdamped) > 0);
unstableModes = [];
if hasPosImag
    for j = 1:length(oscillatoryAndOverdamped)
        if imag(oscillatoryAndOverdamped(j)) > 0
            unstableModes = [unstableModes; oscillatoryAndOverdamped(j)];
        end
    end
    if ~isempty(unstableModes)
        [~, sortIdx] = sort(imag(unstableModes), 'descend');
        unstableModes = unstableModes(sortIdx);
    end
end

% Report instability status
if hasPosImag
    fprintf('\n!!! INSTABILITY DETECTED !!!\n');
    fprintf('%d unstable mode(s) found:\n', length(unstableModes));
    write_progress(sprintf('INSTABILITY_DETECTED:%d', length(unstableModes)));
    for k = 1:length(unstableModes)
        mode_str = sprintf('%.12f%+.12f*I', real(unstableModes(k)), imag(unstableModes(k)));
        fprintf('  Unstable mode %d: %s\n', k, mode_str);
        write_progress(sprintf('UNSTABLE_MODE:%d:%s', k, mode_str));
    end
else
    fprintf('\nNo instabilities detected.\n');
    write_progress('NO_INSTABILITIES');
end

analysis_time = toc(analysis_timer);

% Generate report
reportFile = sprintf('results/report_QNM_s%d_L%d_g%s_%s.txt', s_maple, L_maple, strrep(num2str(g_maple), '/', '_'), datestr(now, 'yyyymmdd_HHMMSS'));

% Use writeToReport function to generate report in old format
generate_detailed_report(reportFile, s_maple, L_maple, g_maple, resolutions, ...
                        commonQNMs, positiveReal, purelyImaginary, oscillatoryAndOverdamped, ...
                        unstableModes, imagThreshold, tol);

fprintf('\nDetailed report saved to: %s\n', reportFile);

% Plotting
plotting_timer = tic;
fprintf('\n========================================\n');
fprintf('Generating Visualization\n');
fprintf('========================================\n');

figure('Position', [100, 100, 800, 600]);
set(gcf, 'Color', 'w');
hold on; grid on;

% Plot QNMs for each resolution
markerStyles = {'o', 's', '^'};
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880]};

for idx = 1:L
    e = qnms{idx};
    plot(real(e), imag(e), markerStyles{idx}, ...
         'MarkerFaceColor', colors{idx}, ...
         'MarkerEdgeColor', 'k', ...
         'MarkerSize', 8, ...
         'LineStyle', 'none');
end

% Highlight common QNMs
if ~isempty(commonQNMs)
    lightBlue = [20, 128, 175]/255.0;
    plot(real(commonQNMs), imag(commonQNMs), 'o', ...
         'MarkerSize', 10, ...
         'LineWidth', 2, ...
         'MarkerEdgeColor', lightBlue, ...
         'MarkerFaceColor', 'none');
end

% Add legend and labels
legend(legs{:}, 'Location', 'best', 'Box', 'off');
xlabel('$\mathrm{Re}\,\omega$', 'Interpreter', 'LaTeX', 'FontSize', 14);
ylabel('$\mathrm{Im}\,\omega$', 'Interpreter', 'LaTeX', 'FontSize', 14);
title_text = sprintf('QNM Spectrum: $s=%d$, $L=%d$', s_maple, L_maple);
title(title_text, 'Interpreter', 'LaTeX', 'FontSize', 14);

hold off;

% Save plots
if ~exist('shots', 'dir')
    mkdir('shots');
end

fprintf('Saving plots...\n');
fprintf('  Saving MATLAB figure (.fig)...');
savefig('shots/chaseQNMs.fig');
fprintf(' done\n');

fprintf('  Saving PDF (vector format)...');
exportgraphics(gcf, 'shots/chaseQNMs.pdf', 'BackgroundColor', 'white', 'ContentType', 'vector');
fprintf(' done\n');

fprintf('  Saving PNG (for viewing)...');
exportgraphics(gcf, 'shots/chaseQNMs.png', 'BackgroundColor', 'white');
fprintf(' done\n');

fprintf('  Saving EPS (vector format)...');
exportgraphics(gcf, 'shots/chaseQNMs.eps', 'BackgroundColor', 'white', 'ContentType', 'vector');
fprintf(' done\n');

plotting_time = toc(plotting_timer);
fprintf('All plots saved to shots/ directory\n');

% Final timing summary
total_time = toc(overall_timer);
fprintf('\n========================================\n');
fprintf('Computation Complete!\n');
fprintf('========================================\n');
fprintf('Resolution computations: %s\n', format_time_matlab(sum(resolution_times)));
fprintf('Analysis phase: %s\n', format_time_matlab(analysis_time));
fprintf('Plotting: %s\n', format_time_matlab(plotting_time));
fprintf('Total time: %s\n', format_time_matlab(total_time));
fprintf('========================================\n');

% Save timing data and mode values for Python to read
fprintf('\nSaving timing and analysis data for Python...\n');
timing_data = struct(...
    'total_time', total_time, ...
    'resolution_times', resolution_times, ...
    'analysis_time', analysis_time, ...
    'plotting_time', plotting_time, ...
    'num_common_qnms', length(commonQNMs), ...
    'num_oscillatory', length(positiveReal), ...
    'num_overdamped', length(purelyImaginary), ...
    'num_unstable_modes', length(unstableModes), ...
    'has_instability', hasPosImag, ...
    'jobName', jobName ...
);

% Save mode values for Python (convert mp to double for compatibility)
if ~isempty(positiveReal)
    timing_data.oscillatory_modes = arrayfun(@(x) complex(double(real(x)), double(imag(x))), positiveReal);
else
    timing_data.oscillatory_modes = [];
end

if ~isempty(purelyImaginary)
    % Save first 20 overdamped modes
    numToSave = min(20, length(purelyImaginary));
    timing_data.overdamped_modes = arrayfun(@(x) complex(double(real(x)), double(imag(x))), purelyImaginary(1:numToSave));
    timing_data.total_overdamped = length(purelyImaginary);
else
    timing_data.overdamped_modes = [];
    timing_data.total_overdamped = 0;
end

if ~isempty(unstableModes)
    timing_data.unstable_modes = arrayfun(@(x) complex(double(real(x)), double(imag(x))), unstableModes);
else
    timing_data.unstable_modes = [];
end

% Save using -v7 format for better compatibility with scipy
fprintf('  Saving to results/timing_data.mat...');
save('results/timing_data.mat', 'timing_data', '-v7');
fprintf(' done\n');

% Write final status for Python
write_progress('ANALYSIS_COMPLETE');

% Success indicator for Python
fprintf('\n========================================\n');
fprintf('ALL COMPUTATIONS COMPLETED SUCCESSFULLY\n');
fprintf('========================================\n');
fprintf('[MATLAB_SUCCESS]\n');

% Helper function to generate detailed report matching old format
function generate_detailed_report(filename, s, L, g, resolutions, commonQNMs, positiveReal, purelyImaginary, oscillatoryAndOverdamped, unstableModes, imagThreshold, tol)
    % Write header to the report
    currentTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    writeToReport(filename, sprintf('QNM Analysis Report - %s', currentTime));
    writeToReport(filename, sprintf('=============================================='));
    writeToReport(filename, sprintf('Parameters:'));
    writeToReport(filename, sprintf('  s = %d, L = %d, g = %s', s, L, num2str(g)));
    writeToReport(filename, sprintf('  resolutions = [%s]', num2str(resolutions)));
    
    jobName = sprintf('QNM_s%d_L%d_g%s', s, L, strrep(num2str(g), '/', '_'));
    writeToReport(filename, sprintf('Job Name: %s', jobName));
    writeToReport(filename, sprintf('==============================================\n'));
    
    % Write QNM findings to the report
    writeToReport(filename, sprintf('QNM Analysis Results:'));
    writeToReport(filename, sprintf('------------------------'));
    
    if isempty(commonQNMs)
        writeToReport(filename, sprintf('No common QNMs found within tolerance %s.', sprintf('%.8g', tol)));
    else
        uniqueModeCount = length(positiveReal) + length(purelyImaginary);
        
        writeToReport(filename, sprintf('Common QNMs (found in all resolutions) within tolerance %s:', sprintf('%.8g', tol)));
        writeToReport(filename, sprintf('Found %d modes (%d oscillatory, %d overdamped)', ...
                     uniqueModeCount, length(positiveReal), length(purelyImaginary)));
        writeToReport(filename, sprintf('Sorted by decreasing imaginary part\n'));
        
        % Keep track of the overall mode index for the report
        report_mode_idx_counter = 0;
        
        if ~isempty(positiveReal)
            writeToReport(filename, 'Oscillatory modes (positive real part):');
            for k_report = 1:length(positiveReal)
                report_mode_idx_counter = report_mode_idx_counter + 1;
                modeText_report = sprintf('  Mode #%-3d: %.7f%+.7f*I', report_mode_idx_counter, real(positiveReal(k_report)), imag(positiveReal(k_report)));
                writeToReport(filename, modeText_report);
            end
        end
        
        % Add empty line if both sections have data
        if ~isempty(positiveReal) && ~isempty(purelyImaginary)
            writeToReport(filename, ''); % Empty line for separation
        end
        
        if ~isempty(purelyImaginary)
            writeToReport(filename, 'Overdamped modes (purely imaginary):');
            % Reset counter for purely imaginary modes to start from 1
            imag_mode_counter = 0;
            for k_report = 1:length(purelyImaginary)
                imag_mode_counter = imag_mode_counter + 1;
                report_mode_idx_counter = report_mode_idx_counter + 1;
                current_imag_report = imag(purelyImaginary(k_report));
                delta_str_report = '      N/A'; % Default for the last mode
                
                if k_report < length(purelyImaginary)
                    next_imag_report = imag(purelyImaginary(k_report+1));
                    delta_val_report = abs(current_imag_report - next_imag_report);
                    delta_str_report = sprintf('%9.7f', delta_val_report);
                end
                
                modeText_report = sprintf('  Mode #%-3d: %+11.7f*I   Delta: %s', imag_mode_counter, current_imag_report, delta_str_report);
                writeToReport(filename, modeText_report);
            end
        end
        writeToReport(filename, sprintf('\nNote: %d total common QNMs found (including negative real part counterparts)', length(commonQNMs)));
    end
    
    % Write instability information
    writeToReport(filename, sprintf('\nInstability Analysis:'));
    writeToReport(filename, sprintf('------------------------'));
    if ~isempty(unstableModes)
        writeToReport(filename, sprintf('ALERT! %d unstable mode(s) detected with Im(ω) >= %g (sorted by decreasing Im part):', length(unstableModes), tol));
        for k_report = 1:length(unstableModes)
            modeText_report = sprintf('  Unstable mode: %.7f%+.7f*I', real(unstableModes(k_report)), imag(unstableModes(k_report)));
            writeToReport(filename, modeText_report);
        end
    else
        writeToReport(filename, sprintf('No instabilities (Im(ω) >= %g) have been detected.', tol));
    end
    
    % Write summary for the report file
    writeToReport(filename, sprintf('\nSummary:'));
    writeToReport(filename, sprintf('------------------------'));
    if ~isempty(unstableModes)
        report_instabilityStatus = sprintf('INSTABILITY DETECTED! (%d modes with Im(ω) >= %g)', length(unstableModes), tol);
    else
        report_instabilityStatus = sprintf('No instabilities (Im(ω) >= %g) found', tol);
    end
    writeToReport(filename, sprintf('- Found %d total QNMs', length(commonQNMs)));
    writeToReport(filename, sprintf('- %d modes (%d oscillatory, %d overdamped)', ...
                 length(oscillatoryAndOverdamped), sum(real(oscillatoryAndOverdamped) > imagThreshold), ...
                 sum(abs(real(oscillatoryAndOverdamped)) < imagThreshold)));
    writeToReport(filename, sprintf('- %s', report_instabilityStatus));
    writeToReport(filename, sprintf('\nEnd of Report'));
    writeToReport(filename, sprintf('=============================================='));
end

% Function to append text to the report file
function writeToReport(filename, text)
    % Ensure the directory exists
    [directory, ~, ~] = fileparts(filename);
    if ~isfolder(directory)
        mkdir(directory);
    end
    
    % Open file for appending (or create if it doesn't exist)
    fid = fopen(filename, 'a');
    if fid == -1
        warning('Unable to open report file for writing: %s', filename);
        return;
    end
    
    % Write the text followed by a newline
    fprintf(fid, '%s\n', text);
    
    % Close the file
    fclose(fid);
end

% Helper function to format time in human-readable format
function time_str = format_time_matlab(seconds)
    if seconds < 60
        time_str = sprintf('%.2f seconds', seconds);
    else
        hours = floor(seconds / 3600);
        remainder = mod(seconds, 3600);
        minutes = floor(remainder / 60);
        secs = mod(remainder, 60);
        
        if hours > 0
            time_str = sprintf('%dh %dm %.1fs', hours, minutes, secs);
        else
            time_str = sprintf('%dm %.1fs', minutes, secs);
        end
    end
end