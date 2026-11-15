#!/usr/bin/env python3
"""
Unified QNM computation orchestrator.
Runs Maple and MATLAB computations across all L* subdirectories.
"""

# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)

import os
import shutil
import subprocess
import sys
import time
import threading
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import scipy.io as sio
import re

from utils import (
    ensure_directories,
    extract_parameters,
    find_L_subdirectories,
    format_time,
    load_config,
    send_telegram_message,
    send_telegram_photo,
    validate_executables,
)


def find_g_subdirectories(root_dir: str) -> list:
    """Find all gX_Y, g0, and g1 subdirectories in the root directory."""
    g_dirs = []
    for item in os.listdir(root_dir):
        if os.path.isdir(os.path.join(root_dir, item)):
            # Match gX_Y format, g0, or g1
            if re.match(r'^g\d+_\d+$', item) or item == 'g0' or item == 'g1':
                g_dirs.append(item)
    return sorted(g_dirs)


def parse_g_value(g_dir_name: str) -> float:
    """Parse the g value from a directory name like 'gX_Y', 'g0', or 'g1'."""
    if g_dir_name == 'g0':
        return 0.0
    elif g_dir_name == 'g1':
        return 1.0
    else:
        # Handle gX_Y format
        parts = g_dir_name.replace('g', '').split('_')
        return float(parts[0]) / float(parts[1])


def parse_s_L_from_path(cwd: str) -> tuple:
    """Parse s and L from the current working directory path."""
    parts = cwd.split(os.sep)
    s = int(parts[-2].replace('s', ''))
    L = int(parts[-1].replace('L', ''))
    return s, L


def run_single_maple_resolution(args: tuple) -> tuple:
    """
    Run Maple for a single resolution in parallel.

    Args:
        args: Tuple of (working_dir, s, L, g, resolution, timeout, l_dir_name, maple_executable)

    Returns:
        Tuple of (success, resolution, time_taken)
    """
    working_dir, s, L, g, resolution, timeout, l_dir_name, maple_executable = args
    start_time = time.time()

    # The Maple script is kept minimal. It reads the main procedure, sets parameters,
    # and calls the procedure. All logic is in matrixassembler.mpl.
    script_content = f"""# Auto-generated Maple script for resolution {resolution}
# Increase verbosity to trace execution
printlevel := 2:

read("{os.path.abspath('matrixassembler.mpl')}"):
s := {s}: L := {L}: g := {g}: N := {resolution}:
path := "{working_dir}":
MatrixAssembler(N, N, s, L, g, path):
quit
"""
    script_path = os.path.join(working_dir, f"maple_run_{resolution}.mpl")
    success = False
    try:
        with open(script_path, "w") as f:
            f.write(script_content)

        process = subprocess.Popen(
            [maple_executable, script_path],  # Removed '-q' for full output
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=working_dir
        )

        print(f"[Maple N={resolution}] Starting computation...")
        for line in iter(process.stdout.readline, ''):
            print(f"[Maple N={resolution}] {line.strip()}")

        process.wait(timeout=timeout)
        success = process.returncode == 0
        if not success:
            stderr_output = process.stderr.read()
            print(f"ERROR: Maple N={resolution} failed. Stderr:\n{stderr_output}")

    except subprocess.TimeoutExpired:
        process.kill()
        print(f"ERROR: Maple process for N={resolution} timed out after {timeout} seconds.")
    except FileNotFoundError:
        print(f"ERROR: Maple executable not found at '{maple_executable}'. Please check config.json.")
    except Exception as e:
        print(f"ERROR: An exception occurred in Maple N={resolution}: {e}")
    finally:
        if os.path.exists(script_path):
            os.remove(script_path)

    time_taken = time.time() - start_time
    status_char = "✓" if success else "✗"
    send_telegram_message(f"{status_char} Maple N={resolution} finished in {time_taken:.1f}s.", prefix=l_dir_name)
    return success, resolution, time_taken


def run_maple_parallel(working_dir: str, s: int, L: int, g: float, resolutions: list, timeout: int, l_dir_name: str, config: dict) -> tuple:
    """
    Run Maple for all resolutions in parallel using a process pool.
    """
    total_start_time = time.time()
    send_telegram_message(f"Maple: Starting {len(resolutions)} parallel jobs...", prefix=l_dir_name)
    print(f"Running {len(resolutions)} Maple processes in parallel...")
    
    maple_executable = config["MAPLE_EXECUTABLE"]
    process_args = [(working_dir, s, L, g, res, timeout, l_dir_name, maple_executable) for res in resolutions]
    all_success = True

    with ProcessPoolExecutor(max_workers=min(len(resolutions), 4)) as executor:
        future_to_res = {executor.submit(run_single_maple_resolution, args): args[4] for args in process_args}
        
        for future in as_completed(future_to_res):
            try:
                success, _, _ = future.result()
                if not success:
                    all_success = False
            except Exception as exc:
                print(f'A Maple process generated an exception: {exc}')
                all_success = False

    total_time = time.time() - total_start_time
    if all_success:
        msg = f"All Maple jobs completed successfully in {format_time(total_time)}."
        print(msg)
        send_telegram_message(f"✅ {msg}", prefix=l_dir_name)
    else:
        msg = f"One or more Maple jobs failed after {format_time(total_time)}."
        print(msg)
        send_telegram_message(f"❌ {msg}", prefix=l_dir_name)
        
    return all_success, total_time


def prepare_matlab_parameters(working_dir: str, s: int, L: int, g: float, resolutions: list):
    """
    Create a parameter file for MATLAB to read.

    Args:
        working_dir: Directory where MATLAB will run
        s: Spin parameter
        L: Angular momentum parameter
        g: Gamma parameter
        resolutions: List of resolutions
    """
    config = load_config()
    advanpix_path = config.get("ADVANPIX_PATH", "")
    use_parallel = config.get("USE_PARALLEL", False)
    
    params = {
        's_value': s,
        'L_value': L,
        'g_value': g,
        'resolutions': resolutions,
        'job_timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'advanpix_path': advanpix_path,
        'use_parallel': use_parallel
    }


    param_file = os.path.join(working_dir, 'matlab_params.mat')
    sio.savemat(param_file, params)


def send_qnm_analysis_telegram(results: dict, l_dir_name: str):
    """Send comprehensive QNM analysis report to Telegram."""
    # Build the detailed analysis message
    msg = "\n📊 **QNM Analysis Results**\n"
    msg += f"{'='*35}\n"
    msg += f"Common modes found: **{results['total_common']}** (in all resolutions)\n"
    msg += f"  • Oscillatory modes: **{results['total_oscillatory']}**\n"
    msg += f"  • Overdamped modes: **{results['total_overdamped']}**\n"
    msg += f"{'='*35}\n"
    
    # Oscillatory modes
    if results['oscillatory_modes']:
        msg += "\n🔴 **Oscillatory modes (positive real part):**\n"
        for mode_num, mode_value in results['oscillatory_modes']:
            msg += f"  Mode #{mode_num:2d}: `{mode_value}`\n"
    
    # Overdamped modes
    if results['overdamped_modes']:
        msg += "\n🔵 **Overdamped modes (purely imaginary):**\n"
        for mode_num, mode_value, delta in results['overdamped_modes']:
            msg += f"  Mode #{mode_num:2d}: `{mode_value}`  Delta: `{delta.strip()}`\n"
        if results['more_overdamped'] > 0:
            msg += f"  ... and {results['more_overdamped']} more overdamped modes_\n"
    
    msg += f"\n_Note: {results['total_common']} total common QNMs found_\n"
    
    # Instability status
    if results['has_instability']:
        msg += "\n⚠️ **INSTABILITY DETECTED!**\n"
        msg += f"{results['num_unstable']} unstable mode(s) found:\n"
        for mode_num, mode_value in results['unstable_modes']:
            msg += f"  Unstable mode {mode_num}: `{mode_value}`\n"
    else:
        msg += "\n✅ **No instabilities detected.**\n"
    
    # Send the message (split if too long)
    if len(msg) > 4000:  # Telegram message limit
        # Send in parts
        parts = msg.split('\n\n')
        current_msg = ""
        for part in parts:
            if len(current_msg) + len(part) + 2 < 4000:
                current_msg += part + "\n\n"
            else:
                if current_msg:
                    send_telegram_message(current_msg.strip(), prefix=l_dir_name)
                current_msg = part + "\n\n"
        if current_msg:
            send_telegram_message(current_msg.strip(), prefix=l_dir_name)
    else:
        send_telegram_message(msg, prefix=l_dir_name)


def monitor_matlab_progress(working_dir: str, l_dir_name: str):
    """Monitor MATLAB progress by reading progress file with enhanced status tracking."""
    progress_file = os.path.join(working_dir, 'matlab_progress.txt')
    
    # Storage for QNM analysis results
    analysis_results = {
        'oscillatory_modes': [],
        'overdamped_modes': [],
        'unstable_modes': [],
        'total_common': 0,
        'total_oscillatory': 0,
        'total_overdamped': 0,
        'more_overdamped': 0,
        'has_instability': False,
        'num_unstable': 0
    }
    
    while True:
        if os.path.exists(progress_file):
            try:
                with open(progress_file, 'r') as f:
                    messages = f.readlines()
                
                # Process all messages in the file
                for message_line in messages:
                    message = message_line.strip()
                    if not message:  # Skip empty lines
                        continue
                    
                    if message.startswith('RESOLUTION_START:'):
                        parts = message.split(':', 2)  # Split only on first 2 colons to preserve timestamp
                        if len(parts) >= 3:
                            resolution = parts[1]
                            timestamp = parts[2]
                            print(f"[{timestamp}] MATLAB: Starting eigenvalue computation for resolution {resolution}")
                            send_telegram_message(f"MATLAB: Starting eigenvalue computation for resolution {resolution}", prefix=l_dir_name)
                            
                    elif message.startswith('RESOLUTION_UPDATE:'):
                        parts = message.split(':', 3)  # Split only on first 3 colons to preserve timestamp
                        if len(parts) >= 4:
                            resolution = parts[1]
                            status = parts[2]
                            timestamp = parts[3]
                            status_msg = status.replace('_', ' ').title()
                            print(f"[{timestamp}] MATLAB: Resolution {resolution} - {status_msg}")
                            # Send status updates to Telegram
                            send_telegram_message(f"MATLAB: Resolution {resolution} - {status_msg}", prefix=l_dir_name)
                            
                    elif message.startswith('RESOLUTION_COMPLETE:'):
                        parts = message.split(':', 3)  # Split only on first 3 colons to preserve timestamp
                        if len(parts) >= 4:
                            resolution = parts[1]
                            time_taken_str = parts[2]
                            timestamp = parts[3]
                            
                            # Convert time to human-readable format if > 60 seconds
                            try:
                                time_taken = float(time_taken_str)
                                if time_taken >= 60:
                                    time_formatted = format_time(time_taken)
                                else:
                                    time_formatted = f"{time_taken:.2f}s"
                            except ValueError:
                                time_formatted = f"{time_taken_str}s"
                            
                            print(f"[{timestamp}] MATLAB: Resolution {resolution} completed in {time_formatted}")
                            send_telegram_message(f"MATLAB: Resolution {resolution} completed ({time_formatted})", prefix=l_dir_name)
                            
                    elif message == 'MATLAB_STARTED':
                        print("MATLAB: Started sequential computation with detailed progress")
                        send_telegram_message("MATLAB: Started sequential computation with detailed progress", prefix=l_dir_name)
                        
                    elif message == 'ANALYSIS_START':
                        print("MATLAB: Starting QNM analysis")
                        send_telegram_message("MATLAB: Starting QNM analysis...", prefix=l_dir_name)
                        
                    elif message.startswith('ANALYSIS_SUMMARY:'):
                        print(f"Debug - Received ANALYSIS_SUMMARY: {message}")
                        parts = message.split(':')
                        if len(parts) >= 4:
                            analysis_results['total_common'] = int(parts[1])
                            analysis_results['total_oscillatory'] = int(parts[2])
                            analysis_results['total_overdamped'] = int(parts[3])
                            print(f"Debug - Parsed summary: common={parts[1]}, oscillatory={parts[2]}, overdamped={parts[3]}")
                            
                    elif message.startswith('OSCILLATORY_MODE:'):
                        parts = message.split(':', 2)
                        if len(parts) >= 3:
                            mode_num = int(parts[1])
                            mode_value = parts[2]
                            analysis_results['oscillatory_modes'].append((mode_num, mode_value))
                        
                    elif message.startswith('OVERDAMPED_MODE:'):
                        parts = message.split(':', 3)
                        if len(parts) >= 4:
                            mode_num = int(parts[1])
                            mode_value = parts[2]
                            delta = parts[3]
                            analysis_results['overdamped_modes'].append((mode_num, mode_value, delta))
                            
                    elif message.startswith('OVERDAMPED_MORE:'):
                        parts = message.split(':')
                        if len(parts) >= 2:
                            analysis_results['more_overdamped'] = int(parts[1])
                            
                    elif message.startswith('INSTABILITY_DETECTED:'):
                        parts = message.split(':')
                        if len(parts) >= 2:
                            analysis_results['has_instability'] = True
                            analysis_results['num_unstable'] = int(parts[1])
                            
                    elif message.startswith('UNSTABLE_MODE:'):
                        parts = message.split(':', 2)
                        if len(parts) >= 3:
                            mode_num = int(parts[1])
                            mode_value = parts[2]
                            analysis_results['unstable_modes'].append((mode_num, mode_value))
                            
                    elif message == 'NO_INSTABILITIES':
                        analysis_results['has_instability'] = False
                    
                    elif message == 'ANALYSIS_COMPLETE':
                        print("MATLAB: Analysis phase completed")
                        # Debug: print analysis results before sending
                        print(f"Debug - Analysis results before sending to Telegram:")
                        print(f"  Total common: {analysis_results['total_common']}")
                        print(f"  Total oscillatory: {analysis_results['total_oscillatory']}")
                        print(f"  Total overdamped: {analysis_results['total_overdamped']}")
                        print(f"  Oscillatory modes: {len(analysis_results['oscillatory_modes'])}")
                        print(f"  Overdamped modes: {len(analysis_results['overdamped_modes'])}")
                        # Send comprehensive analysis report to Telegram
                        send_qnm_analysis_telegram(analysis_results, l_dir_name)
                        # Clear the progress file after analysis complete
                        os.remove(progress_file)
                        return  # Exit the monitoring function
                
                # After processing all messages, clear the file
                os.remove(progress_file)
            except (OSError, IOError):
                pass
        
        time.sleep(0.5)  # Check every 500ms


def run_matlab_script(working_dir: str, timeout: int, g_dir_name: str) -> bool:
    """
    Run the MATLAB script (chaseQNMs.m) with progress monitoring.

    Args:
        working_dir: Directory containing the MATLAB script
        timeout: Timeout in seconds
        g_dir_name: Directory name for Telegram notifications

    Returns:
        True if successful, False otherwise
    """
    config = load_config()
    matlab_executable = config["MATLAB_EXECUTABLE"]

    script_name = "chaseQNMs.m"
    script_path = os.path.join(working_dir, script_name)

    if not os.path.exists(script_path):
        print(f"ERROR: MATLAB script not found: {script_path}")
        return False

    print(f"Running MATLAB in {working_dir}...")
    send_telegram_message("Starting MATLAB computation with detailed progress tracking...", prefix=g_dir_name)

    try:
        # MATLAB command to run the script
        matlab_command = (
            f"try; cd('{working_dir}'); run('{os.path.splitext(script_name)[0]}'); "
            f"catch e; disp(getReport(e)); exit(1); end; exit(0);"
        )

        cmd = [matlab_executable, "-nodisplay", "-nosplash", "-batch", matlab_command]

        # Start MATLAB as a background process
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=working_dir
        )

        # Start monitoring thread
        monitor_thread = threading.Thread(target=monitor_matlab_progress, args=(working_dir, g_dir_name))
        monitor_thread.daemon = True
        monitor_thread.start()

        # Wait for MATLAB to complete
        try:
            stdout, stderr = process.communicate(timeout=timeout)
            
            # Print MATLAB output to console
            if stdout:
                print("MATLAB stdout:")
                print(stdout)
            if stderr:
                print("MATLAB stderr:")
                print(stderr)
                
            # Wait for monitoring thread to finish (it should exit when it sees ANALYSIS_COMPLETE)
            monitor_thread.join(timeout=30)  # Give it 30 seconds to finish
            
            if process.returncode == 0:
                print("MATLAB computation completed successfully.")
                send_telegram_message("✅ MATLAB computation completed successfully!", prefix=g_dir_name)
                return True
            else:
                print(f"ERROR: MATLAB script failed with return code {process.returncode}")
                if stderr:
                    send_telegram_message(f"❌ MATLAB failed: {stderr[:500]}", prefix=g_dir_name)
                return False
                
        except subprocess.TimeoutExpired:
            process.kill()
            print(f"ERROR: MATLAB script timed out after {timeout} seconds")
            send_telegram_message(f"⏰ MATLAB computation timed out after {timeout}s", prefix=g_dir_name)
            return False

    except FileNotFoundError:
        print(f"ERROR: MATLAB executable not found: {matlab_executable}")
        send_telegram_message(f"❌ MATLAB executable not found: {matlab_executable}", prefix=g_dir_name)
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error running MATLAB: {e}")
        send_telegram_message(f"❌ Unexpected MATLAB error: {str(e)}", prefix=g_dir_name)
        return False


def deploy_matlab_script(l_dir_path: str):
    """
    Deploy MATLAB script from template to the L directory.
    
    Args:
        l_dir_path: Path to the L* directory
        
    Returns:
        True if successful, False otherwise
    """
    template_path = os.path.join(os.path.dirname(__file__), 'chaseQNMs_template.m')
    dest_path = os.path.join(l_dir_path, 'chaseQNMs.m')
    
    if not os.path.exists(template_path):
        print(f"ERROR: MATLAB template not found at {template_path}")
        return False
        
    try:
        shutil.copy(template_path, dest_path)
        return True
    except Exception as e:
        print(f"ERROR: Failed to deploy MATLAB script: {e}")
        return False


def process_g_directory(g_dir_path: str, s: int, L: int, g: float, config: dict):
    """
    Process a single g directory: clear old data, deploy scripts, run Maple, run MATLAB.

    Args:
        g_dir_path: Path to the gX_Y directory
        s: Spin parameter
        L: Angular momentum parameter
        g: Gamma parameter
        config: Configuration dictionary

    Returns:
        Tuple of (success, timing_dict) where timing_dict contains timing info
    """
    start_time = time.time()
    g_dir_name = os.path.basename(g_dir_path)

    # --- Setup Phase: Create directories ---
    print(f"Creating subdirectories in {g_dir_name}...")
    for subdir in ["data", "results", "shots"]:
        subdir_path = os.path.join(g_dir_path, subdir)
        os.makedirs(subdir_path, exist_ok=True)

    # Deploy MATLAB script
    if not deploy_matlab_script(g_dir_path):
        send_telegram_message("CRITICAL: Failed to deploy MATLAB script.", prefix=g_dir_name)
        return False, {}

    # --- Maple Computation ---
    maple_timeout = config.get("MAPLE_TIMEOUT", 3600)

    # Run Maple computations in parallel
    maple_success, maple_time = run_maple_parallel(
        g_dir_path, s, L, g, config['RESOLUTIONS'], maple_timeout, g_dir_name, config
    )

    if not maple_success:
        err_msg = f"Maple computation failed after {format_time(maple_time)}."
        print(err_msg)
        send_telegram_message(f"ERROR: {err_msg}", prefix=g_dir_name)
        return False, {"maple": maple_time, "matlab": 0, "total": maple_time}

    # --- MATLAB Computation ---
    matlab_timeout = config.get("MATLAB_TIMEOUT", 3600)

    prepare_matlab_parameters(g_dir_path, s, L, g, config['RESOLUTIONS'])
    
    matlab_start_time = time.time()
    
    matlab_success = run_matlab_script(g_dir_path, matlab_timeout, g_dir_name)
    matlab_time = time.time() - matlab_start_time

    if not matlab_success:
        err_msg = f"MATLAB computation failed after {format_time(matlab_time)}."
        print(err_msg)
        send_telegram_message(f"ERROR: {err_msg}", prefix=g_dir_name)
        return False, {"maple": maple_time, "matlab": matlab_time, "total": time.time() - start_time}

    # --- Analysis and Reporting ---
    # The rest of the function for reporting would go here

    elapsed_time = time.time() - start_time
    completion_msg = f"✅ Completed in {format_time(elapsed_time)} (Maple: {format_time(maple_time)}, MATLAB: {format_time(matlab_time)})"
    send_telegram_message(completion_msg, prefix=g_dir_name)

    return True, {"maple": maple_time, "matlab": matlab_time, "total": elapsed_time}


def clean_all_g_directories(root_dir: str, g_dirs: list) -> None:
    """
    Clean all g* subdirectories at the start of execution.
    
    Args:
        root_dir: Root directory containing g* subdirectories
        g_dirs: List of g directory names to clean
    """
    print("\n--- Initial Cleanup Phase ---")
    print(f"Cleaning {len(g_dirs)} g* directories...")
    
    for g_dir in g_dirs:
        g_dir_path = os.path.join(root_dir, g_dir)
        subdirs_to_clean = ["data", "results", "shots"]
        
        for subdir in subdirs_to_clean:
            subdir_path = os.path.join(g_dir_path, subdir)
            if os.path.exists(subdir_path):
                shutil.rmtree(subdir_path)
                print(f"  Removed {g_dir}/{subdir}/")
        
        # Also remove any stale MATLAB files
        matlab_files = ["chaseQNMs.m", "matlab_params.mat", "matlab_progress.txt"]
        for file in matlab_files:
            file_path = os.path.join(g_dir_path, file)
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"  Removed {g_dir}/{file}")
    
    print("Initial cleanup completed.\n")


def main():
    """Main orchestrator function."""
    total_start_time = time.time()

    # Load configuration
    config = load_config()

    # Validate executables
    if not validate_executables(config):
        sys.exit(1)

    # Get current directory and parse s, L
    root_dir = os.getcwd()
    try:
        s, L = parse_s_L_from_path(root_dir)
    except (IndexError, ValueError) as e:
        print(f"Error: Could not parse s and L from path: {root_dir}")
        print(f"Please run from a directory like '.../sX/LY/'. Error: {e}")
        sys.exit(1)

    print("QNM Computation Orchestrator")
    print(f"Working directory: {root_dir}")
    print(f"Parameters: s = {s}, L = {L}")
    print(f"Resolutions: {config['RESOLUTIONS']}")

    # Find all g* subdirectories
    g_dirs = find_g_subdirectories(root_dir)

    if not g_dirs:
        message = "No gX_Y subdirectories found"
        print(message)
        send_telegram_message(message, prefix="Orchestrator")
        return

    # Clean all g directories before starting
    clean_all_g_directories(root_dir, g_dirs)

    # Send start message
    dir_list = ", ".join(g_dirs)
    start_msg = f"Starting computation for s={s}, L={L} in {len(g_dirs)} directories: {dir_list}"
    print(f"{start_msg}")
    send_telegram_message(start_msg, prefix="Orchestrator", log_errors=True)

    # Process each directory
    successful = 0
    failed = 0
    all_timings = []

    for i, g_dir in enumerate(g_dirs, 1):
        g_dir_path = os.path.join(root_dir, g_dir)
        g_value = parse_g_value(g_dir)

        progress_msg = f"[{i}/{len(g_dirs)}] Processing {g_dir} (g={g_value:.4f})"
        print(f"\n{progress_msg}")
        send_telegram_message(progress_msg, prefix="Orchestrator", log_errors=True)

        success, timing_info = process_g_directory(g_dir_path, s, L, g_value, config)

        if success:
            successful += 1
            timing_info['directory'] = g_dir
            all_timings.append(timing_info)
        else:
            failed += 1

    # Final report with statistics
    total_time = time.time() - total_start_time
    
    # Calculate timing statistics
    if all_timings:
        total_maple_time = sum(t['maple'] for t in all_timings)
        total_matlab_time = sum(t['matlab'] for t in all_timings)
        avg_maple_time = total_maple_time / len(all_timings)
        avg_matlab_time = total_matlab_time / len(all_timings)
    else:
        total_maple_time = total_matlab_time = avg_maple_time = avg_matlab_time = 0

    # Build detailed final report
    print(f"\n{'='*70}")
    print("QNM COMPUTATION SUMMARY")
    print(f"{'='*70}")
    print(f"Parameters: s={s}, L={L}")
    print(f"Total directories processed: {len(g_dirs)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print("\nTiming Statistics:")
    print(f"  Total computation time: {format_time(total_time)}")
    if all_timings:
        print(f"  Total Maple time: {format_time(total_maple_time)}")
        print(f"  Total MATLAB time: {format_time(total_matlab_time)}")
        print("  Average per directory:")
        print(f"    Maple: {format_time(avg_maple_time)}")
        print(f"    MATLAB: {format_time(avg_matlab_time)}")
    
    print(f"{'='*70}")

    # Telegram summary message
    final_msg = (
        f"🎯 QNM Computation Complete for s={s}, L={L}!\n\n"
        f"📊 Summary:\n"
        f"• Processed: {len(g_dirs)} directories\n"
        f"• Successful: {successful} ✅\n"
        f"• Failed: {failed} ❌\n"
        f"• Total time: {format_time(total_time)}\n\n"
        f"⏱️ Timing breakdown:\n"
        f"• Maple total: {format_time(total_maple_time)}\n"
        f"• MATLAB total: {format_time(total_matlab_time)}\n"
    )
    
    send_telegram_message(final_msg, prefix="Orchestrator", log_errors=True)


if __name__ == "__main__":
    # Check for required dependencies
    try:
        import importlib.util
        if importlib.util.find_spec("dotenv") is None:
            raise ImportError("python-dotenv")
        if importlib.util.find_spec("requests") is None:
            raise ImportError("requests")
    except ImportError as e:
        print(f"ERROR: Missing required dependency: {e}")
        print("Please install required packages: pip install python-dotenv requests")
        sys.exit(1)

    main()
