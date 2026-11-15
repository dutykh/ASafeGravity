#!/usr/bin/env python3
"""
Utility functions for QNM computation workflow.
"""

# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)

import os
import re
from typing import List, Optional, Tuple

import requests
from dotenv import load_dotenv

# Locations of Maple and of the Advanpix toolbox, used when the file .env does not set them
DEFAULT_MAPLE_EXECUTABLE = "/path/to/maple"
DEFAULT_ADVANPIX_PATH = "/path/to/advanpix/"


def load_config() -> dict:
    """Load configuration from .env file."""
    load_dotenv()

    config = {
        # Telegram
        "BOT_TOKEN": os.getenv("BOT_TOKEN", ""),
        "CHAT_ID": os.getenv("CHAT_ID", ""),

        # Resolutions
        "RESOLUTIONS": [int(x.strip()) for x in os.getenv("RESOLUTIONS", "400,390,380").split(",")],

        # Executables
        "MAPLE_EXECUTABLE": os.getenv("MAPLE_EXECUTABLE", DEFAULT_MAPLE_EXECUTABLE),
        "MATLAB_EXECUTABLE": os.getenv("MATLAB_EXECUTABLE", "matlab"),

        # Timeouts
        "MAPLE_TIMEOUT": int(os.getenv("MAPLE_TIMEOUT", "3600")),
        "MATLAB_TIMEOUT": int(os.getenv("MATLAB_TIMEOUT", "14400")),
        
        # MATLAB Toolbox Paths
        "ADVANPIX_PATH": os.getenv("ADVANPIX_PATH", DEFAULT_ADVANPIX_PATH),
        
        # Parallel Computation Settings
        "MATLAB_PARALLEL": os.getenv("MATLAB_PARALLEL", "false").lower() == "true",
    }

    return config


import time as time_module

# Rate limiting for Telegram messages
_last_telegram_message_time = 0
_telegram_message_delay = 0.5  # Minimum delay between messages (seconds)

def send_telegram_message(message: str, prefix: Optional[str] = None, retry_count: int = 3, log_errors: bool = False) -> bool:
    """
    Send a message to Telegram bot with rate limiting and error handling.

    Args:
        message: The message to send
        prefix: Optional prefix to add to the message
        retry_count: Number of retry attempts for failed messages
        log_errors: Whether to print error messages to console

    Returns:
        True if message was sent successfully, False otherwise
    """
    global _last_telegram_message_time
    
    config = load_config()
    bot_token = config["BOT_TOKEN"]
    chat_id = config["CHAT_ID"]

    if not bot_token or not chat_id:
        # Silently skip if Telegram is not configured
        return False

    # Rate limiting: ensure minimum delay between messages
    current_time = time_module.time()
    time_since_last = current_time - _last_telegram_message_time
    if time_since_last < _telegram_message_delay:
        time_module.sleep(_telegram_message_delay - time_since_last)
    
    _last_telegram_message_time = time_module.time()

    # Add prefix if provided
    if prefix:
        full_message = f"{prefix}: {message}"
    else:
        full_message = message

    # Truncate very long messages to avoid Telegram limits
    if len(full_message) > 4000:
        full_message = full_message[:3997] + "..."

    url = f"https://api.telegram.org/bot{bot_token}/sendMessage"
    data = {
        "chat_id": chat_id,
        "text": full_message,
        "parse_mode": "Markdown"
    }

    for attempt in range(retry_count):
        try:
            response = requests.post(url, data=data, timeout=10)

            if response.status_code == 200:
                return True
            elif response.status_code == 429:  # Rate limited
                # Wait longer for rate limit
                retry_after = int(response.headers.get('Retry-After', 5))
                if log_errors:
                    print(f"Telegram rate limited, waiting {retry_after}s...")
                time_module.sleep(retry_after)
                continue
            else:
                if log_errors:
                    print(f"Telegram error {response.status_code}: {response.text}")
                if attempt < retry_count - 1:
                    time_module.sleep(2 ** attempt)  # Exponential backoff
                    
        except requests.exceptions.Timeout:
            if log_errors:
                print(f"Telegram timeout (attempt {attempt + 1}/{retry_count})")
            if attempt < retry_count - 1:
                time_module.sleep(2 ** attempt)
        except Exception as e:
            if log_errors:
                print(f"Telegram error (attempt {attempt + 1}/{retry_count}): {e}")
            if attempt < retry_count - 1:
                time_module.sleep(2 ** attempt)

    # All attempts failed
    if log_errors:
        print(f"Failed to send Telegram message after {retry_count} attempts: {full_message[:100]}...")
    return False


def send_telegram_photo(photo_path: str, caption: Optional[str] = None) -> None:
    """
    Send a photo to Telegram bot.

    Args:
        photo_path: Path to the photo file
        caption: Optional caption for the photo
    """
    config = load_config()
    bot_token = config["BOT_TOKEN"]
    chat_id = config["CHAT_ID"]

    if not bot_token or not chat_id:
        return

    if not os.path.exists(photo_path):
        # Silently skip if photo doesn't exist
        return

    try:
        url = f"https://api.telegram.org/bot{bot_token}/sendPhoto"

        with open(photo_path, "rb") as photo:
            files = {"photo": photo}
            data = {"chat_id": chat_id}
            if caption:
                data["caption"] = caption

            response = requests.post(url, files=files, data=data, timeout=30)

            if response.status_code != 200:
                # Silently ignore Telegram errors
                pass
    except Exception:
        # Silently ignore Telegram errors
        pass


def natural_sort(strings: List[str]) -> List[str]:
    """
    Sort strings with numbers in natural order (L0 before L10).

    Args:
        strings: List of strings to sort

    Returns:
        Naturally sorted list
    """
    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(strings, key=alphanum_key)


def find_L_subdirectories(root_dir: str) -> List[str]:
    """
    Find all subdirectories that start with 'L' (e.g., L0, L1, L2).

    Args:
        root_dir: Root directory to search in

    Returns:
        List of L directory names, naturally sorted
    """
    subdirs = []

    for item in os.listdir(root_dir):
        item_path = os.path.join(root_dir, item)

        # Only consider directories that match pattern L<number>
        if os.path.isdir(item_path) and re.match(r"^L(\d+)$", item):
            subdirs.append(item)

    return natural_sort(subdirs)


def extract_parameters(directory: str) -> Tuple[int, int]:
    """
    Extract s and L parameters from directory structure.

    Args:
        directory: Path to L* subdirectory

    Returns:
        Tuple of (s, L) values

    Raises:
        ValueError: If parameters cannot be extracted
    """
    # Get directory name (e.g., "L0")
    dir_name = os.path.basename(directory)

    # Get parent directory name (e.g., "s0")
    parent_dir = os.path.basename(os.path.dirname(directory))

    # Extract s value
    s_match = re.match(r"^s(\d+)$", parent_dir)
    if not s_match:
        raise ValueError(f"Could not extract s value from parent directory: {parent_dir}")
    s = int(s_match.group(1))

    # Extract L value
    l_match = re.match(r"^L(\d+)$", dir_name)
    if not l_match:
        raise ValueError(f"Could not extract L value from directory: {dir_name}")
    L = int(l_match.group(1))

    return s, L


def ensure_directories(base_dir: str) -> None:
    """
    Ensure that data/, results/, and shots/ directories exist.

    Args:
        base_dir: Base directory where subdirectories should be created
    """
    subdirs = ["data", "results", "shots"]

    for subdir in subdirs:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)


def format_time(seconds: float) -> str:
    """
    Format time duration in human-readable format.

    Args:
        seconds: Time duration in seconds

    Returns:
        Formatted string like "2h 15m 30.5s"
    """
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    parts = []
    if hours > 0:
        parts.append(f"{int(hours)}h")
    if minutes > 0 or hours > 0:
        parts.append(f"{int(minutes)}m")
    parts.append(f"{seconds:.1f}s")

    return " ".join(parts)


def validate_executables(config: dict) -> bool:
    """
    Validate that required executables exist.

    Args:
        config: Configuration dictionary

    Returns:
        True if all executables are found, False otherwise
    """
    maple_path = config["MAPLE_EXECUTABLE"]
    matlab_cmd = config["MATLAB_EXECUTABLE"]

    # Check Maple
    if not os.path.exists(maple_path):
        print(f"ERROR: Maple executable not found at: {maple_path}")
        return False

    # Check MATLAB (try to find it in PATH)
    import shutil
    if not shutil.which(matlab_cmd):
        print(f"ERROR: MATLAB executable '{matlab_cmd}' not found in PATH")
        return False

    return True
