#!/bin/bash
# file_utils.sh - Utility functions for file operations with retry logic and validation
# 
# Functions:
#   get_job_id: Generate a unique job ID (SLURM-aware)
#   copy_with_retry: Copy file with md5sum verification and exponential backoff
#   verify_gzipped: Validate that a file is properly gzipped

set -euo pipefail

# Function: get_job_id
# Description: Generate a unique job identifier, using SLURM_JOB_ID if available,
#              otherwise generating a unique hash based on PID, timestamp, and random data
# Arguments: None
# Returns: Prints a unique job ID to stdout
# Usage:
#   JOB_ID=$(get_job_id)
get_job_id() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        # Running under SLURM - use the job ID
        echo "$SLURM_JOB_ID"
    else
        # Not running under SLURM - generate a unique hash
        # Combine PID, timestamp, hostname, and random data for uniqueness
        local UNIQUE_STRING="$$_$(date +%s%N)_$(hostname)_$RANDOM$RANDOM"
        # Generate SHA256 hash and take first 16 characters
        echo -n "$UNIQUE_STRING" | shasum -a 256 | cut -c1-16
    fi
}

# Function: copy_with_retry
# Description: Copy a file with retry logic, md5sum verification, and exponential backoff
# Arguments:
#   --source: Source file path (required)
#   --dest: Destination file path (required)
#   --description: Human-readable description for logging (required)
# Returns:
#   0 on success, 1 on failure after all retries
# Usage:
#   copy_with_retry --source /path/to/source --dest /path/to/dest --description "my file"
copy_with_retry() {
    local SOURCE=""
    local DEST=""
    local DESCRIPTION=""
    
    # Parse command-line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --source)
                SOURCE="$2"
                shift 2
                ;;
            --dest)
                DEST="$2"
                shift 2
                ;;
            --description)
                DESCRIPTION="$2"
                shift 2
                ;;
            *)
                echo "ERROR: Unknown argument: $1" >&2
                echo "Usage: copy_with_retry --source FILE --dest FILE --description DESC" >&2
                return 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "$SOURCE" ]]; then
        echo "ERROR: --source is required" >&2
        return 1
    fi
    if [[ -z "$DEST" ]]; then
        echo "ERROR: --dest is required" >&2
        return 1
    fi
    if [[ -z "$DESCRIPTION" ]]; then
        echo "ERROR: --description is required" >&2
        return 1
    fi
    
    # Verify source file exists
    if [[ ! -f "$SOURCE" ]]; then
        echo "ERROR: Source file does not exist: $SOURCE" >&2
        return 1
    fi
    
    echo "Copying $DESCRIPTION to $DEST with retry logic..."
    echo "  Source: $SOURCE"
    echo "  Destination: $DEST"
    
    local COPY_SUCCESS=0
    for ATTEMPT in {1..10}; do
        echo "[$DESCRIPTION] Copy attempt $ATTEMPT of 10..."
        
        # Copy the file
        if cp "$SOURCE" "$DEST" 2>/dev/null; then
            sync
            
            # Verify with md5sum
            echo "[$DESCRIPTION] Verifying file integrity with md5sum..."
            local SRC_MD5=$(md5sum "$SOURCE" | awk '{print $1}')
            local DST_MD5=$(md5sum "$DEST" | awk '{print $1}')
            
            echo "[$DESCRIPTION]   Source md5: $SRC_MD5"
            echo "[$DESCRIPTION]   Dest md5:   $DST_MD5"
            
            if [[ "$SRC_MD5" == "$DST_MD5" ]]; then
                echo "[$DESCRIPTION] Copy successful on attempt $ATTEMPT (md5: $DST_MD5)"
                COPY_SUCCESS=1
                break
            else
                echo "[$DESCRIPTION] WARNING: md5sum mismatch on attempt $ATTEMPT" >&2
                echo "[$DESCRIPTION]   Source: $SRC_MD5" >&2
                echo "[$DESCRIPTION]   Dest:   $DST_MD5" >&2
                rm -f "$DEST"
            fi
        else
            echo "[$DESCRIPTION] WARNING: Copy failed on attempt $ATTEMPT" >&2
        fi
        
        # Exponential backoff with jitter: 2, 5, 10, 17, 26, 37, 50, 65, 82, 101 seconds
        if [[ $ATTEMPT -lt 10 ]]; then
            local DELAY=$((2 + (ATTEMPT * ATTEMPT)))
            echo "[$DESCRIPTION] Waiting ${DELAY} seconds before retry..."
            sleep $DELAY
        fi
    done
    
    if [[ $COPY_SUCCESS -eq 0 ]]; then
        echo "[$DESCRIPTION] ERROR: Failed to copy file after 10 attempts" >&2
        return 1
    fi
    
    # Verify destination file exists and is readable
    if [[ ! -f "$DEST" ]]; then
        echo "[$DESCRIPTION] ERROR: Destination file not found at $DEST after successful copy" >&2
        return 1
    fi
    
    echo "[$DESCRIPTION] Copy complete and verified"
    return 0
}

# Function: verify_gzipped
# Description: Validate that a file is properly gzipped
# Arguments:
#   --file: File path to validate (required)
#   --description: Human-readable description for logging (required)
# Returns:
#   0 if file is properly gzipped, 1 if not gzipped or validation fails
# Usage:
#   verify_gzipped --file /path/to/file.gz --description "my gzipped file"
verify_gzipped() {
    local FILE=""
    local DESCRIPTION=""
    
    # Parse command-line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --file)
                FILE="$2"
                shift 2
                ;;
            --description)
                DESCRIPTION="$2"
                shift 2
                ;;
            *)
                echo "ERROR: Unknown argument: $1" >&2
                echo "Usage: verify_gzipped --file FILE --description DESC" >&2
                return 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "$FILE" ]]; then
        echo "ERROR: --file is required" >&2
        return 1
    fi
    if [[ -z "$DESCRIPTION" ]]; then
        echo "ERROR: --description is required" >&2
        return 1
    fi
    
    # Verify file exists
    if [[ ! -f "$FILE" ]]; then
        echo "[$DESCRIPTION] ERROR: File does not exist: $FILE" >&2
        return 1
    fi
    
    echo "[$DESCRIPTION] Verifying gzip format for: $FILE"
    
    # Test gzip integrity
    if gzip -t "$FILE" 2>/dev/null; then
        echo "[$DESCRIPTION] ✓ File is properly gzipped"
        return 0
    else
        echo "[$DESCRIPTION] ERROR: File failed gzip integrity test" >&2
        echo "[$DESCRIPTION]   File: $FILE" >&2
        echo "[$DESCRIPTION]   This file must be properly gzipped" >&2
        return 1
    fi
}
