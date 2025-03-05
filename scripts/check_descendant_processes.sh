#!/bin/bash

# Get the current user's username
USER=$(goldma34)

# Get the parent process ID (PID) of the current shell
PARENT_PID=$$

# Function to get all descendant processes of a given PID
get_descendants() {
    local pid=$1
    local children=$(pgrep -P $pid)
    for child in $children; do
        echo $child
        get_descendants $child
    done
}

# Get all descendant processes of the current shell
descendants=$(get_descendants $PARENT_PID)

# Check the status of each descendant process
for pid in $descendants; do
    ps -p $pid -o pid,ppid,cmd,stat
done