#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 [--escaped] <input_file> <output_file>" >&2
    exit 1
}

escaped=0

# Optional flag
if [[ $# -ge 1 && "$1" == "--escaped" ]]; then
    escaped=1
    shift
fi

if [ $# -ne 2 ]; then
    usage
fi

input_file="$1"
output_file="$2"

if [ ! -f "$input_file" ]; then
    echo "Input file not found: $input_file" >&2
    exit 1
fi

# snake_case -> camelCase
to_camel() {
    local word="$1"
    echo "$word" | awk -F'_' '{
        for (i=2; i<=NF; i++) $i = toupper(substr($i,1,1)) substr($i,2);
        for (i=1; i<=NF; i++) printf "%s", $i;
        printf "\n";
    }'
}

# Produce converted content (multi-line, no escaping)
convert_multiline() {
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^[[:space:]]*# || -z "$line" ]]; then
            printf '%s\n' "$line"
        elif [[ "$line" =~ ^[[:space:]]*[a-z_]+[[:space:]]*= ]]; then
            key=$(echo "$line" | sed 's/^[[:space:]]*//' | cut -d'=' -f1 | sed 's/[[:space:]]*$//')
            camel_key=$(to_camel "$key")
            # No leading spaces, normalize to "camelKey = "
            new_line=$(echo "$line" | sed "s/^[[:space:]]*$key[[:space:]]*=/$(printf '%s' "$camel_key") = /")
            printf '%s\n' "$new_line"
        else
            printf '%s\n' "$line"
        fi
    done < "$input_file"
}

if (( escaped )); then
    converted="$(convert_multiline)"
    qout=$(printf '%q\n' "$converted")
    if [[ "$qout" == \$\'*\' ]]; then
        qout=${qout:2}
        qout=${qout::-1}
    fi
    printf '%s\n' "$qout" > "$output_file"
else
    convert_multiline > "$output_file"
fi
