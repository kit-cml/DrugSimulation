#!/bin/bash

watch_progress_percentage() {
  local log_file="$1"
  local output_file="$2"

  mkdir -p "$(dirname "$output_file")"
  touch "$log_file"

  tail -n 0 -F "$log_file" | while IFS= read -r line; do
    if [[ "$line" =~ DrugSimulation[[:space:]]+Global[[:space:]]+progress:[[:space:]]+([0-9]+)% ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}" > "${output_file}.tmp"
      mv "${output_file}.tmp" "$output_file"
    fi
  done
}

emit_progress_percentage() {
  local percent="$1"
  local log_file="$2"

  if [[ -z "$percent" || -z "$log_file" ]]; then
    return 1
  fi

  if (( percent < 0 )); then
    percent=0
  elif (( percent > 100 )); then
    percent=100
  fi

  printf 'DrugSimulation Global progress: %d%%\n' "$percent" >> "$log_file"
}
