#!/usr/bin/env bash
set -euo pipefail

create_escaped_b64_file(){
  local INPUT_FILE=$1
  local INPUT_FILE_ESCAPED=$2
  local INPUT_FILE_ESCAPED_B64=$3

  # Make JSON-safe string from the camelized param file.
  local ESCAPED_PARAM_STRING
  ESCAPED_PARAM_STRING=$(jq -Rs . "${INPUT_FILE}")

  # Remove the unnecessary quotes.
  ESCAPED_PARAM_STRING="${ESCAPED_PARAM_STRING:1:-1}"

  # Write as-is (one line) to the output file.
  printf '%s\n' "${ESCAPED_PARAM_STRING}" > "${INPUT_FILE_ESCAPED}"

  # Convert the param as the base64 file.
  base64 -w 0 "${INPUT_FILE_ESCAPED}" > "${INPUT_FILE_ESCAPED_B64}"
}
