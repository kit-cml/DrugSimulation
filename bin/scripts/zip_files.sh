#!/bin/bash

zip_files() {
  local ROOT_FOLDER="$1"
  local FILE_SUBSTRING="$2"
  local OUTPUT_ZIP="$3"

  find "${ROOT_FOLDER}" -type f -name "${FILE_SUBSTRING}" \
  | zip -r "${OUTPUT_ZIP}" -@ > /dev/null

  echo "Zip file ${OUTPUT_ZIP} consists of ${FILE_SUBSTRING} files has been created!!"
}
