#!/bin/bash
for i in "$@"; do
  CONTENT_LENGTH=`ls -la "$i" | awk '{ print $5}'`
  curl --request PUT --header "Content-Length: "${CONTENT_LENGTH}"" --header "Content-Type: multipart/mixed" --data-binary "@"$i"" "https://api-content.dropbox.com/1/files_put/sandbox/$i?access_token="${ACCESS_TOKEN}""
  printf "\n"
done
