#!/usr/bin/env bash

# download shunit2 in order to run tests:
# curl -L "https://dl.dropboxusercontent.com/u/7916095/shunit2-2.0.3.tgz" | tar zx --overwrite

./.test_suite.sh | tee /dev/stderr | grep -q 'success rate: 100%'
