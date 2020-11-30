#!/usr/bin/env bash
if [[ $OSTYPE == linux-gnu ]]; then
    sed -i -r 's/^case ([0-9]+)/goto case; case \1/g' $1
else
    sed -i -r '.bak' 's/^case ([0-9]+)/goto case; case \1/g' $1
fi
