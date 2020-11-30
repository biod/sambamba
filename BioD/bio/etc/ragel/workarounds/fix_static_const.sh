#!/usr/bin/env bash
if [[ $OSTYPE == linux-gnu ]]; then
	sed -i "s/static const/static/g" $1
else
	sed -i ".bak" "s/static const/static/g" $1
fi
