#!/usr/bin/env bash
if [[ $OSTYPE == linux-gnu ]]; then
	sed -i "s/static const/static/g" sam_alignment.d
else
	sed -i ".bak" "s/static const/static/g" sam_alignment.d
fi
