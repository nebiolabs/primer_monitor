#!/usr/bin/env bash
source ~/.bash_profile
source ~/.zshrc
rubocop_status=$(bundle exec rubocop -F > /dev/null && echo 'success')

if [[ "$rubocop_status" = "success" ]]
then
    exit 0
else
    cat <<EOF
Error: Rubocop detected style conflicts, commit cancelled. For full details, run:
bundle exec rubocop
EOF
    exit 1
fi
