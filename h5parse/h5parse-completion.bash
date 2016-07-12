#!/usr/local/bin/bash

_h5parse() 
{
    local cur prev file opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

	#TODO complete correctly with -o in the middle
	if [[ ${prev} == "-o" ]] ; then
		return 0

	elif [[ COMP_CWORD -eq 1 ]] ; then
		COMPREPLY=( $(compgen -f -X '!*.h5' ${cur}) ); compopt -o filenames -o plusdirs; return 0

	else 
		file="${COMP_WORDS[1]}"
		opts="$(eval "h5dump --contents $file | perl -ne 'if (m~^\s*dataset\s*($cur.*?($|/))~) {print \"\$1\\n\"}'" | uniq)"
        
        #if only one possible argument check to see if dataset or group to decide whether to add extra space
        n=$(echo $opts | wc -w) 
        if [[ $n -eq 1 ]] ; then
            if [[ ${opts: -1} != "/" ]] ; then
                COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) ); return 0
            fi
        fi
		COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) ); compopt -o nospace; return 0
	fi
}
complete -F _h5parse h5parse
