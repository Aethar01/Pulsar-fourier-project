set quiet

default:
    just --list

concatenate OUTPUTFILE:
    tail -n +1 src/* > {{OUTPUTFILE}}
