set quiet

default:
    just --list

quick_graph FILE:
    [ -x "$(command -v gnuplot)" ] || { echo "gnuplot is not installed"; exit 1; }
    [ -x "$(command -v cargo)" ] || { echo "cargo is not installed"; exit 1; }
    cargo r -- analyze -i {{FILE}} | tail -n +8 | head -n -56 | gnuplot -p -e "plot '/dev/stdin' using 1:2 with lines"

synthetic_graph FREQ1 FREQ2:
    [ -x "$(command -v gnuplot)" ] || { echo "gnuplot is not installed"; exit 1; }
    [ -x "$(command -v cargo)" ] || { echo "cargo is not installed"; exit 1; }
    cargo r -- synthetic {{FREQ1}} {{FREQ2}} | tail -n +8 | head -n -57 | gnuplot -p -e "plot '/dev/stdin' using 1:2 with lines"

