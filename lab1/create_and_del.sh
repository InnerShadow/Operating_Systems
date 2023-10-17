#!/usr/bin/env bash

for i in {1..100}; do
    filename="file${i}.txt"
    touch "$filename"
    echo "Created file: $filename"
done

for i in {1..100}; do
    if (( i % 5 == 0 )); then
        filename="file${i}.txt"
        rm -f "$filename"
        echo "Deleted file: $filename"
    fi
done
