#!/usr/bin/env bash

for i in {1..100}; do
    filename="file${i}.txt"
    touch "$filename"
done

for filename in file*0*.txt; do
    new_filename="${filename//0/_}"
    mv "$filename" "$new_filename"
    echo "Renamed: $filename -> $new_filename"
done
