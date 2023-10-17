#!/usr/bin/env bash

charset='qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM!@#$%^&*()-=_+[]{}|;:,.<>?/'

for i in {1..100}; do
    filename="file${i}.txt"
    random_chars=""
    for j in {1..20}; do
        random_index=$((RANDOM % ${#charset}))
        random_chars="${random_chars}${charset:${random_index}:1}"
    done
    echo "$random_chars" > "$filename"
done

symbols_to_check="S1"

for symbol in $(echo "$symbols_to_check" | grep -o .); do
    files_to_delete=$(grep -l "$symbol" file*.txt)
    for file_to_delete in $files_to_delete; do
        rm -f "$file_to_delete"
        echo "Delited file: $file_to_delete"
    done
done
