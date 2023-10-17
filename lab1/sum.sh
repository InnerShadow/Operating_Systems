#!/usr/bin/env bash

N=$1
sum=0

for ((n=1; n<=N; n++)); do
    sum=$(bc -l <<< "$sum + 1 / ($n * $n)")
done

echo "Sum of first N=$sum"
