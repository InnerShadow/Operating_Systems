#!/usr/bin/env bash

input_names="$1"
input_weights="$2"
selected_person="$3"

input_names=$(echo "$input_names" | tr '[:lower:]' '[:upper:]')

names_array=($(echo "$input_names" | awk -F ',' '
{
    for (i=1; i<=NF; i++) 
        print $i
}'))

weights_array=($(echo "$input_weights" | awk -F ',' '
{
    for (i=1; i<=NF; i++) 
        print $i
}'))

if [ "${#names_array[@]}" -eq "0" ]; then
	echo "Empty names list!!!"
    exit 1
fi

if [ "${#weights_array[@]}" -lt "${#names_array[@]}" ]; then
    echo "Weights array size should be at least equal to the names array size!!!"
    exit 2
fi

if [ $((selected_person - 1)) -gt "${#names_array[@]}" ]; then
	echo "Invalid selected person!!!"
	exit 3
fi

if [ $((selected_person)) -eq "0" ]; then
    echo "Select some one!!!"
    exit 3
fi

if [ "$USER" != "masikol" ]; then
	echo "Stop cheating!!!"
	exit 4
fi

declare -A som_name_map

for ((i = 0; i < ${#names_array[@]}; i++)); do
    element="${names_array[$i]}"
    weight="${weights_array[$i]}"
    
    som=${#element}
    for ((j = 0; j < ${#element}; j++)); do
        char="${element:$j:1}"
        ascii_value=$(($(printf "%d" "'$char") - $(printf "%d" "'A") + 1))
        som=$((som + ascii_value))
    done

    result=$((som * weight))
    som_name_map["$result"]=$element
done

sorted_soms=($(
    for som in "${!som_name_map[@]}"; do
        echo "$som"
    done | sort -n))

for som in "${sorted_soms[@]}"; do
    echo "Som: $som, Name: ${som_name_map[$som]}"
done

selected_som="${sorted_soms[$((selected_person - 1))]}"
selected_name="${som_name_map[$selected_som]}"
echo "Selected person #$selected_person - Som: $selected_som, Name: $selected_name"

exit 0

