#!/usr/bin/env bash

length=$1
charset='qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM!@#$%^&*()-=_+[]{}|;:,.<>?/'

password=""

for i in $(seq 1 "$length"); do
    random_char=${charset:$((RANDOM % ${#charset})):1}
    password="${password}${random_char}"
done

echo "Password: $password"
