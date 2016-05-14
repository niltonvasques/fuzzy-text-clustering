#!/bin/bash

cat <(wc -l discover.names | awk '{print $1}') <(wc -l discover.data | awk '{print $1}') discover.names discover.data > in
