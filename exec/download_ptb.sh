#!/bin/sh
# Utility script to download the PTB database from Physionet.
# This will create a directory names 'ptbdb' in the current directory.

wget -r --cut-dirs=2 -np -nH http://www.physionet.org/physiobank/database/ptbdb/
