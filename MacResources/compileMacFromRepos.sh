#!/bin/bash -l
git clone https://github.com/NickKarpowicz/LightwaveExplorer
git clone https://github.com/davisking/dlib
cd LightwaveExplorer
chmod +x makeMacApp.sh
./makeMacApp.sh