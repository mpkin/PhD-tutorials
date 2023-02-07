#! /bin/bash

make clean; make; mkdir -p ./data; cat *.rtparam *.fparam > mg.param; ./mg mg.param

