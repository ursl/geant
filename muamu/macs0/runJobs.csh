#! /bin/csh

cd /Users/ursl/muamu/macs0/_build
./runG4 -n 10000 -bgE 0.005 -f g4run-005kev.root -c vis.mac >&! 005kev.log
./runG4 -n 10000 -bgE 0.010 -f g4run-010kev.root -c vis.mac >&! 010kev.log
./runG4 -n 10000 -bgE 0.020 -f g4run-020kev.root -c vis.mac >&! 020kev.log
./runG4 -n 10000 -bgE 0.050 -f g4run-050kev.root -c vis.mac >&! 050kev.log
./runG4 -n 10000 -bgE 0.100 -f g4run-100kev.root -c vis.mac >&! 100kev.log
./runG4 -n 10000 -bgE 0.200 -f g4run-200kev.root -c vis.mac >&! 200kev.log
./runG4 -n 10000 -bgE 0.500 -f g4run-500kev.root -c vis.mac >&! 500kev.log
