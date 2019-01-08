#!/usr/bin/bash

rm -f boinc_finish_called cambala_depths_out cambala_out chpt out stderr.txt

time ../cambala_boinc_app

echo

#diff -qs cambala_depths_out cambala_depths_out.ref
#diff -qs cambala_out cambala_out.ref
diff -qs out out.ref
