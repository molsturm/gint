#!/bin/bash

FOLDERS=".last_pull krims lazyten rapidcheck sturmint"
read -p "Press Enter to delete the following folders:  $FOLDERS "
rm -rf $FOLDERS
exit $?
