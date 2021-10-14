#!/usr/bin/env bash


old_dir="14/$1/debug/"
new_dir="31/$2/0/"
old_path="/data/chodera/zhangi/perses_benchmark/neq/"
new_path="/data/chodera/zhangi/perses_benchmark/repex/"
rbd_file="4_rbd_ace2_mutant_rbd_tleap_final.pdb"
ace2_file="4_rbd_ace2_mutant_ace2_tleap_final.pdb"

cp "${old_path}${old_dir}${rbd_file}" "${new_path}${new_dir}"0_rbd.pdb
cp "${old_path}${old_dir}${ace2_file}" "${new_path}${new_dir}"0_ace2.pdb
