#!/bin/bash

mkdir -p af2_initial_guess
cd af2_initial_guess

mkdir -p model_weights/params && cd model_weights/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar --extract --verbose --file=alphafold_params_2022-12-06.tar \
&& rm alphafold_params_2022-12-06.tar
