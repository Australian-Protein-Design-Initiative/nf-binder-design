#!/bin/bash

# From: https://github.com/meilerlab/HyperMPNN.git

mkdir -p HyperMPNN/retrained_models

(cd HyperMPNN/retrained_models && \
 wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_002_epoch240_hyper.pt && \
 wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_010_epoch300_hyper.pt && \
 wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_020_epoch300_hyper.pt && \
 wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_030_epoch300_hyper.pt
)
