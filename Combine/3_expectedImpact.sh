#!/bin/bash

cmsenv

mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
# mAs=( 5 15 30)
# mAs=( 1 )

mkdir -p output_impacts


for mA in "${mAs[@]}"; do
    echo "=========================================="
    echo "Processing Text2Workspace for mA = ${mA}"
    echo "=========================================="

    # 1. 生成workspace文件
    python3 RunText2Workspace.py --batch local --queue workday --mA ${mA}
done