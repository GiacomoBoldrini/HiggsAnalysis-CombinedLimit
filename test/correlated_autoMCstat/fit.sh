combine -M MultiDimFit model.root  --algo=grid --points 100  -m 125   -t -1     \
    --redefineSignalPOIs k_cQl1 \
    --freezeParameters r  \
    --setParameters r=1    --setParameterRanges k_cQl1=-5,5     \
    --verbose -1