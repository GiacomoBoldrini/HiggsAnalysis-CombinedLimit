combine -M MultiDimFit model.root  --algo=grid --points 5  -m 125   -t -1     \
    --redefineSignalPOIs k_cQl1 \
    --freezeParameters r  \
    --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-5,5      \
    --verbose 0 -n Morphed


# combine -M MultiDimFit model_EFTNeg.root  --algo=grid --points 1000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cQl1 \
#     --freezeParameters r  \
#     --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-5,5      \
#     --verbose 0 -n EFTNeg