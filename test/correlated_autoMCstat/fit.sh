combine -M MultiDimFit model_1.root  --algo=grid --points 100  -m 125   -t -1     \
    --redefineSignalPOIs k_cQl1 \
    --freezeParameters r  \
    --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-1,1      \
    --verbose 2 -n Morphed_1 --X-rtd MINIMIZER_no_analytic

combine -M MultiDimFit model_2.root  --algo=grid --points 100  -m 125   -t -1     \
    --redefineSignalPOIs k_cQl1 \
    --freezeParameters r  \
    --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-1,1      \
    --verbose 2 -n Morphed_2 --X-rtd MINIMIZER_no_analytic

combine -M MultiDimFit model_5.root  --algo=grid --points 100  -m 125   -t -1     \
    --redefineSignalPOIs k_cQl1 \
    --freezeParameters r  \
    --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-1,1      \
    --verbose 2 -n Morphed_5 --X-rtd MINIMIZER_no_analytic

# combine -M MultiDimFit model_10.root  --algo=grid --points 1  -m 125   -t -1     \
#     --redefineSignalPOIs k_cQl1 \
#     --freezeParameters r  \
#     --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-1,1      \
#     --verbose 1000000 -n Morphed_10 --X-rtd MINIMIZER_no_analytic


#combine -M MultiDimFit model_EFTNeg.root  --algo=grid --points 300  -m 125   -t -1     \
#    --redefineSignalPOIs k_cQl1 \
#    --freezeParameters r  \
#    --setParameters r=1,k_cQl1=0    --setParameterRanges k_cQl1=-5,5      \
#    --verbose 0 -n EFTNeg