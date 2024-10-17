combine -M MultiDimFit model_1.root  --algo=grid --points 10000  -m 125   -t -1     \
    --redefineSignalPOIs k_cje,k_cHj1 \
    --freezeParameters r  \
    --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-0.2,0.2:k_cHj1=-10,10      \
    --verbose 3 -n Morphed_1_bbth_100 --X-rtd MINIMIZER_no_analytic

# combine -M MultiDimFit model_2.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje,k_cHj1 \
#     --freezeParameters r  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_2
# 
# combine -M MultiDimFit model_5.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje,k_cHj1 \
#     --freezeParameters r  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_5

# combine -M MultiDimFit model_EFTNeg.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje,k_cHj1 \
#     --freezeParameters r  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose 2 -n Morphed_EFTNeg

# --X-rtd MINIMIZER_no_analytic
