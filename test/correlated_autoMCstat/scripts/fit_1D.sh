combine -M MultiDimFit model_10.root  --algo=grid --points 100  -m 125   -t -1     \
    --redefineSignalPOIs k_cje \
    --freezeParameters r,k_cHj1  \
    --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
    --verbose 3 -n Morphed_10_cje --X-rtd MINIMIZER_no_analytic





# # cHj1
# 
# combine -M MultiDimFit model_1.root  --algo=grid --points 100  -m 125   -t -1     \
#     --redefineSignalPOIs k_cHj1 \
#     --freezeParameters r,k_cje  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_1_cHj1 
# 
# combine -M MultiDimFit model_2.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cHj1 \
#     --freezeParameters r,k_cje  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_2_cHj1
# 
# combine -M MultiDimFit model_5.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cHj1 \
#     --freezeParameters r,k_cje  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_5_cHj1
# 
# #cje
# 
# combine -M MultiDimFit model_1.root  --algo=grid --points 100  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje \
#     --freezeParameters r,k_cHj1  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_1_cje
# 
# combine -M MultiDimFit model_2.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje \
#     --freezeParameters r,k_cHj1  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_2_cje
# 
# combine -M MultiDimFit model_5.root  --algo=grid --points 10000  -m 125   -t -1     \
#     --redefineSignalPOIs k_cje \
#     --freezeParameters r,k_cHj1  \
#     --setParameters r=1,k_cje=0,k_cHj1=0    --setParameterRanges k_cje=-1,1:k_cHj1=-1,1      \
#     --verbose -1 -n Morphed_5_cje
# 
# # --X-rtd MINIMIZER_no_analytic
