pycbc_coinc_findtrigs \
    --verbose \
    --trigger-files Htrig.hdf5 Ltrig.hdf5 \
    --template-bank Tbank.hdf5 \
    --pivot-ifo H1 \
    --fixed-ifo L1 \
    --ranking-statistic quadsum \
    --sngl-ranking snr \
    --output-file coincs.hdf5 \
    --loudest-keep-values [1:1] 
