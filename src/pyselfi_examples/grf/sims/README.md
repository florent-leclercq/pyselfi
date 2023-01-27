Due to their size, the raw simulation data are not stored in the git repository. They can be regenerated using the notebooks [selfi\_GRF.ipynb](https://github.com/florent-leclercq/pyselfi/blob/master/pyselfi_examples/grf/selfi\_GRF.ipynb) or downloaded from the Aquila Consortium website at https://cloud.aquila-consortium.org/s/H4jeyqdtq2dTMCa.

The files should be placed in the current directory (pyselfi/src/pyselfi\_examples/grf/sims/).

Content:

```lang-sh
    [       4096] .
    ├── [   34229228]  data_input_power.h5
    ├── [     271080]  fiducial.npz
    ├── [     271098]  groundtruth.npz
    ├── [   34092208]  G_sim.h5
    ├── [   33822608]  G_ss.h5
    ├── [        528]  P_0.npy
    ├── [   34229228]  P_exp.h5
    ├── [        300]  phi_obs.npy
    ├── [       4096] pools_GRF_iter0
    │   ├── [     105672]  pool_0.h5
    │   ├── [      72904]  pool_100.h5
    │   ├── [      72904]  pool_10.h5
    │   ├── [      72904]  pool_11.h5
    │   ├── [      72904]  pool_12.h5
    │   ├── [      72904]  pool_13.h5
    │   ├── [      72904]  pool_14.h5
    │   ├── [      72904]  pool_15.h5
    │   ├── [      72904]  pool_16.h5
    │   ├── [      72904]  pool_17.h5
    │   ├── [      72904]  pool_18.h5
    │   ├── [      72904]  pool_19.h5
    │   ├── [      72904]  pool_1.h5
    │   ├── [      72904]  pool_20.h5
    │   ├── [      72904]  pool_21.h5
    │   ├── [      72904]  pool_22.h5
    │   ├── [      72904]  pool_23.h5
    │   ├── [      72904]  pool_24.h5
    │   ├── [      72904]  pool_25.h5
    │   ├── [      72904]  pool_26.h5
    │   ├── [      72904]  pool_27.h5
    │   ├── [      72904]  pool_28.h5
    │   ├── [      72904]  pool_29.h5
    │   ├── [      72904]  pool_2.h5
    │   ├── [      72904]  pool_30.h5
    │   ├── [      72904]  pool_31.h5
    │   ├── [      72904]  pool_32.h5
    │   ├── [      72904]  pool_33.h5
    │   ├── [      72904]  pool_34.h5
    │   ├── [      72904]  pool_35.h5
    │   ├── [      72904]  pool_36.h5
    │   ├── [      72904]  pool_37.h5
    │   ├── [      72904]  pool_38.h5
    │   ├── [      72904]  pool_39.h5
    │   ├── [      72904]  pool_3.h5
    │   ├── [      72904]  pool_40.h5
    │   ├── [      72904]  pool_41.h5
    │   ├── [      72904]  pool_42.h5
    │   ├── [      72904]  pool_43.h5
    │   ├── [      72904]  pool_44.h5
    │   ├── [      72904]  pool_45.h5
    │   ├── [      72904]  pool_46.h5
    │   ├── [      72904]  pool_47.h5
    │   ├── [      72904]  pool_48.h5
    │   ├── [      72904]  pool_49.h5
    │   ├── [      72904]  pool_4.h5
    │   ├── [      72904]  pool_50.h5
    │   ├── [      72904]  pool_51.h5
    │   ├── [      72904]  pool_52.h5
    │   ├── [      72904]  pool_53.h5
    │   ├── [      72904]  pool_54.h5
    │   ├── [      72904]  pool_55.h5
    │   ├── [      72904]  pool_56.h5
    │   ├── [      72904]  pool_57.h5
    │   ├── [      72904]  pool_58.h5
    │   ├── [      72904]  pool_59.h5
    │   ├── [      72904]  pool_5.h5
    │   ├── [      72904]  pool_60.h5
    │   ├── [      72904]  pool_61.h5
    │   ├── [      72904]  pool_62.h5
    │   ├── [      72904]  pool_63.h5
    │   ├── [      72904]  pool_64.h5
    │   ├── [      72904]  pool_65.h5
    │   ├── [      72904]  pool_66.h5
    │   ├── [      72904]  pool_67.h5
    │   ├── [      72904]  pool_68.h5
    │   ├── [      72904]  pool_69.h5
    │   ├── [      72904]  pool_6.h5
    │   ├── [      72904]  pool_70.h5
    │   ├── [      72904]  pool_71.h5
    │   ├── [      72904]  pool_72.h5
    │   ├── [      72904]  pool_73.h5
    │   ├── [      72904]  pool_74.h5
    │   ├── [      72904]  pool_75.h5
    │   ├── [      72904]  pool_76.h5
    │   ├── [      72904]  pool_77.h5
    │   ├── [      72904]  pool_78.h5
    │   ├── [      72904]  pool_79.h5
    │   ├── [      72904]  pool_7.h5
    │   ├── [      72904]  pool_80.h5
    │   ├── [      72904]  pool_81.h5
    │   ├── [      72904]  pool_82.h5
    │   ├── [      72904]  pool_83.h5
    │   ├── [      72904]  pool_84.h5
    │   ├── [      72904]  pool_85.h5
    │   ├── [      72904]  pool_86.h5
    │   ├── [      72904]  pool_87.h5
    │   ├── [      72904]  pool_88.h5
    │   ├── [      72904]  pool_89.h5
    │   ├── [      72904]  pool_8.h5
    │   ├── [      72904]  pool_90.h5
    │   ├── [      72904]  pool_91.h5
    │   ├── [      72904]  pool_92.h5
    │   ├── [      72904]  pool_93.h5
    │   ├── [      72904]  pool_94.h5
    │   ├── [      72904]  pool_95.h5
    │   ├── [      72904]  pool_96.h5
    │   ├── [      72904]  pool_97.h5
    │   ├── [      72904]  pool_98.h5
    │   ├── [      72904]  pool_99.h5
    │   └── [      72904]  pool_9.h5
    ├── [   33824828]  P_ss.h5
    ├── [       5841]  README.md
    └── [    1072584]  selfi_GRF_iter0.h5

    1 directory, 112 files
