# MatLab Implementation for Binaural Decoding from Ambisonics
This script decodes **higher-order ambisonic room impulse responses** (HOA RIRs) stored in HDF5 (.h5) to **Binaural RIRs** (.wav).


## Pipeline:

1. Build a **Spherical Harmonic (SH) basis** on a coarse decode grid

2. Decode the **HOA IRs** to **Directional RIRs (DRIRs)**

3. Optionally **upsample** directions with **SARITA** *(coarse → dense grid)*

4. Pick the **closest HRTF** direction from a SOFA file and **convolve** to get BRIRs

5. **Save** per-direction BRIRs as *azXXX_el±YY.wav*


## Requirements:

* MATLAB (R202x recommended)

* SFS Toolbox (SFS_start)

* Real SH toolbox providing getSH (<https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m>)

* SOFA API for MATLAB (sofaread)

* SARITA Toolbox (<https://github.com/AudioGroupCologne/SARITA>)

* HRTF SOFA file (I am using HRIRs from <https://www.riec.tohoku.ac.jp/pub/hrtf/hrtf_data.html>)

* HOA RIRs in .h5 with datasets at /spatial_ir/<uuid> shaped [T × (N+1)^2]



## Quick start

1. **Edit paths** at the top of the script:

```Matlab
addpath(genpath('c:\...\sfs-matlab'));
addpath(genpath('c:\...\Spherical-Harmonic-Transform'));
SFS_start;

sofa_filename = 'C:/.../RIEC_hrir_subject_001.sofa';
data_root = '\\server\path\to\TrebleSDK_Results\Apartment_45';
```

2. Set **HOA order** and **grids** *(defaults in the script)*

* **HOA order** = 8 → 81 channels

* **Decode (coarse) grid**: az=0:15:345, el=−30:30:60

* **Target (dense) grid** for SARITA: az=0:10:350, el=−30:10:60

3. Run the script. BRIRs are written under a subfolder with the same name as each .h5 file.


## License / credits

* Uses SFS Toolbox (license in their repo)

* Uses SOFA API & HRTF dataset (observe their licenses)

* Uses a real-SH implementation (see its license)
