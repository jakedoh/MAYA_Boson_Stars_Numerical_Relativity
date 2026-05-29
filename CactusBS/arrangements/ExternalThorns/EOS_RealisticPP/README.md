# EOS_RealisticPP

Based on [astro-ph:0812.2163](https://arxiv.org/abs/0812.2163). The parameters for all equations of state are shown in Table III. Additionally, includes an optional thermal component (see [10.1103/PhysRevD.82.084043](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.82.084043) for details) to allow for shock heating.

Registers the equation of state with EOSBase with table name `RealisticPP`. Use with Whisky requries the parameters
```
Whisky::whisky_eos_type = "General"
Whisky::whisky_eos_table = "RealisticPP"
```

For initial data, it is often useful to calculate the internal energy from just the cold part of the equation of state and not the thermal part. For this purpose, this thorn implements the function `EOS_HybridColdEps(CCTK_REAL IN rho)` to be called for initial data. If Lorene data is being loaded, be sure to use `Meudon_Bin_NS_Scotch_Hybrid` instead of `Meudon_Bin_NS_Scotch`.

The only relevant parameters are:
 - `code_unit`: The code unit of mass in solar masses (to convert the EOS parameters from CGS to code units). Default value is 1.
 - `gamma_thermal`: Gamma for the thermal part of the hybrid equation of state. Set to one for no thermal component. Default value is 1.8.
 - `realistic_eos_name`: Name of the EOS from Table III of the above paper to use. Default value is 'SLy'.

For example, to use the `MS1b` equation of state with a code unit of one solar mass and a thermal gamma value of 1.8, the parameters would need to be 
```
EOS_RealisticPP::code_unit = 1.0
EOS_RealisticPP::gamma_thermal = 1.8
EOS_RealisticPP::realistic_eos_name = "MS1b"
```
