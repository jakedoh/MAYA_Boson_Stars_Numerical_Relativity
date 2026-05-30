# Meudon_Bin_NS_Scotch_Hybrid

A minor alteration of the thorn `Meudon_Bin_NS_Scotch` to properly interface with hybrid equations of state. Requires the EOS thorn to implement the function
```
CCTK_REAL FUNCTION EOS_HybridColdEps(CCTK_REAL IN rho)
```
in order to calculte just the cold contribution to the specific internal energy. For an example, note the implementation in the thorn `Scotch/EOS_RealisticPP`.
