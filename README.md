# OGCM Experiments Model Run and Analysis Code for Nagura (2026, JPO)

This repository contains Python/Julia scripts and related utilities used for numerical experiments and the analysis done by Nagura (2026, "Triad resonance between vertically propagating equatorial waves in an idealized OGCM", submitted to J. Phys. Oceanogr.).

The repository includes selected analysis scripts together with the local Python modules required by those scripts. The directory structure is preserved from the original working environment so that relationships between analysis scripts and shared utilities remain clear.

## Structure

* `themes/MidDepth/OFES_IdealExp/` — Analysis and visualization code.
* `common/` — Shared Python utilities used by multiple analysis scripts.
* `OGCM_Exp/` — Files to prepare initial and surface boundary conditions, model grids, and parameters to conduct numerical experiments.

## Python imports

The original code is organized with the source root on `PYTHONPATH`. For example:

```python
from common.some_module import some_function
from themes.MidDepth.OFES_IdealExp.Tools.some_module import some_function
```

Depending on where the repository is cloned, the repository root may need to be added to `PYTHONPATH` before running the scripts.

## Notes

This repository contains selected research code rather than the complete original working directory. Some scripts may depend on datasets or software packages that are not included in the repository.
