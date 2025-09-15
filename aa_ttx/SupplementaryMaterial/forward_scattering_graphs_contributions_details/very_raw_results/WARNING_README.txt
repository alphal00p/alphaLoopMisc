WARNING: the cross-section per forward-scattering graph indicated in these files cannot be used without further processing.
In particular:
1) An overall minus sign must be applied.
2) alphaEW^-1 is set to 132.507 for this generation.
3) The number of massless quark flavour is set to 1.
4) The process considered is photon photon > t tbar.

However, the files in this folder contain a wealth of integration statistics, in particular regarding the number of sample points as well as the fraction of those assigned to each forward-scattering graphs. The chi^2 per d.o.f. of the integration as well as the maximum weight influence are also given.
The file `very_raw_results/scan_singlets_with_multichanneling_over_e_surfaces_mathematica.txt` contains multiple entries for two of the singlet graphs, due to the multichanelling procedure over centers for the E-surface subtraction described in our publication. They can however be added coherently.
If you are interested in the physical contribution of each forward scattering graph, with the same input parameters as given in our publication, refer instead to the file `../raw_data_with_singlet.py`.

The content of the files given here contain detailed result for each point in our scan in the sqrt{s} of the photon-fusion process and are meant to be processed by Mathematica.
