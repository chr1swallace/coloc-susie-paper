Code to evaluate SuSiE in coloc analysis. Everything is run via rake.

``` sh
rake approx:clean    # remove ALL output files
rake approx:run      # run approximation simulations
rake approx:summary  # summarise susie sims to date

rake coloc:clean     # remove ALL output files
rake coloc:run       # run coloc on sim data
rake coloc:summary   # summarise sims to date
```
