1. For calculations of paramagnetic (PM) phase with different U, type

```sh
$ mkdir -p work && cd work
$ python3.7 ../scan_checkboard.py
```

2. For calculations of antiferromagnetic (AFM) phase, type

```sh
$ python3.7 ../scan_checkboard.py -sp
```

3. To compare the results of PM and AFM phase, type

```sh
$ python3.7 ../plot_pmafm.py
```

4. For calculations of AFM phase at Hartree-Fock (HF) mean-field level, type

```sh
$ python3.7 ../scan_checkboard.py -uhf
```

5. To compare the three calculations results, type

```sh
$ python3.7 ../plot_pmafm_gh.py
```

6. For calculations of PM phase at HF level, type

```sh
$ python3.7 ../scan_checkboard.py -rhf
```

7. For comparisons with all calculations, type

```sh
$ python3.7 ../plot_pmafm_gh2.py
```

