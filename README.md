# coulombpy
Simple python library for the Coulomb scattering functions which provides bindings for Ian Thompson's Fortran 77 [coulfg1 library](http://www.fresco.org.uk/functions.htm) using [numpy.f2py](https://numpy.org/doc/stable/f2py/)

## build

```bash
make 
```

To delete the `.pyf` configuration file and make a new one:

```bash
make confclean
make conf
```

## use

```python
import coulfg
print(coulfg.__doc__)
```
produces:

 ```
This module 'coulfg' is auto-generated with f2py (version:1.26.4).
Functions:
    fc,gc,fcp,gcp = coulfg(xx,eta1,xlmin,xlmax,mode1=1,kfn=0)
.
 ```

```python
z = 2.3
eta = 3.1
lmin= 0
lmax = 3
f, g, fp, gp = coulfg.coulfg(z, eta, lmin, lmax)
print(f, g, fp, gp)
```

produces:

```
[0.03388977 0.02178185 0.00964314] [11.88688022 16.69887226 32.28176782] [0.04805471 0.03390169 0.0172948
9] [-12.65217925 -19.91933897 -45.8035791 ]

```

 
[`./coulfg_benchmark.ipynb`](https://github.com/beykyle/coulombpy/blob/main/coulfg_benchmark.ipynb) contains a comparison to the arbitrary percision implementation in [`mpmath.coulombf`](https://mpmath.org/doc/0.18/functions/bessel.html#coulombf). Obviously, Ian's floating point implementation is much faster.
