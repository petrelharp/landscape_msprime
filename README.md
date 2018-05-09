These are some tools for using [msprime](http://msprime.readthedocs.io/en/latest/index.html)
to run coalescent simulations on discretizations of geographic landscapes.

In this repository:

- [landscape_msprime/](landscape_msprime/) A python package; install it by running
    ```
    python3 setup.py install --user
    ```

    * this includes tools for selecting and computing a large, randomly chosen set
      of statistics, in [landscape_msprime/stats.py](landscape_msprime/stats.py) 

- [tortoises/](tortoises/) : An application to desert tortoises.

    * [tortoises/make_migration.R](tortoises/make_migration.R) : an example of using the R package
      [landsim](https://githum.com/petrelharp/landsim) to make a discretization of a migration kernel
      on a landscape defined by a raster
