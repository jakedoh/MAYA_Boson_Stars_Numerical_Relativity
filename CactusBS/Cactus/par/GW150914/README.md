
# Binary black hole GW150914

See [http://einsteintoolkit.org/gallery/bbh/index.html](http://einsteintoolkit.org/gallery/bbh/index.html).

We ask that if you make use of the parameter file or the example data,
then please cite the GW150914 Einstein Toolkit example and data, the
Einstein Toolkit, the Llama multi-block infrastructure, the Carpet
mesh-refinement driver, the apparent horizon finder AHFinderDirect,
the TwoPunctures initial data code, QuasiLocalMeasures, Cactus, and
the McLachlan spacetime evolution code, the Kranc code generation
package, and the Simulation Factory.

An appropriate [bibtex file](etgw150914.bib) is in this directory.

Submission command:

    simfactory/bin/sim create-submit GW150914_28 --define N 28 --parfile par/GW150914.rpar --procs 128 --walltime 24:00:00
The total amount of memory required is about 100 GB, and it runs for about 3 days on 128 Intel Haswell Xeon cores.

