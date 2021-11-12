*******************
Benchmark
*******************

We tested MoP on different datasets at the CRG's HPC where we can run up to 100 jobs in parallel (maximum 8 CPUs each) and using up to 10 GPU cards (GeForce RTX 2080 Ti).

MOP_PREPROCESS
-----------------

.. list-table:: Dataset
   
 * - 
   - Toy sample
   - Flongle
   - MinION
   - gridION
   - PromethION
 * - Fast5 files
   - 10 
   - 20 
   - 100 
   - 500 
   - 2,916 
 * - Reads
   - 8,000
   - 16,000
   - 400,000 
   - 2,000,000
   - 11,600,000
 * - Execution time
   - 10m
   - 18m
   - 1h 4m
   - 4h 19m
   - x

MOP_MOD
-----------------

.. list-table:: Dataset

 * - 
   - Toy sample
   - Flongle
   - MinION
   - gridION
   - PromethION
 * - Fast5 files
   - 10 
   - 20 
   - 100 
   - 500 
   - 2,916 
 * - Reads
   - 8,000
   - 16,000
   - 400,000 
   - 2,000,000
   - 11,600,000
 * - Execution time
   - 45m
   - x
   - x
   - x
   - x

