# Run Distributed Simulations 

While KomaMRI provides built-in support for CPU and GPU parallelization, running simulations distributed across multiple GPUs or compute nodes is not automatically supported. However, it is possible to do so manually with Distributed.jl. The following examples demonstrate how:

## Using Multiple GPUs

To run a simulation using multiple GPUs, the phantom object can be divided using the kfoldperm function. Due to the independent spin property of the system, the signal result of simulating the entire phantom is equal to the sum of results from simulating each subdivision of the phantom. The following script divides the phantom among available CUDA GPUs, having each part taken by a different worker. The worker processes run their part of the simulation on a different CUDA device, after which the results are fetched asynchronously in the main process and recombined to produce the final signal.

```julia
using Distributed
using CUDA

#Add workers based on the number of available devices
addprocs(length(devices()))

#Define inputs on each worker process
@everywhere begin
    using KomaMRI, CUDA
    sys = Scanner()
    seq = PulseDesigner.EPI_example()
    obj = brain_phantom2D()
    #Divide phantom
    parts = kfoldperm(length(obj), nworkers())
end

#Distribute simulation across workers, fetch into array of worker signals
worker_signals = fetch.([ @spawnat pid begin
    KomaMRICore.set_device!(i-1) #Sets device for this worker, note that CUDA devices are indexed from 0
    simulate(obj[parts[i]], seq, sys)
end for (i, pid) in enumerate(workers()) ])

#Final signal
signal = reduce(+, worker_signals)
```

## Using Multiple Nodes in an HPC Cluster

In other cases, it may be useful to run a simulation on multiple compute nodes if the problem is too large to fit into memory for a single computer, or if the number of desired workers is greater than the amount of CPU cores available. The following script uses the package ClusterManagers.jl to initialize worker process on a SLURM cluster based on the number of tasks specified in the #SBATCH --ntasks directive:

```julia
using Distributed
using ClusterManagers

#Add workers based on the number of specified SLURM tasks 
addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])))

#Define inputs on each worker process
@everywhere begin
    using KomaMRI
    sys = Scanner()
    seq = PulseDesigner.EPI_example()
    obj = brain_phantom2D()
    parts = kfoldperm(length(obj), nworkers())
end

#Distribute simulation across workers, fetch into array of worker signals
worker_signals = fetch.([ @spawnat pid begin
    simulate(obj[parts[i]], seq, sys)
end for (i, pid) in enumerate(workers()) ])

#Final Signal
signal = reduce(+, worker_signals)
```
