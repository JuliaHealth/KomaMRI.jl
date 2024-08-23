# Run Distributed Simulations 

While KomaMRI provides built-in support for CPU and GPU parallelization, it is sometimes desirable to distribute simulation work even further across multiple GPUs or compute nodes. This can be done by using Distributed.jl and making use of the independent spin property: each spin in the system is independent from the rest, so the phantom spins can be subdivided into separate simulations and results recombined, as in the diagram below:

```@raw html
<p align="center"><img width="90%" src="../../assets/KomamultiNode.svg"/></p>
```

The following two examples demonstrate how to use Distributed.jl to run a simulation using multiple GPUS, and using multiple nodes in an HPC cluster.

## Using Multiple GPUs

To run a simulation using multiple GPUs, the phantom object can be divided using the kfoldperm function. Distributed.jl can then be used to start one Julia worker process per available device so that each device simulates a different part of the object. The results can then be fetched asynchronously by the main process and combined to produce a final signal. This is shown in the following diagram: 

```@raw html
<p align="center"><img width="90%" src="../../assets/KomamultiGPU.svg"/></p>
```

The SLURM SBATCH script below requests 4 GPUs, all on a single computer. The Julia code then distributes work among each GPU:

!!! details "SLURM Script Requesting Multiple GPUs"

    ```sh
    #!/bin/bash
    #SBATCH --job-name KomaDistributed                 # Job name
    #SBATCH -t 0-00:30                                 # Max runtime for job
    #SBATCH -p batch                                   # Enter partition on which to run the job
    #SBATCH --ntasks=1                                 # 1 task
    #SBATCH --cpus-per-task=1                          # Request 1 CPU
    #SBATCH --gpus=4                                   # Request 4 GPUs
    #SBATCH -o /mnt/workspace/%u/slurm-out/%test.out   # Enter file path to write stdout to
    #SBATCH -e /mnt/workspace/%u/slurm-out/%test.err   # Enter file path to write stderr to

    module load julia/1.10.2
    julia script.jl
    ```

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

#Distribute simulation across workers
raw = Distributed.@distributed (+) for i=1:nworkers()
    KomaMRICore.set_device!(i-1) #Sets device for this worker, note that CUDA devices are indexed from 0
    simulate(obj[parts[i]], seq, sys)
end
```

## Using Multiple Nodes in an HPC Cluster

This example uses the package ClusterManagers.jl to initialize worker processes on a SLURM cluster based on the number of tasks specified in the #SBATCH --ntasks directive. This can be useful to divide simulation work among multiple compute nodes if the problem is too large to fit into memory for a single computer, or if the number of desired workers is greater than the typical number of CPU cores available. An illustration of this is shown below:

```@raw html
<p align="center"><img width="90%" src="../../assets/KomamultiNodeCPU.svg"/></p>
```

This SBATCH script requests 20 separate nodes, with each taking a single task. The Julia code is similar to the example for multiple GPUs, but initializes the processes slightly differently:

!!! details "SLURM Script Requesting Multiple Nodes"

    ```sh
    #!/bin/bash
    #SBATCH --job-name KomaDistributed                 # Job name
    #SBATCH -t 0-00:30                                 # Max runtime for job
    #SBATCH -p batch                                   # Enter partition on which to run the job
    #SBATCH --nodes=20                                 # 20 nodes
    #SBATCH --ntasks=20                                # 20 tasks
    #SBATCH --ntasks-per-node=1                        # 1 task per node
    #SBATCH --cpus-per-task=4                          # 4 CPUs per task
    #SBATCH -o /mnt/workspace/%u/slurm-out/%test.out   # Enter file path to write stdout to
    #SBATCH -e /mnt/workspace/%u/slurm-out/%test.err   # Enter file path to write stderr to

    module load julia/1.10.2
    julia script.jl
    ```

```julia
using Distributed
using ClusterManagers

#Add workers based on the specified number of SLURM tasks
addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])))

#Define inputs on each worker process
@everywhere begin
    using KomaMRI
    sys = Scanner()
    seq = PulseDesigner.EPI_example()
    obj = brain_phantom2D()
    parts = kfoldperm(length(obj), nworkers())
end

#Distribute simulation across workers
raw = Distributed.@distributed (+) for i=1:nworkers()
    simulate(obj[parts[i]], seq, sys)
end
```
