window.BENCHMARK_DATA = {
  "lastUpdate": 1719628652608,
  "repoUrl": "https://github.com/JuliaHealth/KomaMRI.jl",
  "entries": {
    "KomaMRI Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "cncastillo@uc.cl",
            "name": "Carlos Castillo Passi",
            "username": "cncastillo"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0f528f4cfbf49582d4752e1c1b3692c53a9e516d",
          "message": "Changed permissions to github_token to push to gh-pages",
          "timestamp": "2024-06-28T22:14:32-04:00",
          "tree_id": "ad2259aec21b43e43f4788cd12125847bc048a7f",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/0f528f4cfbf49582d4752e1c1b3692c53a9e516d"
        },
        "date": 1719628641300,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1515831906.5,
            "unit": "ns",
            "extra": "gctime=149093509.5\nmemory=3393663024\nallocs=179795\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 873515619,
            "unit": "ns",
            "extra": "gctime=99634718.5\nmemory=3399680848\nallocs=252660\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 952723201,
            "unit": "ns",
            "extra": "gctime=233062771\nmemory=3411067456\nallocs=402349\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3432264757.5,
            "unit": "ns",
            "extra": "gctime=1063611924.5\nmemory=3390558176\nallocs=143859\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 165219104.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50365608\nallocs=767038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 35961041138,
            "unit": "ns",
            "extra": "gctime=0\nmemory=261202344\nallocs=2136402\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3821415521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=121763328\nallocs=1853348\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1837883119,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44732352\nallocs=346676\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3274567292,
            "unit": "ns",
            "extra": "gctime=639881143\nmemory=6675455568\nallocs=179199\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1743775313,
            "unit": "ns",
            "extra": "gctime=192383033\nmemory=6677693848\nallocs=191285\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1338030231.5,
            "unit": "ns",
            "extra": "gctime=280259215.5\nmemory=6682244664\nallocs=220413\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6933313740,
            "unit": "ns",
            "extra": "gctime=1781078786.5\nmemory=6674377184\nallocs=171172\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 268630755,
            "unit": "ns",
            "extra": "gctime=174577599\nmemory=28127160\nallocs=274715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2411448166.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71822808\nallocs=592076\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1439350584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42342688\nallocs=513344\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 838445296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27696296\nallocs=219617\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}