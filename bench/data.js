window.BENCHMARK_DATA = {
  "lastUpdate": 1719967965257,
  "repoUrl": "https://github.com/JuliaHealth/KomaMRI.jl",
  "entries": {
    "KomaMRI Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "ryanakierulf@gmail.com",
            "name": "Ryan Kierulf",
            "username": "rkierulf"
          },
          "committer": {
            "email": "ryanakierulf@gmail.com",
            "name": "Ryan Kierulf",
            "username": "rkierulf"
          },
          "distinct": true,
          "id": "124b2ad0959dcb740959181b487db44ab04924ea",
          "message": "Change benchmark-data-dir-path to \"bench\" from \"bench/dev\"",
          "timestamp": "2024-07-02T19:28:47-05:00",
          "tree_id": "57cdd722037d1e33b4676a8e59090a893fa7aed4",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/124b2ad0959dcb740959181b487db44ab04924ea"
        },
        "date": 1719967954418,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1717249974,
            "unit": "ns",
            "extra": "gctime=304050617\nmemory=3375685936\nallocs=189814\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 903442534,
            "unit": "ns",
            "extra": "gctime=114673776\nmemory=3382669968\nallocs=272208\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 881169783,
            "unit": "ns",
            "extra": "gctime=191852171\nmemory=3396010200\nallocs=441691\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3473828701.5,
            "unit": "ns",
            "extra": "gctime=690062534.5\nmemory=3372192016\nallocs=148020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 154280769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49418608\nallocs=732302\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 27802380580,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239142064\nallocs=1950595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3663731833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116083360\nallocs=1657494\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1790099598,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44053168\nallocs=332737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3758684814,
            "unit": "ns",
            "extra": "gctime=1123620647\nmemory=6666067536\nallocs=181823\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1644532356,
            "unit": "ns",
            "extra": "gctime=160940433\nmemory=6668505080\nallocs=197533\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1366658335,
            "unit": "ns",
            "extra": "gctime=247871189\nmemory=6673524936\nallocs=232654\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6090359461.5,
            "unit": "ns",
            "extra": "gctime=470777534.5\nmemory=6664850912\nallocs=172484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265363308,
            "unit": "ns",
            "extra": "gctime=175080911\nmemory=27800120\nallocs=263092\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2175088901.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64448968\nallocs=529486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1423441708.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40448304\nallocs=448082\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 743231180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27469144\nallocs=215011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}