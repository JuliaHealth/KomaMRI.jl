window.BENCHMARK_DATA = {
  "lastUpdate": 1720031307776,
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
          "id": "ca8948b03a70041c0c44eba77b9455bcb364e2f0",
          "message": "Change other places from bench to benchmarks",
          "timestamp": "2024-07-02T19:39:34-05:00",
          "tree_id": "aed88491d66a9980111ad4a7db19f3856f6f17a9",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/ca8948b03a70041c0c44eba77b9455bcb364e2f0"
        },
        "date": 1719969227194,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1608040346,
            "unit": "ns",
            "extra": "gctime=175566684\nmemory=3375686560\nallocs=189859\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 1222363546,
            "unit": "ns",
            "extra": "gctime=250897084.5\nmemory=3382655880\nallocs=271738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 574514332,
            "unit": "ns",
            "extra": "gctime=61501534\nmemory=3396003248\nallocs=441447\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 2944794969,
            "unit": "ns",
            "extra": "gctime=757533600\nmemory=3372192016\nallocs=148020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 154454932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49388144\nallocs=732302\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 27116431712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239141760\nallocs=1950594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3378098021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116084040\nallocs=1657564\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1719966325,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44052648\nallocs=332877\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3330275254,
            "unit": "ns",
            "extra": "gctime=573419438\nmemory=6666067536\nallocs=181823\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 2321501465,
            "unit": "ns",
            "extra": "gctime=435055267\nmemory=6668518632\nallocs=197961\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1040737268,
            "unit": "ns",
            "extra": "gctime=85232695\nmemory=6673518656\nallocs=232472\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6069086826,
            "unit": "ns",
            "extra": "gctime=721580645\nmemory=6664850912\nallocs=172484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265068560,
            "unit": "ns",
            "extra": "gctime=174306537\nmemory=27797208\nallocs=263092\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2168012395,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64448664\nallocs=529485\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1424272625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40447336\nallocs=448024\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 745579010,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27469792\nallocs=215020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "f29f433186b6e0c90364c5ab52bd45a47d189b0c",
          "message": "Unindent env part",
          "timestamp": "2024-07-03T12:52:02-05:00",
          "tree_id": "2970e8dc6c643b0d6c8f2802204011f9438d777b",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/f29f433186b6e0c90364c5ab52bd45a47d189b0c"
        },
        "date": 1720031295902,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1513355443,
            "unit": "ns",
            "extra": "gctime=137917225.5\nmemory=3375686288\nallocs=189836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 1145395599.5,
            "unit": "ns",
            "extra": "gctime=224797439\nmemory=3382655880\nallocs=271738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 588303070,
            "unit": "ns",
            "extra": "gctime=77493807\nmemory=3396003176\nallocs=441446\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3189930588.5,
            "unit": "ns",
            "extra": "gctime=841311678\nmemory=3372192016\nallocs=148020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 153521272.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49383792\nallocs=732166\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 27439786432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239149064\nallocs=1950831\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3361973604,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116082472\nallocs=1657432\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1750871538.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44051264\nallocs=332859\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3170037917,
            "unit": "ns",
            "extra": "gctime=432985328\nmemory=6666067536\nallocs=181823\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 2192883721.5,
            "unit": "ns",
            "extra": "gctime=377091736.5\nmemory=6668511872\nallocs=197724\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1045831072,
            "unit": "ns",
            "extra": "gctime=99746673\nmemory=6673518656\nallocs=232472\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6094763205,
            "unit": "ns",
            "extra": "gctime=732602487\nmemory=6664850912\nallocs=172484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265646287,
            "unit": "ns",
            "extra": "gctime=175177332\nmemory=27796792\nallocs=263079\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2216409264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64448664\nallocs=529485\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1449848229,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40448248\nallocs=448041\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 723525528,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27469720\nallocs=215018\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}