window.BENCHMARK_DATA = {
  "lastUpdate": 1765232910907,
  "repoUrl": "https://github.com/JuliaHealth/KomaMRI.jl",
  "entries": {
    "KomaMRI Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "ryanakierulf@gmail.com",
            "name": "rkierulf",
            "username": "rkierulf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5b6e9fbbda2311e57f3988185f1d3b13eefe3da1",
          "message": "Change 'main' to 'master' so benchmarking action actually runs",
          "timestamp": "2024-07-03T16:57:33-05:00",
          "tree_id": "ab37319cdb32d8b3cf2ee7562c717325ff0d804d",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/5b6e9fbbda2311e57f3988185f1d3b13eefe3da1"
        },
        "date": 1720045399303,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1533407266.5,
            "unit": "ns",
            "extra": "gctime=150396760\nmemory=3375678928\nallocs=189569\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 945216808,
            "unit": "ns",
            "extra": "gctime=119187213\nmemory=3382669968\nallocs=272208\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 882899126,
            "unit": "ns",
            "extra": "gctime=188518044\nmemory=3396010680\nallocs=441704\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3215848346.5,
            "unit": "ns",
            "extra": "gctime=807588706\nmemory=3372192016\nallocs=148020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 153244287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49383792\nallocs=732166\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 27470305723,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239149368\nallocs=1950842\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3258257667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116099448\nallocs=1657491\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1759079282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44060224\nallocs=332971\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3351480010,
            "unit": "ns",
            "extra": "gctime=527048355\nmemory=6666067536\nallocs=181823\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1780232282,
            "unit": "ns",
            "extra": "gctime=197555148\nmemory=6668505080\nallocs=197533\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1209111139,
            "unit": "ns",
            "extra": "gctime=214285048\nmemory=6673524744\nallocs=232643\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 5198071177,
            "unit": "ns",
            "extra": "gctime=569367311\nmemory=6664850912\nallocs=172484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265341221.5,
            "unit": "ns",
            "extra": "gctime=175018922\nmemory=27796792\nallocs=263079\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2180574102,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64448664\nallocs=529485\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1136992729,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40446752\nallocs=447982\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 727557861,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27469144\nallocs=215011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "eb0512871bb8f681dde735c191f1e3a075152d3b",
          "message": "Merge pull request #422 from JuliaHealth/fix-type-piracy\n\nFixing issues in Julia 1.11 and Julia 1.12",
          "timestamp": "2024-07-05T13:25:04-04:00",
          "tree_id": "f0a00e9f59558ceca807dab163496fb2ef70ba81",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/eb0512871bb8f681dde735c191f1e3a075152d3b"
        },
        "date": 1720203255075,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1905493127.5,
            "unit": "ns",
            "extra": "gctime=334592868.5\nmemory=3375681728\nallocs=189621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 914487178,
            "unit": "ns",
            "extra": "gctime=104862195\nmemory=3382672728\nallocs=272249\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 814946114,
            "unit": "ns",
            "extra": "gctime=163922344.5\nmemory=3396012776\nallocs=441722\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3299300532,
            "unit": "ns",
            "extra": "gctime=609020826\nmemory=3372194776\nallocs=148061\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 155313870,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49392696\nallocs=732231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 27987109817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239144824\nallocs=1950636\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3394142854,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116078920\nallocs=1657272\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1755297484,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44052280\nallocs=332881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3786600886,
            "unit": "ns",
            "extra": "gctime=1087371612\nmemory=6666145944\nallocs=183046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1794912707,
            "unit": "ns",
            "extra": "gctime=194705699\nmemory=6668583488\nallocs=198756\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1220462267,
            "unit": "ns",
            "extra": "gctime=196411304\nmemory=6673602816\nallocs=233844\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 5938874198,
            "unit": "ns",
            "extra": "gctime=422172859\nmemory=6664929320\nallocs=173707\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 270408353,
            "unit": "ns",
            "extra": "gctime=177803549.5\nmemory=27880320\nallocs=264322\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2166473718.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64527072\nallocs=530708\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1392987166,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40525848\nallocs=449227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 757974598,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27546512\nallocs=216221\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "cncastillo@uc.cl",
            "name": "Carlos Castillo Passi",
            "username": "cncastillo"
          },
          "committer": {
            "email": "cncastillo@uc.cl",
            "name": "Carlos Castillo Passi",
            "username": "cncastillo"
          },
          "distinct": true,
          "id": "a2837ae95f2f2a31564ee922b20e85c86d180be7",
          "message": "Bump to 0.9.0-DEV",
          "timestamp": "2024-07-09T23:42:30-04:00",
          "tree_id": "51c72fc79b968611a0f667e63bf3ddb8d1abf121",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/a2837ae95f2f2a31564ee922b20e85c86d180be7"
        },
        "date": 1720586664570,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1518177926,
            "unit": "ns",
            "extra": "gctime=136099989\nmemory=3375696264\nallocs=190090\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 1264921787.5,
            "unit": "ns",
            "extra": "gctime=263300426.5\nmemory=3382659008\nallocs=271788\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 586016460,
            "unit": "ns",
            "extra": "gctime=72552888.5\nmemory=3396006200\nallocs=441486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3944206835,
            "unit": "ns",
            "extra": "gctime=1004839515.5\nmemory=3372195144\nallocs=148070\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 154008923.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49386920\nallocs=732216\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 44650285812.5,
            "unit": "ns",
            "extra": "gctime=13675195\nmemory=239194216\nallocs=1952262\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3202722750,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115139040\nallocs=1621158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1785064601,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44052144\nallocs=332864\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3329807846,
            "unit": "ns",
            "extra": "gctime=524473651\nmemory=6666146312\nallocs=183055\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 2645750395,
            "unit": "ns",
            "extra": "gctime=455163688\nmemory=6668604696\nallocs=199438\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1053370396,
            "unit": "ns",
            "extra": "gctime=97551656\nmemory=6673597432\nallocs=233704\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6503649361.5,
            "unit": "ns",
            "extra": "gctime=731576996\nmemory=6664929688\nallocs=173716\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 267654789,
            "unit": "ns",
            "extra": "gctime=176089514\nmemory=27875568\nallocs=264311\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2930224071,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64575800\nallocs=532289\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1120090167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40427936\nallocs=445511\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 791311651,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27546888\nallocs=216231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "e7bfc1f7edc46fb73a8c6c03798692b31d41f6e8",
          "message": "Merge pull request #444 from JuliaHealth/fix-get-samples\n\nChange `reduce(vcat, itr)` to `reduce(vcat, [itr])`",
          "timestamp": "2024-07-12T18:03:18-04:00",
          "tree_id": "f7dee3d3990ad35dbfcf85a879bce50b6e7bb957",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/e7bfc1f7edc46fb73a8c6c03798692b31d41f6e8"
        },
        "date": 1720823253236,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1597348663.5,
            "unit": "ns",
            "extra": "gctime=213075219.5\nmemory=3365586504\nallocs=188829\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 920684554,
            "unit": "ns",
            "extra": "gctime=109214966\nmemory=3372577512\nallocs=271457\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 828881248,
            "unit": "ns",
            "extra": "gctime=168714278\nmemory=3385917864\nallocs=440932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3226260318,
            "unit": "ns",
            "extra": "gctime=597657156.5\nmemory=3362099560\nallocs=147269\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 154787415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39291336\nallocs=731415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 26694868553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=229058632\nallocs=1950104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3402176271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=105027304\nallocs=1620654\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1790112478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=33957088\nallocs=332096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3662737461,
            "unit": "ns",
            "extra": "gctime=865481298\nmemory=6663220696\nallocs=179871\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1725491180,
            "unit": "ns",
            "extra": "gctime=169277416\nmemory=6665658240\nallocs=195581\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1153779016,
            "unit": "ns",
            "extra": "gctime=186446719\nmemory=6670677392\nallocs=230658\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6320449643.5,
            "unit": "ns",
            "extra": "gctime=469672746.5\nmemory=6662004072\nallocs=170532\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265782897,
            "unit": "ns",
            "extra": "gctime=174350908\nmemory=24949952\nallocs=261127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2107401044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61596688\nallocs=527329\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1410314000,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37443848\nallocs=442337\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 724434160,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24621280\nallocs=213046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "454e0d6bda174da8fc42be959b9c227f1e2fc827",
          "message": "Merge pull request #455 from JuliaHealth/pkg-versions-and-tets-vscode\n\nPkg version handling and choosable test backend for VSCode (again)",
          "timestamp": "2024-07-15T16:23:33-04:00",
          "tree_id": "b210cf3f3bd7a4c6b73f2830e12cf85bd9f8f240",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/454e0d6bda174da8fc42be959b9c227f1e2fc827"
        },
        "date": 1721076291218,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1650175239,
            "unit": "ns",
            "extra": "gctime=169191644\nmemory=3365599424\nallocs=189119\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 929752408,
            "unit": "ns",
            "extra": "gctime=117208199\nmemory=3372583368\nallocs=271502\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 734991030,
            "unit": "ns",
            "extra": "gctime=139332277\nmemory=3385916584\nallocs=440750\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3164337887,
            "unit": "ns",
            "extra": "gctime=726548186.5\nmemory=3362102328\nallocs=147249\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 151580421.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39384728\nallocs=731749\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 26565103907,
            "unit": "ns",
            "extra": "gctime=0\nmemory=229059264\nallocs=1950059\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3280979167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=105035080\nallocs=1620615\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1799927920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=33959768\nallocs=332065\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3301650757.5,
            "unit": "ns",
            "extra": "gctime=507579833.5\nmemory=6663224440\nallocs=179872\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1730789585,
            "unit": "ns",
            "extra": "gctime=181102659\nmemory=6665661984\nallocs=195582\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1104013943.5,
            "unit": "ns",
            "extra": "gctime=166726967.5\nmemory=6670681472\nallocs=230680\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6349771663.5,
            "unit": "ns",
            "extra": "gctime=731478736.5\nmemory=6662006840\nallocs=170512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 264922273,
            "unit": "ns",
            "extra": "gctime=173649948\nmemory=24970128\nallocs=261175\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2111466088.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61597648\nallocs=527288\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1137546146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37451144\nallocs=442350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 757370815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24624136\nallocs=213026\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "1de7627abea287ed1c6127f06a6e2be29bca8dcb",
          "message": "Merge pull request #431 from gsahonero/simulate_docs_change\n\nIncluding details about prints of simulate",
          "timestamp": "2024-07-17T12:34:15-04:00",
          "tree_id": "406949e6605d66eb1d83970a27a83202c089dc26",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/1de7627abea287ed1c6127f06a6e2be29bca8dcb"
        },
        "date": 1721237083175,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 1866110763,
            "unit": "ns",
            "extra": "gctime=327879343\nmemory=3365592720\nallocs=188885\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 945776837,
            "unit": "ns",
            "extra": "gctime=120999970\nmemory=3372569584\nallocs=271033\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 911378432,
            "unit": "ns",
            "extra": "gctime=205815003\nmemory=3385923312\nallocs=440954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 3387409372.5,
            "unit": "ns",
            "extra": "gctime=773777302\nmemory=3362102328\nallocs=147249\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 153551114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39294104\nallocs=731395\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 26913270216,
            "unit": "ns",
            "extra": "gctime=0\nmemory=229061184\nallocs=1950083\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3261239667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=105029384\nallocs=1620437\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1817370709.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=33959936\nallocs=332075\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 3756525851,
            "unit": "ns",
            "extra": "gctime=898978255.5\nmemory=6663224440\nallocs=179872\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 1843011618,
            "unit": "ns",
            "extra": "gctime=221494082\nmemory=6665661984\nallocs=195582\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 1330616070,
            "unit": "ns",
            "extra": "gctime=253096202\nmemory=6670681488\nallocs=230681\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 6715606353,
            "unit": "ns",
            "extra": "gctime=770984273\nmemory=6662006840\nallocs=170512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 265940605.5,
            "unit": "ns",
            "extra": "gctime=174793240.5\nmemory=24952720\nallocs=261107\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 2149205885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61597944\nallocs=527289\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1409920896,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37451312\nallocs=442326\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 712193120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24624152\nallocs=213026\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "e09ba447591d16b8488da1e32dbe8332f67a4d70",
          "message": "Merge pull request #443 from JuliaHealth/optimize-precession\n\nOptimize run_spin_precession! and run_spin_excitation! for CPU",
          "timestamp": "2024-07-19T16:20:40-04:00",
          "tree_id": "890e11958f59e330869bd3aa95fba731e2fedb5b",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/e09ba447591d16b8488da1e32dbe8332f67a4d70"
        },
        "date": 1721421805685,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 244928608.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75067152\nallocs=714858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 140758946,
            "unit": "ns",
            "extra": "gctime=0\nmemory=110637312\nallocs=1327722\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 171370827,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181687440\nallocs=2554308\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 415982190.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57307952\nallocs=408666\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 138386450.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37888408\nallocs=674315\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 18715897786.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=203547904\nallocs=1725816\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 2988417291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92762472\nallocs=1451574\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1806396508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32636920\nallocs=294544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1141348868.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77260576\nallocs=995077\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 622082227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123354696\nallocs=1829990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 462315010,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215627240\nallocs=3501604\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2273113325,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54283224\nallocs=579133\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 255739261,
            "unit": "ns",
            "extra": "gctime=174160274.5\nmemory=24712992\nallocs=251406\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 1709103436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57496488\nallocs=490584\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1203673021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35531320\nallocs=414665\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 688213296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24388040\nallocs=206110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "24b21f4390f254ba53b042e934e91a627368af7e",
          "message": "Merge pull request #457 from JuliaHealth/fix-benchmark-pr\n\nFixing benchmark comments on PRs",
          "timestamp": "2024-07-19T18:41:01-04:00",
          "tree_id": "f82ecbad21a6d70773dd1d00bd7183762e03a66a",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/24b21f4390f254ba53b042e934e91a627368af7e"
        },
        "date": 1721432081975,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 235229966.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75067152\nallocs=714858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 140495905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=110637312\nallocs=1327722\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 169591756.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181687440\nallocs=2554307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 419227547,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57307944\nallocs=408665\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 135837984,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37888408\nallocs=674315\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 18356788557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=203547904\nallocs=1725816\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 2931106125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92756744\nallocs=1451468\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 1750964243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32636920\nallocs=294544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1174040352,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77260576\nallocs=995077\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 622515059.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123354696\nallocs=1829990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 492840880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215627240\nallocs=3501604\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2264093136,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54283224\nallocs=579133\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 257306603,
            "unit": "ns",
            "extra": "gctime=175798595\nmemory=24712992\nallocs=251406\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 1678945735.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57496488\nallocs=490584\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 1129875875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35530576\nallocs=414689\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 679066674,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24388040\nallocs=206110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "ryanakierulf@gmail.com",
            "name": "rkierulf",
            "username": "rkierulf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1457a4c3ae1e3e7cc40819c534d8c36620bccb75",
          "message": "Optimize run_spin_precession! for GPU (#459)\n\nOptimize run_spin_precession! for GPU",
          "timestamp": "2024-07-31T13:51:06-05:00",
          "tree_id": "398ae599b2a908ab6ef06a7daf9bc39a3c41c702",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/1457a4c3ae1e3e7cc40819c534d8c36620bccb75"
        },
        "date": 1722453168241,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 227517325.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624800\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 135033124,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101997264\nallocs=1139789\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 171880824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164652008\nallocs=2178490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 396561930.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54964448\nallocs=361646\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 138134905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37198504\nallocs=651145\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 14155999496.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=192731080\nallocs=1624525\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 3171338479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89581600\nallocs=1531211\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 75482754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31729600\nallocs=272376\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1168211452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71980896\nallocs=883236\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 612565463,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065576\nallocs=1606339\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 495427593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195319240\nallocs=3054333\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2245843835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51508264\nallocs=523197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 108701927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24423344\nallocs=240432\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 776956866,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49979192\nallocs=418547\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 769082459,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32250784\nallocs=381329\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 64232156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23826160\nallocs=190141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "ryanakierulf@gmail.com",
            "name": "rkierulf",
            "username": "rkierulf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1b6c5becc33bcf9d4b95a38b7c15a0f14cab564a",
          "message": "Optimize run_spin_excitation! for GPU (#462)",
          "timestamp": "2024-08-19T14:23:58-05:00",
          "tree_id": "9906d6c481fecd856f0f8afe53c6f23314a4c3c6",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/1b6c5becc33bcf9d4b95a38b7c15a0f14cab564a"
        },
        "date": 1724096743910,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 226897725.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624800\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 134821691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101997264\nallocs=1139789\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 143505624,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164652040\nallocs=2178491\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 343336904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54957664\nallocs=361440\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 56897573.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26574352\nallocs=197923\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 516675827,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53178912\nallocs=404912\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 561507500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32666752\nallocs=343656\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 36448199,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25961768\nallocs=147384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1157773672,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71980896\nallocs=883236\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 621102018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065576\nallocs=1606339\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 383824312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195319240\nallocs=3054333\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 1947221684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51446624\nallocs=521237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101669032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23400080\nallocs=196396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 649875888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36376768\nallocs=302157\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 565000312.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26752720\nallocs=266155\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 60318280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23284848\nallocs=178490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "dae2df143cb4418ceb497c97f4c1fa2ca23e74ef",
          "message": "Merge pull request #465 from gsahonero/fix/brain_phantom\n\nFixing brain phantom values",
          "timestamp": "2024-08-19T23:05:00-04:00",
          "tree_id": "a938a89d0a2f31404226c509e8b91bec0a8b01ef",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/dae2df143cb4418ceb497c97f4c1fa2ca23e74ef"
        },
        "date": 1724124401466,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 225453139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624800\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 135802285.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101997264\nallocs=1139789\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 145800343.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164652040\nallocs=2178491\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 347754103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54957464\nallocs=361428\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 56476716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26574352\nallocs=197923\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 514499972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53179152\nallocs=404915\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 566912396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32666048\nallocs=343635\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 36477480,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25961480\nallocs=147357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1157077065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71980896\nallocs=883236\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 624295507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065576\nallocs=1606339\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 383113892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195312400\nallocs=3054122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 1930369445.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51446624\nallocs=521237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101349123.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23400080\nallocs=196395\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 648132271.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36376768\nallocs=302158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 563862104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26753024\nallocs=266174\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 60560585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23284816\nallocs=178486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4d7519e0c4cc545e931ba86929c579af040f1fe3",
          "message": "Merge pull request #460 from JuliaHealth/compathelper/new_version/2024-08-09-00-40-44-455-03316379457\n\nCompatHelper: bump compat for AMDGPU in [weakdeps] to 1 for package KomaMRICore, (keep existing compat)",
          "timestamp": "2024-08-20T10:55:32-04:00",
          "tree_id": "ad92b4020922527d3cc676cf9866e96f52f3b9d6",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/4d7519e0c4cc545e931ba86929c579af040f1fe3"
        },
        "date": 1724167160007,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 243211040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624808\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 135699416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101983024\nallocs=1139383\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 146055228.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164652032\nallocs=2178491\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 408609577.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54964576\nallocs=361652\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 57044600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26574352\nallocs=197923\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 520674051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53185744\nallocs=405130\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 558646396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32666800\nallocs=343692\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 36907276.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26020936\nallocs=148501\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1018197659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71961128\nallocs=882631\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 623218694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065576\nallocs=1606339\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 382864848.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195319240\nallocs=3054333\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2261868516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51508264\nallocs=523197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 100807488.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23400080\nallocs=196396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 652815075.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36376768\nallocs=302158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 563682791,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26753184\nallocs=266184\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 61205236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23315576\nallocs=179040\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "57f8698e6cb72673ce947bc948fe8ddcbfa8952a",
          "message": "Merge pull request #461 from JuliaHealth/compathelper/new_version/2024-08-14-00-40-24-701-00541212617\n\nCompatHelper: bump compat for MRIReco to 0.9, (keep existing compat)",
          "timestamp": "2024-08-20T11:43:39-04:00",
          "tree_id": "5bc8ef5eb68380a603c49d02e5a702a2f10aa8c0",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/57f8698e6cb72673ce947bc948fe8ddcbfa8952a"
        },
        "date": 1724171650551,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 226618506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624816\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 174536994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101997128\nallocs=1139789\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 146360095.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164652040\nallocs=2178491\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 347644824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54957584\nallocs=361434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 57253633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26575888\nallocs=197971\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 515042255.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53180608\nallocs=404936\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 541353541,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32666376\nallocs=343653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 37619574.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26018296\nallocs=148351\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1024148878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71961128\nallocs=882631\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 580936747,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065480\nallocs=1606321\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 386777586,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195312400\nallocs=3054122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 1925568005.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51446640\nallocs=521237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 100754922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23401360\nallocs=196436\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 654922437.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36376768\nallocs=302163\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 564653500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26753440\nallocs=266188\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 60779232,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23315608\nallocs=179040\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "ryanakierulf@gmail.com",
            "name": "rkierulf",
            "username": "rkierulf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a26d5985c2843a9ab2ced4b623fbb6e75c243041",
          "message": "Add Distributed Examples to Documentation (#468)\n\n* Create 4-run-distributed-simulations.md\r\n\r\nAdd distributed examples in \"how-to\" section\r\n\r\n* Update 4-run-distributed-simulations.md\r\n\r\nCorrect typo\r\n\r\n* Add space\r\n\r\n* Make use case for multi-node simulation clearer\r\n\r\n* Add set_device! override without backend parameter\r\n\r\n* Add function for combining two RawAcquisitionData structs\r\n\r\n* Simplify example scripts\r\n\r\n* Remove sentence\r\n\r\n* Add isapprox function for testing\r\n\r\n* Add RawAcquisitionData test\r\n\r\n* Add images to assets folder\r\n\r\n* Update examples\r\n\r\n* Delete docs/src/assets/KomamultiGPU.svg\r\n\r\n* Delete docs/src/assets/KomamultiNode.svg\r\n\r\n* Delete docs/src/assets/KomamultiNodeCPU.svg\r\n\r\n* Add files via upload\r\n\r\n* Fix broken image\r\n\r\n* Use correct image for multiGPU\r\n\r\n* Update to use distributed macro",
          "timestamp": "2024-08-23T12:01:27-05:00",
          "tree_id": "636479d01ce73a3098391ec30e9fb9b20575dd00",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/a26d5985c2843a9ab2ced4b623fbb6e75c243041"
        },
        "date": 1724433849937,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 224991845.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70624808\nallocs=620867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 174863887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=101997120\nallocs=1139788\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 90572834,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164637904\nallocs=2178085\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 347400444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54957568\nallocs=361434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 57092571.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26599952\nallocs=198023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 522866366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53185760\nallocs=405130\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 568303458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32741992\nallocs=343064\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 36981171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26020472\nallocs=148469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1151416902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71980896\nallocs=883236\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 580233533,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113065288\nallocs=1606321\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 341164837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195312136\nallocs=3054104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 1930414196.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51446624\nallocs=521237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101438582.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23401104\nallocs=196399\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 636188156.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36376880\nallocs=302152\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 565478250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26781704\nallocs=265921\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 60929383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23315672\nallocs=179039\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "927de27202e5d99568948be68c97b4dc4da1c89a",
          "message": "Merge pull request #442 from JuliaHealth/new-motion\n\nNew Motion",
          "timestamp": "2024-09-19T14:49:25-03:00",
          "tree_id": "091bd6a1064dc7cb04c142720f90cf4ca4a21a43",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/927de27202e5d99568948be68c97b4dc4da1c89a"
        },
        "date": 1726769733211,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 224153374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70907904\nallocs=630309\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 135808380,
            "unit": "ns",
            "extra": "gctime=0\nmemory=102563472\nallocs=1158673\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 149398614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=165784448\nallocs=2216259\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 413189478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55106088\nallocs=366372\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 55925086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26301920\nallocs=184757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 496536737.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53136880\nallocs=404203\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 553973500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32984144\nallocs=360091\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 34739822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25980280\nallocs=147857\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1032666827.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73956480\nallocs=946967\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 624017416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=117042696\nallocs=1734605\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 388957396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=203273264\nallocs=3310848\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2245739282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52502520\nallocs=555262\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101367770.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23287088\nallocs=190819\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 636447326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36364616\nallocs=301972\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 554868916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26892760\nallocs=273406\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 58876104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23303768\nallocs=178854\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "5df34cbcda06cd53c6775b73932298213f3750c1",
          "message": "Merge pull request #483 from JuliaHealth/fix-extra-allocations\n\nFix extra allocations when benchmarking with no motion",
          "timestamp": "2024-09-20T12:52:04-03:00",
          "tree_id": "d35c242a6621c4c061216f21afb1bc2a2368a60a",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/5df34cbcda06cd53c6775b73932298213f3750c1"
        },
        "date": 1726849074107,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 223767376.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69134152\nallocs=585965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 130477217,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99001712\nallocs=1069579\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 138275342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158689408\nallocs=2038883\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 347398917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54212232\nallocs=343983\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 56227105,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26301920\nallocs=184757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 496678721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53136880\nallocs=404203\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 559269812.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32979728\nallocs=360111\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 34678068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25980296\nallocs=147858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1144007818,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70433440\nallocs=854570\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 609258944.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=109970664\nallocs=1549007\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 378398584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189122480\nallocs=2939452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 1928663456,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50672904\nallocs=506904\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101329308.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23287088\nallocs=190819\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 635643734,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36364616\nallocs=301972\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 554670167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26886992\nallocs=273298\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 58892774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23303784\nallocs=178855\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "pablo.villacorta@uva.es",
            "name": "Pablo Villacorta Aylagas",
            "username": "pvillacorta"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "356c48af40517719200f6518aadb6241a02852a8",
          "message": "Merge pull request #486 from JuliaHealth/fix-rotation\n\nFix sum of `Grad`s",
          "timestamp": "2024-09-26T19:09:10+02:00",
          "tree_id": "273ea2a88810046bd7baf90feeb325bfad39a9e3",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/356c48af40517719200f6518aadb6241a02852a8"
        },
        "date": 1727372140980,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 231918438,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69134160\nallocs=585965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 171093596,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99015816\nallocs=1069984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 87065540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158675280\nallocs=2038477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 404792264.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54219152\nallocs=344195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 55451684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26301920\nallocs=184757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 505814210,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53137072\nallocs=404215\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 539627125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32984272\nallocs=360175\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 34748766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25979288\nallocs=147795\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1007165039.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70413672\nallocs=853965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 576725538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=109970376\nallocs=1548989\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 331123481,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189122312\nallocs=2939440\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2261304965.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50734536\nallocs=508864\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 100941879.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23287088\nallocs=190819\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 649245415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36364712\nallocs=301978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 552172687.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26887216\nallocs=273313\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 59069032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23303608\nallocs=178844\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "pablo.villacorta@uva.es",
            "name": "Pablo Villacorta Aylagas",
            "username": "pvillacorta"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "48a5bbe6bf7a99cf33fa3698ee0a1a9069b3aeb8",
          "message": "Merge pull request #488 from JuliaHealth/fix-spinspan\n\nFix bugs related with `SpinRange` and flow",
          "timestamp": "2024-09-27T16:30:02+02:00",
          "tree_id": "e9ee8a0848df8ae0b27978d2e37925405a53d1be",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/48a5bbe6bf7a99cf33fa3698ee0a1a9069b3aeb8"
        },
        "date": 1727449742531,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 230260983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69134144\nallocs=585965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 130246287,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99001712\nallocs=1069579\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 141633208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158689416\nallocs=2038883\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 400677945.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54219152\nallocs=344195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 55535980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26301920\nallocs=184757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 495029548,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53136960\nallocs=404204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 541225208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32984176\nallocs=360204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 34691883.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25980664\nallocs=147881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1153307006.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70433440\nallocs=854570\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 605676119.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=109970664\nallocs=1549007\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 387857453,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189129416\nallocs=2939669\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2232808553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50734536\nallocs=508864\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101291917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23287088\nallocs=190819\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 632996422,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36364616\nallocs=301972\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 552117729,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26886320\nallocs=273207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 59154320.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23303736\nallocs=178852\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "ca742f69d7e60c1515396ceae739ad8cabd25891",
          "message": "Merge pull request #492 from JuliaHealth/koma-v0.9\n\nKomaMRI v0.9",
          "timestamp": "2024-10-01T13:03:25-03:00",
          "tree_id": "40512f85d948666799add85c8057234429ad7ef5",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/ca742f69d7e60c1515396ceae739ad8cabd25891"
        },
        "date": 1727800271279,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 223115374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69133552\nallocs=585949\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 132097936.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99001112\nallocs=1069563\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 146470163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158688792\nallocs=2038867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 370624128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54211632\nallocs=343967\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 55534616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26301304\nallocs=184741\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 498364290,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53136376\nallocs=404193\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 547774125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32983832\nallocs=360198\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 34679784,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25973728\nallocs=147470\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1155465921,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70432840\nallocs=854554\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 615942454.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=109970064\nallocs=1548991\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 418827220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189128816\nallocs=2939653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2098515633.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50672320\nallocs=506889\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 101206570.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23286440\nallocs=190800\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 637714919,
            "unit": "ns",
            "extra": "gctime=0\nmemory=36364016\nallocs=301952\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 553824458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26885240\nallocs=273235\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 59199152.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23302864\nallocs=178833\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4a1d50a7a4e45e6b105f153bdabcc5183a04f55b",
          "message": "Merge pull request #587 from JuliaHealth/fix-oneapi-indexing\n\nFix oneAPI indexing issues",
          "timestamp": "2025-07-11T18:21:58-07:00",
          "tree_id": "75899a646e11ec001be35e471f3730622b9c3687",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/4a1d50a7a4e45e6b105f153bdabcc5183a04f55b"
        },
        "date": 1752285557697,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 359457603,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837272\nallocs=625243\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 205463135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96363232\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 130918492.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=151663456\nallocs=2089761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 648595605,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944584\nallocs=381821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22507306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563704\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 77739627,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704328\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 246340063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27065728\nallocs=247066\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22682784,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22573184\nallocs=158144\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1827544273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019560\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 958711681,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996976\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 521711928.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191462760\nallocs=2976990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3592433967.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 33198120.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051552\nallocs=202623\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120296218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350160\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 142967208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23902216\nallocs=235150\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31094864.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4ee43463c1476245789e627bdddc8afe826691c0",
          "message": "Merge pull request #588 from JuliaHealth/codecov-buildkite\n\nAdd tags to Buildkite codecov process",
          "timestamp": "2025-07-13T22:17:15-07:00",
          "tree_id": "87d531e53a050c4ea7b16b73b8bea86d2d973fb3",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/4ee43463c1476245789e627bdddc8afe826691c0"
        },
        "date": 1752484782802,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 354807825.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68838488\nallocs=625231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 205584513.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96362912\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 147214074,
            "unit": "ns",
            "extra": "gctime=17359755\nmemory=151661872\nallocs=2089761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 644042817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944488\nallocs=381815\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22773783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563768\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76529189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704104\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 258768500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27065728\nallocs=247066\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22702779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22572896\nallocs=158143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1828102378.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019720\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 956719648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996912\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 523435881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191461288\nallocs=2976990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3589905726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32423524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051456\nallocs=202623\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120197654.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350192\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 146164291.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23902152\nallocs=235150\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31037952,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "6583dec24d40ad51a1293ac02cea0e3f404483b7",
          "message": "Merge pull request #567 from JuliaHealth/compathelper/new_version/2025-05-16-00-54-58-192-04056216492\n\nCompatHelper: bump compat for Interpolations to 0.16 for package KomaMRIBase, (keep existing compat)",
          "timestamp": "2025-07-14T13:57:49-07:00",
          "tree_id": "691a8ef3e35e37e2b2f66d0573f2c65a69b56bca",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/6583dec24d40ad51a1293ac02cea0e3f404483b7"
        },
        "date": 1752532987219,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 356274474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837144\nallocs=625231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 205857034,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96361408\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 146882281,
            "unit": "ns",
            "extra": "gctime=16584771\nmemory=151658080\nallocs=2089765\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 644386281,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944488\nallocs=381815\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 23114847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563640\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76761088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704488\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 244550125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27026400\nallocs=245633\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 23381530,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22569984\nallocs=158100\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1827796208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019528\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 957176780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996464\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 538082166,
            "unit": "ns",
            "extra": "gctime=17652739\nmemory=191459928\nallocs=2976994\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3584231528,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 33608659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051440\nallocs=202622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 119909948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350224\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 147911375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24010936\nallocs=239380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31811289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "ff2c2be0f37f6f14630a32e150a7d4f2b6b8172b",
          "message": "Merge pull request #589 from JuliaHealth/codecov-builkite-change\n\nFlags added to official julia-coverage",
          "timestamp": "2025-07-16T13:03:42-07:00",
          "tree_id": "f8c83369ba3319688c675af8e666b8e943fc115c",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/ff2c2be0f37f6f14630a32e150a7d4f2b6b8172b"
        },
        "date": 1752698593339,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 355010524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837272\nallocs=625231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 206532685.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96362496\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 132315816,
            "unit": "ns",
            "extra": "gctime=0\nmemory=151655296\nallocs=2089761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 643700677.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944584\nallocs=381816\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22271830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563448\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 78502763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704872\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 273508916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27023328\nallocs=245621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22809761.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22573184\nallocs=158144\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1828998254.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019624\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 958242582.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996496\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 522380322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191458504\nallocs=2976990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3575610820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32663759.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051424\nallocs=202621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 119952601,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350288\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 141700292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24008152\nallocs=239370\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31498067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "024d87d67a39209c937b7d99faa0775a95098b12",
          "message": "Merge pull request #594 from JuliaHealth/simplify-heart-phantom\n\nReduce the number of points (N) in `heart_phantom`  from 21 to 10",
          "timestamp": "2025-07-29T18:03:38-07:00",
          "tree_id": "f282e9f9fc3bdeb2c70a65fadf1a05e6b248cb3d",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/024d87d67a39209c937b7d99faa0775a95098b12"
        },
        "date": 1753899988534,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 359990543,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837240\nallocs=625237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 205836148,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96363040\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 148852177,
            "unit": "ns",
            "extra": "gctime=17722056\nmemory=151663552\nallocs=2089765\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 653574554,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944584\nallocs=381821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22490861,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563640\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76595674,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704616\nallocs=250234\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96807458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26273600\nallocs=225209\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22763178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22575808\nallocs=158186\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1829167772,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019528\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 958795872.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996752\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 540296550,
            "unit": "ns",
            "extra": "gctime=18779385\nmemory=191461784\nallocs=2976994\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3582546657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32748607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051456\nallocs=202623\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120058664,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350160\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 114030958.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23664688\nallocs=230332\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31456911,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "d68d7d572030796367779a2bf705fd9630a28cd9",
          "message": "Merge pull request #593 from JuliaHealth/fix-rotation-center\n\nChange rotation center from geometric to center of gravity",
          "timestamp": "2025-08-04T10:07:32-07:00",
          "tree_id": "70528a7d88ab2446aa7d634d65c554465ed81ef2",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/d68d7d572030796367779a2bf705fd9630a28cd9"
        },
        "date": 1754432522928,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 359548396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837112\nallocs=625243\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 206093932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96362720\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 146894502,
            "unit": "ns",
            "extra": "gctime=16692319\nmemory=151653504\nallocs=2089765\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 651826277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944584\nallocs=381821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22515131,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563576\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76695663.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28704776\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 95336250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26270528\nallocs=225197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22675457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22575808\nallocs=158186\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1828565461,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019464\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 959094384,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115997104\nallocs=1584836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 537487147.5,
            "unit": "ns",
            "extra": "gctime=18578313.5\nmemory=191459480\nallocs=2976988\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3587326583,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32711926,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051520\nallocs=202623\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120246642,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350288\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 111782521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23661904\nallocs=230321\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31222080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "f87059bc262d838c8faa3df0b9bf154d8c4acfef",
          "message": "Merge pull request #592 from JuliaHealth/compathelper/new_version/2025-07-19-00-57-51-575-00327139571\n\nCompatHelper: bump compat for AMDGPU in [weakdeps] to 2 for package KomaMRICore, (keep existing compat)",
          "timestamp": "2025-08-05T20:15:54-07:00",
          "tree_id": "8ad9c688488f53b300cde837d2fc62e86d1c32fd",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/f87059bc262d838c8faa3df0b9bf154d8c4acfef"
        },
        "date": 1754452432046,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 356743374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68837464\nallocs=625237\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 206193126,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96362656\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 148115920,
            "unit": "ns",
            "extra": "gctime=17221462\nmemory=151656064\nallocs=2089765\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 645904121,
            "unit": "ns",
            "extra": "gctime=0\nmemory=57944584\nallocs=381821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22901517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563576\nallocs=166158\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 78406187,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28707208\nallocs=250232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 97127709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26273600\nallocs=225209\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22461985,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22561536\nallocs=157972\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1828061790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019752\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 959221031,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115995800\nallocs=1584832\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 539744578.5,
            "unit": "ns",
            "extra": "gctime=18180384.5\nmemory=192603744\nallocs=2977217\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3579575351,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123424\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32706449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051504\nallocs=202622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 119677988,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24350288\nallocs=244736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 113714813,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23665040\nallocs=230334\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31119335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "cf3718bfbf524fd75be9eeadff7474447544ece6",
          "message": "Merge pull request #602 from JuliaHealth/bump-version\n\nbump v0.9.1 to v0.9.2",
          "timestamp": "2025-08-27T17:02:25-07:00",
          "tree_id": "1915fc07977a0068b81b7a196ab9be00b4a085f7",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/cf3718bfbf524fd75be9eeadff7474447544ece6"
        },
        "date": 1756341901517,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 325132652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68836640\nallocs=625214\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 230868276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96361024\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 207313117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153651648\nallocs=2090230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 551166951.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56797136\nallocs=381548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22623461,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563416\nallocs=166162\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 74790664,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29769832\nallocs=257961\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96525021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26273600\nallocs=225209\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22430419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22576096\nallocs=158187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1553231608.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019544\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 840490227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115995552\nallocs=1584776\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 531658882,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191458392\nallocs=2976954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3017924449.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123408\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32816501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051552\nallocs=202627\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120603090.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25091152\nallocs=250112\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 113699895.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23664720\nallocs=230333\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31337762.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "9c416d13dfa434113890515c887da5d5584bafda",
          "message": "Merge pull request #605 from JuliaHealth/fix-motion-constructors\n\nChange `Motion` function names to avoid conflicts with `Action` structure names",
          "timestamp": "2025-08-29T08:38:00-07:00",
          "tree_id": "31ee09ac4257734d24ce10e5a96b1e0b82a7d91a",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/9c416d13dfa434113890515c887da5d5584bafda"
        },
        "date": 1756489454949,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 328441559,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68836928\nallocs=625213\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 218426667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96361856\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 255465467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153654504\nallocs=2090231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 548446867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56797136\nallocs=381548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 23067343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563608\nallocs=166162\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 73564443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29769320\nallocs=257960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 93898395.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26270528\nallocs=225197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22773122,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22576096\nallocs=158187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1555127061.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019608\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 841536135.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115995704\nallocs=1584776\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 511540086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191459160\nallocs=2976954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3008203465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123408\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 33307487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051392\nallocs=202626\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120150115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25091312\nallocs=250112\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 110470250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23661904\nallocs=230321\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31494449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "cncastillo@uc.cl",
            "name": "Carlos Castillo Passi",
            "username": "cncastillo"
          },
          "committer": {
            "email": "cncastillo@uc.cl",
            "name": "Carlos Castillo Passi",
            "username": "cncastillo"
          },
          "distinct": true,
          "id": "f6fabaeed92fb14a0d84801f7f4a99ea6d96f6b4",
          "message": "Bumped version",
          "timestamp": "2025-08-29T13:26:31-07:00",
          "tree_id": "a42eb74dbca72d11a7a51b3770ae480b8d991ba9",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/f6fabaeed92fb14a0d84801f7f4a99ea6d96f6b4"
        },
        "date": 1756504565744,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 325414319.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68836656\nallocs=625219\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 208540970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96362888\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 279084255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153651816\nallocs=2090230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 549151693,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56796880\nallocs=381548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 23070514,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22563768\nallocs=166162\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 73495647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29769576\nallocs=257960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96813875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26270528\nallocs=225197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22815788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22534752\nallocs=157583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1547690094.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019448\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 838069253,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115996128\nallocs=1584776\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 508217138.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191458264\nallocs=2976954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3001947386,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123408\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32838450.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21053472\nallocs=202627\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 120576823,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25091056\nallocs=250112\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 114058083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23661616\nallocs=230320\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31626507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "3e8d11d752c411c204f9ad720d91623b6c5f3dd2",
          "message": "Merge pull request #606 from JuliaHealth/dependabot/github_actions/actions/checkout-5\n\nBump actions/checkout from 4 to 5",
          "timestamp": "2025-09-01T15:36:30-07:00",
          "tree_id": "c1e068fde4e17575b6062f2134e5270c35ffd043",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/3e8d11d752c411c204f9ad720d91623b6c5f3dd2"
        },
        "date": 1756768596324,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 328282505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68836656\nallocs=625213\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 229151622.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96360648\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 256871711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153650976\nallocs=2090230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 546203153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56797392\nallocs=381548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22515996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22564152\nallocs=166162\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 78733342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29773384\nallocs=257960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 93634584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26270816\nallocs=225198\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22578298,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22534752\nallocs=157583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1558054826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019480\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 847425190,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115995392\nallocs=1584776\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 516251799,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191458872\nallocs=2976954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3003783230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123408\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 33505273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21051584\nallocs=202627\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 124117518,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25091344\nallocs=250113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 110128896,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23661616\nallocs=230320\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31360909,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21145624\nallocs=201089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "8b62195432077812f3eff37e1b9cf66b206120bd",
          "message": "Merge pull request #617 from JuliaHealth/fix-metal-1_8\n\nChange order of thread reduction so it is done on the CPU",
          "timestamp": "2025-09-18T14:50:23-07:00",
          "tree_id": "4df5cba8b4bf169cfcd1e151942ed55e6745328d",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/8b62195432077812f3eff37e1b9cf66b206120bd"
        },
        "date": 1758234510731,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 328508190,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68836728\nallocs=625214\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 225624482,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96361152\nallocs=1113357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 163355495,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153656112\nallocs=2090230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 548280455,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56797552\nallocs=381548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22724945,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22576768\nallocs=166128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 77230670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30986152\nallocs=260049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96602916.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26230080\nallocs=225289\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22289157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22583464\nallocs=158099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1558872240,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80019512\nallocs=889127\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 859609769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115995112\nallocs=1584776\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 522488232.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191459608\nallocs=2976954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3010268970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72123440\nallocs=543430\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32524226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21132072\nallocs=202592\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 123943367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25989744\nallocs=251388\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 116308437.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23703280\nallocs=230272\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31118611,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21223536\nallocs=201043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4ae2dbdd21286a11e0490d7979a34c9bc57acde9",
          "message": "Merge pull request #558 from aTrotier/extension_label\n\nImplement labels",
          "timestamp": "2025-09-29T11:47:57-07:00",
          "tree_id": "a13cc73b6a0d3c9ece4a7d1149f2a3c5a5b7f203",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/4ae2dbdd21286a11e0490d7979a34c9bc57acde9"
        },
        "date": 1759173896104,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 319052626.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68897568\nallocs=625735\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 218059671,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96428880\nallocs=1113879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 211203146.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153711600\nallocs=2090752\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 548521174.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56857904\nallocs=382070\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22403464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22641568\nallocs=166722\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76763095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31047528\nallocs=260571\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96699896,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26291040\nallocs=225811\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22546962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22605912\nallocs=158059\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1567310315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80299312\nallocs=890971\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 855897469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116278880\nallocs=1586621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 540102520,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191738320\nallocs=2978798\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3016874304,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72403240\nallocs=545274\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 33213097.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21461120\nallocs=205345\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 123795074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26269800\nallocs=253233\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 113197521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23983552\nallocs=232116\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 30871794.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21503368\nallocs=202887\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "34935af3b5c0d61c30c6b6645169c9a3dd4850d7",
          "message": "Merge pull request #621 from JuliaHealth/label-support-bump\n\nLabel support version bump",
          "timestamp": "2025-10-06T21:41:20-07:00",
          "tree_id": "cf68d81d837d164c2f62f71f829673d0e045020a",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/34935af3b5c0d61c30c6b6645169c9a3dd4850d7"
        },
        "date": 1759814537336,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 325159954,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68897824\nallocs=625735\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 208094287,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96422464\nallocs=1113879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 263077165,
            "unit": "ns",
            "extra": "gctime=0\nmemory=153713064\nallocs=2090753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 551195809,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56857872\nallocs=382070\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 22727765,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22641568\nallocs=166722\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 78201803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31046824\nallocs=260571\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96494625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26291040\nallocs=225811\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 23157576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22644424\nallocs=158621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1567874474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80299472\nallocs=890971\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 842065822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116275768\nallocs=1586620\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 522460610,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191738352\nallocs=2978798\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 2999810477,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72403240\nallocs=545274\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32229083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21460672\nallocs=205346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 123823465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26269960\nallocs=253233\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 113251812.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23983552\nallocs=232116\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 31324354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21503368\nallocs=202887\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "04d8576700c24bd19b86461dd46b698bf695ea83",
          "message": "Merge pull request #623 from JuliaHealth/asym-periodic\n\nSupport fully asymmetric `TimeCurve`s with `Periodic` constructor",
          "timestamp": "2025-10-14T08:06:23-07:00",
          "tree_id": "1e307109b30585a4f2bd703ecf5a763b92e2e07f",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/04d8576700c24bd19b86461dd46b698bf695ea83"
        },
        "date": 1760459298592,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 339142516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=65769608\nallocs=624857\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 275457804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94390104\nallocs=1113143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 229902425.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152416744\nallocs=2133447\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 548960033,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51474144\nallocs=380908\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21260701.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22519360\nallocs=159643\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 78948742,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30894248\nallocs=256066\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 93939709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25738640\nallocs=187516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 21976880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22704240\nallocs=163708\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1581887656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69161912\nallocs=888935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 879708727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108586712\nallocs=1585011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 558294982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191433624\nallocs=3222656\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3008922484,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49493776\nallocs=541818\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32541458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21479856\nallocs=206006\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 124783584.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26251992\nallocs=252000\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 111314250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23777216\nallocs=214872\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 28652857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21583920\nallocs=207807\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "424582b967f7a19c3f52e501ac68ab3218db0458",
          "message": "Merge pull request #636 from JuliaHealth/fix-times\n\nFix bug in `times` function when `periods` is not a scalar",
          "timestamp": "2025-10-16T09:34:54-07:00",
          "tree_id": "71306a5d41b0aff2a9f5ff326efe207bfea84770",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/424582b967f7a19c3f52e501ac68ab3218db0458"
        },
        "date": 1760635161929,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 339871015.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=65769832\nallocs=624858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 274858100.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94392520\nallocs=1113143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 211362275.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=151940696\nallocs=2095270\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 549528040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51474384\nallocs=380908\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21524941.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22512080\nallocs=159542\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 79512564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30894984\nallocs=256066\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 95261562.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25704304\nallocs=186437\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 22078766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22761728\nallocs=166428\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1575008589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69162232\nallocs=888935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 883136236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108587608\nallocs=1585017\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 562279327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187535912\nallocs=2978511\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3008745113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49494032\nallocs=541818\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32417441,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21474832\nallocs=205982\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 124873710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26252152\nallocs=252000\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 112302083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23762768\nallocs=214501\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 28728289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21588096\nallocs=207932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "b5c63a9277de82f46335d2d59577cd3c12320e52",
          "message": "Benchmarks: Update rocmgpu setting to allow all GPUs (#653)",
          "timestamp": "2025-11-23T22:12:15-08:00",
          "tree_id": "b7d0a5a6791bc12fba2ab7da4686ce7dc8f26651",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/b5c63a9277de82f46335d2d59577cd3c12320e52"
        },
        "date": 1763966861444,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 340960290,
            "unit": "ns",
            "extra": "gctime=0\nmemory=66424856\nallocs=625347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 272489189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95049160\nallocs=1114095\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 212102146.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152457672\nallocs=2094772\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 553123722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52132336\nallocs=381197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21705981,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23174880\nallocs=159961\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 79470326.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31857080\nallocs=258911\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 95738917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25934384\nallocs=186742\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 26055313,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23427040\nallocs=168196\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1588926486.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68796928\nallocs=889490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 886294756.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108223536\nallocs=1585732\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 561392306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187150720\nallocs=2979241\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3027587378,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49129064\nallocs=542307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32637155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21112696\nallocs=206484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 124598595.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26098304\nallocs=254258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 112846167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23323072\nallocs=215003\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 34082841,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21206744\nallocs=208434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3ca859f23596a8b7e01d1ab49ae9b92aaa123853",
          "message": "Bump actions/checkout from 5 to 6 (#656)\n\nBumps [actions/checkout](https://github.com/actions/checkout) from 5 to 6.\n- [Release notes](https://github.com/actions/checkout/releases)\n- [Changelog](https://github.com/actions/checkout/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/actions/checkout/compare/v5...v6)\n\n---\nupdated-dependencies:\n- dependency-name: actions/checkout\n  dependency-version: '6'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2025-12-05T10:17:19-08:00",
          "tree_id": "04966037f7eeda2c0adb7aff143d025439236581",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/3ca859f23596a8b7e01d1ab49ae9b92aaa123853"
        },
        "date": 1764964141181,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 338417435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=66425864\nallocs=625347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 280440621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95105816\nallocs=1115247\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 275342094,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152561624\nallocs=2097077\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 563429680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52132240\nallocs=381196\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21768197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23171760\nallocs=159792\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76269687,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31854584\nallocs=258911\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96584416.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25939552\nallocs=186940\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 25191780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23434800\nallocs=168452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1592072215,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68797312\nallocs=889490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 885767416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108234640\nallocs=1586075\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 561506775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187169920\nallocs=2979941\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3085344184.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49129064\nallocs=542307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32118358,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21112520\nallocs=206473\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 121160780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26097600\nallocs=254258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 112343520.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23256384\nallocs=212708\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 32825642,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21206856\nallocs=208451\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ee3bc2fec1c9e8133c3bad86083bba195589dcf3",
          "message": "CompatHelper: bump compat for MAT to 0.11 for package KomaMRIBase, (keep existing compat) (#649)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>",
          "timestamp": "2025-12-05T13:08:53-08:00",
          "tree_id": "2e8ac2b02d799274f1f8be0f0fe7b242ccf35fc9",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/ee3bc2fec1c9e8133c3bad86083bba195589dcf3"
        },
        "date": 1764974728313,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 341921413,
            "unit": "ns",
            "extra": "gctime=0\nmemory=66425112\nallocs=625347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 275323882,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95047928\nallocs=1114159\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 285513733,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152420280\nallocs=2094772\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 552505595.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52132624\nallocs=381200\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21193443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23173088\nallocs=159895\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76632822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31855000\nallocs=258911\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 96112375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25939472\nallocs=186939\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 28621139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23403344\nallocs=167146\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1622705599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68796800\nallocs=889490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 891208860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108228416\nallocs=1585977\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 571183580.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187146624\nallocs=2979247\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3037417759.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49129064\nallocs=542307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32085949,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21112328\nallocs=206470\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 121193856.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26097824\nallocs=254258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 112065187.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23256672\nallocs=212709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 36595885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21206200\nallocs=208407\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "aadb879f6951c8718aa64883432b6946d43aff1e",
          "message": "CompatHelper: bump compat for MAT to 0.11, (keep existing compat) (#648)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>",
          "timestamp": "2025-12-05T13:10:06-08:00",
          "tree_id": "b0d843ca409733c570fdebd896ee987b906ac403",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/aadb879f6951c8718aa64883432b6946d43aff1e"
        },
        "date": 1765222489594,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 335355680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=66425544\nallocs=625352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 279185962.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95047608\nallocs=1114095\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 240103346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152498424\nallocs=2096748\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 554125362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52132672\nallocs=381196\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21223767,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23172784\nallocs=159860\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 76121455.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31854680\nallocs=258911\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 91941833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25939648\nallocs=186940\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 24418140.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23428496\nallocs=168181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1605518511,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68797184\nallocs=889490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 892280712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108223072\nallocs=1585726\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 557147082.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187151744\nallocs=2979767\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3068786029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49129128\nallocs=542307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32603299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21112472\nallocs=206478\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 121109825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26097600\nallocs=254258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 111506063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23256848\nallocs=212710\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 32527276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21206888\nallocs=208454\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "pablo.villacorta@uva.es",
            "name": "Pablo Villacorta Aylagas",
            "username": "pvillacorta"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "11d3180b2fc8717e5d695d9edd80ff180d3a1ad8",
          "message": "Add key time points in `DiscreteSequence` when evaluating periodic and flow-related motions (#638)\n\n* First commit\n\n* Remove debug messages and update docstrings\n\n* Solve bug in `export_2_mat_raw` with non-ASCII keys\n\n* Code coverage and method dispatch for motion comparison\n\n* Address requested changes:\n- Avoid type pyracy (new `CenterOfMass` struct and related `≈` methods)\n- Remove `sampling_params[\"Δt_rise\"]` parameter\n- Use `MIN_RISE_TIME` KomaMRIBase constant\n- Clarify tests\n- Improve coverage\n\n* Commit of the following:\n- Improve `==` and `≈` definitions: remove redundant methods and dispatch basing on the abstract type\n- Solve bug when reading writing the rotation center into a phantom file\n\n* KomaMRIFiles coverage\n\n* Address requested changes:\n- Fix type pyracy\n- Fix and test `RotateX` functions\n- Fix type stability in `displacement!` (new `get_center` function)\n- Change variable name `per` -> `periods`\n- Improve test clarity\n- Add tests: `times` and non-periodic `add_key_time_points!`\n\n* Revert type restrictions (since `displacement!` functions are not exported)\n\n* Test `RotateX` and displaced rotations\n\n* Update file format svg files with rotation center information\n\n* Test timeout in docs.\n\n* Optimize `add_key_time_points!` for periodic cases\n\n* Bump base, core and files\n\n* Update Project.toml\n\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>\n\n* Update Project.toml\n\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>\n\n* Remove unnecesary line from CI.yml\n\n* Bump KomaBase version in KomaPlots/Project.toml\n\n* Undo KomaBase bump in KomaPlots\n\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>\n\n---------\n\nCo-authored-by: Carlos Castillo Passi <cncastillo@uc.cl>",
          "timestamp": "2025-12-08T13:49:21-08:00",
          "tree_id": "3d68f90fd2a44c5e2dcecdfbe98a243d0c38b483",
          "url": "https://github.com/JuliaHealth/KomaMRI.jl/commit/11d3180b2fc8717e5d695d9edd80ff180d3a1ad8"
        },
        "date": 1765232862749,
        "tool": "julia",
        "benches": [
          {
            "name": "MRI Lab/Bloch/CPU/2 thread(s)",
            "value": 337633469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=66424872\nallocs=625348\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/4 thread(s)",
            "value": 276282606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95046184\nallocs=1114095\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/8 thread(s)",
            "value": 209852165,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152431128\nallocs=2094572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/CPU/1 thread(s)",
            "value": 555065232,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52132400\nallocs=381200\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/CUDA",
            "value": 21231134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23174320\nallocs=159976\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/oneAPI",
            "value": 77053690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31854552\nallocs=258911\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/Metal",
            "value": 95540333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25939552\nallocs=186940\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "MRI Lab/Bloch/GPU/AMDGPU",
            "value": 24763685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23438656\nallocs=168481\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/2 thread(s)",
            "value": 1592087059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68796896\nallocs=889490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/4 thread(s)",
            "value": 889539516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=108222896\nallocs=1585737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/8 thread(s)",
            "value": 565401528.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187145888\nallocs=2979565\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/CPU/1 thread(s)",
            "value": 3029269071,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49129096\nallocs=542307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/CUDA",
            "value": 32639703,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21109672\nallocs=206287\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/oneAPI",
            "value": 121489006.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26097664\nallocs=254258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/Metal",
            "value": 111152750,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23256416\nallocs=212709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Slice Selection 3D/Bloch/GPU/AMDGPU",
            "value": 32636098,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21207016\nallocs=208464\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":true,\"seconds\":120,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}