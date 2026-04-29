# KomaMRI.jl

## Git
- Use a new branch for new features. Base it on `master`.
- Never push to `master`.
- If the current worktree has unrelated uncommitted changes and you need a clean branch, prefer a separate worktree instead of stashing or rewriting the user's changes.
- If the branch should start from the latest upstream `master`, fetch first and branch from `origin/master` rather than a potentially stale local `master`.
- Keep PRs scoped to one issue or feature when possible.

## Repo Shape
- KomaMRI is a monorepo. The umbrella package lives at the repo root.
- Subpackages: `KomaMRIBase`, `KomaMRICore`, `KomaMRIFiles`, `KomaMRIPlots`.
- The root `KomaMRI` package directly depends on `KomaMRICore`, `KomaMRIFiles`, and `KomaMRIPlots`.
- `KomaMRICore`, `KomaMRIFiles`, and `KomaMRIPlots` depend on `KomaMRIBase`.
- Each package has its own `Project.toml`. Test envs live under `*/test/Project.toml`. Docs use `docs/Project.toml`.
- Do not edit `Manifest.toml` directly. Let `Pkg` own generated manifest state.
- `Manifest.toml` churn is usually not commit-worthy here. Manifests are mostly ignored in this repo, so do not commit manifest changes unless explicitly asked.

## Julia Workflow
- Use idiomatic Julia, not Python habits.
- Treat the relevant project directory as the environment. Use `julia --project=<path>` and `Pkg.activate(...)`.
- Prefer project-local environments over the global default environment.
- Do not create alternate env roots, `env/` folders, or `DEPOT` / `HOME` hacks unless explicitly asked.
- Reuse a persistent Julia REPL with Revise. Do not start a fresh REPL for convenience, precompilation, or one-off commands.
- Only restart the REPL if none exists, it crashed, or the user asks. Say why first.
- Keep the REPL open until explicitly told to close it.
- Keep diffs tight. Do not reformat or "clean up" unrelated code.
- Read `Project.toml` when needed, but change dependencies and compat via `Pkg` APIs, not by hand-editing dependency state casually.

## Local Package Wiring
- This repo defines a root `[workspace]` with `test`, `docs`, `benchmarks`, `KomaMRIBase`, `KomaMRICore`, `KomaMRIFiles`, and `KomaMRIPlots`.
- On modern Julia/Pkg workspaces, prefer the workspace-first path instead of blindly `Pkg.develop`-ing every local subpackage:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
- If you want a child project active, activate that child directly; it still participates in the shared workspace manifest:
```julia
using Pkg
Pkg.activate("KomaMRICore")  # or docs / benchmarks / another workspace child
Pkg.instantiate()
```
- This also applies to docs:
```julia
using Pkg
Pkg.activate("docs")
Pkg.instantiate()
```
- Use explicit `Pkg.develop(path=...)` only when you need path-tracked local packages outside the normal workspace flow, or when matching older/compatibility-oriented setup.
- If you need to mirror the current CI/docs wiring exactly, use the explicit path setup from `.github/workflows/CI.yml`, e.g. for docs:
```julia
using Pkg
Pkg.activate("docs")
Pkg.develop([
    PackageSpec(path=pwd(), subdir="."),
    PackageSpec(path=pwd(), subdir="KomaMRIBase"),
    PackageSpec(path=pwd(), subdir="KomaMRICore"),
    PackageSpec(path=pwd(), subdir="KomaMRIFiles"),
    PackageSpec(path=pwd(), subdir="KomaMRIPlots"),
])
Pkg.instantiate()
```
- CI still uses explicit `Pkg.develop(...)` wiring today because the matrix includes Julia `1.10`, so use that path when reproducing CI exactly.

## Julia Code
- Do not be verbose. Prefer the shortest clear correct Julia solution that matches the surrounding style.
- Prioritize elegance and concision. Do not lie or take shortcuts.
- Prefer the simplest correct design. Do not add helpers, wrappers, aliases, or abstractions unless they materially improve clarity, reuse, or dispatch.
- Avoid review-hostile churn: unnecessary boilerplate, defensive overengineering, speculative abstractions, and large refactors when a small local change is enough.
- Do not overtype function arguments; add method types only when dispatch needs them.
- Do not use abstract types in structs.
- Follow existing repo style and BlueStyle where it matches the surrounding code.
- Add comments or docstrings only when they carry real information.
- Watch type stability. Use `@code_warntype` when warranted.
- For performance work, profile first, then benchmark with interpolation, then optimize the actual hotspot.

## Julia Testing
- Use the correct project context. Wrong env selection causes misleading failures.
- No need to open a new REPL for testing; use `Pkg.test`.
- Prefer `Pkg.test(...)` over directly running `test/runtests.jl`. Only run `runtests.jl` directly when debugging the test runner itself.
- Run the narrowest relevant tests first. If the change crosses package boundaries, test each affected package and then the umbrella package if needed.
- From the repo root, the usual entry points are:
```julia
import Pkg
Pkg.test("KomaMRIBase")
Pkg.test("KomaMRICore")
Pkg.test("KomaMRIFiles")
Pkg.test("KomaMRIPlots")
Pkg.test()  # KomaMRI umbrella package
```
- `Pkg.test()` at the repo root exercises the `KomaMRI` package, including UI-oriented paths. On headless Linux, mirror CI and use `xvfb-run` for `KomaMRIPlots` and root `KomaMRI` tests.
- `KomaMRICore` defaults to CPU via `KomaMRICore/test/Project.toml`. To change backend, set `[preferences.KomaMRICore].test_backend` to `CPU`, `CUDA`, `AMDGPU`, `Metal`, or `oneAPI`, or use `Pkg.test("KomaMRICore"; test_args=\`CUDA\`)`.
- To control CPU thread count for `KomaMRICore`, use `Pkg.test("KomaMRICore"; julia_args=\`--threads=4\`)`.
- GPU backend packages must already be available locally before running non-CPU `KomaMRICore` tests. `oneAPI` is experimental.
- Do not commit backend-specific GPU packages such as `CUDA`, `AMDGPU`, `Metal`, or `oneAPI` to shared test projects like `KomaMRICore/test/Project.toml`. Add them only in the backend-specific CI/local test command or a temporary environment.
- For small exploratory tests, use a temporary environment. For longer isolated tests, use a temp folder with an activated environment. Reuse the same REPL.
- If you add behavior, add or update tests in the owning package's `runtests.jl`.
- Use `@testitem` tags correctly:
  - `KomaMRIBase`: `:base`
  - `KomaMRICore`: `:core` plus exactly one of `:motion` or `:nomotion`
  - `KomaMRIFiles`: `:files`
  - `KomaMRIPlots`: `:plots`
  - `KomaMRI`: `:koma`
- If the user asks for benchmarking, run benchmarks sequentially, never in parallel.
- If the user asks to test a PR implementation, add that version with `Pkg` using `url` / `rev`; do not hand-copy files.

## Docs
- For docs changes, use the `docs` environment.
- KomaMRI docs use Documenter.jl, DocumenterVitepress, Literate.jl, and PlutoSliderServer.
- Treat the source files as canonical: edit the hand-written docs and the `lit-*.jl` / `pluto-*.jl` sources, not the generated tutorial artifacts.
- Keep `DocumenterVitepress.deploydocs(...)` in `docs/make.jl`; do not replace it with `Documenter.deploydocs(...)`.
- Run doctests with:
```bash
julia --project=docs -e 'using Documenter: doctest; using KomaMRI; doctest(KomaMRI)'
```
- Build docs locally with:
```bash
julia --project=docs docs/make.jl
```
- For faster local docs iteration, you can temporarily set `build_vitepress = false` in `docs/make.jl`, regenerate docs, and preview with:
```bash
julia --project=docs -e 'using DocumenterVitepress; DocumenterVitepress.dev_docs("docs/build")'
```
- `dev_docs` expects the build directory path `docs/build`, not `docs/build/.documenter`.
- `dev_docs` does not rerun `makedocs` for you. After content changes, rerun `julia --project=docs docs/make.jl`.
- Restore the normal Vitepress build setting before committing.
- Do not commit files generated by `docs/make.jl`, `Literate.markdown/script/notebook`, or `PlutoSliderServer.export_directory` unless explicitly asked.
- In practice, generated docs artifacts to avoid committing include `docs/src/tutorial/*`, `docs/src/tutorial/pluto/**`, and generated tutorial/download assets under `docs/src/public/`.

## PRs
- Base PRs on `master`.
- Run the relevant package tests before opening or updating the PR.
- Docs PRs can use the `documentation` label to trigger preview docs deployment (`docs/make.jl push_preview` in CI).
- GPU-sensitive PRs should get `run-gpu-ci` at PR creation when the first CI run matters. This is for `KomaMRICore/ext/`, kernels, or backend-specific changes. `oneAPI` is experimental and excluded from the default GPU CI and benchmark path.
- `pre-release` triggers Julia pre-release CI; use it only when intentionally checking upcoming Julia compatibility.
- Give reviewers the missing context: what changed, why, and what you tested.

## Canonical References
- Contributor workflow: `docs/src/how-to/5-contribute-to-koma.md`
- CI truth: `.github/workflows/CI.yml` and `.github/workflows/CIPreRelease.yml`
