# KomaMRI.jl

## Style
- Prioritize non-verbosity and information density.
- Prioritize elegance and simplicity: concise idiomatic Julia, simple data expressions over stateful loops, no abstract fields, no speculative abstractions.
- Type function arguments only when dispatch needs them.
- Do not use `isa`/`typeof` type checks for behavior that can be expressed with dispatch or small trait functions.
- Use `!` only for functions that mutate an argument. Internal mutating helpers should return `nothing`; value-returning helpers should not use `!`.
- Use low-level constructors like `new{typeof(...)}` only in the minimal inner-constructor boundary that actually needs them.
- Add helper functions only when reused or semantically justified.
- Reuse existing APIs before adding private wrappers.
- Add docstrings only for public-facing functions.
- Prefer domain names over implementation-mechanic names and explanatory comments. Add comments only when they add information.
- Make semantic phases visible in code structure.
- Be terse. Assume expert user.
- Prefer commands/diffs over long explanations.
- Keep diffs tight. Do not reformat or clean up unrelated code.
- Reduce lines by simplifying design, not by compressing code.

## Git
- If the user gives an existing branch/worktree, use it after verifying it is the intended one.
- If starting feature/PR work yourself, create a unique branch from latest `origin/master` in a separate worktree; keep the main checkout on `master`.
- If the active checkout is detached, create a unique local branch there before editing.
- Fresh worktree pattern: `git fetch origin && git worktree prune && git worktree add -b codex/<task> ../KomaMRI.jl-<task> origin/master`.
- If a branch/path exists, choose a new name. Do not force-remove unless asked.
- Never push unless the latest user message explicitly asks. Never push to `master`.
- When creating pull requests, open them ready for review unless explicitly asked for a draft PR.
- Do not prefix PR titles with `[codex]` or other author/tool tags unless explicitly requested.
- No destructive Git commands without explicit approval.
- Keep generated scratch artifacts under `.tmp/`; never stage or push `.tmp/` contents.

## Repo
- Monorepo packages: root `KomaMRI`, plus `KomaMRIBase`, `KomaMRICore`, `KomaMRIFiles`, `KomaMRIPlots`.
- `KomaMRICore`, `KomaMRIFiles`, and `KomaMRIPlots` depend on `KomaMRIBase`; the root package depends on those three.
- Each package has its own `Project.toml`; tests use `*/test/Project.toml`; docs use `docs/Project.toml`.
- Root `[workspace]` includes `test`, `docs`, `benchmarks`, and all subpackages.
- Do not edit `Manifest.toml` directly or commit manifest churn unless explicitly asked.

## Julia
- Use the relevant project: `julia --project=<path>` and `Pkg.activate(...)`; never use the global env by accident.
- First action for Julia work: use or start the persistent Julia REPL; do not run Julia through `exec_command` unless isolation is technically required and stated first.
- Use one persistent Julia REPL/session started with `--threads=auto` and Revise for all Julia work, including diagnostics, plotting, scratch scripts, examples, and tests. Do not run `julia -e`, `julia script.jl`, or package tests in fresh shell processes while a threaded REPL is available; send commands to the session. If code is too large to paste safely, write it under `.tmp/` and run `include("...")` from the existing REPL. Use a fresh process only when isolation is technically required or the user asks. Restart only if absent/crashed, incorrectly threaded, or the user asks. Keep it open.
- On Julia 1.12+, Revise can handle struct redefinitions; do not restart only because a struct changed.
- Prefer workspace setup: activate root or the child project and `Pkg.instantiate()`. Use explicit `Pkg.develop(path=...)` only to reproduce CI or older Julia 1.10 wiring.
- Change dependencies/compat with Pkg APIs, not casual `Project.toml` edits.
- Prefer command-line arguments for script options. Use environment variables only for actual environment/CI/backend semantics such as CPU/GPU/CUDA/Metal selection.
- For performance work: profile first, benchmark with interpolation, then optimize the actual hotspot.
- Never write raw `seq += ...` in examples or generated code; use `@addblock` or `@addblocks`.
- Never call `build_*` only to extract events or duration and then rebuild the same block. Use `make_*` for custom blocks, copy a built block when preserving its block semantics, or append the built sequence/block directly.

## Python
- Use `uv` for reproducible Python environments. Do not use bare `pip`.

## Testing
- In the persistent Julia session, run tests from the correct project with `Pkg.test(...)`; use the narrowest package first.
- Usual package tests: `Pkg.test("KomaMRIBase")`, `Pkg.test("KomaMRICore")`, `Pkg.test("KomaMRIFiles")`, `Pkg.test("KomaMRIPlots")`, root `Pkg.test()`.
- Root and `KomaMRIPlots` tests may need `xvfb-run` on headless Linux.
- `KomaMRICore` backend tests use `test_args` (`CPU`, `CUDA`, `AMDGPU`, `Metal`, `oneAPI`) or test preferences; do not commit backend GPU deps to shared test projects.
- Use `@testitem` tags: `:base`, `:files`, `:plots`, `:koma`, and for core `:core` plus exactly one of `:motion` or `:nomotion`.
- If testing a PR implementation, add it with Pkg `url`/`rev`; do not hand-copy files.
- Use representative semantic variable names or formulas for expected values; avoid unexplained magic numbers in assertions.
- For groups of semantically similar tests, add one compact comment above the group when intent is not obvious. Do not comment every assertion.
- Do not add tests that merely check a function equals its own definition or reimplement the same logic in the test. Test behavioral contracts, regressions, edge cases, and cross-implementation parity; do not add random assertions only to increase coverage.

## Docs
- Use the `docs` environment.
- Edit source docs, `lit-*.jl`, and `pluto-*.jl`; do not edit generated tutorial artifacts.
- In the persistent Julia session activated for `docs`, build with `include("docs/make.jl")`; doctest with `using Documenter: doctest; using KomaMRI; doctest(KomaMRI)`.
- With this repo's `MarkdownVitepress` defaults, `docs/make.jl` runs the full VitePress build; the final local site is `docs/build/1`. Preview final docs by serving that directory, as DocumenterVitepress recommends. Use one server, reuse it when possible, and stop it when done.
- The npm scripts in `docs/package.json` run VitePress on `docs/build/.documenter`, which is the generated VitePress source. Use `npm run docs:dev`, `npm run docs:build`, or `npm run docs:preview` only when intentionally using the VitePress workflow; do not treat `.documenter` as the final static site and do not serve `docs/build/1` through VitePress.
- If VitePress components such as tabs are involved, verify the final `docs/build/1` page is interactive and the referenced `/assets/app.*.js` returns 200.
- For Literate docs, hide implementation-only setup and plotting plumbing with `#hide`; show the result the reader should learn from. Avoid leaving a bare final variable such as `p` visible when it only exists to render a plot.
- Prefer examples that generate their figures from source data/code. Keep complex plotting code hidden unless the plotting code itself is the lesson.
- For interactive docs figures, prefer KomaMRI/KomaMRIPlots APIs and rendered interactive HTML. Do not replace them with static assets or low-level Plotly trace construction unless the user asks or there is no higher-level API.
- Keep `DocumenterVitepress.deploydocs(...)` in `docs/make.jl`.
- Do not commit files generated by `docs/make.jl`, Literate, or PlutoSliderServer unless explicitly asked.

## PRs And Releases
- PRs target `master`, stay scoped, and include what changed, why, and what was tested.
- Add trigger labels at PR creation when the first CI run matters: `documentation`, `run-gpu-ci`, or `pre-release`. Adding them after PR creation is too late for the first run.
- Manual fallback release tags are annotated tags on current `origin/master`, not feature branches: `git fetch origin master --tags`, `git tag -a <tag> -F notes.md origin/master`, `git push origin <tag>`, `gh release create <tag> --title <tag> --notes-file notes.md`.
- Register package versions only after the release PR is merged and `origin/master` contains the target package/version. Do not trigger JuliaRegistrator from an open PR.
- For subpackage registration, comment on the merged `origin/master` commit, not the PR: `@JuliaRegistrator register subdir=<PackageName>`. PR comments do not trigger Registrator in this repo.
- Registration comes before tagging/releasing: wait for the corresponding `JuliaRegistries/General` PR to merge.
- JuliaTagBot creates tags and GitHub releases for the root package and subpackages after General registration merges; do not create manual tags/releases unless TagBot fails or the user explicitly asks.
- Registrator notes must mention breaking status: include `## Breaking changes` or state `No breaking changes`.

## References
- Contributor workflow: `docs/src/how-to/5-contribute-to-koma.md`
- Pulseq MATLAB to Koma translation: `docs/src/how-to/pulseq-matlab-to-koma.md`
- CI truth: `.github/workflows/CI.yml` and `.github/workflows/CIPreRelease.yml`
