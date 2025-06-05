# Contribute to Koma

If you're interested in contributing to Koma, this document will guide you through configuring everything you need to get started. By contributing, you help enhance the functionality, usability, and performance of the Koma ecosystem. Your efforts are welcomed because it help us to advancing the project. Before you begin, it's necessary to install and configure a few essential tools on your machine to ensure a smooth development experience:

- Git
- GitHub
- VScode
- VScode Plugins: 
    - Julia 
    - GitHub Pull Requests

## Installing KomaMRI as a developer
### 1. Clone KomaMRI repository

To install the dev version of Koma, we will use the Julia REPL:
```julia
pkg> dev KomaMRI
``` 
This command will clone KomaMRI.jl's repository (`dev` version) to your `~/.julia/dev/KomaMRI/` directory if you are in a MacOS or Linux operative system, or `C:\Users\<user-name>\.julia\dev\KomaMRI\` if you are using Windows, where `<user-name>` should be replaced with your Windows user.

### 2. Create your fork of KomaMRI

If you try to commit or generate a pull request at this point, you will get an `Access denied` error. This is because you need to create a fork before you can contribute to this repository directly (unless you are included as a collaborator!).

To create this fork, go to the official [KomaMRI repository](https://github.com/JuliaHealth/KomaMRI.jl) and follow the steps below:

![](../assets/create-fork-step1.png)
![](../assets/create-fork-step2.png)

### 3. Access your GitHub account in VSCode

Now, you need to ensure that your GitHub account is connected to VSCode. This allows you to clone repositories, create branches, and manage pull request directrly within VSCode.

- Open VSCode.
- Go to the **Source Control** tab.
- Sign in to your GitHub account if you're not already signed in.

>💡You can also check if your `git` credentials are correctlly added to your machine by writing in the VScode terminal:
>```shell
>git config --global user.name
>git config --global user.email
>```

### 4. Open your forked repository in VSCode

In VSCode, click on **File** -> **Open Folder...** and select your `~/.julia/dev/KomaMRI/` directory (`C:\Users\<user-name>\.julia\dev\KomaMRI\` if you are using Windows).

Now add the fork URL by clicking **Source Control** -> **...** -> **Remote** -> **Add Remote...**

```@raw html
    <img width="80%" src="../../assets/add-remote.png">
```
This will create the option to provide a repository URL. Here is where you will paste your fork URL and give it the name `my-fork`.

![](../assets/create-remote-step1.png)
![](../assets/create-remote-step2.png)

>💡Press `Yes` when prompted to constantly fetch in the future.

The Julia extension should automatically detect the `KomaMRI` environment. To check this, look at the status bar (bottom) end you should see `Julia env: KomaMRI`. If this is not the case, click the option in the menu bar and select KomaMRI.jl.

### 5. KomaMRI monorepo setup

As KomaMRI.jl contains multiple packages in one GitHub repository, you need to specify that you want to use your local copies (instead of the ones available on the Julia registries) and using the `instantiate` command to install all the required packages (specified in `Project.toml`) with the following script:

```julia
using Pkg  
# Koma sub-packages dev setup  
koma_subpkgs = ["KomaMRICore", "KomaMRIFiles", "KomaMRIPlots"]  
for pkg in koma_subpkgs  
    Pkg.activate(pkg)  
    Pkg.develop(path = "./KomaMRIBase")  
end  
# Main package (KomaMRI) dev setup  
Pkg.activate(".")  
for pkg in koma_subpkgs  
    Pkg.develop(path = "./$pkg")  
end
Pkg.instantiate()
```
In case you want to contribute specifically in documentation, you will need to use the `docs` enviroment with the following script:

```julia
Pkg.activate("docs")
Pkg.develop(path = ".")
Pkg.instantiate()
```

This will also include all the specific package versions into the `Manifest.toml`. The `Manifest.toml` should not be updated to the repo when making a commit or pull request. Thus, it is present in the `.gitignore`.

### 6. Create a new branch for your feature

If you did correctly follow the previous steps you will have correctly created your fork connected to the original Koma repository. Now, if you want to create your own changes, you will need to create a new branch from your fork.

To create this new branch, go to **Source Control** -> **...** -> **Branch** -> **Create Branch form...**

```@raw html
    <img width="80%" src="../../assets/add-branch.png">
```
This will open a menu to select an starting point for your branch. Select `my-fork/master` as your starting point, and give it the name `my-new-feature`.

![](../assets/create-branch-step1.png)
![](../assets/create-branch-step2.png)

>💡In your VScode terminal use `git status` to check if your branch is correctly created. Your branch should be listed at the top of the output.

## How to commit

If you have already created your first modifications in your local version of the repository, you will want to commit your changes in your public branch.

To do this, in VScode go to the Source Control panel in the Activity Bar.

Assuming you are currently in your `my-new-feature` branch, the Source Control panel should show your changes to the project and the option to create a commit message.

```@raw html
    <img width="40%" src="../../assets/how-to-commit.png">
```
If you hove over the `Changes` tab, it should show a `+` icon. Press it to stage all changes in the project.

Write down a message that describes your changes you are stageing to the project, and press the Commit button.

Press Sync Changes to push your commit into your branch.

>💡 If you want to make sure if the commit was correctly done, check your GitHub repository and see if the changes you commited are present.

## How to Test Your Contributions

In Koma, we have tests to verify the correctness of the current state of the code. Depending on the package where you made your changes, `KomaMRIBase`, `KomaMRICore`, `KomaMRIFiles`, `KomaMRIPlots` or `KomaMRI`, follows the instructions below to test them correctly:

### Test `KomaMRI`:

In the Julia REPL open the Julia package manager mode by pressing `]`, then run the following script:

```julia
activate .
test
```

This should open the UI and all buttons will be clicked to test their functionality.

For all of these packages tests look the same. There are two options to run the tests:


- **Using VSCode**: On the activity bar, open the `Testing` extension, expand the available tests, and select the "▶" icon next to the respective package to run the test. The results will be displayed in the `Test Results` panel.

- **Using the Julia REPL**: Open the Julia package manager mode by pressing `]`, then run the following script, with `[package]` being the selected package to be tested (`KomaMRIBase`, `KomaMRIPlots`, `KomaMRIFiles`):

```julia
activate [package]
test [package]
```

    
### Test `KomaMRICore`:
In this package, you may want to run tests using the CPU or a GPU. By default the tests will run on the CPU, with the number of threads set to `Threads.nthreads()`. You can run them using the Julia REPL (`] activate ; test KomaMRICore`) or VSCode.  To run KomaMRICore's tests, on the activity bar, open the `Testing` extension, expand the available tests, and select the "▶" icon next to the word `KomaMRICore` to run the test. The results will be displayed in the `Test Results` panel.


Using the Julia REPL, by default, tests are run on the CPU with the number of threads set to `Threads.nthreads()`. To run on a specific GPU backend, add the name of the backend package ("AMDGPU", "CUDA", "Metal", or "oneAPI") to the `test/Project.toml` file in `KomaMRICore` and pass the name as a test argument.

Example:
```julia
import Pkg
Pkg.test("KomaMRICore"; test_args=["CUDA"])
```
To run on the CPU with a specific number of threads, pass the number of threads as a Julia argument.

Example:
```julia
import Pkg
Pkg.test("KomaMRICore"; julia_args=`--threads=4`)
```    
To change the default backend used for testing, modify the `[preferences.KomaMRICore]` section in the test/Project.toml file:

```julia
[preferences.KomaMRICore]
test_backend = "CPU"
```

For the backend preference to take effect, you need to:

- **REPL Testing**: No action needed. `] test` should pick up the preference immediately.
- **VSCode Testing**: You need to restart VSCode.

>Sadly, `LocalPreferences.toml` files are not picked up by VSCode (they could be `.gitignore`'d), so we put them into the `test/Project.toml` file instead.

If your contributions do not affect the correct execution of the code, the tests will return a message indicating that your changes have successfully passed.

### Adding a new test

In case your contribution generates a new function that is not currently tested, the code coverage will decrease when the pull request is analyzed, creating an automatic comment indicating the problem. To include this contribution in the tests, you must incorporate your function into the runtest.jl file corresponding to the package where you made your contribution.

To add this new test, you must follow these steps:

1. In VSCode, press `Ctrl + p` to open the search bar and type `runtest.jl` to select the desired file according to the package.
2. Once the file is selected, you will see the tests distributed in `@testitem` that contain `@testset`. Check if the contribution to test fulfills the conditions for an existing `@testitem`. If not, create a new `@testitem`.
3. Open a `@testset` inside the `@testitem`, and create a small proof of concept that validates the contribution (use the surrounding `@testset` as a guide).
4. To finish your test, add a `@test` macro followed by the verification of the proof of concept. This macro will pass the test if the boolean value of the result is true.

## How to create a pull request

If you want to send your commited new version of the repository, you can create a pull request that will be reviewed by a Koma certified developer.

To create this pull request, in VScode, go to the `GitHub Pull Request` panel in the Activity Bar and hove over the `Pull request` tab. This should show a Create pull request icon to press.

```@raw html
    <img width="50%" src="../../assets/create-pull-request.png">
```

In the `Create` tab that appears, select `JuliaHealth/master` as the base and the branch you are working with to merge.

To finish your pull request, give it a name with a clear mention of the  subject of the contribution you made, and a description that explains the issue or feature you are addresing in your branch, and press the Create button.

```@raw html
    <img width="50%" src="../../assets/fill-pull-request.png">
```
>💡 **Tips for a successful Pull Request:**
>   - Try to address one issue or feature per pull request to make it easier for reviewers.
>   - Provide all the context necesary, including all the information of the related issue or added feature.
>   - Respond to feedback and suggestions to make adjustments based on the reviewers comments.
