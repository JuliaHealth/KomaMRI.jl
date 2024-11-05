# Contribute to Koma

If you're interested in contributing to Koma, this document will guide you through configuring everything you need to get started. By contributing, you help enhance the functionality, usability, and performance of the Koma ecosystem. Whether you're fixing bugs, adding features, or improving documentation, your efforts are welcomed because it help us to advancing the project. Before you begin, it's necessary to install and configure a few essential tools on your machine to ensure a smooth development experience:

- **Juliaup:** This is the manager for installing different versions of Julia.
- **Julia:** This is the programming language. It is advisable to install it with Juliaup.
- **Git:** This is a version control system handy for coding.
- **GitHub:** This is a cloud-based Git repository handy for managing the KomaMRI project. You need to create an account.
- **VScode:** This is a code editor with support for development operations.
- **VScode Plugins:** They enable VSCode to have more handy features. We recommend you install the following: 
    - "Julia" for Julia development
    - "GitHub Pull Requests" for collaboration with GitHub.

## Installing KomaMRI as a developer
### Clone KomaMRI repository

To install the dev version of Koma, we will use the Julia REPL:
```julia-repl
pkg> dev KomaMRI
``` 
This command will clone KomaMRI.jl's repository (`dev` version) to your `~/.julia/dev/KomaMRI/` directory.

### Create your Fork of KomaMRI

You need to create a fork of the KomaMRI repository in GitHub. This allows you to have your own copy of the repository to work on without modifying the original repository.

Go to the official [KomaMRI repository](https://github.com/JuliaHealth/KomaMRI.jl) and follow the steps below:

![](../assets/create-fork-step1.png)
![](../assets/create-fork-step2.png)

### Access Your GitHub Account in VSCode

Now, you need to ensure that your GitHub account is connected to VSCode. This allows you to clone repositories, create branches, and manage pull request directrly within VSCode.

- Open VSCode.
- Go to the **Source Control** tab.
- Sign in to your GitHub account if you're not already signed in.

>💡You can also check if your `git` credentials are correctlly added to your machine by writing in the VScode terminal:
```shell
git config --global user.name
git config --global user.email
```

### Open Your Forked Repository in VSCode

In VSCode, click on **File** -> **Open Folder...** and select your `~/.julia/dev/KomaMRI/` directory.

![](../assets/open-folder.png)

Now add the fork URL by clicking **Source Control** -> **...** -> **Remote** -> **Add Remote...**

![](../assets/add-remote.png)

This will create the option to provide a repository URL. Here is where you will paste your fork URL and give it the name `my-fork`.

![](../assets/create-remote-step1.png)
![](../assets/create-remote-step2.png)

>💡Press `Yes` when prompted to constantly fetch in the future.

The Julia extension should automatically detect the `KomaMRI` environment. To check this, look at the status bar (bottom) end you should see `Julia env: KomaMRI`. If this is not the case, click the option in the menu bar and select KomaMRI.jl.

#### KomaMRI monorepo setup

As KomaMRI.jl contains multiple packages in one GitHub repository, you need to specify that you want to use your local copies (instead of the ones available on the Julia registries) with:
```shell
(KomaMRI) pkg> dev ./KomaMRIBase ./KomaMRICore ./KomaMRIFiles ./KomaMRIPlots
```
Finally, use the `instantiate` command to install all the required packages (specified in the `Project.toml`):
```shell
(KomaMRI) pkg> instantiate
```

This will also include all the specific package versions into the `Manifest.toml`. The `Manifest.toml` should not be updated to the repo when making a commit or pull request. Thus, it is present in the `.gitignore`.

### Create a New Branch for your Feature

If you did correctly follow the previous steps you will have correctly created your fork connected to the original Koma repository. Now, if you want to create your own changes, you will need to create a new branch from your fork.

To create this new branch, go to **Source Control** -> **...** -> **Branch** -> **Create Branch form...**

![](../assets/add-branch.png)

This will open a menu to select an starting point for your branch. Select `my-fork/master` as your starting point, and give it the name `my-new-feature`.

![](../assets/create-branch-step1.png)
![](../assets/create-branch-step2.png)

>💡In your VScode terminal use `git status` to check if your branch is correctly created. Your branch should be listed at the top of the output.

## How to Commit

If you have already created your first modifications in your local version of the repository, you will want to commit your changes in your public branch.

To do this, in VScode go to the Source Control panel in the Activity Bar.

Assuming you are currently in your `my-new-feature` branch, the Source Control panel should show your changes to the project and the option to create a commit message.

![](../assets/how-to-commit.png)

If you hove over the `Changes` tab, it should show a `+` icon. Press it to stage all changes in the project.

Write down a message that describes your changes you are stageing to the project, and press the Commit button.

Press Sync Changes to push your commit into your branch.

>💡 If you want to make sure if the commit was correctly done, check your GitHub repository and see if the changes you commited are present.

## How to create a pull request

If you want to send your commited new version of the repository, you can create a pull request that will be reviewed by a Koma certified developer.

To create this pull request, in VScode, go to the `GitHub Pull Request` panel in the Activity Bar and hove over the `Pull request` tab. This should show a Create pull request icon to press.

![](../assets/create-pull-request.png)

In the `Create` tab that appears, select `JuliaHealth/master` as the base and the branch you are working with to merge.

To finish your pull request, give your pull request a description that explains the issue or feature you are addresing in your branch, and press the Create button.

![](../assets/fill-pull-request.png)

>💡 **Tips for a successful Pull Request:**
>   - Try to address one issue or feature per pull request to make it easier for reviewers.
>   - Provide all the context necesary, including all the information of the related issue or added feature.
>   - Respond to feedback and suggestions to make adjustments based on the reviewers comments.