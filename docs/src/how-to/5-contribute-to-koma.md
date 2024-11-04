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

*"create fork steps" IMAGES*
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

*"OPEN FOLDER" IMAGE*

Now add the fork URL by clicking **Source Control** -> **...** -> **Remote** -> **Add Remote...**

*"Add Remote..." IMAGE*

This will create the option to provide a repository URL. Here is where you will paste your fork URL and give it the name `my-fork`.

*"Fork url and name" IMAGES*

>💡Press `Yes` when prompted to constantly fetch in the future.

### Create a New Branch for your Feature

If you did correctly follow the previous steps you will have correctly created your fork connected to the original Koma repository. Now, if you want to create your own changes, you will need to create a new branch from your fork.

To create this new branch, go to **Source Control** -> **...** -> **Branch** -> **Create Branch form...**

*"Branch from..." IMAGE*

This will open a menu to select an starting point for your branch. Select `my-fork/master` as your starting point, and give it the name `my-new-feature`.

*"select my-fork/master and name" IMAGES*

>💡In your VScode terminal use `git status` to check if your branch is correctly created. Your branch should be listed at the top of the output.

## How to Commit

## How to create a pull request