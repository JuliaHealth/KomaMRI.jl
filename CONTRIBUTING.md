# Contributing Guidelines

Thanks for taking the time to contribute to KomaMRI! ❤️

We appreciate your interest in improving this project. Before you start contributing, please take a moment to review the following guidelines.


## Reporting Issues

First off, we assume that you have read the available [Documentation](https://cncastillo.github.io/KomaMRI.jl).

Before you report an issue, it is best to search for existing [Issues](https://github.com/cncastillo/KomaMRI.jl/issues) that might help you. In case you have found a suitable issue and still need clarification, you can write your question in this issue.

If you then still feel the need to report an issue, we recommend the following:

- Open an [Issue](https://github.com/cncastillo/KomaMRI.jl/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions, depending on what seems relevant.

We will then take care of the issue as soon as possible.


## How to Contribute

1. **Fork KomaMRI Repository**: Fork the [KomaMRI.jl repository](https://github.com/cncastillo/KomaMRI.jl) to your GitHub account to create a personal copy.

2. **Clone Fork to Local Machine**: Clone your fork of the repository to your local machine using the following Git command:

   ```bash
   git clone git@github.com:<your-username>/KomaMRI.jl.git
   ```

3. **Create Contribution Branch**: Create a new branch dedicated to your contribution within the cloned repository.

   ```bash
   git checkout -b feature/your-feature
   ```

4. **Implement Code Changes**: Introduce your modifications, ensuring adherence to the [Julia Blue Style Guidelines](https://github.com/invenia/BlueStyle). For new features, include informative comments, docstrings, and consider enriching the documentation with relevant examples.

5. **Validate Changes with Tests**: Execute existing tests to verify the compatibility of your alterations with the current functionality. If applicable, incorporate additional tests to validate your new contributions.

6. **Commit Code Changes**: Commit your modifications with a precise and descriptive commit message using the following Git command:

   ```bash
   git commit -m "Add a concise summary of your changes"
   ```

7. **Push Changes to Your Fork**: Push your committed changes to the corresponding branch on your forked repository:

   ```bash
   git push origin feature/your-feature
   ```

8. **Initiate Pull Request (PR)**: Propose your changes by creating a pull request against the `master` branch of the KomaMRI repository. Provide a clear title and a detailed description of your modifications. Ensure that your pull request is linked to an issue by including its #ID in the description.

9. **Undergo Code Review**: Subject your code to review by maintainers, who will provide feedback and may request further adjustments before merging.


## License

By contributing to KomaMRI, you agree that your contributions will be licensed under the [MIT License](https://github.com/cncastillo/KomaMRI.jl/blob/master/LICENSE).

Thank you for contributing to KomaMRI! 🌟
