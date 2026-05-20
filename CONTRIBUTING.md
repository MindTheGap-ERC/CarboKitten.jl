# Contributing to this repository

Welcome, and thank you for contributing to this codebase! This document specifies contribution guidelines.

## General guidelines

### Reporting bugs

Suspect or found a bug? To report it, use the GitHub issues [here](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues), and tag your issue with the "bug" label. Please describe the bug in as much detail as you can, including (1) a description of the unexpected behavior you observed (2) what behavior you expected and (3) (if possible) a minimum running example. The more detailed your bug report is, the easier it is for us to fix, and the faster we will be able to fix it.

### Request features or enhancements

Do you think the codebase is lacking, could use a cool new feature, or should better integrate with existing codebases or packages? Then submit your feature/enhancement request in the [GitHub issues](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues) using the "enhancement" label. Please describe in detail what the new feature should do, and how it should integrate with the existing codebase or other packages. We will review each enhancement request and decide on a case-by-case basis if we will implement it. Our decision is guided both by the usefulness of the request and available development time.

### Improving documentation

Do you think the documentation is lacking? We're always happy to improve! If you want to improve documentation, fork the repository, make your changes, and submit a pull request with the "documentation" label. We will include your improvements into the code base after a review.  
If you think you found a mistake in the documentation, or are unsure whether your modifications are suitable, please open an [issues](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues) with the "documentation" label and describe what you found.

### Contributing code and features

Would you like to contribute a new feature to the code base, or improve existing code? Then fork the repository, add your features, and submit a pull request using the "enhancement" label. We will decide on a case-by-case basis if the pull request is accepted or not. This will be discussed with all active authors and contributers to the code base. Criteria for inclusion of new code are code quality, clarity and utility of the features, and whether the new feature enhances the original idea of the code base. If you are unsure if your feature would fit into the codebase, please use GitHub issues to discuss your idea (using the "enhancement" label) so we can give you feedback.

While this is not a criterion for code inclusion, we strongly encourage that each new feature includes tests that ensure that the feature works as intended and integrates with the existing codebase seamlessly.

To summarize:

1. **issue first** Before contributing a new feature, it is best to always open an issue before opening a PR. This way we can discuss viability, API choices and even implementation details before you put in too much work up front. Don't be afraid to contact us early in the process.
2. **be consistent** While we don't have strict rules on code formatting, it is important to be consistent with the existing style and architecture of CarboKitten. Please read the architecture documentation to learn about some of the choices we made in making CarboKitten modular. This also applies if your contribution involves new data sets. If you are uncertain about coding style and formatting, it may be helpful to use an automatic code formatter like [`JuliaFormatter.jl`](https://juliaeditorsupport.github.io/JuliaFormatter.jl/stable/). Be careful though not to reformat unrelated code, as that would infringe on the next point.
3. **be atomic** CarboKitten is designed to be modular to its core. Make sure that both the issue and the related PR describe and implement a single feature. If your work can be split up in multiple features, please do. You can create sub-issues to indicate that a feature is related to a larger goal. There is no lower limit to the amount of code that can change in a PR.
4. **be clean** In the process of development you may leave pieces of code that are no longer useful, especially when the code has been modified by coding agents. When you submit code for review, take care that the code is readable and free of unwanted distractions.
5. **document** Make sure your new feature is documented including working examples.
6. **test** Include automated tests to check for code correctness (CarboKitten does not have 100% test coverage, but this is something we want to improve, gardening style).

We understand that most of our users and contributors are scientists and not software engineers. We are always willing to guide you in adhering to these criteria.

## Code quality checklist

Though many of these things are checked in CI, please pay attention to the following guide lines:

- Unit tests should pass: run `make test` before pushing to an active PR.
- Documentation should build: run `make docs` before pushing to an active PR.
- Code and documentation should be synchronized: run `entangled sync` before commiting.

## Entangled

This project is written and documented with [Entangled](https://entangled.github.io). See the README for more information.

## Authors and contributors

For a list of authors and contributors please see the README file.
