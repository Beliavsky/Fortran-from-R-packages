# Contributing to creditr

Thanks for taking the time to contribute! Below our some guidelines that might be helpful. 

## Have a question?

A good place to start is going through the [vignette](https://github.com/yanrong-stacy-song/creditr/blob/main/vignettes/creditr.pdf). 
We have taken care to document several use cases for our package in detail. If your question is still not answered, please open an issue.

## Have a bug report or feature request?

1. Search current open GitHub [issues](https://github.com/yanrong-stacy-song/creditr/issues)
   to check if the bug/feature has already been reported/requested.
2. Additionally, for bug reports, make sure your local copy is up-to-date, and verify the bug still
   exists in the HEAD of the main branch.
3. Open an issue with a minimal, reproducible example of the bug you wish to report. If you are requesting a feature, please detail a relevant use case.

## Want to submit a pull request?

1. Unless the change is extremely trivial (e.g. typos), please
   create an issue and wait for
   our feedback *before* you start work on a pull request. That will avoid the
   possibility you spend time on a patch that won't be merged.
2. Create a branch for the feature/bug fix reported in the issue. Please use a
   short and descriptive branch name that starts with the issue number (e.g.
   123_custom_function). Use that branch as the base for your pull request.
3. Please write detailed commit messages.
4. It would be much appreciated if you also add tests that cover your changes.

## Commit Messages

Follow the [The Seven Rules](http://chris.beams.io/posts/git-commit/#seven-rules)
of [How to Write a Git Commit Message](http://chris.beams.io/posts/git-commit/).
Pay particular attention to [rule 7: Use the body to explain what and why
versus how](http://chris.beams.io/posts/git-commit/#why-not-how). The body
should also include the motivation for the change and how it compares to prior
behavior.

If the commit is to fix a bug or add a feature, the commit message should
contain enough information to understand the bug/feature without having to
reference an external tracker (e.g. GitHub issues). But please *do reference
the GitHub issue* on the last line of your commit message body.