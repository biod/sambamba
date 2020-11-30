---
name: Bug report 🐞
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''
---
**Only bug reports!**

The D version of Sambamba is in maintenance mode. Use the github issue
tracker to report bugs *only*. For comments, questions and features,
please use the google group mailing list as stated on the
[README](https://github.com/biod/sambamba/blob/master/README.md)!

**Describe the bug**

A clear and concise description of what the bug is.

**To Reproduce**

Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error
**Expected behavior**
A clear and concise description of what you expected to happen.
**Screenshots**
If applicable, add screenshots to help explain your problem.
**Desktop (please complete the following information):**
- Browser [e.g. chrome, safari]
- Version [e.g. 22]
**Additional context**
Add any other context about the problem here.

Include a set of BAM/BED files to reproduce the issue

+ bonus points if you try to minimize the test case yourself, as issues are often localized:
  - try to use sambamba slice to first extract the reference where the error occurs
  - if that succeeds (the error is still reproducible), continue to crop the file in binary-search fashion
