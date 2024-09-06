# Release workflow

The release workflow is mostly automated using continuous integration builds,
but some actions and triggering a release needs to be done manually by the developer.

To create a release:

1. Create a branch and edit the `CHANGELOG.md` and `CITATION.cff` files to update it with a new version number.
2. Create a new PR from the branch. The PR must have `RELEASE` in the title.
3. This will trigger a build with multiple versions of python.
4. Review the branch, check that all tests for all python versions pass and if so merge.
5. Once merged, the CI should create a github release in "Draft" mode.
6. Check that the release page is correct (has the wheel and `mltbx` files and the release notes are ok).
7. Check that the wheel and `mltbx` toolbox can be installed and work.
8. Then manually trigger the `Publish to PyPI` action to upload the wheel to PyPI.


In particular, in step 1:

* in `CHANGELOG.md` the first title line must have the form:

```
# [<newver>](https://github.com/spinw/spinw/compare/<oldver>...<newver>)
```

* in `CITATION.cff` the `version` field must be updated with a new version

If the version string in these two files do not match, or if the version string matches an existing git tag,
then the CI build will fail.

Also note that in step 8, after uploading to PyPI the release cannot be changed on PyPI (only deleted).
If a release is deleted, you have to then create a new release version (PyPI does not allow overwriting previous releases).
