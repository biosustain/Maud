Contributing to Maud
====================

This page explains how you can contribute to Maud.


## Reporting bugs and making feature requests

If you find a bug in Maud or would like to request a new feature, consider
adding to the [issues list](https://github.com/biosustain/Maud/issues). Bugs
should get the red `bug` label and features the turquoise
`enhancement` label.


## Contributing code changes

If you want to add a new feature or fix a bug yourself, the first step is to
clone the code. If you have ssh access to the `Maud github repository
<https://github.com/biosustain/Maud>`_, you can clone it by running


```sh
git clone git@github.com:biosustain/Maud.git
```

Next, check out a new branch

```sh
git checkout -b descriptive_branch_name
```

When you are happy with your changes, commit them to your new branch and then
push it to github


```sh
git commit -m "Short description of the changes I've made"
git push
```

Finally, use the online github interface to make a new pull request from your
branch, complete with a longer description.

At least one approving review is required before a pull request can be
merged. You can make your reviewers' lives a lot easier by adding as much
detailed commentary as possible, referring to existing issues where appropriate
and keeping your pull requests small and logically distinct.


## Running tests locally

To run the tests locally, run the command

```sh
tox
```

## Releasing Maud

To release a new version of Maud, edit the field `version` in the file
`pyproject.toml`, then make a pull request with this change.

Once the changes are merged into the `origin/master` branch, add a tag
with the new version number prefixed by "v", e.g. `v0.2.1`. You can do
this using the github ui, or from the command line like this:


```sh
git tag v0.2.1
git push origin "v0.2.1"
```

Now go to Maud's [releases](https://github.com/biosustain/Maud/releases)
page and click "draft a new release". Choose your tag from the menu, then add
a description and click "publish release".
