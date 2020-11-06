# Workflow

## Create a branch

All your development should happen in branches. Do not develop in master or devel, you will bring yourself into troubles.

To create a branch:

1. Make sure you do not have any local changes: `git status`
1. Branch off of devel: `git checkout -b <branch name> devel`

   Follow these rules for __naming your branch__:

   1. Prefix the branch name with your user name, so the form is like `<username>/<topic>`. This will prevent name collisions with other users.
   1. Use lower-case letters.
   1. Dashes (`-`) work better than underscores (`_`).


+Useful commands:+

* To get a list of branches: `git branch -a`
* To see your active branch: `git branch`

## Creating a commit

In git, you have to stage your changes first. Then you check in those changes locally and it will form a commit (or a patch):

* To add (stage) a new file: `git add <file>`
* To add (stage) modification on a file: `git add <file>`
* To add all new files and local changes: `git add -A`
* To add only modified files: `git add -u`

* To move a file: `git mv <source> <destination>`
* To remove a file: `git rm <file>`

* To see the status (staged, upstages changes and non-tracked files): `git status`

* To commit locally: `git commit`

    The commit message should have this format:
    ```
Short description (up to 70 chars)
<empty line>
Detailed description. List as much useful and relevant information as possible.
This can have multiple lines and even paragraphs.
<empty line>
{Refs | Closes} #XYZ.
    ```

    The detailed description will help a) the reviewer, b) in future when we forget what we did and will need to look it up. Even though the detailed description is optional, it is *highly* recommended that you fill it in.

+Rule of thumb:+ You want your commits to be as small as possible, so they are easier (and therefore faster) to review.  Stay on one topic per commit. Try +not+ to do multiple things in one commit. Do rather a series of patches than one big patch.

+Useful commands:+

* To see what files have changed: `git status`
* To look at patches: `git log`
* To see the code changes: `git diff`

## Sending a pull/merge request

When you think you are finished with your branch, send a pull request for review.

1. `git push upstream <branch name>`
1. Go to [branches](https://hpcgitlab.inl.gov/idaholab/nitrogen/branches).
1. Find your branch and click `+ Merge request`.
1. It will show you an overview of your merge request. If you are happy with it, click on `Submit Merge Request`.

At this point, someone will review your branch and possibly leave comments on what to fix. Otherwise your branch will be merged. You get a notification in both cases. If your branch was merged, you should go update your `devel` branch by pulling from the upstream
```
git checkout devel
git pull upstream devel
```
and you can delete your local branch with:
```
git branch -d <branch name>
```

## Fixing your branch

The easiest way how to fix your branch is by using fixup commits (see [fun with autosquash](http://technosorcery.net/blog/2010/02/07/fun-with-the-upcoming-1-7-release-of-git-rebase---interactive---autosquash/)).

1. Go patch by patch and fix the code according to the reviewer comments.
1. Stage your changes. Typically `git add -u` will do.
1. Do `git commit`

     Type `fixup! <first line of the patch you are fixing>`

     Alternatively, you can use the aliases mentioned in the article above which will make this a lot easier.

1. Continue until all is fixed.
1. Now, do `git rebase --interactive --autosquash master`. This will open an editor, just save and quit the file you see and let git do the hard work.
1. If all went right, you have your fixed branch. You can look at it by: `git log -p`.
1. Push your fixed up branch: `git push -f upstream <branch name>`

## Updating your branch

If you wish to use the latest `devel` branch (assuming it changed since you started your branch). You have to rebase it on top of it, i.e.:
```
<update devel branch>
git checkout <my-branch-name>
git rebase devel
```

Note: You will most likely run into conflicts that you will have to resolve, thus do this only if you know how to do it.


## Reviewing a merge request

* Go to [Merge request list](https://hpcgitlab.inl.gov/idaholab/nitrogen/merge_requests)
* Find an open merge request that is not assigned to any one
* Assign the merge request to yourself so others know you are working on it
* You can see the commits and on the `Changes` tab you can see the complete diff.
* If you want to grab the branch and try it out on your machine
  * Fetch the branch: `git fetch upstream <branch name>:<branch name>`, where `<branch name>` is shown on the merge request page
  * Switch to the branch: `git checkout <branch name>`
  * Now you can try to compile the changes, and do other work with the branch.
* If all is good, go back to the merge request page.
  * If you do not have merge rights, tag the merge request with tag `Approved`.
  * If you have merge rights, click `Accept Merge Request`.
* If you think there should be changes made to the MR, add comments to each line where you think the change should be made. When you are done reviewing, tag the merge request with `Address comments`.

## Collaborating on a shared branch

It is often useful to allow others access to your branch and vice versa. The following section will describe how to setup remote branches and the associated workflow. The branches will be shared in the `upstream` repository.

To **share a branch**, you just need to push it into upstream like you would do for a merge request:
```
git push -u upstream <branch-name>
```
Notice the `-u` which will set your branch to track the remote branch you just created by the push command.


To **get someone else's branch**:
```
git fetch upstream
git branch <somebody's/branch> upstream/<somebody's/branch>
git checkout <somebody's/branch>
```

You and the other users can now push and pull to this branch in the normal fashion (e.g. if you are currently in the branch, execute `git pull` to pull the upstream commits the other users have pushed).

+IMPORTANT+ When working on shared branch, +do not+ change the history on the branch, i.e. do not `rebase` the branch, do not `--amend` on published commits, etc. If you do so, git will not let you push into the branch. Doing a force-push is a big NO-NO, it will hose up your collaborators!


Occasionally it is necessary to incorporate features from the most current version of ```devel``` in your remote branch. To avoid mucking up the history and other collaborators on the branch, one should merge ```devel``` in the following manner:

```
git fetch upstream
git checkout <branch>
git merge --no-ff upstream/devel
```
