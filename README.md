# GebPy

GebPy is a Python-based tool for the generation of geophysical and geochemical 
data of the subsurface. The main assumption of GepPy is that all rock properties 
are determined by the mineral assemblage besides structural features.

### Getting Started

We use GitHub for the cooperative work among the members of our working group. 
By using this tool, we can keep all codes under version control and can easily 
share a stable version of GebPy.
Before a new feature will be added to the stable version of the master branch,
all new features have to be developed and tested in at least one additional 
branch first, before they will be committed to the stable version.  

More detailed information about GitHub with focus on the functionality and the 
workflow can be found [here](https://github.com/features).

### Prerequisites/Installing

Please install the latest version of 
[Git](https://github.com/git-guides/install-git) or keep your
already existing version up to date, before you clone our shared repository by:

```
git clone git@github.com:MABeeskow/GebPy.git
```

The path where you carry out this command will now include the git repository as 
a new folder. You need to navigate into this new folder to execute any git 
commands. 

### Creating a new branch

As a member of this repository, a new branch can easiliy created within the 
browser environment. Otherwise, a new branch can also be created by using the
following command:

```
git checkout -b name_of_your_branch
```

You can use your new branch as a personal version-control when new features will 
be developed by you. Please always make sure that you are working on this 
particular branch before any changes will be added to the master branch! New 
features have to be tested sufficiently to ensure a stable working version of 
GebPy on the master branch.

### Workflow
#### Updating your branch

It can be useful to update your branch with the stable version of GebPy on the
master branch, but you can of course update your branch with any other branch
of the reposity. Before you do that, you should be sure that you have committed 
all changes of your branch to the local version. You can check this by using:

```
git status
```

Now you can update your branch with the master branch by:

```
git pull origin master
```

This command will merge any changes between your local branch and the latest 
version of the development branch.

#### Adding new features to your branch

If you have made any changes within your branch, for example the creation of a 
new file, you have to add this file to the index by using the following command:

```
git add .
```

By using `git status`, all new or modified files are now green colored. Before 
these files were added, they appeared red colored.  
The next step is to commit these changes to your local version of the current
branch. This will be done by using:

```
git commit -m "additional information"
```

It can be useful, especially if conflicts have to be solved, if the comment 
function was used.  
The last step is to update the version of your branch with the local changes. 
For this purpose, the following command has to be used:

```
git push origin name_of_your_branch
```

#### Adding new features to the master branch

This step will take place if new features that were sufficiently tested, can
now be added to the stable version of GebPy on the master branch. The workflow
is similar to the chapter before but it ends now with the (additional) command

```
git push origin master
```
