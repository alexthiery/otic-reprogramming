---
title: Linux
category: Cheat Sheets
order: 2
---

# Linux 101

## Git
- `git add <file1>` or `git add <file1> <file2> <file3>` - Adds one or several files to the staging area
- `git add .` - Adds all the files inside the project folder to the staging area
- `git branch` - Visualize the branches of the repository and which one is active
- `git branch -d <branch-name>` - Deletes the specified branch
- `git checkout <branch-name>` - Switches to this branch. Future commits will be made on the checked-out branch
- `git checkout -b <new-branch-name>` - Creates new branch and change to it at the same time
- `git clone <URL>` - Locally clone a repository from github with the repo URL
- `git commit -m <'commit message'>` - Commits changes to the local branch
- `git init` - Initialize a git repository for the working directory
- `git log` - Displays all the commits that were made for the current project
- `git merge <branch-name>` - Merges current branch into the specified branch.
- `git push` - Pushes commits to the remote branch
- `git submodule update --remote` - Updates a submodule to the latest commit


## Screen
- `screen -S nextflow -m bash -c 'sh run.sh; exec sh'` - Run named screen session with wait command at end so window doesnt auto-close
- `screen -ls` - List screens
- `screen -r 344074` - Attach to running screen with number before decimal point
- `screen -X -S [session # you want to kill] quit` - kill session
- Screen `ctrl+a+d` detatch

## macos
- `sudo nano /etc/paths` - Edit path on mac 

## Nextflow
- `nextflow run <REPO>` - Runs pipeline from github
- `nextflow run <REPO> -with-docker <CONTAINER>:<TAG>` - Runs pipeline from github with specified docker image
- `nextflow pull <REPO>` - Updates repo from github
- `nextflow run <SCRIPT> -w <WORKDIR>` - Specify work dir
- `nextflow run <SCRIPT> -ansi-log false` - Change log to each task (good for debugging)
- `singularity pull  --name <NAME>.img docker://<PATH>:<TAG> > /<TAG>/null` - Pull custom docker image
- `rm -rf ~/.nextflow/assets/<REPO-PATH>` - Removes downloaded repo if switching branches etc
- `rm -rf .nextflow* results work` - Clean a working directory if having problems
