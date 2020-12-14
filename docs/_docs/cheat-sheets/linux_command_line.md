---
title: Linux command line
category: Cheat Sheets
order: 1
---

**Note:** Don't write Bash scripts, only one-liners! (Use Python for scripting and Nextflow for running tools.)

See the full Bash manual [here](https://www.gnu.org/software/bash/manual/bash.html).

## General
- `echo $?` - Print **exit status** of the last command: 0 if the last command exited successfully and a non-zero value otherwise 
- `command1 | command2` - Execute `command1` and redirect its `stdout` into `stdin` of `command2`
- `command1 && command2` - Execute `command1` and then execute `command2` only if `command1` exited successfully
- `command1; command2` - Execute `command1` and then `command2` independently of the success of `command1`
- `command > file` - Redirect `stdout` of `command` into `file`; the `file` gets overwritten
- `command >> file` - Add the output of `command` to the end of `file`
- `command 2> file` - Redirect `stderr` to `file`
- `command 2>&1` - Redirect `stderr` to `stdout`
- `command 2> /dev/null` - Remove `stderr`
- `<(command2)` - **Process substitution**: treat `stdout` of `command2` as a file. **Note:** There should be no space between `<` and `(`. 
- `diff <(ls dir1) <(ls dir2)` - Example of process substitution: find differences between lists of files in `dir1` and `dir2` (`diff` needs two files as arguments)
- `command1 | tee file | command2` - `tee` splits `stdout` of `command1` and writes one copy into `file` while feeding the other copy into `stdin` of `command2`
- `{ command1; command2; ...; commandN; }` - **Command grouping**: `stdout` streams of `command1`, `command2`, ..., `commandN` are collected and can be redirected together. The semicolon after the last command and spaces after `{` and before `}` are mandatory.
- `{ echo -n "a"; echo "b"; } > f.txt` - Example of command grouping: print `ab` into `f.txt`. 
- `PATH=PATH:/path/to/add` - Add `/path/to/add` to the value of the `PATH` variable

## Misc
- `whoami` - Print your current username
- `watch command` - Execute `command` every 2 seconds (it is useful if you need to watch how the output changes)
- `watch -n 4 squeue -u yourlogin` - Execute `squeue -u yourlogin` every 4 seconds to see how your processes scheduled by slurm are getting on

## Operations with directories and files
- `pwd` - Show current path
- `cd dir` - Change directory to `dir`
- `cd ..` - Change to parent directory
- `cd` - Change to the home directory
- `rm file1 [file2 ... fileN]` - Remove `file1`, `file2`, ..., `fileN`
- `rm -r dir` - Remove directory `dir` with all its contents
- `rm -f file1 [file2 ... fileN]` - Remove `file1`, `file2`, ..., `fileN` ignoring any nonexisting files and suppressing prompts from `rm`
- `ls` - Lists contents of the current directory
- `ls -a` - List _all_ contents of the current directory, including files/directories whose names begin with `.` (for example, `.bashrc`)
- `ls -l` - List contents of the current directory in a long format (include permissions, owner, modification date, etc.)
- `ls -lh` - List contents of the current directory in a long format representing file sizes in KBs, MBs, GBs, instead of bytes
- `ls -lt` - List contents of the current directory in a long format and put most recently modified files first
- `ls -alth | head -5` - List 5 most recently modified files in a long format with human-readable file sizes (in KBs, MBs or GBs)
- `tree` - Show the directory hierarchy as a tree structure (useful only for small hierarchies)
- `mv path/to/source path/to/target` - Renames or moves a file from source to target
- `touch file` - Creates a new empty file called `file`
- `cat > file` - Creates a file called `file` by asking user to type its contents. Press `Ctrl+D` when you are finished
- `nano file` - Edit a text file called `file` with a simple text editor called Nano. It is useful for editing small files. You can use [Visual Studio Code](../../Coding/VS_Code) for editing big files, like scripts or TSV/CSV tables
- `vim file` Edit a text file called `file` with an advanced text editor called Vim (see the [Vim manual](https://www.vim.org/docs.php) for details)
- `mkdir newdir` - Create directory `newdir`
- `mkdir -p path/to/newdir` - Create directory `newdir`, along with its parent directories (`path/to/`) if they do not yet exist
- `ln -s file linkname` - Create a symbolic link called `linkname` to the `file`
- `du -h -d1` - List top level directory space usage
- `cp path/to/file another/path` - Copy `file` from `path/to` to `another/path`
- `cp -r path/to/dir another/path` - Copy directory `dir` from `path/to` to `another path`
- `cp -a path/to/dir path/to/file another/path` - Copy `dir` and `file` from `path/to` to `another/path` and preserve all links
- `cp -s path/to/file another/path` - Instead of copying, make a symbolic link `another/path/file` to `path/to/file`
- `cp -u path/to/file another/path` - Copy `file` from `path/to` to `another/path` only if there is no `file` in `another/path` or `file` in `path/to` is newer
- `find ./myfolder -mindepth 1 ! -regex '^./myfolder/test2\(/.*\)?' -delete` - this will delete all folders inside ./myfolder except that ./myfolder/test2 and all its contents will be preserved
- stat

## Viewing and comparing files
- `tail -f` - Monitor a log file
- head
- less
- more
- wc
- grep
- diff
- comm

## Editing and summarising files
- awk (Note: AWK a language, not just a command, but please don't write scripts in AWK, only short one-liners)
- sed (Note: sed a language, not just a command, but please don't write scripts in sed, only short one-liners)
- tr
- sort
- uniq
- cut
- paste
- shuffle
