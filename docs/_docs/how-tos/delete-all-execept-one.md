---
title: Delete All files/folders except one folder in a directory
category: How-to's
order: 1
---

This will delete all folders inside ./myfolder except that ./myfolder/test2 and all its contents will be preserved:

```
find ./myfolder -mindepth 1 ! -regex '^./myfolder/test2\(/.*\)?' -delete
```

## How it works
- `find` starts a find command.
- `./myfolder` tells find to start with the directory `./myfolder` and its contents.
- `-mindepth 1` not to match `./myfolder` itself, just the files and directories under it.
- `! -regex '^./myfolder/test2\(/.*\)?'` tells find to exclude (!) any file or directory matching the regular expression `^./myfolder/test2\(/.*\)?.` `^` matches the start of the path name. The expression `(/.*\)?` matches either (a) a slash followed by anything or (b) nothing at all.
- `delete` tells find to delete the matching (that is, non-excluded) files.

Consider a directory structure that looks like;

```
$ find ./myfolder
./myfolder
./myfolder/test1
./myfolder/test1/dir1
./myfolder/test1/dir1/test2
./myfolder/test1/dir1/test2/file4
./myfolder/test1/file1
./myfolder/test3
./myfolder/test3/file3
./myfolder/test2
./myfolder/test2/file2
./myfolder/test2/dir2
```

We can run the find command (without -delete) to see what it matches:

```
$ find ./myfolder -mindepth 1 ! -regex '^./myfolder/test2\(/.*\)?'
./myfolder/test1
./myfolder/test1/dir1
./myfolder/test1/dir1/test2
./myfolder/test1/dir1/test2/file4
./myfolder/test1/file1
./myfolder/test3
./myfolder/test3/file3
```

We can verify that this worked by looking at the files which remain:

```
$ find ./myfolder
./myfolder
./myfolder/test2
./myfolder/test2/file2
./myfolder/test2/dir2
```