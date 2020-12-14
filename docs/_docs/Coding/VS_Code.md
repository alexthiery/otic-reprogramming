---
title: Visual Studio Code
category: Coding
order: 2
---

[Visual Studio Code (or VS Code)](https://code.visualstudio.com/) is an open-source text editor released by Microsoft and 
[becoming increasingly popular](https://insights.stackoverflow.com/survey/2019#development-environments-and-tools) among software developers. 
Thanks to a huge number of handy extensions, VS Code is a great option for anything, from Python or Nextflow scripts 
to Markdown to Dockerfiles to tables in CSV/TSV format. Also, VS Code is tightly integrated with [Git](../../reproducibility/git.md) and allows a developer to use multiple command line terminals. If you would like to try VS Code (which we really recommend!), then please read further.

**Note.** For R we recommend using [RStudio](https://rstudio.com/), for data analysis in Python - [Jupyter Notebook or Jupyter Lab](https://jupyter.org/), 
for a big software development project in Python [PyCharm](https://www.jetbrains.com/pycharm/) may be the best environment; 
finally, for coding in Java [IntelliJ IDEA](https://www.jetbrains.com/idea/) is a standard _de facto_. 

To begin to use VS Code:

1. Install VS Code from [its official website](https://code.visualstudio.com/).

2. Choose a color theme you like using the Color theme section of the Welcome page.

3. If you are moving to VS Code from another text editor in which you have worked for quite a while, you may want to use familiar keybindings. If so, check out the Settings and keybindings section of the Welcome page. If you want to use native VS Code keybindings, please check out its cheatsheet, available from the Welcome page as a PDF (see Help, Printable keyboard cheatsheet). You can also easily add your own keybindings (see paragraph 7 below).

4. Checkout the Interactive Playground from the Welcome page to learn some useful ways of efficient code editing in VS Code.

5. Install the following useful extensions to benefit from the start:

- [Nextflow](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) (Nextflow syntax support)

- [Python](https://marketplace.visualstudio.com/items?itemName=ms-python.python) (Python syntax support)

- [Markdown Shortcuts](https://marketplace.visualstudio.com/items?itemName=mdickin.markdown-shortcuts) (Markdown support)

- [markdownlint](https://marketplace.visualstudio.com/items?itemName=DavidAnson.vscode-markdownlint) (Markdown style checking)

- [Docker](https://marketplace.visualstudio.com/items?itemName=ms-azuretools.vscode-docker) (Dockerfile syntax support and more)

- [Gitconfig Syntax](https://marketplace.visualstudio.com/items?itemName=sidneys1.gitconfig) (Syntax support for git config files)

- [gitignore](https://marketplace.visualstudio.com/items?itemName=codezombiech.gitignore) (Syntax support for .gitignore)

- [Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh) (Support for coding on a remote machine)

- [Edit csv](https://marketplace.visualstudio.com/items?itemName=janisdd.vscode-edit-csv) (Edit CSV files as an Excel-style spreadsheet)

- [WordCounter](https://marketplace.visualstudio.com/items?itemName=kirozen.wordcounter) (Count words as you write; it is useful when writing abstracts)

6. Install and try some useful commands:

- `code` allows you to open a text file from the command line in VS Code. For example, `code file.txt`. See [this post on Stackoverflow](https://stackoverflow.com/questions/29955500/code-not-working-in-command-line-for-visual-studio-code-on-osx-mac) on how to install this command.

- `open` allows you to open a viewer for files in different formats (e. g., PDF or PNG) from the command line in VS Code. For example, `open file.pdf`. You don't need to install the command - it should be available by default.

7. It is convenient to be able to open a terminal next to your code, instead of having a terminal at the bottom of the screen, and to smoothly move from code to the terminal and back, as well as between multiple terminals. To set the necessary shortcuts, open the Keyboard Shortcut editor (`Command + K + S`) and then click on "Open Keyboard Shortcuts (JSON)" button (a page with a bent arrow).

**Note**: Be careful if you switched to keybindings of another editor (see paragraph 3 above): the shortcuts proposed below may overlap with the ones you switched to.

Paste the following shortcut descriptions into the opened JSON file:

```
[
    { "key": "ctrl+alt+`", "command": "workbench.action.positionPanelRight" },
    { "key": "ctrl+`", "command": "workbench.action.terminal.focus"},
    { "key": "ctrl+`", "command": "workbench.action.focusActiveEditorGroup", "when": "terminalFocus" },
    { "key": "cmd+shift+j", "command": "workbench.action.terminal.focusNext" },
    { "key": "cmd+shift+k", "command": "workbench.action.terminal.focusPrevious" },
    { "key": "cmd+shift+w", "command": "workbench.action.terminal.kill" }
]
```

Save the changes and close the file and the Keyboard Shortcut editor.

Now you can mix these shortcuts with the default ones in the following way:

- Create a new file (and hence, a new tab) with `Command + N` (a default shortcut).

- Create a new terminal next to your editing area with `` Command + Option + ` ``.

- Create more terminals with `` Shift + Ctrl + ` `` (a default shortcut).

- Move to the next terminal with `Command + Shift + J` (when you reach the last terminal, you jump to the first one).

- Move to the previous terminal with `Command + Shift + K` (when you reach the first terminal, you jump to the last one).

- Remove the current terminal with `Command + Shift + W`.

- Move back and forth between the current tab and the current terminal with `` Ctrl + ` ``.

- Cycle through the editor tabs with `Command + Option + RightArrow` (a default shortcut).

- Close the current tab with `Command + W` (a default shortcut).

You can learn more about setting VS Code shortcuts from [the VS Code documentation](https://code.visualstudio.com/docs/getstarted/keybindings). Also, have a look a the cheatsheet of the default VS Code keybindings (see paragraph 3 above).

8. You can work with Git straight from VS Code. Please see [VS Code documentation about using Git](https://code.visualstudio.com/docs/editor/versioncontrol). Of course, instead, you can still use Git from the command line or a dedicated graphical software like [GitHub Desktop](https://desktop.github.com/) or [GitKraken](https://www.gitkraken.com/). For details, please see our pages about [Git](../../reproducibility/git) and [GitHub](../../reproducibility/luslab-github).

You are now well equipped to use VS Code! Happy coding!
