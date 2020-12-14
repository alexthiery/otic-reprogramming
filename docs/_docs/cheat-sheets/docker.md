---
title: Docker
category: Cheat Sheets
order: 4
---

## Resources

All official documentation: [https://docs.docker.com/reference/](https://docs.docker.com/reference/)

Command line documentation: [https://docs.docker.com/engine/reference/commandline/cli/](https://docs.docker.com/engine/reference/commandline/cli/)

Dockerfile format documentation: [https://docs.docker.com/engine/reference/builder/](https://docs.docker.com/engine/reference/builder/)

Unofficial Dockerfile cheatsheet: [https://kapeli.com/cheat_sheets/Dockerfile.docset/Contents/Resources/Documents/index](https://kapeli.com/cheat_sheets/Dockerfile.docset/Contents/Resources/Documents/index)

## Build

- `docker build -t <CONTAINER>:<TAG> <FOLDER CONTAINING DOCKERFILE>` -- build new image from docker file

## Run

- `docker run -it --rm <CONTAINER>:<TAG> /bin/bash` -- Run a container with an interactive shell

## Share

- `docker pull <CONTAINER>:<TAG>` -- Pulls image from docker hub to local repo
- `docker push <CONTAINER>:<TAG>` -- Pushes local image to docker hub

## Maintain

- `docker images` -- View local images
- `docker system prune` -- Docker provides a single command that will clean up any resources — images, containers, volumes, and networks — that are dangling (not associated with a container)
- `docker image prune` -- By default, docker image prune only cleans up dangling images. A dangling image is one that is not tagged and is not referenced by any container.
- `docker image prune -a` -- To remove all images which are not used by existing containers, use the -a flag:
- `docker rmi $(docker images -a -q)` - Remove all docker images
- `docker rm $(docker ps -a -f status=exited -f status=created -q)` - Remove all exited containers
