## Basic usage instructions


### Get:

```
docker pull immcantation/lab:devel
```

### Run:

Start RStudio at http://localhost:8787/ by executing the command below. You
must set a password with `-e`. You can give the image a name with `--name`.


```
docker run --network=host \
	--name MyImmcantationProject \
	--detach \
	-e PASSWORD=MyPassword \
	--rm -p 8787:8787 immcantation/lab:devel
```

The user will be `magus`and the password the one you provided (`MyPassword` in the example).

## Save:

If you made changes and want to save the new container status, follow these steps. 

With the container running, in a new terminal, find out the container id with:

```
docker ps
``` 

Then commit the changes to immcantation/lab:<tag>:

```
# <tag> can be devel if you want to overwrite the devel image
docker commit <container id> immcantation/lab:<tag>

# list images
docker image ls
``` 


## Stop:

You can stop a named container with this command:

```
docker stop <name>
```

## Work with notebooks

Bind the local notebooks' dir to the containers' notebooks dir and make `notebooks` the working directory. Open RStudio in the browser and start working with the notebook.

```
# Example. How to edit the Immcantation training notebooks.
docker run --network=host \
	--name MyImmcantationProject \
	--detach \
	-v <Path-to-immcantation-training>:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	-e PASSWORD=MyPassword \
	--rm -p 8787:8787 immcantation/lab:devel
```

## Render a notebook

```
docker run --network=host \
	--name MyImmcantationProject \
	-v <Path-to-immcantation-training>:/home/magus/notebooks:z \
	-e PASSWORD=MyPassword \
	--rm -p 8787:8787 immcantation/lab:devel Rscript -e "rmarkdown::render('notebook.Rmd')"
```