## Basic usage instructions


### Get:

```
docker pull immcantation/studio:devel
```

### Run:

Start RStudio at http://localhost:8787/ by executing the command below. You
must set a password with `-e`. You can give the imaga a name with `--name`.


```
docker run --network=host \
	--name MyImmcantationProject \
	--detach \
	-e PASSWORD=MyPassword \
	--rm -p 8787:8787 immcantation/studio:devel
```

The user will be `magus`and the password the one you provided (`MyPassword`).

## Save:

If you made changes and want to save new container status. With the container running,
in a terminal find out the container id with:

```
docker ps
``` 

Then commit the changes to immcantation/studio:<tag>:

```
# <tag> can be devel if you want to overwrite the devel image
docker commit <container id> immcantation/studio:<tag>

# list images
docker image ls
``` 


## Stop:

```
docker stop MyImmcantationProject
```


