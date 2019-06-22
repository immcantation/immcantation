# Immcantation training materials

### Introduction to B cell repertoire analysis 

Get an global overview of how the different tools in the Immcantation framework work together with a jupyter notebook based on the materials presented in the [webinar](https://immcantation.eventbrite.com)

To use with [Docker](https://www.docker.com/) locally:

1. Pull the Immcantation Lab container image:

    ```
    # Example: pull release version 2.7.0-lab
    docker pull kleinstein/immcantation:2.7.0-lab
    ```
    
1. Run the container:

    ```
    docker run --network=host -it --rm -p 8888:8888 kleinstein/immcantation:2.7.0-lab
    ```

    Or, if you want to save the results in your computer:
    
    ```
    # Note: change my-out-dir for the full path to the local directory where 
    # you want to have the results saved to
    docker run --network=host -it --rm -v my-out-dir:/home/magus/notebooks/results:z -p 8888:8888 kleinstein/immcantation:2.7.0-lab
    ```
    
    Once the container is running, You will see a message asking you to visit a url like `http://<hostname>:8888/?token=<token>`

1. Open your internet browser and visit that url

    # Example: http://localhost:8888/?token=18303237b2521e72f00685e4fdf754f955f82a958a8e57ec

1. Use CTRL+Enter to execute the commands inside the code cells
