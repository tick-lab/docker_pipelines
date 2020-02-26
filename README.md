# Docker Pipelines

Once you have docker installed on your machine you can run a docker container by running two commands (the first of which only needs to be run once):

```
sudo docker build -t mann/dada2
```

This installs all of the necessary packages/libraries/etc contained in the Dockerfile and corresponding script. To run the pipeline:

```
sudo docker run -it --rm -v ~/lab_members/mann/docker_test/01_data:/01_data -v ~/lab_members/mann/docker_test/03_output:/03_output mann/dada2
```
