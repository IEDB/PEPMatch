# Nmer Match

## Docker container

To run the script from the docker container, you'll need to give it access to the directory where you have your input files with the -v option, e.g.:
```bash
docker run -v $PWD:/scratch -w /scratch IMAGE_ID -a build -l 15 -s test_data/human.fasta -c catalogs/humb1
```