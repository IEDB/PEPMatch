FROM perl:5.28

# Build from this directory with:
# docker build -f Dockerfile ../
# or from the root directory with:
# docker build -f docker/Dockerfile .

# copy the code
ADD . /app/nmer-match
WORKDIR /app/nmer-match

# install perl dependencies
# the dependencies where we skip tests need to be installed
# first
RUN cat requirements-notest.txt | xargs -I {} cpanm -n {}
RUN cat requirements.txt | xargs -I {} cpanm {}

# copy & define the entrypoint
COPY ./docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]