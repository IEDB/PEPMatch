FROM perl:5.30.1

# install perl dependencies
# the dependencies where we skip tests need to be installed
# first
COPY requirements-notest.txt .
COPY requirements.txt .
RUN cat requirements-notest.txt | xargs -I {} cpanm -n {}
RUN cat requirements.txt | xargs -I {} cpanm {}

# copy the code
ADD . /app/nmer-match
WORKDIR /app/nmer-match

# copy & define the entrypoint
COPY ./docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]