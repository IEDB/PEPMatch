FROM perl:5.30.1

# install sqlite
RUN apt update && apt -y install sqlite3

# install perl dependencies
# the dependencies where we skip tests need to be installed
# first
COPY requirements.txt .
RUN cat requirements.txt | xargs -I {} cpanm {}

# copy the code
ADD . /app/nmer-match
WORKDIR /app/nmer-match

# copy & define the entrypoint
COPY ./docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]