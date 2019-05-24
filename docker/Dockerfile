FROM perl:5.28

# install perl dependencies
RUN cpanm JSON::XS Inline::C Bio::Perl Test::More Data::Compare

# copy the code
ADD . /app/nmer-match

# copy & define the entrypoint
COPY ./docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]

# define default arguments
#CMD []