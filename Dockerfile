FROM perl:5.28

# install perl dependencies
RUN cpanm JSON::XS Inline::C Bio::Perl

# copy the code
ADD . /app/nmer-match

# define the entrypoint
ENTRYPOINT perl /app/nmer-match/bin/run_nmer_match.pl