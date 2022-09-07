FROM gitpod/workspace-full

RUN brew install R
RUN brew install gcc@5
# For XML install as part of MATT
RUN brew install libxml2
RUN brew install xmlsec1
# For ragg install as part of devtools
RUN brew install libtiff
