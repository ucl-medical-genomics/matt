FROM gitpod/workspace-full

RUN brew install R
RUN brew install gcc@5
RUN brew install libxml2
RUN brew install xmlsec1
RUN brew install libtiff
