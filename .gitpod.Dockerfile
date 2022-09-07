FROM gitpod/workspace-full

RUN brew install R
RUN brew install gcc@5

ENV R_LIBS_USER "/home/linuxbrew/.linuxbrew/R/%v/site-library"
RUN sudo apt-get install -y r-cran-xml
