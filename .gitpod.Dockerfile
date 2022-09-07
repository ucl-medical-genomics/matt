FROM gitpod/workspace-full

RUN brew install R
RUN brew install gcc@5
RUN brew install libxml2
RUN brew install xmlsec1

#ENV R_LIBS_USER "/home/linuxbrew/.linuxbrew/lib/R/%v/site-library/"
#RUN sudo apt-get install -y libxml2-dev
#RUN sudo apt-get install -y r-cran-xml
#ENV PATH="/home/linuxbrew/.linuxbrew/Cellar/libxml2/2.10.2/bin:${PATH}"