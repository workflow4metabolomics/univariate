FROM r-base:3.3.1

MAINTAINER Etienne Thevenot (etienne.thevenot@cea.fr)

# Setup package repos
RUN echo "deb http://mirrors.ebi.ac.uk/CRAN/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Update system
RUN apt-get update && apt-get -y install git

# Install R packages needed
RUN R -e "install.packages('batch', dependencies = TRUE)"

# Clone tool
RUN git clone -b docker https://github.com/phnmnl/univariate /files/univariate

# Clean up
RUN apt-get clean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*

# Make wrapper executable
RUN cp /files/univariate/*.R /usr/local/bin/
RUN chmod a+rx /usr/local/bin/univariate_wrapper.R

# Define Entry point script
ENTRYPOINT ["/files/univariate/univariate_wrapper.R"]
