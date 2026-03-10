# Use the official RStudio image as a base
# download pre-build rocker image that comes with R studio server installed.  
# I am using a fixed version of R (4.5.1), instead of latest, for reproducibility reasons (i.e. to build always same version).
FROM rocker/rstudio:4.5.1  

# Install system dependencies for python
# the -y option gives yes automatically to any response, otherwise it defaults to no and installation fails.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    libgsl-dev \
    cmake \
    build-essential \
    libpq-dev \
    libgeos-dev \
    libproj-dev \
    libgdal-dev \
    libudunits2-dev \
    libabsl-dev \
    tcl-dev \
    tk-dev \
    libglpk40 \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R packages, other possible packages: tidyverse
RUN R -e "install.packages(c('caret','ggplot2','stringr','parallel','randomForest','concaveman','sp','sf','stats','grDevices','doParallel','nnet','sm','pracma','patchwork','ggExtra','RColorBrewer','dplyr','data.tree','visNetwork','klaR','remotes','devtools','tinytex'))"

RUN R -e "devtools::install_github('semontante/flowMagic')"

RUN Rscript -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('flowCore', 'flowWorkspace','Biobase', 'CytoML', 'sva'), update = TRUE, ask = FALSE)"

RUN R -e "tinytex::install_tinytex()"


# Use legacy auth behaviour to avoid timeouts and keep session alive
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf && \
    echo "auth-stay-signed-in-days=30" >> /etc/rstudio/rserver.conf

# Set up a working directory
WORKDIR /home/rstudio

# Set username
ENV USER=rstudio
# Set password
ENV PASSWORD=1234

# Expose the RStudio port
EXPOSE 8787

# Set the entrypoint for RStudio
# /init is a initiation script included in the rocker image
CMD ["/init"]



