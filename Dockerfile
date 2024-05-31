# Rocker base image
FROM rocker/rstudio

# Set work dir
WORKDIR /code

COPY . .

RUN chmod +x /code/scripts/batch_correct.R

RUN apt clean && apt-get update


# Install required system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libz-dev\
    libfontconfig1-dev\
    libgit2-dev\
    libharfbuzz-dev\ 
    libfribidi-dev\
    libtiff-dev\ 
    libfreetype6-dev\ 
    libpng-dev\ 
    libjpeg-dev\ 
    libbz2-dev\
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN Rscript -e "install.packages(c('data.table','cowplot','tidyverse','configr','optparse','emdist','systemfonts','BiocManager','ggridges','outliers','uwot', 'plyr','dplyr','dotenv' ), repos = 'https://cloud.r-project.org/',Ncpus=16)"

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site
RUN R -e "install.packages('devtools', dependencies = TRUE)"
RUN R -e "setRepositories(ind = c(1:6, 8))"
RUN R -e "BiocManager::install(c('flowCore', 'sva', 'Biobase'), ask = FALSE)"
RUN R -e "devtools::install_github('biosurf/cyCombine')"

ENTRYPOINT [ "Rscript", "/code/scripts/batch_correct.R" ]