# Rocker base image
FROM rocker/rstudio

# Set work dir
WORKDIR /app

RUN apt clean && apt-get update

# Copy files into image
COPY . .