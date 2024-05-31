# Use the existing image as the base image
FROM pennsieveci/normalize_it:latest

WORKDIR /code

# Copy files into image
COPY . .
RUN ls /code
RUN ls -la /code/scripts

ENTRYPOINT [ "Rscript", "/code/scripts/batch_correct.R" ]
