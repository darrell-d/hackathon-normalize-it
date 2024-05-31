# Use the existing image as the base image
FROM pennsieveci/normalize_it:latest

# Copy files into image
COPY . .

ENTRYPOINT [ "Rscript", "/code/scripts/batch_correct.R" ]