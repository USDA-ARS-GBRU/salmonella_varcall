FROM alpine:3.16.0

RUN apk update
RUN apk add wget openssl python3


# create working directories
RUN mkdir /home/SentericaVariantPipeline/
RUN mkdir /home/programs/

# set the working directory
WORKDIR /home/SentericaVariantPipeline/
