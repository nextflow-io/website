FROM openjdk:8-alpine

RUN apk update && apk upgrade && apk add bash curl
RUN curl -fsSL https://dl.bintray.com/jbake/binary/jbake-2.3.1-bin.zip > jbake-2.3.1-bin.zip \
 && unzip jbake-2.3.1-bin.zip \
 && rm jbake-2.3.1-bin.zip
 
ENV PATH=/jbake-2.3.1/bin/:$PATH 