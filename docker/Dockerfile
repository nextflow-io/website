FROM openjdk:8-jre

RUN curl -fsSL https://github.com/jbake-org/jbake/releases/download/v2.3.1/jbake-2.3.1-bin.zip > jbake-2.3.1-bin.zip \
 && unzip jbake-2.3.1-bin.zip \
 && rm jbake-2.3.1-bin.zip
 
ENV PATH=/jbake-2.3.1/bin/:$PATH 