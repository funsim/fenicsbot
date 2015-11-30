# Temporarily use minrk/fenics-dev while waiting for fenicsproject/dev to update
# FROM fenicsproject/dev
FROM minrk/fenics-dev

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install python-matplotlib python-pip
RUN pip install python-twitter pytest
ADD . /srv/fenicsbot
RUN chown -R fenics /srv/fenicsbot
WORKDIR /srv/fenicsbot
RUN bash -c 'source /home/fenics/fenics.conf && py.test tests'
CMD ["/srv/fenicsbot/docker.sh"]
