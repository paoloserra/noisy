FROM kernsuite/base:3
MAINTAINER bhugo@ska.ac.za

RUN apt-get install -y software-properties-common
RUN apt-add-repository -s ppa:kernsuite/kern-2
RUN apt-add-repository multiverse
RUN apt-add-repository restricted
RUN apt-get update
RUN apt-get install -y python-casacore
RUN apt-get install -y python-pip
RUN pip install pip setuptools wheel -U

ADD noisy /src/noisy
ADD MANIFEST.in /src/MANIFEST.in
ADD requirements.txt /src/requirements.txt
ADD setup.py /src/setup.py
ADD setup.cfg /src/setup.cfg
ADD README.md /src/README.md
ADD bin /src/bin
RUN pip install /src/ -U

ENTRYPOINT ["noisy_predictrms.py"]
CMD ["--help"]
