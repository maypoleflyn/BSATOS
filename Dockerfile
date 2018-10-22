#debian distro
FROM debian
RUN echo "deb http://cran.rstudio.com/bin/linux/debian stretch-cran34/" >>  /etc/apt/sources.list
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

#update it
#RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get update
#RUN gpg --keyserver pgp.mit.edu --recv-key --recv-keys FCAE2A0E115C3D8A
#RUN gpg -a --export FCAE2A0E115C3D8A | apt-key add -

#RUN gpg --keyserver keyring.gnupg.net --recv-keys FCAE2A0E115C3D8A
RUN apt-get -y install apt-utils
#RUN gpg --keyserver pgp.mit.edu --recv-key --recv-keys FCAE2A0E115C3D8A
#RUN gpg -a --export FCAE2A0E115C3D8A | apt-key add -


#RUN apt-get -y install apt-utils
#RUN gpg --keyserver keys.gnupg.net --recv-key FCAE2A0E115C3D8A

#install R
RUN apt-get -y --allow-unauthenticated install r-base
#install python
RUN apt-get -y install python
#install git
RUN apt-get -y install git-core
#install autoheader
RUN apt-get -y install autoconf
#install samtools
RUN cd /root ; git clone https://github.com/samtools/htslib.git; cd /root/htslib;  autoheader; autoconf ; ./configure ; make ; make install
RUN cd /root; git clone https://github.com/samtools/samtools.git ; cd /root/samtools/; autoheader ; autoconf -Wno-syntax;  ./configure; make ; make install
#install bedtools 
RUN cd /root;  git clone https://github.com/arq5x/bedtools2.git;  cd /root/bedtools2 ; make    
#install bcftools
RUN cd /root ; git clone https://github.com/samtools/bcftools.git; cd /root/bcftools ;  autoheader ; autoconf -Wno-syntax;  ./configure; make ; make install

#install bwa
RUN cd /root ; git clone https://github.com/lh3/bwa.git; cd /root/bwa ; make

#install HAPCUT2
RUN cd /root ; git clone https://github.com/vibansal/HapCUT2.git; cd /root/HapCUT2/ ; make ; make install-hairs; make install-hapcut2

#install gffread
RUN cd /root; git clone https://github.com/gpertea/gclib; git clone https://github.com/gpertea/gffread ; cd gffread ; make

#install delly
RUN apt-get -y install libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev
RUN cd /root; git clone --recursive https://github.com/dellytools/delly.git ; cd delly ; make all; make install;
#install ANNOVAR
COPY ./src/annovar.latest.tar.gz /root/
RUN cd /root ; tar xvf annovar.latest.tar.gz;


RUN cp /root/bwa/bwa /usr/local/bin/bwa
RUN cp /root/bcftools/bcftools /usr/local/bin/bcftools
RUN cp /root/bedtools2/bin/* /usr/local/bin/
RUN cp /root/samtools/samtools /usr/local/bin/samtools
RUN cp /root/HapCUT2/build/HAPCUT2 /usr/local/bin/ ; cp /root/HapCUT2/build/extractFOSMID /usr/local/bin/ ; cp /root/HapCUT2/build/extractHAIRS /usr/local/bin/
RUN cp /root/gffread/gffread /usr/local/bin/
RUN cp /root/delly/bin/delly /usr/local/bin/
RUN cp /root/annovar/*.pl /usr/local/bin/


#Finally, get BSATOS
COPY ./src/BSATOS.tar.gz /root/
RUN cd /root; tar zxvf BSATOS.tar.gz; cd BSATOS; chmod +x *; chmod -R +x scripts/*.pl
RUN cp -r /root/BSATOS/* /usr/local/bin

#install modeest library 
RUN R CMD INSTALL /root/BSATOS/modeest_2.1.tar.gz 
RUN rm -r /root/*
VOLUME /BSATOS
