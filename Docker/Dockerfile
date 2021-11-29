FROM debian
FROM conda/miniconda3


RUN apt-get update && apt-get install -y git wget

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda && rm ~/miniconda.sh && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && echo "conda activate base" >> ~/.bashrc && find /opt/conda/ -follow -type f -name '.a' -delete && find /opt/conda/ -follow -type f -name '.js.map' -delete && /opt/conda/bin/conda clean -afy
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin



#VOLUME /BSATOS
RUN  git clone https://github.com/maypoleflyn/BSATOS.git; 
#RUN cd /root/BSATOS ; chmod +x *; chmod -R +x scripts/*.pl

#RUN cd /BSATOS; conda config --set channel_priority strict; 
RUN cd /BSATOS/Docker; conda env create -f bsatos.yml 
#RUN conda init bash 
#RUN conda activate bsatos

RUN apt-get -y --allow-unauthenticated install r-base
RUN R CMD INSTALL /BSATOS/dependancy/modeest_2.1.tar.gz

RUN echo "source /opt/conda/etc/profile.d/conda.sh ; conda init bash; conda activate bsatos" > ~/.bashrc

RUN chmod 777 /BSATOS/scripts/*

RUN chmod 777 /BSATOS/bsatos

RUN cp -r /BSATOS/* /usr/local/bin


