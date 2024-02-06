FROM admiralenola/niph_tb_pipeline
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

USER root
RUN echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
		r-base \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c( 'seqinr','ggplot2',  'reshape2',  'readxl','writexl',  'data.table'),repos = 'http://cran.us.r-project.org')"
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
RUN rm /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline.py
RUN rm /opt/niph_tb_pipeline/setup.py
COPY Tomnt/Reference /mnt/Reference
COPY Tomnt/Latex_template /mnt/Latex_template
COPY Tomnt/FigTree_v1.4.3 /mnt/FigTree_v1.4.3
COPY Scripts /home/tbuser/Code
COPY WHO-TB_Catalog-2021.7-eng.xlsx /home/tbuser/WHO-TB_Catalog-2021.7-eng.xlsx
COPY niph_tb_pipeline/niph_tb_pipeline.py /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline.py
COPY niph_tb_pipeline/niph_tb_pipeline_part1.py /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline_part1.py
COPY niph_tb_pipeline/niph_tb_pipeline_part2.py /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline_part2.py
COPY niph_tb_pipeline/setup.py /opt/niph_tb_pipeline/setup.py

RUN chmod -R +rx /mnt/Latex_template/*
RUN chmod -R +rx /mnt/FigTree_v1.4.3/*
RUN chmod -R +rx /home/tbuser/Code/*
RUN chmod +r /home/tbuser/WHO-TB_Catalog-2021.7-eng.xlsx
RUN chmod +x /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline.py
RUN chmod +x /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline_part1.py
RUN chmod +x /opt/niph_tb_pipeline/niph_tb_pipeline/niph_tb_pipeline_part2.py
RUN chmod +x /opt/niph_tb_pipeline/setup.py

RUN cd /opt/niph_tb_pipeline && pip install --upgrade .

USER tbuser
WORKDIR /data
CMD ["bash", "-c", "niph_tb_pipeline"]
