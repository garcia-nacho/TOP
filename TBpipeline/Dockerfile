FROM ghcr.io/garcia-nacho/top_tbpipeline
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

USER root

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
